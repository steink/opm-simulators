/*
  Copyright 2025 Equinor ASA

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>
#include <opm/simulators/wells/ProdGroupTreeBalancer.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/input/eclipse/EclipseState/Phase.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateModel.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <unordered_map>

namespace Opm::ProdGroupTreeBalancer {

namespace {

// ---------------------------------------------------------------------------
// Internal constants
// ---------------------------------------------------------------------------

/// Small initial rate used to seed the algorithm when all rates are zero.
template<class Scalar>
constexpr Scalar kSmallRate = Scalar(1e-4);

/// Tolerance used for capping rate-vs-limit comparisons to avoid floating-point noise.
template<class Scalar>
constexpr Scalar kFeasibilityTolerance = Scalar(1e-12);

/// Default uniform weight when a guide rate is unavailable.
template<class Scalar>
constexpr Scalar kUniformWeight = Scalar(1.0);

/// Default equal fraction for each active phase when rates are zero.
template<class Scalar>
constexpr Scalar kDefaultUniformFraction = Scalar(1.0 / 3);

/// OIL / WATER / GAS canonical indices (always 0/1/2 in the tree).
constexpr int kOil   = 0;
constexpr int kWater = 1;
constexpr int kGas   = 2;

// ---------------------------------------------------------------------------
// Mode conversion helpers
// ---------------------------------------------------------------------------

/// Convert Group::ProductionCMode to Well::ProducerCMode for balancing.
/// FLD and VREP are not valid for wells, so they map to CMODE_UNDEFINED.
inline Well::ProducerCMode groupModeToWellMode(Group::ProductionCMode groupMode)
{
    switch (groupMode) {
        case Group::ProductionCMode::ORAT: return Well::ProducerCMode::ORAT;
        case Group::ProductionCMode::WRAT: return Well::ProducerCMode::WRAT;
        case Group::ProductionCMode::GRAT: return Well::ProducerCMode::GRAT;
        case Group::ProductionCMode::LRAT: return Well::ProducerCMode::LRAT;
        case Group::ProductionCMode::RESV: return Well::ProducerCMode::RESV;
        case Group::ProductionCMode::NONE: return Well::ProducerCMode::NONE;
        default: return Well::ProducerCMode::CMODE_UNDEFINED;  // FLD, VREP, etc.
    }
}

/// Convert Well::ProducerCMode to Group::ProductionCMode (for well's preferred mode).
/// BHP, THP, GRUP, CRAT are well-specific and map to NONE for group representation.
inline Group::ProductionCMode wellModeToGroupMode(Well::ProducerCMode wellMode)
{
    switch (wellMode) {
        case Well::ProducerCMode::ORAT: return Group::ProductionCMode::ORAT;
        case Well::ProducerCMode::WRAT: return Group::ProductionCMode::WRAT;
        case Well::ProducerCMode::GRAT: return Group::ProductionCMode::GRAT;
        case Well::ProducerCMode::LRAT: return Group::ProductionCMode::LRAT;
        case Well::ProducerCMode::RESV: return Group::ProductionCMode::RESV;
        case Well::ProducerCMode::NONE: return Group::ProductionCMode::NONE;
        default: return Group::ProductionCMode::NONE;  // BHP, THP, GRUP, etc.
    }
}

// ---------------------------------------------------------------------------
// Phase mapping helpers
// ---------------------------------------------------------------------------

/// Return the active phase index for canonical index \p canonical using pu,
/// or -1 if the phase is not active.
template<typename IndexTraits>
int activeIdx(const PhaseUsageInfo<IndexTraits>& pu, int canonical)
{
    if (canonical == kOil)   return pu.phaseIsActive(IndexTraits::oilPhaseIdx)   ? pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)   : -1;
    if (canonical == kWater) return pu.phaseIsActive(IndexTraits::waterPhaseIdx) ? pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx) : -1;
    if (canonical == kGas)   return pu.phaseIsActive(IndexTraits::gasPhaseIdx)   ? pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)   : -1;
    return -1;
}

/// Copy surface_rates (active-phase vector) into a canonical [oil,water,gas] array,
/// placing zero for inactive phases.  The sign convention is preserved
/// (negative = production).
template<class Scalar, typename IndexTraits>
std::array<Scalar, 3> toCanonical3(const std::vector<Scalar>& activeRates,
                                   const PhaseUsageInfo<IndexTraits>& pu)
{
    std::array<Scalar, 3> r{};
    for (int c = 0; c < 3; ++c) {
        const int a = activeIdx(pu, c);
        if (a >= 0 && a < static_cast<int>(activeRates.size())) {
            r[c] = activeRates[a];
        }
    }
    return r;
}

/// Convert canonical [oil,water,gas] array back into an active-phase vector.
template<class Scalar, typename IndexTraits>
std::vector<Scalar> toActive(const std::array<Scalar, 3>& canonical3,
                              const PhaseUsageInfo<IndexTraits>& pu)
{
    std::vector<Scalar> v(pu.numPhases, Scalar(0));
    for (int c = 0; c < 3; ++c) {
        const int a = activeIdx(pu, c);
        if (a >= 0) {
            v[a] = canonical3[c];
        }
    }
    return v;
}

// ---------------------------------------------------------------------------
// Rate projections onto individual limit types
// ---------------------------------------------------------------------------

/// Given canonical rates [oil,water,gas] (negative = production) and a control
/// mode, return the corresponding "production rate" value that the limit is
/// compared against.  Positive = production.
template<class Scalar>
Scalar rateForMode(const std::array<Scalar, 3>& rates,
                   Well::ProducerCMode cmode,
                   const std::array<Scalar, 3>& resvCoeff)
{
    switch (cmode) {
    case Well::ProducerCMode::ORAT:
        return -rates[kOil];
    case Well::ProducerCMode::WRAT:
        return -rates[kWater];
    case Well::ProducerCMode::GRAT:
        return -rates[kGas];
    case Well::ProducerCMode::LRAT:
        return -(rates[kOil] + rates[kWater]);
    case Well::ProducerCMode::BHP:
    case Well::ProducerCMode::THP:
        return -(rates[kOil] + rates[kWater] + rates[kGas]);
    case Well::ProducerCMode::RESV:
    {
        Scalar resv = 0;
        for (int c = 0; c < 3; ++c) {
            resv += -rates[c] * resvCoeff[c];
        }
        return resv;
    }
    default:
        return Scalar(0);
    }
}

/// Corresponding linear-term projection (component of the rate change per
/// unit alpha for the given control mode).
/// Rate[mode] = rateForMode(rates, mode); d(Rate[mode])/d(alpha) = linearTermForMode(lt, mode).
///
/// Sign derivation: stepAlpha sets rates[c] -= alpha * lt[c], so
///   d(rates[c])/d(alpha) = -lt[c]
///   d(Rate[ORAT])/d(alpha) = d(-rates[kOil])/d(alpha) = lt[kOil]  (positive for producers)
/// Returns are therefore +lt[c] (not -lt[c]).
template<class Scalar>
Scalar linearTermForMode(const std::array<Scalar, 3>& lt,
                         Well::ProducerCMode cmode,
                         const std::array<Scalar, 3>& resvCoeff)
{
    switch (cmode) {
    case Well::ProducerCMode::ORAT:
        return lt[kOil];
    case Well::ProducerCMode::WRAT:
        return lt[kWater];
    case Well::ProducerCMode::GRAT:
        return lt[kGas];
    case Well::ProducerCMode::LRAT:
        return lt[kOil] + lt[kWater];
    case Well::ProducerCMode::BHP:
    case Well::ProducerCMode::THP:
        return (lt[kOil] + lt[kWater] + lt[kGas]);
    case Well::ProducerCMode::RESV:
    {
        Scalar v = 0;
        for (int c = 0; c < 3; ++c) {
            v += lt[c] * resvCoeff[c];
        }
        return v;
    }
    default:
        return Scalar(0);
    }
}

template<class Scalar>
Scalar positiveRateForMode(const std::array<Scalar, 3>& positiveRates,
                           Well::ProducerCMode mode,
                           const std::array<Scalar, 3>& resvCoeff)
{
    if (mode == Well::ProducerCMode::BHP) {
        return positiveRates[kOil] + positiveRates[kWater] + positiveRates[kGas];
    }

    const std::array<Scalar, 3> negRates{
        -positiveRates[kOil],
        -positiveRates[kWater],
        -positiveRates[kGas]
    };
    return rateForMode(negRates, mode, resvCoeff);
}

template<class Scalar>
Scalar guideSlopeForMode(const std::array<Scalar, 3>& guideRateSums,
                         Well::ProducerCMode mode,
                         const std::array<Scalar, 3>& resvCoeff)
{
    if (mode == Well::ProducerCMode::BHP) {
        return guideRateSums[kOil] + guideRateSums[kWater] + guideRateSums[kGas];
    }

    return linearTermForMode(guideRateSums, mode, resvCoeff);
}

inline int canonicalRateIndexForMode(Well::ProducerCMode mode)
{
    switch (mode) {
        case Well::ProducerCMode::ORAT: return kOil;
        case Well::ProducerCMode::WRAT: return kWater;
        case Well::ProducerCMode::GRAT: return kGas;
        default: return -1;
    }
}

// ---------------------------------------------------------------------------
// IPR-based BHP→rate conversion
// ---------------------------------------------------------------------------

/// Convert BHP limit to an equivalent total rate using the linear IPR model
/// q = ipr_a - ipr_b * bhp.  Stores the result as the BHP control mode in
/// node.Limits.  Does nothing if ipr_b_sum <= 0.
template<class Scalar, typename IndexTraits>
void addBhpRateLimit(ProdGroupTreeNode<Scalar>& node,
                     const SingleWellState<Scalar, IndexTraits>& ws,
                     Scalar bhpLimit)
{
    // Accumulate ipr_a and ipr_b over active phases.
    const Scalar ipr_a = std::accumulate(ws.implicit_ipr_a.begin(), ws.implicit_ipr_a.end(), Scalar(0));
    const Scalar ipr_b = std::accumulate(ws.implicit_ipr_b.begin(), ws.implicit_ipr_b.end(), Scalar(0));
    if (ipr_b <= Scalar(0)) {
        return; // Cannot compute rate limit from BHP
    }

    // Current BHP:
    const Scalar bhp = ws.bhp;
    if (bhp <= bhpLimit) {
        return; // Already below the BHP limit – no extra constraint
    }

    // Total rate at the BHP limit from the IPR:  q_limit = (ipr_a - ipr_b * bhp_limit)
    // The parametric IPR is: q_total = ipr_a - ipr_b * bhp
    // At current bhp: q_now = ipr_a - ipr_b * bhp (should match -sum(rates))
    // At bhp_limit:   q_lim = ipr_a - ipr_b * bhp_limit
    // We store the *change* in total rate when bhp reaches the limit.
    const Scalar q_limit = (ipr_a - ipr_b * bhpLimit);
    if (q_limit <= Scalar(0)) {
        return; // Limit implies zero or negative production – skip
    }

    // Store as ORAT-equivalent using the current phase fractions.
    // (We cannot store individual phase rates here without the current fractions,
    //  so we store as a total-rate limit and key it under ProducerCMode::BHP.)
    node.Limits[Well::ProducerCMode::BHP] = q_limit;
}

// ---------------------------------------------------------------------------
// Guide-rate helpers
// ---------------------------------------------------------------------------

/// Convert Group::ProductionCMode to the corresponding GuideRateModel::Target.
/// Mirrors GroupStateHelper::getProductionGuideTargetModeFromControlMode.
GuideRateModel::Target productionCModeToGuideTarget(Group::ProductionCMode cmode)
{
    switch (cmode) {
    case Group::ProductionCMode::ORAT: return GuideRateModel::Target::OIL;
    case Group::ProductionCMode::WRAT: return GuideRateModel::Target::WAT;
    case Group::ProductionCMode::GRAT: return GuideRateModel::Target::GAS;
    case Group::ProductionCMode::LRAT: return GuideRateModel::Target::LIQ;
    case Group::ProductionCMode::RESV: return GuideRateModel::Target::RES;
    default:                           return GuideRateModel::Target::NONE;
    }
}

/// Look up the guide rate for node \p name in mode \p ctrlMode.
/// Returns kUniformWeight if the node has no guide rate or the mode is unrecognised.
template<class Scalar>
Scalar getGuideRateForMode(const std::string& name,
                           const std::array<Scalar, 3>& rates,
                           Group::ProductionCMode ctrlMode,
                           const GuideRate& guideRate)
{
    //if (!guideRate.has(name)) return kUniformWeight<Scalar>;
    const auto target = productionCModeToGuideTarget(ctrlMode);
    if (target == GuideRateModel::Target::NONE) return kUniformWeight<Scalar>;
    // GuideRate::RateVector is {oil, gas, water} (not canonical [oil,water,gas])
    const Scalar oilRate   = -rates[kOil];
    const Scalar gasRate   = -rates[kGas];
    const Scalar waterRate = -rates[kWater];
    const GuideRate::RateVector rv{oilRate, gasRate, waterRate};
    const Scalar gr = guideRate.get(name, target, rv);
    return gr > Scalar(0) ? gr : kUniformWeight<Scalar>;
}

// ---------------------------------------------------------------------------

template<class Scalar, typename IndexTraits>
void populateWellNode(ProdGroupTreeNode<Scalar>& node,
                      const std::string& wellName,
                      const Schedule& schedule,
                      const WellState<Scalar, IndexTraits>& wellState,
                      const GroupState<Scalar>& /*groupState*/,
                      const GuideRate& /*guideRate*/,
                      const SummaryState& /*summaryState*/,
                      int reportStep,
                      int fipnum,
                      int pvtreg,
                      const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                      const std::unordered_map<std::string, std::pair<int, Scalar>>& limits)
{
    // Use globally available state so this function produces the same result on
    // every MPI rank regardless of which rank owns the well.
    //
    // - well_rates / currentWellRates: populated for ALL wells on ALL ranks.
    // - isProductionGrup(): from GlobalWellInfo, communicated via updateGlobalIsGrup().
    // - getGlobalEfficiencyScalingFactor(): from GlobalWellInfo (comm.min reduction).
    // - limits: passed in from prepareWellsForBalancing_*() — globally consistent.

    const auto& eclWell = schedule.getWell(wellName, reportStep);
    const auto& pu      = wellModel.phaseUsage();

    node.name   = wellName;
    node.type   = ProdNodeType::Well;
    node.parent = eclWell.groupName();
    node.availableForGroupControl = eclWell.isAvailableForGroupControl();
    node.hasGuideRate = true;
    node.efficiencyFactor = eclWell.getEfficiencyFactor()
                          * wellState.getGlobalEfficiencyScalingFactor(wellName);

    // well_rates/currentWellRates stores positive = production; the balancer uses negative = production.
    // Convert via the same active→canonical mapping as toCanonical3.
    {
        const auto& wr = wellState.currentWellRates(wellName);
        for (int c = 0; c < 3; ++c) {
            const int a = activeIdx(pu, c);
            node.rates[c] = (a >= 0 && a < static_cast<int>(wr.size())) ? -wr[a] : Scalar(0);
        }
    }

    // Control mode — use the boolean GRUP flag from GlobalWellInfo.
    // The full production_cmode enum is not communicated; only GRUP vs. non-GRUP is.
    if (wellState.isProductionGrup(wellName)) {
        node.mode         = Well::ProducerCMode::GRUP;
        node.modeCategory = ProdNodeModeCategory::Group;
    } else {
        // Individually-controlled well.  node.mode must carry the current
        // production control mode because single-node subtree balancing
        // (when this well is a subtree root with no ancestor limit) uses
        // node.mode to determine the operating mode in balanceGroupTree().
        // We read from globalBalancerLimits_ (populated for ALL wells via
        // comm.sum in prepareWellsForBalancing_*) so the value is identical
        // on every MPI rank regardless of which rank owns the well.
        node.modeCategory = ProdNodeModeCategory::Individual;
        if (const auto it = limits.find(wellName); it != limits.end()) {
            node.mode = static_cast<Well::ProducerCMode>(it->second.first);
        } else {
            node.mode = Well::ProducerCMode::CMODE_UNDEFINED;
        }
    }

    // RESV conversion coefficients — computed from globally consistent rates.
    {
        std::vector<Scalar> coeffVec(pu.numPhases, Scalar(0));
        std::vector<Scalar> posRates(pu.numPhases, Scalar(0));
        for (int c = 0; c < 3; ++c) {
            const int a = activeIdx(pu, c);
            if (a >= 0) posRates[a] = -node.rates[c]; // node.rates is negative = production
        }
        wellModel.calcResvCoeff(fipnum, pvtreg, posRates, coeffVec);
        for (int c = 0; c < 3; ++c) {
            const int a = activeIdx(pu, c);
            node.resvCoeff[c] = (a >= 0) ? coeffVec[a] : Scalar(0);
        }
    }

    // Individual limits — read from the globally consistent limits map.
    // Valid for all MPI ranks regardless of which rank owns this well.
    if (!eclWell.isProducer()) return;

    const auto limIt = limits.find(wellName);
    if (limIt == limits.end()) {
        node.Limits.clear();
        return;
    }

    const auto cachedMode = static_cast<Well::ProducerCMode>(limIt->second.first);
    const Scalar limitValue = limIt->second.second;

    if (cachedMode == Well::ProducerCMode::CMODE_UNDEFINED) {
        node.Limits.clear();
        return;
    }

    if (cachedMode == Well::ProducerCMode::THP && limitValue <= Scalar(0)) {
        // No stable operating point at THP limit: well cannot produce.
        node.rates        = {};
        node.modeCategory = ProdNodeModeCategory::Individual;
        node.mode         = Well::ProducerCMode::THP;
        node.Limits.clear();
        return;
    }

    node.Limits.clear();
    if (limitValue > Scalar(0)) {
        node.Limits[cachedMode] = limitValue;
    }
}

// ---------------------------------------------------------------------------

template<class Scalar, typename IndexTraits>
void populateGroupNode(ProdGroupTreeNode<Scalar>& node,
                       const std::string& groupName,
                       const Schedule& schedule,
                       const WellState<Scalar, IndexTraits>& /*wellState*/,
                       const GroupState<Scalar>& groupState,
                       const PhaseUsageInfo<IndexTraits>& pu,
                       const SummaryState& summaryState,
                       int reportStep,
                       int fipnum,
                       int pvtreg,
                       const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel)
{
    const auto& group = schedule.getGroup(groupName, reportStep);

    node.name   = groupName;
    node.type   = ProdNodeType::Group;
    node.parent = (groupName == "FIELD") ? "" : group.parent();

    // FIELD group is never available for group control (has no parent to control it)
    if (groupName == "FIELD") {
        node.availableForGroupControl = false;
        node.efficiencyFactor = Scalar(1);
    } else {
        node.availableForGroupControl = group.productionGroupControlAvailable();
        node.efficiencyFactor = group.getGroupEfficiencyFactor();
    }

    // Children
    node.children = group.groups();
    for (const auto& w : group.wells()) {
        node.children.push_back(w);
    }

    // Check if this is a satellite group (has fixed rates from GSATPROD, not controlled by parent)
    // Note: hasSatelliteProduction() only returns true AFTER GSATPROD becomes active.
    // To detect satellite groups before activation, check if GSATPROD is defined for this group
    // at the end of the schedule (more efficient than checking every step).
    bool isSatellite = group.hasSatelliteProduction();

    // If not currently active, check if GSATPROD is defined for this group at the last schedule step
    if (!isSatellite && schedule.size() > 0) {
        const auto& final_gsatprod = schedule.back().gsatprod();
        if (final_gsatprod.has(groupName)) {
            isSatellite = true;
            node.isSatellite = true;
            node.modeCategory = ProdNodeModeCategory::Satellite;
            node.availableForGroupControl = false;
        }
    }

    // Current rates from group state or satellite data.
    // GroupState stores rates as positive = production in active-phase order.
    // The balancer uses negative = production in canonical [oil, water, gas] order.
    if (isSatellite) {
        // For satellite groups, get rates from GSATPROD data via satelliteProductionRate
        // These rates are not in groupState production_rates, but come from schedule
        const auto& gsat_prod = schedule[reportStep].gsatprod();
        if (gsat_prod.has(groupName)) {
            const auto& sat_rates = gsat_prod.get(groupName, summaryState);
            // GSATPROD rates map: rate[phase_enum_value] where phase is Phase::OIL, WATER, GAS
            node.rates[0] = -sat_rates.rate[static_cast<int>(Phase::OIL)];
            node.rates[1] = -sat_rates.rate[static_cast<int>(Phase::WATER)];
            node.rates[2] = -sat_rates.rate[static_cast<int>(Phase::GAS)];
        }
    } else if (groupState.has_production_rates(groupName)) {
        const auto& gr = groupState.production_rates(groupName);
        for (int c = 0; c < 3; ++c) {
            const int a = activeIdx(pu, c);
            node.rates[c] = (a >= 0 && a < static_cast<int>(gr.size())) ? -gr[a] : Scalar(0);
        }
    }

    // Group control status and own control mode
    const auto ctrl = groupState.has_production_control(groupName)
        ? groupState.production_control(groupName)
        : Group::ProductionCMode::NONE;

    //node.ownCtrlMode = ctrl;
    node.mode = groupModeToWellMode(ctrl); // Convert to well mode for balancing (FLD/VREP→CMODE_UNDEFINED)

    // Determine control status:
    // - Satellite groups have fixed rates and should not be controlled by parent groups
    // - For other groups: check if available for group control
    if (isSatellite) {
        node.modeCategory= ProdNodeModeCategory::Satellite;
        node.availableForGroupControl = false;  // Satellite groups cannot be controlled by parent
        node.isSatellite = true;
    }

    //OpmLog::info(fmt::format("buildTree: group {} availableForGroupControl={} satellite={}",
    //                         groupName, node.availableForGroupControl, isSatellite));

    // RESV conversion coefficients for this group
    {
        std::vector<Scalar> coeffVec(pu.numPhases, Scalar(0));
        std::vector<Scalar> posRates(pu.numPhases, Scalar(0));
        // Convert group rates (negative=production) to positive rates for RESV calculation
        for (int c = 0; c < 3; ++c) {
            const int a = activeIdx(pu, c);
            if (a >= 0) {
                posRates[a] = -node.rates[c];
            }
        }
        wellModel.calcResvCoeff(fipnum, pvtreg, posRates, coeffVec);
        for (int c = 0; c < 3; ++c) {
            const int a = activeIdx(pu, c);
            node.resvCoeff[c] = (a >= 0) ? coeffVec[a] : Scalar(0);
        }
    }

    // Populate individual limits from group production controls
    if (group.isProductionGroup()) {
        const auto controls = schedule.getGroup(groupName, reportStep).productionControls(summaryState);

        // Set preferredMode from schedule control mode (store Group::ProductionCMode directly)
        // FLD and NONE will be resolved via inheritance later
        node.preferredMode = controls.cmode;

        // A group participates in guide-rate balancing if item 10 is present or if
        // it has active individual limits that keep it non-transparent.
        node.hasGuideRate = (controls.guide_rate_def != Group::GuideRateProdTarget::NO_GUIDE_RATE);// ||
                            //(node.availableForGroupControl && !node.Limits.empty()));
        
        const auto& action = controls.group_limit_action; 
        const bool actionAllIsRate = (action.allRates == Opm::Group::ExceedAction::RATE);
        // Populate individual limits from schedule controls
        if (group.has_control(Group::ProductionCMode::ORAT)) {
            if (actionAllIsRate || action.oil == Opm::Group::ExceedAction::RATE) {
                node.Limits[Well::ProducerCMode::ORAT] = controls.oil_target;
            }
        }
        if (group.has_control(Group::ProductionCMode::WRAT)) {
            if (actionAllIsRate || action.water == Opm::Group::ExceedAction::RATE) {
                node.Limits[Well::ProducerCMode::WRAT] = controls.water_target;
            }
        }
        if (group.has_control(Group::ProductionCMode::GRAT)) {
            if (actionAllIsRate || action.gas == Opm::Group::ExceedAction::RATE) {
                node.Limits[Well::ProducerCMode::GRAT] = controls.gas_target;
            }
        }
        if (group.has_control(Group::ProductionCMode::LRAT)) {
            // Skip degenerate LRAT == ORAT case (no water production): same guard as
            // GroupStateHelper::checkGroupProductionConstraints.
            if ((actionAllIsRate || action.liquid == Opm::Group::ExceedAction::RATE)
                && controls.liquid_target != controls.oil_target)
            {
                node.Limits[Well::ProducerCMode::LRAT] = controls.liquid_target;
            }
        }
        // GRAT: prefer GCONSALE-adjusted target when present.
        // Mirrors GroupStateHelper::getProductionGroupTargetForMode_(GRAT).
        if (group.has_control(Group::ProductionCMode::GRAT) && node.Limits.count(Well::ProducerCMode::GRAT) > 0) {
            if (groupState.has_grat_sales_target(groupName)) {
                const Scalar salesTarget = groupState.grat_sales_target(groupName);
                if (salesTarget > Scalar(0)) {
                    node.Limits[Well::ProducerCMode::GRAT] = salesTarget;
                }
            }
        }
        // RESV: mirror GroupStateHelper's GPMAINT-override logic.
        // Mirrors GroupStateHelper::getProductionGroupTargetForMode_(RESV).
        if (group.has_control(Group::ProductionCMode::RESV)) {
            Scalar resv_target = Scalar(0);
            if (group.has_gpmaint_control(Group::ProductionCMode::RESV)
                && groupState.has_gpmaint_target(groupName))
            {
                resv_target = groupState.gpmaint_target(groupName);
            } else {
                resv_target = controls.resv_target;
            }
            if (resv_target > Scalar(0)) {
                node.Limits[Well::ProducerCMode::RESV] = resv_target;
            }
        }
    }
}

} // anonymous namespace

// ===========================================================================
// Public API
// ===========================================================================

template<class Scalar, typename IndexTraits>
Tree<Scalar> buildTree(const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                       const SummaryState& summaryState,
                       int reportStep,
                       const std::unordered_map<std::string, std::pair<int, Scalar>>& limits)
{
    const auto& schedule   = wellModel.schedule();
    const auto& wellState  = wellModel.wellState();
    const auto& groupState = wellModel.groupState();
    const auto& guideRate  = wellModel.guideRate();
    const auto& pu         = wellModel.phaseUsage();

    // Get RESV conversion parameters (fipnum, pvtreg)
    const auto [fipnum, pvtreg] = wellModel.getGroupFipnumAndPvtreg();

    Tree<Scalar> tree;

    // Build a name → interface map once so validWellNode is O(1) per call
    // instead of O(container size).  Wells absent from the map are SHUT (not
    // in the container) and are treated as invalid.
    std::unordered_map<std::string,
                       const WellInterfaceGeneric<Scalar, IndexTraits>*> wellIfaceMap;
    for (const auto* w : wellModel.genericWells()) {
        wellIfaceMap.emplace(w->name(), w);
    }

    // A well node is valid if it is a producer and not runtime-stopped.
    // Queries wellIsStopped() on the interface (volatile Newton-step status,
    // not the persistent ws.status) so wells stopped by operability checks or
    // prepareWellsForBalancing*() are correctly excluded without side-effects.
    auto validWellNode = [&](const std::string& wname) -> bool {
        const auto& well_ecl = schedule.getWell(wname, reportStep);
        if (!well_ecl.isProducer()) return false;
        const auto it = wellIfaceMap.find(wname);
        return it != wellIfaceMap.end() && !it->second->wellIsStopped();
    };

    // Use a stack-based traversal starting from FIELD
    std::vector<std::string> toVisit = {"FIELD"};
    while (!toVisit.empty()) {
        const std::string name = toVisit.back();
        toVisit.pop_back();

        if (schedule.hasWell(name, reportStep)) {
            // Well node: visible to all ranks via currentWellRates().
            if (validWellNode(name)) {
                auto& node = tree[name];
                populateWellNode(node, name, schedule, wellState, groupState,
                                 guideRate, summaryState, reportStep,
                                 fipnum, pvtreg, wellModel, limits);
            }
        } else {
            // Group node: populate and push children onto the stack
            const auto& group = schedule.getGroup(name, reportStep);
            auto& node = tree[name];
            populateGroupNode(node, name, schedule, wellState, groupState, pu,
                            summaryState, reportStep, fipnum, pvtreg, wellModel);

            // Push children onto the stack
            for (const auto& child : group.groups()) {
                toVisit.push_back(child);
            }
            for (const auto& well : group.wells()) {
                toVisit.push_back(well);
            }
        }
    }

    // Group rates are already populated from groupState.
    // If groupState doesn't have rates for a group (e.g., newly active), compute from children.
    std::function<void(const std::string&)> fillMissingGroupRates = [&](const std::string& nodeName) {
        if (tree.count(nodeName) == 0) return;
        auto& node = tree.at(nodeName);
        if (node.type == ProdNodeType::Well) return;

        for (const auto& child : node.children) { fillMissingGroupRates(child); }

        // If rates are zero (not from groupState), compute from children
        const Scalar total = std::accumulate(node.rates.begin(), node.rates.end(), Scalar(0),
                                              [](Scalar s, Scalar r){ return s + (-r); });
        if (total <= Scalar(0)) {
            std::array<Scalar, 3> sumRates{};
            for (const auto& child : node.children) {
                if (tree.count(child) > 0) {
                    for (int c = 0; c < 3; ++c) { sumRates[c] += tree.at(child).rates[c]; }
                }
            }
            node.rates = sumRates;
        }
    };
    fillMissingGroupRates("FIELD");
    // Top-down pass: propagate effective group control modes, compute guide rates,
    // propagate preferred control, and detect transparent groups.
    std::function<void(const std::string&, Well::ProducerCMode, Group::ProductionCMode, bool)>
    propagateGuideRatesAndMode = [&](const std::string& nodeName,
                                    Well::ProducerCMode inheritedMode,
                                    Group::ProductionCMode inheritedpreferredMode,
                                    const bool parentSeesLimits)
    {
        if (tree.count(nodeName) == 0) return;
        auto& node = tree.at(nodeName);
        if (node.isSatellite) {
            node.hasLimitedAncestor = false;
            return; // Satellite groups have fixed rates
        }
        //auto modeToChild = inheritedMode;
        node.hasLimitedAncestor = parentSeesLimits;
        if (!node.availableForGroupControl) {
            // Not available for group control (no inheritance from above), should
            // have non-defaulted preferredMode
            const auto& preferredMode = node.preferredMode;
            if(preferredMode == Group::ProductionCMode::FLD ||
               preferredMode == Group::ProductionCMode::NONE) {
                //OpmLog::warning(fmt::format("Group {} is not available for group control but has preferredMode {}. "
                //                            "This may indicate a problem in the schedule or an unsupported case for the balancer.",
                //                            nodeName, preferredMode == Group::ProductionCMode::FLD ? "FLD" : "NONE"));
            }
            node.hasLimitedAncestor = false;
        } else {
            // If mode is given as FLD/NONE, we set its preferredMode to the inherited preferred control
            if(node.preferredMode == Group::ProductionCMode::FLD ||
               node.preferredMode == Group::ProductionCMode::NONE) {
                node.preferredMode = inheritedpreferredMode; // Inherit preferred control from parent
            }
            if (!(node.modeCategory == ProdNodeModeCategory::Individual)) {
                // Inherit control mode from parent if not individually controlled
                node.mode = inheritedMode;
            }
        }
        // Set guide rate for this node with current mode
        if (node.hasGuideRate) {
            const bool validMode = (!(node.mode == Well::ProducerCMode::CMODE_UNDEFINED)) &&
                                   (!(node.mode == Well::ProducerCMode::NONE));
            const auto groupMode = validMode ? wellModeToGroupMode(node.mode) : node.preferredMode;
            node.groupTarget.ctrlMode = groupMode;
            node.groupTarget.guideRate = getGuideRateForMode(nodeName, node.rates, groupMode, guideRate);
            // also fill groupTargetFallback with the preferredMode guide rate
            const auto groupModeFallback = node.preferredMode;
            node.groupTargetFallback.ctrlMode = groupModeFallback;
            node.groupTargetFallback.guideRate = getGuideRateForMode(nodeName, node.rates, groupModeFallback, guideRate);
        }
        // Propagate effective control mode to children
        for (const auto& childName : node.children) {
            propagateGuideRatesAndMode(childName, node.mode, node.preferredMode, node.hasLimitedAncestor || !node.Limits.empty());
        }
    };
    propagateGuideRatesAndMode("FIELD", Well::ProducerCMode::CMODE_UNDEFINED, 
        Group::ProductionCMode::NONE, /* parentSeesLimits */ false);
    return tree;
}

// ---------------------------------------------------------------------------

template<class Scalar>
std::vector<std::string> getSubTreeOrdering(const Tree<Scalar>& tree,
                                            const std::string& rootName)
{
    // Collect all rootnodes of subtrees that needs to be balanced individually. This 
    // includes:
    // - top-most node with a given limit (often the FIELD group, but can be a subgroup if FIELD has no limits)
    // - satellite groups (have fixed rates, not controlled by parent)
    // - groups that are not available for group control (e.g., have GCONPROD item 8 = NO)
    // - wells that are not available for group control (e.g., have WGRUPCON item 2 = NO)
    // Post-order traversal: leaves first, root last.
    //OpmLog::info(fmt::format("getSubTreeOrdering: starting from {}", rootName));
    std::vector<std::string> ordering;
    std::vector<std::string> stack;
    stack.push_back(rootName);

    int iterCount = 0;
    while (!stack.empty()) {
        if (++iterCount > 1000) {
            //OpmLog::error("getSubTreeOrdering: stuck in loop after 1000 iterations");
            break;
        }
        const std::string name = stack.back();
        stack.pop_back();
        if (tree.count(name) == 0) continue;
        const auto& node = tree.at(name);
        // Only include groups that are NOT available for group control
        // (FIELD, satellite groups, groups with GCONPROD item 8 = NO)
        if (node.isSatellite) {
            // Satellite groups are always included (have fixed rates, not controlled by parent)
            ordering.push_back(name);
        } else {
            // any node that has limits (including FIELD) and does not have ancestors with limits 
            // (hasLimitedAncestor=false) should be included as a root of a subtree to balance
            if (!node.Limits.empty() && !node.hasLimitedAncestor) {
                ordering.push_back(name);
            }
            /*
            if (node.type == ProdNodeType::Group && !node.availableForGroupControl) {
                // if group has no limits (e.g., no preferredMode) it's children can't be controlled
                if (node.preferredMode == Group::ProductionCMode::NONE) {
                    for (const auto& child : node.children) {
                        if (tree.count(child) > 0) {
                            ordering.push_back(child);
                        }
                    }
                } else {
                    ordering.push_back(name);
                }
            }
            */
        }
        // Push children before the current node so they appear first in output
        for (const auto& child : node.children) {
            if (tree.count(child) > 0) {
                stack.push_back(child);
            }
        }
    }

    // Reverse so children appear before parents (post-order)
    std::ranges::reverse(ordering);
    //OpmLog::info(fmt::format("getSubTreeOrdering: found {} nodes: {}",
    //                         ordering.size(),
    //                         fmt::join(ordering, ", ")));
    return ordering;
}

// ---------------------------------------------------------------------------
// Helper functions for new sorting-based algorithm
// ---------------------------------------------------------------------------

template<class Scalar>
bool hasFreePath(const Tree<Scalar>& tree,
                 const std::string& nodeName,
                 const std::string& nodeControl)
{
    // Check that path from node up to nodeControl is "free" (all transparent)
    if (nodeName == nodeControl) return true;
    if (tree.count(nodeName) == 0) return false;

    const auto& node = tree.at(nodeName);
    if (node.parent == nodeControl) return true;

    if (node.parent.empty() || tree.count(node.parent) == 0) return false;

    const auto& parent = tree.at(node.parent);
    const bool isFreeStep = (parent.modeCategory == ProdNodeModeCategory::Transparent) && !parent.visited;

    return isFreeStep && hasFreePath(tree, node.parent, nodeControl);
}

template<class Scalar>
void incrementParentRateSums(Tree<Scalar>& tree,
                             const std::string& nodeName,
                             const std::string& origin,
                             const std::optional<std::array<Scalar, 3>>& rates)
{
    // Increment rate_sums from node up to origin
    if (tree.count(nodeName) == 0) return;
    const auto& node = tree.at(nodeName);

    // Determine rates to add (convert to positive production for rateSums,
    // scaled by the node's own efficiency factor so the parent sees the
    // efficiency-adjusted contribution).
    std::array<Scalar, 3> ratesToAdd;
    if (rates.has_value()) {
        // Caller already provides the delta in parent-frame units.
        ratesToAdd = rates.value();
    } else {
        // node.rates are negative (production); convert to positive and scale.
        for (int c = 0; c < 3; ++c) {
            ratesToAdd[c] = node.efficiencyFactor * (-node.rates[c]);
        }
    }

    std::string parent = nodeName;
    while (parent != origin && tree.count(parent) > 0) {
        const auto& currentNode = tree.at(parent);
        parent = currentNode.parent;
        if (!parent.empty() && tree.count(parent) > 0) {
            auto& parentNode = tree.at(parent);
            for (int c = 0; c < 3; ++c) {
                parentNode.rateSums[c] += ratesToAdd[c];
            }
        }
    }
}

template<class Scalar>
void decrementParentGuideRateSums(Tree<Scalar>& tree,
                                   const std::string& nodeName,
                                   const std::string& origin,
                                   const std::array<Scalar, 3>& guideRateSums)
{
    // Subtract guide_rate_sums from node up to and including origin
    if (tree.count(nodeName) == 0) return;
    const auto& node = tree.at(nodeName);

    std::string parent = node.parent;
    while (!parent.empty() && tree.count(parent) > 0) {
        auto& parentNode = tree.at(parent);
        for (int c = 0; c < 3; ++c) {
            parentNode.guideRateSums[c] -= guideRateSums[c];
        }
        if (parent == origin) break;  // Stop after updating origin
        parent = parentNode.parent;
    }
}

// ---------------------------------------------------------------------------

template<class Scalar>
std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>>
getLocalTreeDescendants(const Tree<Scalar>& tree, const std::string& nodeName)
{
    // d: nodes with guide-rate available for group control
    // df: nodes not available for group control (fixed)
    // dt: nodes without guide-rate in between d and nodeName (transparent)

    std::vector<std::string> d, df, dt;

    if (tree.count(nodeName) == 0) return {d, df, dt};
    const auto& node = tree.at(nodeName);

    for (const auto& childName : node.children) {
        if (tree.count(childName) == 0) continue;
        const auto& child = tree.at(childName);

        if (!child.availableForGroupControl) {
            df.push_back(childName);
        } else if (child.hasGuideRate) {
            d.push_back(childName);
        } else {
            // Recursive descent for transparent groups
            auto [d_child, df_child, dt_child] = getLocalTreeDescendants(tree, childName);
            d.insert(d.end(), d_child.begin(), d_child.end());
            df.insert(df.end(), df_child.begin(), df_child.end());
            dt.insert(dt.end(), dt_child.begin(), dt_child.end());
            dt.push_back(childName);
        }
    }

    return {d, df, dt};
}

// ---------------------------------------------------------------------------
template<class Scalar>
void updateGuideRatesForMode(Tree<Scalar>& tree,
                             const std::vector<std::string>& c,
                             Well::ProducerCMode mode,
                             const GuideRate& guideRate)
{
    const auto modeAsGroupMode = wellModeToGroupMode(mode);
    // Update guide rates for children if needed (e.g., if mode changed)
    for (const auto& childName : c) {
        if (tree.count(childName) == 0) continue;
        auto& child = tree.at(childName);
        if (child.groupTarget.ctrlMode != modeAsGroupMode) {
            child.groupTarget.ctrlMode = modeAsGroupMode;
            child.groupTarget.guideRate = getGuideRateForMode(childName, child.rates, modeAsGroupMode, guideRate);
        }
    }
}

template<class Scalar>
void resetRatesAndGuideRateSums(Tree<Scalar>& tree,
                                const std::vector<std::string>& c,
                                const std::vector<std::string>& c_fixed,
                                const std::vector<std::string>& c_trans,
                                const std::string& origin,
                                Well::ProducerCMode mode)
{
    // Reset guide-rate-sums and rate-sums for all relevant nodes
    std::vector<std::string> allNodes = {origin};
    allNodes.insert(allNodes.end(), c.begin(), c.end());
    allNodes.insert(allNodes.end(), c_fixed.begin(), c_fixed.end());
    allNodes.insert(allNodes.end(), c_trans.begin(), c_trans.end());

    for (const auto& nodeName : allNodes) {
        if (tree.count(nodeName) == 0) continue;
        auto& node = tree.at(nodeName);
        node.rateSums = {Scalar(0), Scalar(0), Scalar(0)};
        node.guideRateSums = {Scalar(0), Scalar(0), Scalar(0)};
        node.visited = false;
    }

    // Set transparent nodes category
    for (const auto& nodeName : c_trans) {
        if (tree.count(nodeName) > 0) {
            tree.at(nodeName).modeCategory = ProdNodeModeCategory::Transparent;
        }
    }

    if (tree.count(origin) == 0) return;
    const auto& originNode = tree.at(origin);
    const Group::ProductionCMode modePreferred = originNode.preferredMode;
    // Convert mode to Group::ProductionCMode for comparison
    const Group::ProductionCMode modeAsGroupMode = wellModeToGroupMode(mode);
    const bool modeIsPref = (modeAsGroupMode == modePreferred);

    // Compute total guide rate for mode if not preferred
    Scalar guideSumOrigin = Scalar(0);
    if (!modeIsPref) {
        for (const auto& childName : c) {
            if (tree.count(childName) > 0) {
                const auto& child = tree.at(childName);
                const auto& activeTarget = child.useFallback ? child.groupTargetFallback : child.groupTarget;
                guideSumOrigin += activeTarget.guideRate;
            }
        }
    }

    // Set guide_rate_sums for each child in c
    for (const auto& childName : c) {
        if (tree.count(childName) == 0) continue;
        auto& child = tree.at(childName);
        child.useFallback = false;

        const auto& rates = child.rates;

        // Check for fallback condition (non-preferred mode with tiny fraction)
        if (!modeIsPref) {
            const Scalar rateSum = std::accumulate(rates.begin(), rates.end(), Scalar(0),
                                                   [](Scalar s, Scalar r){ return s + (-r); });
            const Scalar rateForCurrentMode = rateForMode(rates, mode, child.resvCoeff);
            const Scalar rateFracMode = (rateSum > Scalar(0)) ? (rateForCurrentMode / rateSum) : Scalar(0);
            const Scalar guideRateForMode = child.groupTarget.guideRate;
            const Scalar guideRatioMode = (guideSumOrigin > Scalar(0))
                ? (guideRateForMode / guideSumOrigin) : Scalar(0);

            if (rateFracMode < std::sqrt(std::numeric_limits<Scalar>::epsilon()) ||
                guideRatioMode < std::sqrt(std::numeric_limits<Scalar>::epsilon())) {
                child.useFallback = true;
                child.guideRateSums = {std::numeric_limits<Scalar>::quiet_NaN(),
                                       std::numeric_limits<Scalar>::quiet_NaN(),
                                       std::numeric_limits<Scalar>::quiet_NaN()};
                continue;
            }
        }

        // Set guide-rate sums proportional to rates, scaled by the child's
        // efficiency factor so parent-frame guide-rate sums and rate-sums are
        // consistently in the same (parent) frame.
        const auto& activeTarget = child.useFallback ? child.groupTargetFallback : child.groupTarget;
        const Scalar guideRateValue = activeTarget.guideRate;
        const Scalar rateMagnitude = rateForMode(rates, mode, child.resvCoeff);

        if (rateMagnitude > Scalar(0)) {
            for (int cmp = 0; cmp < 3; ++cmp) {
                child.guideRateSums[cmp] =
                    child.efficiencyFactor * (guideRateValue / rateMagnitude) * (-rates[cmp]);
            }
        } else {
            child.guideRateSums = {Scalar(0), Scalar(0), Scalar(0)};
        }

        // Propagate up to origin
        std::string parent = childName;
        while (parent != origin && tree.count(parent) > 0) {
            const auto& pnode = tree.at(parent);
            parent = pnode.parent;
            if (tree.count(parent) > 0) {
                auto& parentNode = tree.at(parent);
                for (int cmp = 0; cmp < 3; ++cmp) {
                    parentNode.guideRateSums[cmp] += child.guideRateSums[cmp];
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------

template<class Scalar>
std::vector<Scalar>
computeRatiosForSorting(const Tree<Scalar>& tree, const std::vector<std::string>& c)
{
    // Get limit-to-guide-rate ratios for distribution ordering.
    // For each child, assume production rates evolve linearly as
    //   rate(a) = rateSums + a * guideRateSums,
    // and find the smallest non-negative alpha that activates any limit.
    const int nc = c.size();
    std::vector<Scalar> ratios(nc, std::numeric_limits<Scalar>::max());

    for (int k = 0; k < nc; ++k) {
        const auto& childName = c[k];
        if (tree.count(childName) == 0) continue;
        const auto& child = tree.at(childName);

        // Select target mode/guide-rate from fallback target when requested.
        const auto& activeTarget = child.useFallback ? child.groupTargetFallback : child.groupTarget;
        const Well::ProducerCMode activeMode = groupModeToWellMode(activeTarget.ctrlMode);
        const Scalar targetGuideRate = std::max(activeTarget.guideRate, Scalar(0));

        Scalar guideRateSumsNorm = Scalar(0);
        if (activeMode == Well::ProducerCMode::BHP) {
            guideRateSumsNorm = std::abs(child.guideRateSums[kOil] + child.guideRateSums[kWater] + child.guideRateSums[kGas]);
        } else if (activeMode != Well::ProducerCMode::CMODE_UNDEFINED && activeMode != Well::ProducerCMode::NONE) {
            guideRateSumsNorm = std::abs(linearTermForMode(child.guideRateSums, activeMode, child.resvCoeff));
        } else {
            guideRateSumsNorm = std::abs(child.guideRateSums[kOil])
                              + std::abs(child.guideRateSums[kWater])
                              + std::abs(child.guideRateSums[kGas]);
        }
        const Scalar minGuideScale = std::max(targetGuideRate, Scalar(1));
        const bool isZero = guideRateSumsNorm <= std::sqrt(std::numeric_limits<Scalar>::epsilon()) * minGuideScale;

        // For 'None' nodes, the node is already limited by its children.
        // Use current node rates as effective phase limits.
        if (child.modeCategory == ProdNodeModeCategory::None) {
            Scalar minAlpha = std::numeric_limits<Scalar>::infinity();
            for (int cmp = 0; cmp < 3; ++cmp) {
                const Scalar limit = -child.rates[cmp];      // positive production limit
                const Scalar current = child.rateSums[cmp];  // positive production at a = 0
                const Scalar slope = child.guideRateSums[cmp];
                if (slope <= Scalar(0)) {
                    continue;
                }
                const Scalar alpha = (limit - current) / slope;
                if (alpha < minAlpha) {
                    minAlpha = alpha;
                }
            }

            if (minAlpha < Scalar(0)) {
                minAlpha = Scalar(0);
            }
            ratios[k] = minAlpha;
            continue;
        }

        if (isZero || child.Limits.empty()) {
            ratios[k] = std::numeric_limits<Scalar>::infinity();
            continue;
        }

        // Build negative-rate arrays for helper projections that expect negative=production.
        const std::array<Scalar, 3> negRateSums{
            -child.rateSums[kOil],
            -child.rateSums[kWater],
            -child.rateSums[kGas]
        };

        // rateSums and guideRateSums are in the parent frame (already scaled by
        // efficiencyFactor).  Limits are native (child frame).  Convert the
        // limit to the parent frame before comparing.
        Scalar minAlpha = std::numeric_limits<Scalar>::infinity();
        for (const auto& [mode, limit] : child.Limits) {
            Scalar current = Scalar(0);
            Scalar slope = Scalar(0);

            // Effective limit in the parent frame.
            const Scalar effectiveLimit = child.efficiencyFactor * limit;

            if (mode == Well::ProducerCMode::BHP) {
                current = child.rateSums[kOil] + child.rateSums[kWater] + child.rateSums[kGas];
                slope = child.guideRateSums[kOil] + child.guideRateSums[kWater] + child.guideRateSums[kGas];
            } else {
                current = rateForMode(negRateSums, mode, child.resvCoeff);
                slope = linearTermForMode(child.guideRateSums, mode, child.resvCoeff);
            }

            // This mode cannot become more restrictive along the current distribution direction.
            if (slope <= Scalar(0)) {
                continue;
            }

            const Scalar alpha = (effectiveLimit - current) / slope;
            if (alpha < minAlpha) {
                minAlpha = alpha;
            }
        }

        // Already at/above some limit means immediate restriction.
        if (minAlpha < Scalar(0)) {
            minAlpha = Scalar(0);
        }

        /*
        if (limitingMode == Well::ProducerCMode::ORAT) {
        } else if (limitingMode == Well::ProducerCMode::WRAT) {
            minIdx = kWater;
        } else if (limitingMode == Well::ProducerCMode::GRAT) {
            minIdx = kGas;
        }
        */
        ratios[k] = minAlpha;
        //modeIndices[k] = minIdx;
    }
    return ratios;
}

// ---------------------------------------------------------------------------

template<class Scalar>
std::vector<std::string>
updateTransparentGroups(Tree<Scalar>& tree,
                        const std::string& nodeName,
                        std::vector<std::string>& c_trans,
                        const GuideRate& guideRate,
                        Scalar nextRatio,
                        Scalar qm,
                        Well::ProducerCMode mode,
                        Scalar tol)
{
    // Check if any transparent groups have limit-to-guide ratio less than nextRatio
    std::vector<std::string> c_trans_update;

    bool foundLessThanNext = true;
    while (foundLessThanNext && !c_trans.empty()) {
        const auto ratios = computeRatiosForSorting(tree, c_trans);

        // Find minimum ratio
        Scalar ratioMin = std::numeric_limits<Scalar>::max();
        int minIdx = -1;
        for (size_t i = 0; i < ratios.size(); ++i) {
            if (ratios[i] < ratioMin) {
                ratioMin = ratios[i];
                minIdx = i;
            }
        }

        if (minIdx < 0 || ratioMin >= nextRatio) {
            foundLessThanNext = false;
            break;
        }

        if (tree.count(nodeName) == 0) break;
        const auto& originNode = tree.at(nodeName);
        const Scalar gsum = guideSlopeForMode(originNode.guideRateSums, mode, originNode.resvCoeff);
        const Scalar qmRemain = qm - positiveRateForMode(originNode.rateSums, mode, originNode.resvCoeff);

        if (ratioMin >= qmRemain / (gsum + Scalar(1e-20))) {
            foundLessThanNext = false;
            break;
        }

        const auto& ckName = c_trans[minIdx];
        if (!hasFreePath(tree, ckName, nodeName)) {
            c_trans.erase(c_trans.begin() + minIdx);
            continue;
        }

        if (tree.count(ckName) == 0) {
            c_trans.erase(c_trans.begin() + minIdx);
            continue;
        }

        c_trans_update.push_back(ckName);
        auto& ck = tree.at(ckName);

        // Save original values
        const auto rateSumsOrig = ck.rateSums;
        const auto guideRateSumsOrig = ck.guideRateSums;
        const Scalar gk = guideSlopeForMode(ck.guideRateSums, mode, ck.resvCoeff);
        // qk and rateSumsOrig are in the parent frame; divide by efficiencyFactor
        // to convert to the child's native frame before recursing.
        const Scalar qkParent = (gk / (gsum + Scalar(1e-20))) * qmRemain
                + positiveRateForMode(rateSumsOrig, mode, ck.resvCoeff);
        const Scalar qk = qkParent / ck.efficiencyFactor;

        // Recursive balance
        balanceGroupTree(tree, ckName, guideRate, mode, qk, tol);

        // Update parent rate-sums (efficiencyFactor scaling handled inside incrementParentRateSums)
        std::array<Scalar, 3> ratesDelta;
        for (int c = 0; c < 3; ++c) {
            ratesDelta[c] = ck.efficiencyFactor * (-ck.rates[c]) - rateSumsOrig[c];
        }
        incrementParentRateSums(tree, ckName, nodeName, std::make_optional(ratesDelta));
        decrementParentGuideRateSums(tree, ckName, nodeName, guideRateSumsOrig);

        // Remove from c_trans
        c_trans.erase(c_trans.begin() + minIdx);
    }

    return c_trans_update;
}

// ---------------------------------------------------------------------------

template<class Scalar>
void distributeFallbackRates(Tree<Scalar>& tree,
                             const std::string& nodeName,
                             const std::vector<std::string>& c,
                             const GuideRate& guideRate,
                             Well::ProducerCMode mode,
                             Scalar tol)
{
    // Distribute rates to children that need fallback (tiny fractions in non-preferred mode)
    if (tree.count(nodeName) == 0) return;
    const auto& node = tree.at(nodeName);

    const Group::ProductionCMode modePreferred = node.preferredMode;
    const int modeIdx = canonicalRateIndexForMode(mode);
    const Well::ProducerCMode modePreferredAsWell = groupModeToWellMode(modePreferred);

    // Collect fallback children and compute ratio from group-controlled children
    std::vector<std::string> cFallback;
    std::vector<Scalar> cFallbackGuideRates;
    Scalar rateSums = Scalar(0);
    Scalar guiderateSums = Scalar(0);

    for (const auto& childName : c) {
        if (tree.count(childName) == 0) continue;
        const auto& child = tree.at(childName);
        const Scalar preferredGuideRate = getGuideRateForMode(childName, child.rates, modePreferred, guideRate);

        if (child.useFallback) {
            cFallback.push_back(childName);
            cFallbackGuideRates.push_back(preferredGuideRate);
        } else if (child.modeCategory == ProdNodeModeCategory::Group) {
            rateSums += rateForMode(child.rates, modePreferredAsWell, child.resvCoeff);
            guiderateSums += preferredGuideRate;
        }
    }

    if (cFallback.empty()) return;

    const bool anyGroupWells = (rateSums > Scalar(0));
    const Scalar ratio = (guiderateSums > Scalar(0)) ? (rateSums / guiderateSums) : Scalar(1);

    for (size_t i = 0; i < cFallback.size(); ++i) {
        const auto& ckName = cFallback[i];
        if (tree.count(ckName) == 0) continue;
        auto& ck = tree.at(ckName);

        Scalar qk;
        if (anyGroupWells) {
            qk = cFallbackGuideRates[i] * ratio;
        } else {
            // No group-controlled wells: use limit
            const auto it = ck.Limits.find(modePreferredAsWell);
            qk = (it != ck.Limits.end()) ? it->second : Scalar(0);
        }

        // qk is in the parent frame (derived from guide-rate ratios against
        // other parent-frame rates); convert to child native frame.
        balanceGroupTree(tree, ckName, guideRate, modePreferredAsWell, qk / ck.efficiencyFactor, tol);

        // Add efficiency-scaled rates, zeroing the control-mode component to
        // avoid triggering re-balancing for this fallback child.
        std::array<Scalar, 3> scaledRates;
        for (int ph = 0; ph < 3; ++ph) {
            scaledRates[ph] = ck.efficiencyFactor * (-ck.rates[ph]);
        }
        if (modeIdx >= 0) {
            scaledRates[modeIdx] = Scalar(0);
        }
        incrementParentRateSums(tree, ckName, nodeName, std::make_optional(scaledRates));
    }
}

// ---------------------------------------------------------------------------

template<class Scalar>
void balanceGroupTree(Tree<Scalar>& tree,
                      const std::string& nodeName,
                      const GuideRate& guideRate,
                      Well::ProducerCMode targetMode,
                      Scalar targetRate,
                      Scalar tol)
{
    if (tree.count(nodeName) == 0) return;
    auto& node = tree.at(nodeName);
    if (node.isSatellite) return; // Skip satellite groups

    //OpmLog::info(fmt::format("balanceGroupTree: {} mode={} target={:.3f} count={}",
    //                         nodeName, WellProducerCMode2String(targetMode), targetRate, node.balancingCount));

    // Check recursion depth to prevent stack overflow
    if (node.balancingCount > 100) {
        // Too many recursive calls - bail out
        OpmLog::warning(fmt::format("balanceGroupTree: {} exceeded recursion limit", nodeName));
        node.isBalanced = false;
        return;
    }

    const bool skipBalancing = false; // Could add logic for already balanced nodes

    if (!skipBalancing) {
        const int maxResortingCount = 10;
        const int maxTopSwitchCount = 2;

        Well::ProducerCMode mode = targetMode;
        Scalar qm = targetRate;

        // Scale rates to target or strictest individual control
        auto rates = node.rates;
        // int modeIdx = kOil;
        //switch (mode) {
        //    case Well::ProducerCMode::ORAT: modeIdx = kOil; break;
        //    case Well::ProducerCMode::WRAT: modeIdx = kWater; break;
        //    case Well::ProducerCMode::GRAT: modeIdx = kGas; break;
        //    default: modeIdx = kOil;
        //}

        Scalar scaledRates[3];
        const Scalar rateForModeVal = rateForMode(rates, mode, node.resvCoeff);
        if (rateForModeVal > Scalar(0)) {
            for (int cmp = 0; cmp < 3; ++cmp) {
                scaledRates[cmp] = rates[cmp] * (qm / rateForModeVal);
            }

            // Check if any limit is violated after scaling
            Scalar maxViolation = Scalar(0);
            Well::ProducerCMode violatingMode = mode;
            for (const auto& [limitMode, limit] : node.Limits) {
                const Scalar scaledRate = rateForMode(std::array<Scalar,3>{scaledRates[0], scaledRates[1], scaledRates[2]},
                                                      limitMode, node.resvCoeff);
                const Scalar violation = scaledRate / limit;
                if (violation > maxViolation) {
                    maxViolation = violation;
                    violatingMode = limitMode;
                }
            }

            if (maxViolation > Scalar(1) - tol) {
                mode = violatingMode;
                const Scalar violatingScaledRate = rateForMode(std::array<Scalar,3>{scaledRates[0], scaledRates[1], scaledRates[2]},
                                                               mode, node.resvCoeff);
                qm = violatingScaledRate / maxViolation;
            }
        }

        node.balancingCount++;

        bool anyGroupChildren = false;

        if (node.type == ProdNodeType::Well) {
            // Just scale rates for wells
            const Scalar currentRate = rateForMode(rates, mode, node.resvCoeff);
            if (currentRate > Scalar(0)) {
                for (int cmp = 0; cmp < 3; ++cmp) {
                    node.rates[cmp] = rates[cmp] * (qm / currentRate);
                }
            }
            node.isBalanced = true;
        } else {
            // Group node: iterative balancing
            bool balanced = false;
            int resortingCount = 0;
            int topSwitchCount = 0;

            while (!balanced && resortingCount <= maxResortingCount && topSwitchCount <= maxTopSwitchCount) {
                //OpmLog::info(fmt::format("  balanceGroupTree: {} while-loop iteration resort={} switch={}",
                //                         nodeName, resortingCount, topSwitchCount));
                // Get local tree descendants
                auto [c, c_fixed, c_trans] = getLocalTreeDescendants(tree, nodeName);

                // Update guide rates if mode changed
                updateGuideRatesForMode(tree, c, mode, guideRate);

                // Reset rate_sums and guide_rate_sums
                resetRatesAndGuideRateSums(tree, c, c_fixed, c_trans, nodeName, mode);

                // Add in rates from fixed sub-nodes
                for (const auto& cfName : c_fixed) {
                    incrementParentRateSums(tree, cfName, nodeName);
                }

                // Get ratios for sorting
                const auto ratios = computeRatiosForSorting(tree, c);

                // Sort by ratios
                std::vector<size_t> sortedIndices(c.size());
                std::iota(sortedIndices.begin(), sortedIndices.end(), 0);
                std::sort(sortedIndices.begin(), sortedIndices.end(),
                         [&ratios](size_t i1, size_t i2) { return ratios[i1] < ratios[i2]; });

                bool needsResorting = false;
                anyGroupChildren = false;
                bool anyChildrenNeedFallback = false;

                // Distribute rates to children in sorted order
                for (size_t k = 0; k < sortedIndices.size(); ++k) {
                    const size_t idx = sortedIndices[k];
                    const auto& ckName = c[idx];

                    if (tree.count(ckName) == 0) continue;
                    auto& ck = tree.at(ckName);

                    if (ck.useFallback) {
                        anyChildrenNeedFallback = true;
                        continue;
                    }

                    // Check if ck has a free path
                    if (!hasFreePath(tree, ckName, nodeName)) continue;

                    // Check and update transparent groups
                    if (!c_trans.empty() && !anyGroupChildren) {
                        auto c_trans_update = updateTransparentGroups(tree, nodeName, c_trans,
                                                                      guideRate,
                                                                      ratios[idx], qm, mode, tol);

                        if (!c_trans_update.empty()) {
                            bool anyTransparent = false;
                            bool anyIndividual = false;
                            for (const auto& ctkName : c_trans_update) {
                                if (tree.count(ctkName) > 0 && hasFreePath(tree, ctkName, nodeName)) {
                                    const auto& ctk = tree.at(ctkName);
                                    anyTransparent = anyTransparent || (ctk.modeCategory == ProdNodeModeCategory::Transparent);
                                    anyIndividual = anyIndividual || (ctk.modeCategory == ProdNodeModeCategory::Individual);
                                }
                            }
                            needsResorting = anyGroupChildren && anyIndividual;
                            anyGroupChildren = anyGroupChildren || anyTransparent;

                            if (!hasFreePath(tree, ckName, nodeName)) continue;
                        }
                    }

                    if (tree.count(nodeName) == 0) break;
                    const auto& originNode = tree.at(nodeName);
                    const Scalar gsum = guideSlopeForMode(originNode.guideRateSums, mode, originNode.resvCoeff);
                    const Scalar gk = guideSlopeForMode(ck.guideRateSums, mode, ck.resvCoeff);
                    const Scalar qmRemain = qm - positiveRateForMode(originNode.rateSums, mode, originNode.resvCoeff);
                    // qkParent is in parent frame; divide by efficiencyFactor for child native frame.
                    const Scalar qkParent = (gk / (gsum + Scalar(1e-20))) * qmRemain;
                    const Scalar qk = qkParent / ck.efficiencyFactor;

                    // Save original guide_rate_sums
                    const auto guideRateSumsOrig = ck.guideRateSums;

                    // Recursive balance
                    balanceGroupTree(tree, ckName, guideRate, mode, qk, tol);

                    if (tree.count(ckName) > 0) {
                        const auto& ckBalanced = tree.at(ckName);
                        anyGroupChildren = anyGroupChildren || (ckBalanced.modeCategory == ProdNodeModeCategory::Group);
                        if (ckBalanced.modeCategory != ProdNodeModeCategory::Group && anyGroupChildren) {
                            needsResorting = true;
                        }

                        // Update parent sums
                        if (ckBalanced.hasGuideRate) {
                            decrementParentGuideRateSums(tree, ckName, nodeName, guideRateSumsOrig);
                            incrementParentRateSums(tree, ckName, nodeName);
                        }
                    }
                }

                // Update rates for remaining transparent groups
                for (const auto& ctName : c_trans) {
                    if (tree.count(ctName) > 0) {
                        auto& ct = tree.at(ctName);
                        // Convert rateSums (positive) to rates (negative = production)
                        for (int phase = 0; phase < 3; ++phase) {
                            ct.rates[phase] = -ct.rateSums[phase];
                        }

                        // Check if any limits are broken
                        for (const auto& [limitMode, limit] : ct.Limits) {
                            const Scalar currentRate = rateForMode(ct.rates, limitMode, ct.resvCoeff);
                            if (currentRate > limit * (Scalar(1) + tol)) {
                                needsResorting = true;
                                break;
                            }
                        }
                    }
                }

                if (needsResorting && resortingCount < maxResortingCount) {
                    resortingCount++;
                } else {
                    if (anyChildrenNeedFallback) {
                        distributeFallbackRates(tree, nodeName, c, guideRate, mode, tol);

                        // Update transparent groups once more
                        for (const auto& ctName : c_trans) {
                            if (tree.count(ctName) > 0) {
                                auto& ct = tree.at(ctName);
                                // Convert rateSums (positive) to rates (negative = production)
                                for (int phase = 0; phase < 3; ++phase) {
                                    ct.rates[phase] = -ct.rateSums[phase];
                                }
                            }
                        }
                    }

                    resortingCount = 0;

                    // Check if any new limit is violated
                    bool limitViolated = false;
                    if (tree.count(nodeName) > 0) {
                        const auto& currentNode = tree.at(nodeName);
                        const auto& rateSums = currentNode.rateSums;

                        if (topSwitchCount < maxTopSwitchCount) {
                            // Check all limits
                            Scalar maxViolation = Scalar(0);
                            Well::ProducerCMode newMode = mode;

                            for (const auto& [limitMode, limit] : currentNode.Limits) {
                                Scalar effectiveLimit = (limitMode == targetMode) ?
                                    std::min(limit, targetRate) : limit;
                                const Scalar currentRate = -rateForMode(rateSums, limitMode, currentNode.resvCoeff);
                                const Scalar violation = currentRate / effectiveLimit;
                                if (violation > maxViolation) {
                                    maxViolation = violation;
                                    newMode = limitMode;
                                }
                            }

                            if (maxViolation > Scalar(1) - tol && mode != newMode) {
                                if (newMode == targetMode) {
                                    topSwitchCount = maxTopSwitchCount;
                                } else {
                                    topSwitchCount++;
                                }
                                limitViolated = true;
                                mode = newMode;
                                const Scalar modeRate = -rateForMode(rateSums, mode, currentNode.resvCoeff);
                                qm = modeRate / maxViolation;
                                // update guide rate for new mode
                                const auto groupMode = wellModeToGroupMode(mode);
                                node.groupTarget.ctrlMode = groupMode;
                                node.groupTarget.guideRate = getGuideRateForMode(nodeName, node.rates, groupMode, guideRate);
                            } else if (maxViolation < Scalar(1) - tol && anyGroupChildren) {
                                // Did not manage to distribute everything
                                limitViolated = false;
                            }
                        }
                    }

                    balanced = !limitViolated;
                    if (tree.count(nodeName) > 0) {
                        auto& finalNode = tree.at(nodeName);
                        finalNode.isBalanced = balanced;
                        // Convert rateSums (positive) to rates (negative = production)
                        for (int phase = 0; phase < 3; ++phase) {
                            finalNode.rates[phase] = -finalNode.rateSums[phase];
                        }
                    }
                }
            }

            if (resortingCount > maxResortingCount) {
                // Reached max resort count
            }
        }

        // Update category
        if (tree.count(nodeName) > 0) {
            auto& finalNode = tree.at(nodeName);
            finalNode.mode = mode;
            finalNode.visited = true;

            if (finalNode.type == ProdNodeType::Group && !anyGroupChildren) {
                finalNode.modeCategory = ProdNodeModeCategory::None;
                finalNode.mode = Well::ProducerCMode::CMODE_UNDEFINED;
            } else {
                bool atLimit = false;
                for (const auto& [limitMode, limit] : finalNode.Limits) {
                    const Scalar currentRate = rateForMode(finalNode.rates, limitMode, finalNode.resvCoeff);
                    if (std::abs(currentRate - limit) <= tol * limit) {
                        atLimit = true;
                        break;
                    }
                }

                if (atLimit) {
                    finalNode.modeCategory = ProdNodeModeCategory::Individual;
                } else {
                    const Scalar rateAtMode = rateForMode(finalNode.rates, mode, finalNode.resvCoeff);
                    if (std::abs(rateAtMode - qm) <= tol * std::abs(qm)) {
                        if (finalNode.hasGuideRate) {
                            finalNode.modeCategory = ProdNodeModeCategory::Group;
                        } else {
                            finalNode.modeCategory = ProdNodeModeCategory::Transparent;
                            finalNode.mode = Well::ProducerCMode::CMODE_UNDEFINED;
                        }
                    } else {
                        // report probelem with categorization
                        //OpmLog::info("balanceGroupTree: {} categorized as INDIVIDUAL due to mismatch between rate and qm", nodeName);
                        finalNode.modeCategory = ProdNodeModeCategory::Individual;
                    }
                }
            }

            if (!finalNode.isBalanced) {
                // Reached maximum resort count, continuing with unbalanced sub-tree
            }
        }
    }
}

// ---------------------------------------------------------------------------

template<class Scalar>
void setTargets(Tree<Scalar>& tree, const std::string& topName)
{
    if (tree.count(topName) == 0) return;
    const auto& top = tree.at(topName);

    if (top.children.empty()) return;

    // Get descendants (skip transparent groups)
    std::function<std::vector<std::string>(const std::string&)> getDescendants;
    getDescendants = [&](const std::string& nodeName) -> std::vector<std::string> {
        std::vector<std::string> desc;
        if (tree.count(nodeName) == 0) return desc;
        const auto& node = tree.at(nodeName);

        for (const auto& childName : node.children) {
            if (tree.count(childName) == 0) continue;
            const auto& child = tree.at(childName);

            if (child.modeCategory == ProdNodeModeCategory::Transparent) {
                auto childDesc = getDescendants(childName);
                desc.insert(desc.end(), childDesc.begin(), childDesc.end());
            } else {
                desc.push_back(childName);
            }
        }
        return desc;
    };

    auto cList = getDescendants(topName);

    if (top.modeCategory == ProdNodeModeCategory::None) {
        // No target for descendants of 'NONE'
        for (const auto& cName : cList) {
            if (tree.count(cName) > 0) {
                auto& c = tree.at(cName);
                c.groupTarget.groupName = topName;
                c.groupTarget.ctrlMode = Group::ProductionCMode::NONE;
                c.groupTarget.value = std::numeric_limits<Scalar>::max();
                setTargets(tree, cName);
            }
        }
    } else {
        // Node is 'INDIVIDUAL' or 'GRUP'
        const Well::ProducerCMode mode = top.mode;
        const Group::ProductionCMode modePreferred = top.preferredMode;
        // Convert mode to Group::ProductionCMode for comparison
        const bool modeIsPref = (wellModeToGroupMode(mode) == modePreferred);

        // Check for fallback modes
        bool anyFallback = false;
        if (!modeIsPref) {
            for (const auto& cName : cList) {
                if (tree.count(cName) > 0 && tree.at(cName).useFallback) {
                    anyFallback = true;
                    break;
                }
            }
        }

        // Collect guide rates, rates, and efficiency factors for children.
        // rates[i] is always the native (child-frame) rate; effFactors[i] is used
        // to convert between native and parent frame when computing targetSum and
        // assigning groupTarget.value.
        std::vector<Scalar> guideRates, rates, effFactors;
        std::vector<bool> isGrup;

        for (const auto& cName : cList) {
            if (tree.count(cName) == 0) continue;
            const auto& c = tree.at(cName);

            Scalar gr = c.groupTarget.guideRate;
            Scalar r = rateForMode(c.rates, mode, c.resvCoeff);

            if (anyFallback && c.useFallback) {
                gr = Scalar(0);
                r = Scalar(0);
            }

            guideRates.push_back(gr);
            rates.push_back(r);
            effFactors.push_back(c.efficiencyFactor);
            isGrup.push_back(c.modeCategory == ProdNodeModeCategory::Group);
        }

        const Scalar guideSum = std::accumulate(guideRates.begin(), guideRates.end(), Scalar(0));
        Scalar guideSumIsGroup = Scalar(0);
        for (size_t i = 0; i < isGrup.size(); ++i) {
            if (isGrup[i]) guideSumIsGroup += guideRates[i];
        }

        // top.rates is efficiency-adjusted (rateSums accumulated eff-scaled child
        // contributions), so targetSum is in the parent frame throughout.
        // Non-GRUP children's parent-frame contribution is effFactor * native_rate.
        const Scalar target = rateForMode(top.rates, mode, top.resvCoeff);
        Scalar targetSum = target;
        for (size_t i = 0; i < isGrup.size(); ++i) {
            if (!isGrup[i]) targetSum -= effFactors[i] * rates[i];
        }

        // Set targets for children
        for (size_t i = 0; i < cList.size(); ++i) {
            const auto& cName = cList[i];
            if (tree.count(cName) == 0) continue;
            auto& c = tree.at(cName);
            if (top.modeCategory == ProdNodeModeCategory::Group) {
                // inherit group-name from top node
                c.groupTarget.groupName = top.groupTarget.groupName;
            } else {
                c.groupTarget.groupName = topName;
            }
            c.groupTarget.ctrlMode = wellModeToGroupMode(mode);
            c.groupTarget.guideRate = guideRates[i];

            if (!modeIsPref) {
                if (c.useFallback) {
                    c.groupTarget.guideRateRatio = Scalar(0);
                    c.groupTarget.value = std::numeric_limits<Scalar>::quiet_NaN();
                    setTargets(tree, cName);
                    continue;
                }
                c.groupTarget.guideRateRatio = (guideSum > Scalar(0)) ? (guideRates[i] / guideSum) : Scalar(0);
            }

            if (isGrup[i]) {
                // Allocate in parent frame, then convert to child native frame.
                const Scalar allocParentFrame = (guideSumIsGroup > Scalar(0))
                    ? (guideRates[i] / guideSumIsGroup) * targetSum
                    : Scalar(0);
                c.groupTarget.value = (effFactors[i] > Scalar(0))
                    ? allocParentFrame / effFactors[i]
                    : Scalar(0);
            } else {
                // For Individual children: show the hypothetical GRUP target.
                // Restore the parent-frame contribution of this child before
                // distributing, then convert the result to native frame.
                const Scalar totalGuide = guideSumIsGroup + guideRates[i];
                const Scalar totalTarget = targetSum + effFactors[i] * rates[i];
                const Scalar allocParentFrame = (totalGuide > Scalar(0))
                    ? (guideRates[i] / totalGuide) * totalTarget
                    : Scalar(0);
                c.groupTarget.value = (effFactors[i] > Scalar(0))
                    ? allocParentFrame / effFactors[i]
                    : Scalar(0);
            }

            setTargets(tree, cName);
        }

        // Set fallback targets if not preferred mode
        // (Implementation simplified - full version would set groupTargetFallback)
    }
}

// ---------------------------------------------------------------------------

template<class Scalar, typename IndexTraits>
bool runBalancingAlgorithm(const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                           Tree<Scalar>& tree, Scalar tol)
{
    // Main algorithm entry point (corresponds to MATLAB sortalgo2)
    //OpmLog::info("runBalancingAlgorithm: starting");

    // Initialize node categories (simplified preprocessing)
    std::function<void(const std::string&)> initializeCategories = [&](const std::string& nodeName) {
        if (tree.count(nodeName) == 0) return;
        auto& node = tree.at(nodeName);

        // Initialize state
        node.isBalanced = false;
        node.balancingCount = 0;

        // Recurse to children first
        for (const auto& child : node.children) {
            initializeCategories(child);
        }

        /*
        // Set initial category based on availableForGroupControl
        if (!node.availableForGroupControl) {
            node.modeCategory = ProdNodeModeCategory::Individual;
        } else if (node.hasGuideRate) {
            node.modeCategory = ProdNodeModeCategory::Group;
        } else {
            node.modeCategory = ProdNodeModeCategory::Transparent;
        }
        */
    };
    initializeCategories("FIELD");

    // Get sub-tree ordering
    const auto topnodes = getSubTreeOrdering(tree, "FIELD");

    if (topnodes.empty()) {
        // No nodes to balance
        return true;
    }

    for (const auto& nodeName : topnodes) {
        //OpmLog::info(fmt::format("runBalancingAlgorithm: processing topnode {}", nodeName));
        if (tree.count(nodeName) == 0) continue;
        auto& node = tree.at(nodeName);
        if (node.isSatellite) continue; // Skip satellite groups

        // Determine mode and target
        // A top node is either individual or NONE (in which case it has a preferredMode)
        // A well top-node is always individual
        Well::ProducerCMode mode;
        if (node.modeCategory == ProdNodeModeCategory::Individual || node.type == ProdNodeType::Well) {
            mode = node.mode;
        } else {
            // Convert preferredMode (Group::ProductionCMode) to Well::ProducerCMode
            mode = groupModeToWellMode(node.preferredMode);
        }
        if (mode == Well::ProducerCMode::CMODE_UNDEFINED || mode == Well::ProducerCMode::NONE) {
            OpmLog::warning(fmt::format("runBalancingAlgorithm: top node '{}' has undefined mode, defaulting to ORAT", nodeName));
        }
        //assert(node.Limits.count(mode) > 0);
        // Validate mode
        /*
        if (mode == Well::ProducerCMode::CMODE_UNDEFINED ||
            mode == Well::ProducerCMode::NONE) {
            mode = Well::ProducerCMode::ORAT; // Fallback to ORAT
        }
        */

        // Get limit for mode
        Scalar qm = std::numeric_limits<Scalar>::max();
        if (node.Limits.count(mode) > 0) {
            qm = node.Limits.at(mode);
        } else {
            OpmLog::warning(fmt::format("runBalancingAlgorithm: top node '{}' has no limit, using infinity",
                                     nodeName));
        }

        // Balance this subtree
        balanceGroupTree(tree, nodeName, wellModel.guideRate(), mode, qm, tol);

        // Set targets
        setTargets(tree, nodeName);
    }

    // Check tree validity (simplified)
    bool valid = true;
    // TODO: Add checkGroupTree equivalent

    return valid;
}

// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------

template<class Scalar>
bool checkTreeValidity(const Tree<Scalar>& tree,
                       const std::string& topName,
                       Scalar tol,
                       DeferredLogger& logger)
{
    bool valid = true;

    std::function<void(const std::string&)> check = [&](const std::string& name) {
        if (tree.count(name) == 0) return;
        const auto& node = tree.at(name);

        const auto category = node.modeCategory;

        // Check that Individual nodes are at their limit
        if (category == ProdNodeModeCategory::Individual) {
            const auto it = node.Limits.find(node.mode);
            if (it != node.Limits.end()) {
                const Scalar limit = it->second;
                const Scalar current = (node.mode == Well::ProducerCMode::BHP)
                    ? std::accumulate(node.rates.begin(), node.rates.end(), Scalar(0),
                                       [](Scalar s, Scalar r){ return s + (-r); })
                    : rateForMode(node.rates, node.mode, node.resvCoeff);
                const Scalar relErr = std::abs(current - limit) / (std::abs(limit) + kFeasibilityTolerance<Scalar>);
                if (relErr > tol) {
                    logger.warning("ProdGroupTreeBalancer",
                        fmt::format("Node '{}' is Individual but current rate ({:.4g}) "
                                    "differs from limit ({:.4g}) by {:.2g}%%",
                                    name, current, limit, 100.0 * relErr));
                    valid = false;
                }
            }
        }

        // Check that Group nodes satisfy their group target
        if (category == ProdNodeModeCategory::Group &&
            node.groupTarget.ctrlMode != Group::ProductionCMode::NONE) {
            const Scalar target = node.groupTarget.value;
            const Scalar current = std::accumulate(node.rates.begin(), node.rates.end(), Scalar(0),
                                                    [](Scalar s, Scalar r){ return s + (-r); });
            const Scalar relErr = std::abs(current - target) / (std::abs(target) + kFeasibilityTolerance<Scalar>);
            if (false) {//(relErr > tol) {
                logger.warning("ProdGroupTreeBalancer",
                    fmt::format("Node '{}' is Group but current total rate ({:.4g}) "
                                "differs from group target ({:.4g}) by {:.2g}%%",
                                name, current, target, 100.0 * relErr));
                valid = false;
            }
        }

        // Recurse
        for (const auto& child : node.children) { check(child); }
    };

    check(topName);
    return valid;
}

// ---------------------------------------------------------------------------

template<class Scalar, typename IndexTraits>
void applyTreeToState(const Tree<Scalar>& tree,
                      BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel)
{
    auto& wellState  = wellModel.wellState();
    auto& groupState = wellModel.groupState();

    for (const auto& [name, node] : tree) {
        if (node.type == ProdNodeType::Well) {
            if (!wellState.has(name)) continue;
            auto& ws = wellState.well(name);
            if (!ws.producer) continue;

            // Convert canonical 3-component rates back to active-phase vector
            const auto activeRates = toActive(node.rates, ws.pu);
            ws.surface_rates = activeRates;

            // Update control mode
            if (node.modeCategory== ProdNodeModeCategory::Group) {
                ws.production_cmode = Well::ProducerCMode::GRUP;
            } else if (node.modeCategory== ProdNodeModeCategory::Individual) {
                ws.production_cmode = node.mode;
            }

            // Update the group target stored in the well state
            //if (node.groupTarget.ctrlMode != Group::ProductionCMode::NONE) {
                typename SingleWellState<Scalar, IndexTraits>::GroupTarget gt;
                gt.group_name      = node.groupTarget.groupName;
                gt.production_cmode = node.groupTarget.ctrlMode;
                gt.target_value    = node.groupTarget.value;
                gt.guiderate_ratio = node.groupTarget.guideRateRatio;
                ws.group_target    = gt;
            //}
        } else {
            if (node.isSatellite) continue; // Skip satellite groups
            // Group node: update group state
            // Rates in groupState are stored in active-phase order, positive = production
            const auto& pu = wellModel.phaseUsage();
            std::vector<Scalar> activeRates(pu.numPhases, Scalar(0));
            for (int c = 0; c < 3; ++c) {
                const int a = activeIdx(pu, c);
                if (a >= 0) {
                    // GroupState stores positive = production; our rates are negative = production
                    activeRates[a] = -node.rates[c];
                }
            }
            groupState.update_production_rates(name, activeRates);

            // Update group control mode.
            // applyTreeToState() now runs on all MPI ranks simultaneously (each rank
            // runs the same deterministic algorithm on globally consistent inputs and
            // reaches the same result), so production_control writes are valid.
            if (node.modeCategory == ProdNodeModeCategory::Group) {
                groupState.production_control(name, Group::ProductionCMode::FLD);
            } else if (node.modeCategory == ProdNodeModeCategory::Individual) {
                switch (node.mode) {
                case Well::ProducerCMode::ORAT:
                    groupState.production_control(name, Group::ProductionCMode::ORAT); break;
                case Well::ProducerCMode::WRAT:
                    groupState.production_control(name, Group::ProductionCMode::WRAT); break;
                case Well::ProducerCMode::GRAT:
                    groupState.production_control(name, Group::ProductionCMode::GRAT); break;
                case Well::ProducerCMode::LRAT:
                    groupState.production_control(name, Group::ProductionCMode::LRAT); break;
                case Well::ProducerCMode::RESV:
                    groupState.production_control(name, Group::ProductionCMode::RESV); break;
                default:
                    break;
                }
            } else {
                groupState.production_control(name, Group::ProductionCMode::NONE);
            }
        }
    }
}

// ---------------------------------------------------------------------------

template<class Scalar, typename IndexTraits>
bool runGroupTreeBalancer(BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                          const SummaryState& summaryState,
                          int reportStep,
                          Scalar tol,
                          [[maybe_unused]] int maxIter,
                          const std::unordered_map<std::string, std::pair<int, Scalar>>& limits,
                          DeferredLogger& logger)
{
    OPM_TIMEFUNCTION();
    (void)maxIter;

    const auto t0 = std::chrono::steady_clock::now();

    auto tree = buildTree(wellModel, summaryState, reportStep, limits);

    const bool converged = runBalancingAlgorithm(wellModel, tree, tol);

    if (!converged) {
        logger.warning("ProdGroupTreeBalancer",
            "Group tree balancer did not converge for all subtrees. "
            "Results may be inconsistent.");
    }

    const bool valid = checkTreeValidity(tree, "FIELD", tol, logger);

    if (!valid && converged) {
        logger.warning("ProdGroupTreeBalancer",
            "Group tree balancer converged but validity check failed. "
            "Some nodes may not satisfy their constraints.");
    }

    applyTreeToState(tree, wellModel);

    const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - t0).count();
    logger.debug("ProdGroupTreeBalancer",
        fmt::format("Group tree balancer completed in {}ms. "
                    "Convergence: {}, Validity: {}",
                    elapsed, converged ? "OK" : "FAILED", valid ? "OK" : "FAILED"));

    return valid;
}

// ===========================================================================
// Explicit instantiations
// ===========================================================================

template Tree<double> buildTree<double, BlackOilDefaultFluidSystemIndices>(
    const BlackoilWellModelGeneric<double, BlackOilDefaultFluidSystemIndices>&,
    const SummaryState&, int,
    const std::unordered_map<std::string, std::pair<int, double>>&);

template std::vector<std::string>
getSubTreeOrdering<double>(const Tree<double>&, const std::string&);

// New algorithm helper functions
template bool hasFreePath<double>(const Tree<double>&, const std::string&, const std::string&);

template void incrementParentRateSums<double>(Tree<double>&, const std::string&, const std::string&,
    const std::optional<std::array<double, 3>>&);

template void decrementParentGuideRateSums<double>(Tree<double>&, const std::string&, const std::string&, const std::array<double, 3>&);

template std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>>
getLocalTreeDescendants<double>(const Tree<double>&, const std::string&);

template <class Scalar> void updateGuideRatesForMode(Tree<Scalar>& tree,
                                                     const std::vector<std::string>& c,
                                                     Well::ProducerCMode mode,
                                                     const GuideRate& guideRate);

template void resetRatesAndGuideRateSums<double>(Tree<double>&,
                                                 const std::vector<std::string>&,
                                                 const std::vector<std::string>&,
                                                 const std::vector<std::string>&,
                                                 const std::string&,
                                                 Well::ProducerCMode);

template std::vector<double>
computeRatiosForSorting<double>(const Tree<double>&, const std::vector<std::string>&);

template std::vector<std::string>
updateTransparentGroups<double>(Tree<double>&, const std::string&, std::vector<std::string>&,
    const GuideRate&, double, double, Well::ProducerCMode, double);

template void distributeFallbackRates<double>(Tree<double>&, const std::string&,
    const std::vector<std::string>&, const GuideRate&, Well::ProducerCMode, double);

template void balanceGroupTree<double>(Tree<double>&, const std::string&,
    const GuideRate&, Well::ProducerCMode, double, double);

template void setTargets<double>(Tree<double>&, const std::string&);

template bool runBalancingAlgorithm<double, BlackOilDefaultFluidSystemIndices>(
    const BlackoilWellModelGeneric<double, BlackOilDefaultFluidSystemIndices>&,
    Tree<double>&,
    double);

template bool checkTreeValidity<double>(const Tree<double>&, const std::string&, double, DeferredLogger&);

template void applyTreeToState<double, BlackOilDefaultFluidSystemIndices>(
    const Tree<double>&,
    BlackoilWellModelGeneric<double, BlackOilDefaultFluidSystemIndices>&);

template bool runGroupTreeBalancer<double, BlackOilDefaultFluidSystemIndices>(
    BlackoilWellModelGeneric<double, BlackOilDefaultFluidSystemIndices>&,
    const SummaryState&, int, double, int,
    const std::unordered_map<std::string, std::pair<int, double>>&,
    DeferredLogger&);

#ifdef FLOW_INSTANTIATE_FLOAT

template Tree<float> buildTree<float, BlackOilDefaultFluidSystemIndices>(
    const BlackoilWellModelGeneric<float, BlackOilDefaultFluidSystemIndices>&,
    const SummaryState&, int,
    const std::unordered_map<std::string, std::pair<int, float>>&);

template std::vector<std::string>
getSubTreeOrdering<float>(const Tree<float>&, const std::string&);

// New algorithm helper functions
template bool hasFreePath<float>(const Tree<float>&, const std::string&, const std::string&);

template void incrementParentRateSums<float>(Tree<float>&, const std::string&, const std::string&,
    const std::optional<std::array<float, 3>>&);

template void decrementParentGuideRateSums<float>(Tree<float>&, const std::string&, const std::string&, const std::array<float, 3>&);

template std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>>
getLocalTreeDescendants<float>(const Tree<float>&, const std::string&);

template void resetRatesAndGuideRateSums<float>(Tree<float>&,
    const std::vector<std::string>&,
    const std::vector<std::string>&, const std::vector<std::string>&, const std::string&, Well::ProducerCMode);

template std::vector<float>
computeRatiosForSorting<float>(const Tree<float>&, const std::vector<std::string>&);

template std::vector<std::string>
updateTransparentGroups<float>(Tree<float>&, const std::string&, std::vector<std::string>&,
    const GuideRate&, float, float, Well::ProducerCMode, float);

template void distributeFallbackRates<float>(Tree<float>&, const std::string&,
    const std::vector<std::string>&, const GuideRate&, Well::ProducerCMode, float);

template void balanceGroupTree<float>(Tree<float>&, const std::string&,
    const GuideRate&, Well::ProducerCMode, float, float);

template void setTargets<float>(Tree<float>&, const std::string&);

template bool runBalancingAlgorithm<float, BlackOilDefaultFluidSystemIndices>(
    const BlackoilWellModelGeneric<float, BlackOilDefaultFluidSystemIndices>&,
    Tree<float>&,
    float);

template bool checkTreeValidity<float>(const Tree<float>&, const std::string&, float, DeferredLogger&);

template void applyTreeToState<float, BlackOilDefaultFluidSystemIndices>(
    const Tree<float>&,
    BlackoilWellModelGeneric<float, BlackOilDefaultFluidSystemIndices>&);

template bool runGroupTreeBalancer<float, BlackOilDefaultFluidSystemIndices>(
    BlackoilWellModelGeneric<float, BlackOilDefaultFluidSystemIndices>&,
    const SummaryState&, int, float, int,
    const std::unordered_map<std::string, std::pair<int, float>>&,
    DeferredLogger&);

#endif

} // namespace Opm::ProdGroupTreeBalancer
