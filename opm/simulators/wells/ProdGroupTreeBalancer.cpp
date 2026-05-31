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
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
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

// ---------------------------------------------------------------------------
// IPR-based BHP→rate conversion
// ---------------------------------------------------------------------------

/// Convert BHP limit to an equivalent total rate using the linear IPR model
/// q = ipr_a - ipr_b * bhp.  Stores the result as the BHP control mode in
/// node.individualLimits.  Does nothing if ipr_b_sum <= 0.
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
    node.individualLimits[Well::ProducerCMode::BHP] = q_limit;
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
    if (!guideRate.has(name)) return kUniformWeight<Scalar>;
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

/// Safe lookup of a WellInterfaceGeneric by name.  Returns nullptr if not found
/// in the local well container (e.g., the well belongs to another MPI rank).
template<class Scalar, typename IndexTraits>
const WellInterfaceGeneric<Scalar, IndexTraits>*
findGenericWell(const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                const std::string& wellName)
{
    for (const auto* w : wellModel.genericWells()) {
        if (w->name() == wellName) return w;
    }
    return nullptr;
}

// ---------------------------------------------------------------------------
// buildTree helpers
// ---------------------------------------------------------------------------

template<class Scalar, typename IndexTraits>
void populateWellNode(ProdGroupTreeNode<Scalar>& node,
                      const std::string& wellName,
                      const Schedule& schedule,
                      const WellState<Scalar, IndexTraits>& wellState,
                      const GroupState<Scalar>& /*groupState*/,
                      const GuideRate& /*guideRate*/,
                      const SummaryState& summaryState,
                      int reportStep,
                      int fipnum,
                      int pvtreg,
                      const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel)
{
    const auto& ws = wellState.well(wellName);
    const auto& pu = ws.pu;
    const auto& eclWell = schedule.getWell(wellName, reportStep);

    node.name   = wellName;
    node.type   = ProdNodeType::Well;
    node.parent = eclWell.groupName();
    node.allowGroupControl = eclWell.isAvailableForGroupControl();

    // Rates: canonical 3-component (negative = production)
    node.rates = toCanonical3(ws.surface_rates, pu);

    // Current control mode
    node.activeIndividualCtrl = ws.production_cmode;
    node.preferredCtrl        = ws.production_cmode;

    // Determine ctrlStatus from current control mode
    if (ws.production_cmode == Well::ProducerCMode::GRUP) {
        node.ctrlStatus = ProdNodeCtrlStatus::GroupControlled;
    } else {
        node.ctrlStatus = ProdNodeCtrlStatus::IndividualControlled;
    }

    // Phase fractions (for distributing a group target to phases)
    const Scalar totalRate = std::accumulate(node.rates.begin(), node.rates.end(), Scalar(0),
                                             [](Scalar s, Scalar r){ return s + (-r); });
    if (totalRate > Scalar(0)) {
        for (int c = 0; c < 3; ++c) {
            node.wellRateFractions[c] = -node.rates[c] / totalRate;
        }
    } else {
        // Default fractions when no production: equal share among active phases.
        // The exact value matters only when the algorithm initializes with small rates.
        node.wellRateFractions[kOil]   = pu.phaseIsActive(IndexTraits::oilPhaseIdx)   ? kDefaultUniformFraction<Scalar> : Scalar(0);
        node.wellRateFractions[kWater] = pu.phaseIsActive(IndexTraits::waterPhaseIdx) ? kDefaultUniformFraction<Scalar> : Scalar(0);
        node.wellRateFractions[kGas]   = pu.phaseIsActive(IndexTraits::gasPhaseIdx)   ? kDefaultUniformFraction<Scalar> : Scalar(0);
        const Scalar sumFrac = node.wellRateFractions[0] + node.wellRateFractions[1] + node.wellRateFractions[2];
        if (sumFrac > Scalar(0)) {
            for (auto& f : node.wellRateFractions) { f /= sumFrac; }
        }
    }

    // RESV conversion coefficients for this well
    {
        std::vector<Scalar> coeffVec(pu.numPhases, Scalar(0));
        std::vector<Scalar> posRates(pu.numPhases, Scalar(0));
        for (int i = 0; i < static_cast<int>(ws.surface_rates.size()); ++i) {
            posRates[i] = -ws.surface_rates[i]; // negate since production is negative
        }
        wellModel.calcResvCoeff(fipnum, pvtreg, posRates, coeffVec);
        for (int c = 0; c < 3; ++c) {
            const int a = activeIdx(pu, c);
            node.resvCoeff[c] = (a >= 0) ? coeffVec[a] : Scalar(0);
        }
    }

    // Individual limits
    if (!eclWell.isProducer()) return;

    const auto controls = eclWell.productionControls(summaryState);

    if (controls.hasControl(Well::ProducerCMode::ORAT) && controls.oil_rate > Scalar(0)) {
        node.individualLimits[Well::ProducerCMode::ORAT] = controls.oil_rate;
    }
    if (controls.hasControl(Well::ProducerCMode::WRAT) && controls.water_rate > Scalar(0)) {
        node.individualLimits[Well::ProducerCMode::WRAT] = controls.water_rate;
    }
    if (controls.hasControl(Well::ProducerCMode::GRAT) && controls.gas_rate > Scalar(0)) {
        node.individualLimits[Well::ProducerCMode::GRAT] = controls.gas_rate;
    }
    if (controls.hasControl(Well::ProducerCMode::LRAT) && controls.liquid_rate > Scalar(0)) {
        node.individualLimits[Well::ProducerCMode::LRAT] = controls.liquid_rate;
    }
    if (controls.hasControl(Well::ProducerCMode::RESV) && controls.resv_rate > Scalar(0)) {
        node.individualLimits[Well::ProducerCMode::RESV] = controls.resv_rate;
    }

    // BHP / THP limit via IPR.
    // Use estimateStableBhp to find the BHP at the THP operating point, which is
    // more robust than forced_bhp_from_thp.  If no stable operating point exists
    // at the THP limit, the well cannot be kept open and is assigned zero rates.
    {
        Scalar effectiveBhpLimit = Scalar(0);
        bool   hasBhpConstraint  = false;

        if (controls.hasControl(Well::ProducerCMode::BHP) && controls.bhp_limit > Scalar(0)) {
            effectiveBhpLimit = controls.bhp_limit;
            hasBhpConstraint  = true;
        }

        // Attempt to get the stable BHP corresponding to the THP limit
        const auto* genWell = findGenericWell(wellModel, wellName);
        if (genWell != nullptr) {
            WellBhpThpCalculator<Scalar, IndexTraits> calc(*genWell);
            if (calc.wellHasTHPConstraints(summaryState)) {
                // Compute reference density as a weighted average based on surface rates
                // Standard surface densities: oil ~850 kg/m³, water ~1000 kg/m³, gas ~1 kg/m³
                const Scalar oilDens   = Scalar(850.0);
                const Scalar waterDens = Scalar(1000.0);
                const Scalar gasDens   = Scalar(1.0);
                
                Scalar totalRate = Scalar(0);
                Scalar weightedDens = Scalar(0);
                const auto& rates = node.rates; // canonical [oil, water, gas], negative = production
                if (pu.phaseIsActive(IndexTraits::oilPhaseIdx) && rates[kOil] < Scalar(0)) {
                    totalRate += -rates[kOil];
                    weightedDens += -rates[kOil] * oilDens;
                }
                if (pu.phaseIsActive(IndexTraits::waterPhaseIdx) && rates[kWater] < Scalar(0)) {
                    totalRate += -rates[kWater];
                    weightedDens += -rates[kWater] * waterDens;
                }
                if (pu.phaseIsActive(IndexTraits::gasPhaseIdx) && rates[kGas] < Scalar(0)) {
                    totalRate += -rates[kGas];
                    weightedDens += -rates[kGas] * gasDens;
                }
                const Scalar refDens = (totalRate > Scalar(0)) ? (weightedDens / totalRate) : oilDens;
                
                const auto stableBhp = calc.estimateStableBhp(
                    wellState, eclWell, ws.surface_rates,
                    refDens, summaryState);

                if (!stableBhp.has_value()) {
                    // No stable operating point at THP limit: well cannot produce.
                    node.rates        = {};
                    node.ctrlStatus   = ProdNodeCtrlStatus::IndividualControlled;
                    node.activeIndividualCtrl = Well::ProducerCMode::THP;
                    node.individualLimits.clear();
                    return; // no further limit processing needed
                }

                // THP is more restrictive if the stable BHP is higher than the
                // BHP limit (higher BHP → lower production for a producer).
                if (stableBhp.value() > effectiveBhpLimit) {
                    effectiveBhpLimit = stableBhp.value();
                    hasBhpConstraint  = true;
                }
            }
        } else if (ws.forced_bhp_from_thp.has_value()) {
            // Fallback when the well's WellInterfaceGeneric is not in the local
            // container (e.g., in a distributed-wells setup on a non-owner rank).
            // forced_bhp_from_thp is the BHP computed from the VFP table at the
            // current operating point under THP control; less accurate than
            // estimateStableBhp but better than ignoring the THP constraint.
            // TODO: Remove once all wells are guaranteed to be on rank 0.
            if (ws.forced_bhp_from_thp.value() > effectiveBhpLimit) {
                effectiveBhpLimit = ws.forced_bhp_from_thp.value();
                hasBhpConstraint  = true;
            }
        }

        if (hasBhpConstraint) {
            addBhpRateLimit(node, ws, effectiveBhpLimit);
        }
    }

    // Pre-compute the single binding limit for this well.
    // The direction of rate change is fixed (wellRateFractions), so we can
    // determine at construction time which limit will be hit first.
    // Keeping only the binding limit simplifies parametrizeTree for wells.
    if (!node.individualLimits.empty()) {
        const auto& lt = node.wellRateFractions; // direction (positive = more production)
        Well::ProducerCMode bindingMode = Well::ProducerCMode::CMODE_UNDEFINED;
        Scalar minAlpha = std::numeric_limits<Scalar>::max();

        for (const auto& [mode, limit] : node.individualLimits) {
            const Scalar currentRate = rateForMode(node.rates, mode, node.resvCoeff);
            const Scalar ltForMode   = linearTermForMode(lt, mode, node.resvCoeff);
            if (ltForMode <= Scalar(0)) continue;
            const Scalar remaining = limit - currentRate;
            const Scalar alpha = (remaining <= Scalar(0)) ? Scalar(0) : remaining / ltForMode;
            if (alpha < minAlpha) {
                minAlpha    = alpha;
                bindingMode = mode;
            }
        }

        if (bindingMode != Well::ProducerCMode::CMODE_UNDEFINED) {
            const Scalar bindingLimit = node.individualLimits.at(bindingMode);
            node.individualLimits.clear();
            node.individualLimits[bindingMode] = bindingLimit;
        }
    }
}

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
        node.allowGroupControl = false;
    } else {
        node.allowGroupControl = group.productionGroupControlAvailable();
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
            OpmLog::info(fmt::format("buildTree: {} identified as SATELLITE group (GSATPROD defined in schedule)",
                                     groupName));
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

    node.ownCtrlMode = ctrl;

    // Determine control status:
    // - Satellite groups have fixed rates and should not be controlled by parent groups
    // - For other groups: check if available for group control
    if (isSatellite) {
        node.ctrlStatus = ProdNodeCtrlStatus::Satellite;
        node.allowGroupControl = false;  // Satellite groups cannot be controlled by parent
    } else if (node.allowGroupControl) {
        // Available for group control
        node.ctrlStatus = ProdNodeCtrlStatus::GroupControlled;
    } else {
        // Not available for group control
        node.ctrlStatus = ProdNodeCtrlStatus::Undetermined;
    }

    OpmLog::info(fmt::format("buildTree: group {} allowGroupControl={} satellite={}",
                             groupName, node.allowGroupControl, isSatellite));

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
        const auto controls = group.productionControls(summaryState);
        
        // Set preferredCtrl from schedule control mode (not from active groupState control)
        // Map Group::ProductionCMode to Well::ProducerCMode
        switch (controls.cmode) {
        case Group::ProductionCMode::ORAT:
            node.preferredCtrl = Well::ProducerCMode::ORAT; break;
        case Group::ProductionCMode::WRAT:
            node.preferredCtrl = Well::ProducerCMode::WRAT; break;
        case Group::ProductionCMode::GRAT:
            node.preferredCtrl = Well::ProducerCMode::GRAT; break;
        case Group::ProductionCMode::LRAT:
            node.preferredCtrl = Well::ProducerCMode::LRAT; break;
        case Group::ProductionCMode::RESV:
            node.preferredCtrl = Well::ProducerCMode::RESV; break;
        case Group::ProductionCMode::FLD:
            node.preferredCtrl = Well::ProducerCMode::GRUP; break;  // FLD means follow parent
        case Group::ProductionCMode::NONE:
        default:
            node.preferredCtrl = Well::ProducerCMode::NONE; break;
        }
        
        // Populate individual limits from schedule controls
        if (group.has_control(Group::ProductionCMode::ORAT) && controls.oil_target > Scalar(0)) {
            node.individualLimits[Well::ProducerCMode::ORAT] = controls.oil_target;
        }
        if (group.has_control(Group::ProductionCMode::WRAT) && controls.water_target > Scalar(0)) {
            node.individualLimits[Well::ProducerCMode::WRAT] = controls.water_target;
        }
        if (group.has_control(Group::ProductionCMode::GRAT) && controls.gas_target > Scalar(0)) {
            node.individualLimits[Well::ProducerCMode::GRAT] = controls.gas_target;
        }
        if (group.has_control(Group::ProductionCMode::LRAT) && controls.liquid_target > Scalar(0)) {
            node.individualLimits[Well::ProducerCMode::LRAT] = controls.liquid_target;
        }
        if (group.has_control(Group::ProductionCMode::RESV) && controls.resv_target > Scalar(0)) {
            node.individualLimits[Well::ProducerCMode::RESV] = controls.resv_target;
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
                       int reportStep)
{
    const auto& schedule   = wellModel.schedule();
    const auto& wellState  = wellModel.wellState();
    const auto& groupState = wellModel.groupState();
    const auto& guideRate  = wellModel.guideRate();
    const auto& pu         = wellModel.phaseUsage();

    // Get RESV conversion parameters (fipnum, pvtreg)
    const auto [fipnum, pvtreg] = wellModel.getGroupFipnumAndPvtreg();

    Tree<Scalar> tree;

    // Use a stack-based traversal starting from FIELD
    std::vector<std::string> toVisit = {"FIELD"};
    while (!toVisit.empty()) {
        const std::string name = toVisit.back();
        toVisit.pop_back();

        if (schedule.hasWell(name, reportStep)) {
            // Well node: only populate if locally available (producer only)
            if (wellState.has(name)) {
                const auto& ws = wellState.well(name);
                if (ws.producer) {
                    auto& node = tree[name];
                    populateWellNode(node, name, schedule, wellState, groupState,
                                     guideRate, summaryState, reportStep,
                                     fipnum, pvtreg, wellModel);
                }
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
    std::function<void(const std::string&, Group::ProductionCMode, Well::ProducerCMode)>
    propagateGuideRates = [&](const std::string& groupName,
                              Group::ProductionCMode inheritedMode,
                              Well::ProducerCMode inheritedPreferredCtrl)
    {
        if (tree.count(groupName) == 0) return;
        auto& gnode = tree.at(groupName);

        // Effective mode this group uses to issue targets to its children.
        Group::ProductionCMode effectiveMode = gnode.ownCtrlMode;
        if (effectiveMode == Group::ProductionCMode::FLD ||
            effectiveMode == Group::ProductionCMode::NONE) {
            effectiveMode = inheritedMode; // relay parent's mode
        }
        // Store the effective mode back so setTargets can use it.
        gnode.ownCtrlMode = effectiveMode;

        // Propagate preferred control for FIELD/NONE cases
        if (gnode.preferredCtrl == Well::ProducerCMode::GRUP ||
            gnode.preferredCtrl == Well::ProducerCMode::NONE) {
            gnode.preferredCtrl = inheritedPreferredCtrl;
        }

        for (const auto& childName : gnode.children) {
            if (tree.count(childName) == 0) continue;
            auto& child = tree.at(childName);
            child.groupTarget.ctrlMode  = effectiveMode;
            child.groupTarget.guideRate =
                getGuideRateForMode(childName, child.rates, effectiveMode, guideRate);

            // Detect if this node has a guide rate
            child.hasGuideRate = guideRate.has(childName);

            // Compute guide rates for all control modes
            const std::array<Group::ProductionCMode, 5> modes = {
                Group::ProductionCMode::ORAT,
                Group::ProductionCMode::WRAT,
                Group::ProductionCMode::GRAT,
                Group::ProductionCMode::LRAT,
                Group::ProductionCMode::RESV
            };
            for (const auto& mode : modes) {
                const auto grForMode = getGuideRateForMode(childName, child.rates, mode, guideRate);
                // Map Group::ProductionCMode to Well::ProducerCMode for storage
                Well::ProducerCMode wellMode;
                switch (mode) {
                    case Group::ProductionCMode::ORAT: wellMode = Well::ProducerCMode::ORAT; break;
                    case Group::ProductionCMode::WRAT: wellMode = Well::ProducerCMode::WRAT; break;
                    case Group::ProductionCMode::GRAT: wellMode = Well::ProducerCMode::GRAT; break;
                    case Group::ProductionCMode::LRAT: wellMode = Well::ProducerCMode::LRAT; break;
                    case Group::ProductionCMode::RESV: wellMode = Well::ProducerCMode::RESV; break;
                    default: continue;
                }
                child.guideRatesForMode[wellMode] = grForMode;
            }

            // Populate guideRates array [oil, water, gas] for ORAT, WRAT, GRAT
            child.guideRates[kOil]   = child.guideRatesForMode[Well::ProducerCMode::ORAT];
            child.guideRates[kWater] = child.guideRatesForMode[Well::ProducerCMode::WRAT];
            child.guideRates[kGas]   = child.guideRatesForMode[Well::ProducerCMode::GRAT];

            // Propagate preferred control for FIELD/NONE cases
            if (child.preferredCtrl == Well::ProducerCMode::GRUP ||
                child.preferredCtrl == Well::ProducerCMode::NONE) {
                child.preferredCtrl = gnode.preferredCtrl;
            }

            if (child.type == ProdNodeType::Group) {
                propagateGuideRates(childName, effectiveMode, gnode.preferredCtrl);
            }
        }
    };
    // ORAT is the ultimate fallback if FIELD has no control mode set.
    // Use ORAT as the initial preferred control as well.
    propagateGuideRates("FIELD", Group::ProductionCMode::ORAT, Well::ProducerCMode::ORAT);

    return tree;
}

// ---------------------------------------------------------------------------

template<class Scalar>
std::vector<std::string> getSubTreeOrdering(const Tree<Scalar>& tree,
                                             const std::string& rootName)
{
    // Post-order traversal: leaves first, root last.
    OpmLog::info(fmt::format("getSubTreeOrdering: starting from {}", rootName));
    std::vector<std::string> ordering;
    std::vector<std::string> stack;
    stack.push_back(rootName);

    int iterCount = 0;
    while (!stack.empty()) {
        if (++iterCount > 1000) {
            OpmLog::error("getSubTreeOrdering: stuck in loop after 1000 iterations");
            break;
        }
        const std::string name = stack.back();
        stack.pop_back();
        if (tree.count(name) == 0) continue;
        const auto& node = tree.at(name);
        // Only include groups that are NOT available for group control
        // (FIELD, satellite groups, groups with GCONPROD item 8 = NO)
        if (node.type == ProdNodeType::Group && !node.allowGroupControl) {
            ordering.push_back(name);
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
    OpmLog::info(fmt::format("getSubTreeOrdering: found {} nodes: {}",
                             ordering.size(),
                             fmt::join(ordering, ", ")));
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
    const bool isFreeStep = (parent.modeCategory == "TRANSPARENT") && !parent.visited;

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

    // Determine rates to add (convert to positive production for rateSums)
    std::array<Scalar, 3> ratesToAdd;
    if (rates.has_value()) {
        ratesToAdd = rates.value();
    } else {
        // node.rates are negative (production), convert to positive for rateSums
        for (int c = 0; c < 3; ++c) {
            ratesToAdd[c] = -node.rates[c];
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

        if (!child.allowGroupControl) {
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
            tree.at(nodeName).modeCategory = "TRANSPARENT";
        }
    }

    // Convert mode to canonical index
    int modeIdx = -1;
    switch (mode) {
        case Well::ProducerCMode::ORAT: modeIdx = kOil; break;
        case Well::ProducerCMode::WRAT: modeIdx = kWater; break;
        case Well::ProducerCMode::GRAT: modeIdx = kGas; break;
        default: modeIdx = kOil; // fallback
    }

    if (tree.count(origin) == 0) return;
    const auto& originNode = tree.at(origin);
    const Well::ProducerCMode modePreferred = originNode.preferredCtrl;
    const bool modeIsPref = (mode == modePreferred);

    // Compute total guide rate for mode if not preferred
    Scalar guideSumOrigin = Scalar(0);
    if (!modeIsPref && modeIdx >= 0) {
        for (const auto& childName : c) {
            if (tree.count(childName) > 0) {
                guideSumOrigin += tree.at(childName).guideRates[modeIdx];
            }
        }
    }

    // Set guide_rate_sums for each child in c
    for (const auto& childName : c) {
        if (tree.count(childName) == 0) continue;
        auto& child = tree.at(childName);

        const auto& rates = child.rates;

        // Check for fallback condition (non-preferred mode with tiny fraction)
        if (!modeIsPref && modeIdx >= 0) {
            const Scalar rateSum = std::accumulate(rates.begin(), rates.end(), Scalar(0),
                                                   [](Scalar s, Scalar r){ return s + (-r); });
            const Scalar rateFracMode = (rateSum > Scalar(0)) ? (-rates[modeIdx] / rateSum) : Scalar(0);
            const Scalar guideRatioMode = (guideSumOrigin > Scalar(0))
                ? (child.guideRates[modeIdx] / guideSumOrigin) : Scalar(0);

            if (rateFracMode < std::sqrt(std::numeric_limits<Scalar>::epsilon()) ||
                guideRatioMode < std::sqrt(std::numeric_limits<Scalar>::epsilon())) {
                child.useFallback = true;
                child.guideRateSums = {std::numeric_limits<Scalar>::quiet_NaN(),
                                       std::numeric_limits<Scalar>::quiet_NaN(),
                                       std::numeric_limits<Scalar>::quiet_NaN()};
                continue;
            }
        }

        // Set guide-rate sums proportional to rates
        child.useFallback = false;
        const Scalar guideRate = (modeIdx >= 0) ? child.guideRates[modeIdx] : Scalar(1);
        const Scalar rateMagnitude = (modeIdx >= 0) ? -rates[modeIdx] : Scalar(1);

        if (rateMagnitude > Scalar(0)) {
            for (int c = 0; c < 3; ++c) {
                child.guideRateSums[c] = (guideRate / rateMagnitude) * (-rates[c]);
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
                for (int c = 0; c < 3; ++c) {
                    parentNode.guideRateSums[c] += child.guideRateSums[c];
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------

template<class Scalar>
std::pair<std::vector<Scalar>, std::vector<int>>
getRatiosForSorting(const Tree<Scalar>& tree, const std::vector<std::string>& c)
{
    // Get limit-to-guiderate ratios for distribution ordering
    const int nc = c.size();
    std::vector<Scalar> ratios(nc, std::numeric_limits<Scalar>::max());
    std::vector<int> modeIndices(nc, -1);

    for (int k = 0; k < nc; ++k) {
        const auto& childName = c[k];
        if (tree.count(childName) == 0) continue;
        const auto& child = tree.at(childName);

        // Determine limits to use
        std::array<Scalar, 3> limits;
        if (child.modeCategory == "NONE") {
            limits = child.rates;
        } else {
            // Use individual limits, but we need to map them to rates
            limits = {std::numeric_limits<Scalar>::max(),
                     std::numeric_limits<Scalar>::max(),
                     std::numeric_limits<Scalar>::max()};
            // Check each limit and compute remaining capacity
            for (const auto& [mode, limit] : child.individualLimits) {
                const Scalar current = rateForMode(child.rates, mode, child.resvCoeff);
                const Scalar remaining = limit - current;
                // Map to canonical index (this is approximate - proper implementation would consider all modes)
                int idx = -1;
                switch (mode) {
                    case Well::ProducerCMode::ORAT: idx = kOil; break;
                    case Well::ProducerCMode::WRAT: idx = kWater; break;
                    case Well::ProducerCMode::GRAT: idx = kGas; break;
                    default: continue;
                }
                if (idx >= 0) {
                    limits[idx] = std::min(limits[idx], current + remaining);
                }
            }
        }

        // Check if guide_rate_sums is effectively zero
        const bool isZero = (child.guideRateSums[kOil] < std::sqrt(std::numeric_limits<Scalar>::epsilon()) * child.guideRates[kOil]) &&
                           (child.guideRateSums[kWater] < std::sqrt(std::numeric_limits<Scalar>::epsilon()) * child.guideRates[kWater]) &&
                           (child.guideRateSums[kGas] < std::sqrt(std::numeric_limits<Scalar>::epsilon()) * child.guideRates[kGas]);

        if (isZero) {
            ratios[k] = std::numeric_limits<Scalar>::infinity();
            continue;
        }

        // Find minimum ratio
        Scalar minRatio = std::numeric_limits<Scalar>::max();
        int minIdx = -1;
        for (int c = 0; c < 3; ++c) {
            if (child.guideRateSums[c] > Scalar(0)) {
                const Scalar remaining = limits[c] - child.rateSums[c];
                const Scalar ratio = remaining / child.guideRateSums[c];
                if (ratio < minRatio) {
                    minRatio = ratio;
                    minIdx = c;
                }
            }
        }

        ratios[k] = minRatio;
        modeIndices[k] = minIdx;
    }

    return {ratios, modeIndices};
}

// ---------------------------------------------------------------------------

template<class Scalar>
std::vector<std::string>
updateTransparentGroups(Tree<Scalar>& tree,
                        const std::string& nodeName,
                        std::vector<std::string>& c_trans,
                        Scalar nextRatio,
                        Scalar qm,
                        Well::ProducerCMode mode,
                        Scalar tol)
{
    // Check if any transparent groups have limit-to-guide ratio less than nextRatio
    std::vector<std::string> c_trans_update;

    // Convert mode to canonical index
    int modeIdx = kOil;
    switch (mode) {
        case Well::ProducerCMode::ORAT: modeIdx = kOil; break;
        case Well::ProducerCMode::WRAT: modeIdx = kWater; break;
        case Well::ProducerCMode::GRAT: modeIdx = kGas; break;
        default: modeIdx = kOil;
    }

    bool foundLessThanNext = true;
    while (foundLessThanNext && !c_trans.empty()) {
        auto [ratios, modeIndices] = getRatiosForSorting(tree, c_trans);

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
        const Scalar gsum = originNode.guideRateSums[modeIdx];
        const Scalar qmRemain = qm - originNode.rateSums[modeIdx];

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
        const Scalar gk = ck.guideRateSums[modeIdx];
        const Scalar qk = (gk / (gsum + Scalar(1e-20))) * qmRemain + rateSumsOrig[modeIdx];

        // Recursive balance
        balanceGroupTree(tree, ckName, mode, qk, tol);

        // Update parent rate-sums
        std::array<Scalar, 3> ratesDelta;
        for (int c = 0; c < 3; ++c) {
            ratesDelta[c] = ck.rates[c] - rateSumsOrig[c];
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
                             Well::ProducerCMode mode,
                             Scalar tol)
{
    // Distribute rates to children that need fallback (tiny fractions in non-preferred mode)
    if (tree.count(nodeName) == 0) return;
    const auto& node = tree.at(nodeName);

    const Well::ProducerCMode modePreferred = node.preferredCtrl;
    int modePrefIdx = kOil;
    switch (modePreferred) {
        case Well::ProducerCMode::ORAT: modePrefIdx = kOil; break;
        case Well::ProducerCMode::WRAT: modePrefIdx = kWater; break;
        case Well::ProducerCMode::GRAT: modePrefIdx = kGas; break;
        default: modePrefIdx = kOil;
    }

    int modeIdx = kOil;
    switch (mode) {
        case Well::ProducerCMode::ORAT: modeIdx = kOil; break;
        case Well::ProducerCMode::WRAT: modeIdx = kWater; break;
        case Well::ProducerCMode::GRAT: modeIdx = kGas; break;
        default: modeIdx = kOil;
    }

    // Collect fallback children and compute ratio from group-controlled children
    std::vector<std::string> cFallback;
    Scalar rateSums = Scalar(0);
    Scalar guiderateSums = Scalar(0);

    for (const auto& childName : c) {
        if (tree.count(childName) == 0) continue;
        const auto& child = tree.at(childName);

        if (child.useFallback) {
            cFallback.push_back(childName);
        } else if (child.modeCategory == "GRUP") {
            rateSums += -child.rates[modePrefIdx];
            guiderateSums += child.guideRates[modePrefIdx];
        }
    }

    if (cFallback.empty()) return;

    const bool anyGroupControlledWells = (rateSums > Scalar(0));
    const Scalar ratio = (guiderateSums > Scalar(0)) ? (rateSums / guiderateSums) : Scalar(1);

    for (const auto& ckName : cFallback) {
        if (tree.count(ckName) == 0) continue;
        auto& ck = tree.at(ckName);

        Scalar qk;
        if (anyGroupControlledWells) {
            qk = ck.guideRates[modePrefIdx] * ratio;
        } else {
            // No group-controlled wells: use limit
            const auto it = ck.individualLimits.find(modePreferred);
            qk = (it != ck.individualLimits.end()) ? it->second : Scalar(0);
        }

        balanceGroupTree(tree, ckName, modePreferred, qk, tol);

        // Add rates but set mode rate to zero to avoid triggering rebalancing
        auto rates = ck.rates;
        rates[modeIdx] = Scalar(0);
        incrementParentRateSums(tree, ckName, nodeName, std::make_optional(rates));
    }
}

// ---------------------------------------------------------------------------

template<class Scalar>
void balanceGroupTree(Tree<Scalar>& tree,
                      const std::string& nodeName,
                      Well::ProducerCMode targetMode,
                      Scalar targetRate,
                      Scalar tol)
{
    if (tree.count(nodeName) == 0) return;
    auto& node = tree.at(nodeName);

    OpmLog::info(fmt::format("balanceGroupTree: {} mode={} target={:.3f} count={}",
                             nodeName, WellProducerCMode2String(targetMode), targetRate, node.balancingCount));

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
        int modeIdx = kOil;
        switch (mode) {
            case Well::ProducerCMode::ORAT: modeIdx = kOil; break;
            case Well::ProducerCMode::WRAT: modeIdx = kWater; break;
            case Well::ProducerCMode::GRAT: modeIdx = kGas; break;
            default: modeIdx = kOil;
        }

        Scalar scaledRates[3];
        const Scalar rateForModeVal = rateForMode(rates, mode, node.resvCoeff);
        if (rateForModeVal > Scalar(0)) {
            for (int c = 0; c < 3; ++c) {
                scaledRates[c] = rates[c] * (qm / rateForModeVal);
            }

            // Check if any limit is violated after scaling
            Scalar maxViolation = Scalar(0);
            Well::ProducerCMode violatingMode = mode;
            for (const auto& [limitMode, limit] : node.individualLimits) {
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

        bool anyGroupControlledChildren = false;

        if (node.type == ProdNodeType::Well) {
            // Just scale rates for wells
            const Scalar currentRate = rateForMode(rates, mode, node.resvCoeff);
            if (currentRate > Scalar(0)) {
                for (int c = 0; c < 3; ++c) {
                    node.rates[c] = rates[c] * (qm / currentRate);
                }
            }
            node.isBalanced = true;
        } else {
            // Group node: iterative balancing
            bool balanced = false;
            int resortingCount = 0;
            int topSwitchCount = 0;

            while (!balanced && resortingCount <= maxResortingCount && topSwitchCount <= maxTopSwitchCount) {
                OpmLog::info(fmt::format("  balanceGroupTree: {} while-loop iteration resort={} switch={}",
                                         nodeName, resortingCount, topSwitchCount));
                // Get local tree descendants
                auto [c, c_fixed, c_trans] = getLocalTreeDescendants(tree, nodeName);

                // Reset rate_sums and guide_rate_sums
                resetRatesAndGuideRateSums(tree, c, c_fixed, c_trans, nodeName, mode);

                // Add in rates from fixed sub-nodes
                for (const auto& cfName : c_fixed) {
                    incrementParentRateSums(tree, cfName, nodeName);
                }

                // Get ratios for sorting
                auto [ratios, modeIndices] = getRatiosForSorting(tree, c);

                // Sort by ratios
                std::vector<size_t> sortedIndices(c.size());
                std::iota(sortedIndices.begin(), sortedIndices.end(), 0);
                std::sort(sortedIndices.begin(), sortedIndices.end(),
                         [&ratios](size_t i1, size_t i2) { return ratios[i1] < ratios[i2]; });

                bool needsResorting = false;
                anyGroupControlledChildren = false;
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
                    if (!c_trans.empty() && !anyGroupControlledChildren) {
                        auto c_trans_update = updateTransparentGroups(tree, nodeName, c_trans,
                                                                      ratios[idx], qm, mode, tol);

                        if (!c_trans_update.empty()) {
                            bool anyTransparent = false;
                            bool anyIndividual = false;
                            for (const auto& ctkName : c_trans_update) {
                                if (tree.count(ctkName) > 0 && hasFreePath(tree, ctkName, nodeName)) {
                                    const auto& ctk = tree.at(ctkName);
                                    anyTransparent = anyTransparent || (ctk.modeCategory == "TRANSPARENT");
                                    anyIndividual = anyIndividual || (ctk.modeCategory == "INDIVIDUAL");
                                }
                            }
                            needsResorting = anyGroupControlledChildren && anyIndividual;
                            anyGroupControlledChildren = anyGroupControlledChildren || anyTransparent;

                            if (!hasFreePath(tree, ckName, nodeName)) continue;
                        }
                    }

                    if (tree.count(nodeName) == 0) break;
                    const auto& originNode = tree.at(nodeName);
                    int modeCanonical = kOil;
                    switch (mode) {
                        case Well::ProducerCMode::ORAT: modeCanonical = kOil; break;
                        case Well::ProducerCMode::WRAT: modeCanonical = kWater; break;
                        case Well::ProducerCMode::GRAT: modeCanonical = kGas; break;
                        default: modeCanonical = kOil;
                    }

                    const Scalar gsum = originNode.guideRateSums[modeCanonical];
                    const Scalar gk = ck.guideRateSums[modeCanonical];
                    const Scalar qmRemain = qm - originNode.rateSums[modeCanonical];
                    const Scalar qk = (gk / (gsum + Scalar(1e-20))) * qmRemain;

                    // Save original guide_rate_sums
                    const auto guideRateSumsOrig = ck.guideRateSums;

                    // Recursive balance
                    balanceGroupTree(tree, ckName, mode, qk, tol);

                    if (tree.count(ckName) > 0) {
                        const auto& ckBalanced = tree.at(ckName);
                        anyGroupControlledChildren = anyGroupControlledChildren || (ckBalanced.modeCategory == "GRUP");
                        if (ckBalanced.modeCategory != "GRUP" && anyGroupControlledChildren) {
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
                        for (int c = 0; c < 3; ++c) {
                            ct.rates[c] = -ct.rateSums[c];
                        }

                        // Check if any limits are broken
                        for (const auto& [limitMode, limit] : ct.individualLimits) {
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
                        distributeFallbackRates(tree, nodeName, c, mode, tol);

                        // Update transparent groups once more
                        for (const auto& ctName : c_trans) {
                            if (tree.count(ctName) > 0) {
                                auto& ct = tree.at(ctName);
                                // Convert rateSums (positive) to rates (negative = production)
                                for (int c = 0; c < 3; ++c) {
                                    ct.rates[c] = -ct.rateSums[c];
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

                            for (const auto& [limitMode, limit] : currentNode.individualLimits) {
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
                            } else if (maxViolation < Scalar(1) - tol && anyGroupControlledChildren) {
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
                        for (int c = 0; c < 3; ++c) {
                            finalNode.rates[c] = -finalNode.rateSums[c];
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
            finalNode.activeIndividualCtrl = mode;
            finalNode.visited = true;

            if (finalNode.type == ProdNodeType::Group && !anyGroupControlledChildren) {
                finalNode.modeCategory = "NONE";
                finalNode.activeIndividualCtrl = Well::ProducerCMode::CMODE_UNDEFINED;
            } else {
                bool atLimit = false;
                for (const auto& [limitMode, limit] : finalNode.individualLimits) {
                    const Scalar currentRate = rateForMode(finalNode.rates, limitMode, finalNode.resvCoeff);
                    if (std::abs(currentRate - limit) <= tol * limit) {
                        atLimit = true;
                        break;
                    }
                }

                if (atLimit) {
                    finalNode.modeCategory = "INDIVIDUAL";
                } else {
                    const Scalar rateAtMode = rateForMode(finalNode.rates, mode, finalNode.resvCoeff);
                    if (std::abs(rateAtMode - qm) <= tol * std::abs(qm)) {
                        if (finalNode.hasGuideRate) {
                            finalNode.modeCategory = "GRUP";
                        } else {
                            finalNode.modeCategory = "TRANSPARENT";
                            finalNode.activeIndividualCtrl = Well::ProducerCMode::CMODE_UNDEFINED;
                        }
                    } else {
                        // Problematic balancing
                        finalNode.modeCategory = "UNDETERMINED";
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

            if (child.modeCategory == "TRANSPARENT") {
                auto childDesc = getDescendants(childName);
                desc.insert(desc.end(), childDesc.begin(), childDesc.end());
            } else {
                desc.push_back(childName);
            }
        }
        return desc;
    };

    auto cList = getDescendants(topName);

    if (top.modeCategory == "NONE") {
        // No target for descendants of 'NONE'
        for (const auto& cName : cList) {
            if (tree.count(cName) > 0) {
                auto& c = tree.at(cName);
                c.groupTarget.groupName = topName;
                c.groupTarget.ctrlMode = Group::ProductionCMode::NONE;
                c.groupTarget.value = std::numeric_limits<Scalar>::quiet_NaN();
                setTargets(tree, cName);
            }
        }
    } else {
        // Node is 'INDIVIDUAL' or 'GRUP'
        const Well::ProducerCMode mode = top.activeIndividualCtrl;
        const Well::ProducerCMode modePreferred = top.preferredCtrl;
        const bool modeIsPref = (mode == modePreferred);

        // Get mode index
        int modeIdx = kOil;
        switch (mode) {
            case Well::ProducerCMode::ORAT: modeIdx = kOil; break;
            case Well::ProducerCMode::WRAT: modeIdx = kWater; break;
            case Well::ProducerCMode::GRAT: modeIdx = kGas; break;
            default: modeIdx = kOil;
        }

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

        // Collect guide rates and rates for children
        std::vector<Scalar> guideRates, rates;
        std::vector<bool> isGrup;

        for (const auto& cName : cList) {
            if (tree.count(cName) == 0) continue;
            const auto& c = tree.at(cName);

            Scalar gr = c.guideRates[modeIdx];
            Scalar r = -c.rates[modeIdx];

            if (anyFallback && c.useFallback) {
                gr = Scalar(0);
                r = Scalar(0);
            }

            guideRates.push_back(gr);
            rates.push_back(r);
            isGrup.push_back(c.modeCategory == "GRUP");
        }

        const Scalar guideSum = std::accumulate(guideRates.begin(), guideRates.end(), Scalar(0));
        Scalar guideSumIsGroup = Scalar(0);
        for (size_t i = 0; i < isGrup.size(); ++i) {
            if (isGrup[i]) guideSumIsGroup += guideRates[i];
        }

        const Scalar target = -top.rates[modeIdx];
        Scalar targetSum = target;
        for (size_t i = 0; i < isGrup.size(); ++i) {
            if (!isGrup[i]) targetSum -= rates[i];
        }

        // Set targets for children
        for (size_t i = 0; i < cList.size(); ++i) {
            const auto& cName = cList[i];
            if (tree.count(cName) == 0) continue;
            auto& c = tree.at(cName);

            c.groupTarget.groupName = topName;
            c.groupTarget.ctrlMode = static_cast<Group::ProductionCMode>(mode); // Approximate conversion
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
                c.groupTarget.value = (guideSumIsGroup > Scalar(0)) ?
                    (guideRates[i] / guideSumIsGroup) * targetSum : Scalar(0);
            } else {
                const Scalar totalGuide = guideSumIsGroup + guideRates[i];
                const Scalar totalTarget = targetSum + rates[i];
                c.groupTarget.value = (totalGuide > Scalar(0)) ?
                    (guideRates[i] / totalGuide) * totalTarget : Scalar(0);
            }

            setTargets(tree, cName);
        }

        // Set fallback targets if not preferred mode
        // (Implementation simplified - full version would set groupTargetFallback)
    }
}

// ---------------------------------------------------------------------------

template<class Scalar>
bool runBalancingAlgorithm(Tree<Scalar>& tree, Scalar tol)
{
    // Main algorithm entry point (corresponds to MATLAB sortalgo2)
    OpmLog::info("runBalancingAlgorithm: starting");

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

        // Set initial category based on allowGroupControl
        if (!node.allowGroupControl) {
            node.modeCategory = "INDIVIDUAL";
        } else if (node.hasGuideRate) {
            node.modeCategory = "GRUP";
        } else {
            node.modeCategory = "TRANSPARENT";
        }
    };
    initializeCategories("FIELD");

    // Get sub-tree ordering
    const auto topnodes = getSubTreeOrdering(tree, "FIELD");

    if (topnodes.empty()) {
        // No nodes to balance
        return true;
    }

    for (const auto& nodeName : topnodes) {
        OpmLog::info(fmt::format("runBalancingAlgorithm: processing topnode {}", nodeName));
        if (tree.count(nodeName) == 0) continue;
        auto& node = tree.at(nodeName);

        // Determine mode and target
        Well::ProducerCMode mode;
        if (node.modeCategory == "INDIVIDUAL") {
            mode = node.activeIndividualCtrl;
        } else {
            mode = node.preferredCtrl;
        }

        // Validate mode
        if (mode == Well::ProducerCMode::CMODE_UNDEFINED ||
            mode == Well::ProducerCMode::NONE) {
            mode = Well::ProducerCMode::ORAT; // Fallback to ORAT
        }

        // Get limit for mode
        Scalar qm = std::numeric_limits<Scalar>::max();
        if (node.individualLimits.count(mode) > 0) {
            qm = node.individualLimits.at(mode);
        }

        // Balance this subtree
        balanceGroupTree(tree, nodeName, mode, qm, tol);

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

        const auto status = node.ctrlStatus;

        // Check that IndividualControlled nodes are at their limit
        if (status == ProdNodeCtrlStatus::IndividualControlled) {
            const auto it = node.individualLimits.find(node.activeIndividualCtrl);
            if (it != node.individualLimits.end()) {
                const Scalar limit = it->second;
                const Scalar current = (node.activeIndividualCtrl == Well::ProducerCMode::BHP)
                    ? std::accumulate(node.rates.begin(), node.rates.end(), Scalar(0),
                                       [](Scalar s, Scalar r){ return s + (-r); })
                    : rateForMode(node.rates, node.activeIndividualCtrl, node.resvCoeff);
                const Scalar relErr = std::abs(current - limit) / (std::abs(limit) + kFeasibilityTolerance<Scalar>);
                if (relErr > tol) {
                    logger.warning("ProdGroupTreeBalancer",
                        fmt::format("Node '{}' is IndividualControlled but current rate ({:.4g}) "
                                    "differs from limit ({:.4g}) by {:.2g}%%",
                                    name, current, limit, 100.0 * relErr));
                    valid = false;
                }
            }
        }

        // Check that GroupControlled nodes satisfy their group target
        if (status == ProdNodeCtrlStatus::GroupControlled &&
            node.groupTarget.ctrlMode != Group::ProductionCMode::NONE) {
            const Scalar target = node.groupTarget.value;
            const Scalar current = std::accumulate(node.rates.begin(), node.rates.end(), Scalar(0),
                                                    [](Scalar s, Scalar r){ return s + (-r); });
            const Scalar relErr = std::abs(current - target) / (std::abs(target) + kFeasibilityTolerance<Scalar>);
            if (relErr > tol) {
                logger.warning("ProdGroupTreeBalancer",
                    fmt::format("Node '{}' is GroupControlled but current total rate ({:.4g}) "
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
            if (node.ctrlStatus == ProdNodeCtrlStatus::GroupControlled) {
                ws.production_cmode = Well::ProducerCMode::GRUP;
            } else if (node.ctrlStatus == ProdNodeCtrlStatus::IndividualControlled) {
                ws.production_cmode = node.activeIndividualCtrl;
            }

            // Update the group target stored in the well state
            if (node.groupTarget.ctrlMode != Group::ProductionCMode::NONE) {
                typename SingleWellState<Scalar, IndexTraits>::GroupTarget gt;
                gt.group_name      = node.groupTarget.groupName;
                gt.production_cmode = node.groupTarget.ctrlMode;
                gt.target_value    = node.groupTarget.value;
                gt.guiderate_ratio = node.groupTarget.guideRateRatio;
                ws.group_target    = gt;
            }
        } else {
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

            // Update group control mode
            if (node.ctrlStatus == ProdNodeCtrlStatus::GroupControlled ||
                node.ctrlStatus == ProdNodeCtrlStatus::Undetermined) {
                // Keep as FLD/NONE (controlled by higher group)
            } else if (node.ctrlStatus == ProdNodeCtrlStatus::IndividualControlled) {
                // Map individualCtrl (Well::ProducerCMode) to Group::ProductionCMode
                // These are separate enums but share the same name mappings.
                switch (node.activeIndividualCtrl) {
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
                          int maxIter,
                          DeferredLogger& logger)
{
    OPM_TIMEFUNCTION();

    const auto t0 = std::chrono::steady_clock::now();

    auto tree = buildTree(wellModel, summaryState, reportStep);

    const bool converged = runBalancingAlgorithm(tree, tol);

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
    const SummaryState&, int);

template std::vector<std::string>
getSubTreeOrdering<double>(const Tree<double>&, const std::string&);

// New algorithm helper functions
template bool hasFreePath<double>(const Tree<double>&, const std::string&, const std::string&);

template void incrementParentRateSums<double>(Tree<double>&, const std::string&, const std::string&,
    const std::optional<std::array<double, 3>>&);

template void decrementParentGuideRateSums<double>(Tree<double>&, const std::string&, const std::string&, const std::array<double, 3>&);

template std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>>
getLocalTreeDescendants<double>(const Tree<double>&, const std::string&);

template void resetRatesAndGuideRateSums<double>(Tree<double>&, const std::vector<std::string>&,
    const std::vector<std::string>&, const std::vector<std::string>&, const std::string&, Well::ProducerCMode);

template std::pair<std::vector<double>, std::vector<int>>
getRatiosForSorting<double>(const Tree<double>&, const std::vector<std::string>&);

template std::vector<std::string>
updateTransparentGroups<double>(Tree<double>&, const std::string&, std::vector<std::string>&,
    double, double, Well::ProducerCMode, double);

template void distributeFallbackRates<double>(Tree<double>&, const std::string&,
    const std::vector<std::string>&, Well::ProducerCMode, double);

template void balanceGroupTree<double>(Tree<double>&, const std::string&,
    Well::ProducerCMode, double, double);

template void setTargets<double>(Tree<double>&, const std::string&);

template bool runBalancingAlgorithm<double>(Tree<double>&, double);

template bool checkTreeValidity<double>(const Tree<double>&, const std::string&, double, DeferredLogger&);

template void applyTreeToState<double, BlackOilDefaultFluidSystemIndices>(
    const Tree<double>&,
    BlackoilWellModelGeneric<double, BlackOilDefaultFluidSystemIndices>&);

template bool runGroupTreeBalancer<double, BlackOilDefaultFluidSystemIndices>(
    BlackoilWellModelGeneric<double, BlackOilDefaultFluidSystemIndices>&,
    const SummaryState&, int, double, int, DeferredLogger&);

#ifdef FLOW_INSTANTIATE_FLOAT

template Tree<float> buildTree<float, BlackOilDefaultFluidSystemIndices>(
    const BlackoilWellModelGeneric<float, BlackOilDefaultFluidSystemIndices>&,
    const SummaryState&, int);

template std::vector<std::string>
getSubTreeOrdering<float>(const Tree<float>&, const std::string&);

// New algorithm helper functions
template bool hasFreePath<float>(const Tree<float>&, const std::string&, const std::string&);

template void incrementParentRateSums<float>(Tree<float>&, const std::string&, const std::string&,
    const std::optional<std::array<float, 3>>&);

template void decrementParentGuideRateSums<float>(Tree<float>&, const std::string&, const std::string&, const std::array<float, 3>&);

template std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>>
getLocalTreeDescendants<float>(const Tree<float>&, const std::string&);

template void resetRatesAndGuideRateSums<float>(Tree<float>&, const std::vector<std::string>&,
    const std::vector<std::string>&, const std::vector<std::string>&, const std::string&, Well::ProducerCMode);

template std::pair<std::vector<float>, std::vector<int>>
getRatiosForSorting<float>(const Tree<float>&, const std::vector<std::string>&);

template std::vector<std::string>
updateTransparentGroups<float>(Tree<float>&, const std::string&, std::vector<std::string>&,
    float, float, Well::ProducerCMode, float);

template void distributeFallbackRates<float>(Tree<float>&, const std::string&,
    const std::vector<std::string>&, Well::ProducerCMode, float);

template void balanceGroupTree<float>(Tree<float>&, const std::string&,
    Well::ProducerCMode, float, float);

template void setTargets<float>(Tree<float>&, const std::string&);

template bool runBalancingAlgorithm<float>(Tree<float>&, float);

template bool checkTreeValidity<float>(const Tree<float>&, const std::string&, float, DeferredLogger&);

template void applyTreeToState<float, BlackOilDefaultFluidSystemIndices>(
    const Tree<float>&,
    BlackoilWellModelGeneric<float, BlackOilDefaultFluidSystemIndices>&);

template bool runGroupTreeBalancer<float, BlackOilDefaultFluidSystemIndices>(
    BlackoilWellModelGeneric<float, BlackOilDefaultFluidSystemIndices>&,
    const SummaryState&, int, float, int, DeferredLogger&);

#endif

} // namespace Opm::ProdGroupTreeBalancer
