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

#include <opm/common/TimingMacros.hpp>

#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateModel.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/GroupState.hpp>
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
                const auto stableBhp = calc.estimateStableBhp(
                    wellState, eclWell, ws.surface_rates,
                    genWell->getRefDensity(), summaryState);

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
            // Fallback when the well interface is not available locally (e.g., on a
            // different MPI rank): use the previously-computed forced_bhp_from_thp.
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
                       int reportStep)
{
    const auto& group = schedule.getGroup(groupName, reportStep);

    node.name   = groupName;
    node.type   = ProdNodeType::Group;
    node.parent = (groupName == "FIELD") ? "" : group.parent();
    node.allowGroupControl = group.productionGroupControlAvailable();

    // Children
    node.children = group.groups();
    for (const auto& w : group.wells()) {
        node.children.push_back(w);
    }

    // Current rates from group state.
    // GroupState stores rates as positive = production in active-phase order.
    // The balancer uses negative = production in canonical [oil, water, gas] order.
    if (groupState.has_production_rates(groupName)) {
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

    if (ctrl == Group::ProductionCMode::FLD || ctrl == Group::ProductionCMode::NONE) {
        node.ctrlStatus = ProdNodeCtrlStatus::GroupControlled;
    } else {
        node.ctrlStatus = ProdNodeCtrlStatus::Undetermined;
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
            populateGroupNode(node, name, schedule, wellState, groupState, pu, reportStep);

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

    // Top-down pass: propagate effective group control modes and compute guide rates.
    // Each group distributes targets to its children in its own production control mode
    // (or its parent's mode if it is FLD-controlled).  The guide rate for each child
    // is looked up for that specific mode so it matches the target being distributed.
    std::function<void(const std::string&, Group::ProductionCMode)>
    propagateGuideRates = [&](const std::string& groupName,
                              Group::ProductionCMode inheritedMode)
    {
        if (tree.count(groupName) == 0) return;
        auto& gnode = tree.at(groupName);

        // Effective mode this group uses to issue targets to its children.
        Group::ProductionCMode effectiveMode = gnode.ownCtrlMode;
        if (effectiveMode == Group::ProductionCMode::FLD ||
            effectiveMode == Group::ProductionCMode::NONE) {
            effectiveMode = inheritedMode; // relay parent's mode
        }
        // Store the effective mode back so setFinalTargets can use it.
        gnode.ownCtrlMode = effectiveMode;

        for (const auto& childName : gnode.children) {
            if (tree.count(childName) == 0) continue;
            auto& child = tree.at(childName);
            child.groupTarget.ctrlMode  = effectiveMode;
            child.groupTarget.guideRate =
                getGuideRateForMode(childName, child.rates, effectiveMode, guideRate);
            if (child.type == ProdNodeType::Group) {
                propagateGuideRates(childName, effectiveMode);
            }
        }
    };
    // ORAT is the ultimate fallback if FIELD has no control mode set.
    propagateGuideRates("FIELD", Group::ProductionCMode::ORAT);

    return tree;
}

// ---------------------------------------------------------------------------

template<class Scalar>
std::vector<std::string> getSubTreeOrdering(const Tree<Scalar>& tree,
                                             const std::string& rootName)
{
    // Post-order traversal: leaves first, root last.
    std::vector<std::string> ordering;
    std::vector<std::string> stack;
    stack.push_back(rootName);

    while (!stack.empty()) {
        const std::string name = stack.back();
        stack.pop_back();
        if (tree.count(name) == 0) continue;
        const auto& node = tree.at(name);
        // Only include nodes that still need balancing (groups with Undetermined status)
        if (node.type == ProdNodeType::Group &&
            node.ctrlStatus == ProdNodeCtrlStatus::Undetermined) {
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
    return ordering;
}

// ---------------------------------------------------------------------------

template<class Scalar>
void parametrizeTree(Tree<Scalar>& tree, const std::string& nodeName)
{
    if (tree.count(nodeName) == 0) return;
    auto& node = tree.at(nodeName);

    if (node.type == ProdNodeType::Well) {
        // For a well: the linearTerm is the phase-fraction direction.
        for (int c = 0; c < 3; ++c) {
            node.linearTerm[c] = node.wellRateFractions[c]; // positive = production direction
        }

        // Since populateWellNode pre-computed the single binding limit,
        // individualLimits has at most one entry.  Compute alphaToNextLimit
        // for that limit (and handle BHP as a total-rate limit).
        node.alphaToNextLimit = std::numeric_limits<Scalar>::max();
        node.nextLimitCtrl    = Well::ProducerCMode::CMODE_UNDEFINED;

        for (const auto& [mode, limit] : node.individualLimits) {
            Scalar currentRate = Scalar(0);
            Scalar lt          = Scalar(0);

            if (mode == Well::ProducerCMode::BHP) {
                // BHP limit stored as a total-rate equivalent
                for (int c = 0; c < 3; ++c) { currentRate += -node.rates[c]; }
                for (int c = 0; c < 3; ++c) { lt          += node.linearTerm[c]; }
            } else {
                currentRate = rateForMode(node.rates, mode, node.resvCoeff);
                lt          = linearTermForMode(node.linearTerm, mode, node.resvCoeff);
            }

            if (lt <= Scalar(0)) continue; // no production increase in this direction

            const Scalar remaining = limit - currentRate;
            if (remaining <= Scalar(0)) {
                node.alphaToNextLimit = Scalar(0);
                node.nextLimitCtrl    = mode;
            } else {
                const Scalar alpha = remaining / lt;
                if (alpha < node.alphaToNextLimit) {
                    node.alphaToNextLimit = alpha;
                    node.nextLimitCtrl    = mode;
                }
            }
        }

        return;
    }

    // Group node: aggregate children's linear terms weighted by their guide rates.
    node.linearTerm = {};

    // Sum up GroupControlled children
    Scalar totalGuideRate = Scalar(0);
    for (const auto& childName : node.children) {
        if (tree.count(childName) == 0) continue;
        auto& child = tree.at(childName);
        if (child.ctrlStatus != ProdNodeCtrlStatus::GroupControlled) continue;

        // First recurse
        parametrizeTree(tree, childName);

        // Weight by the child's pre-computed guide rate (matches the parent's ctrl mode)
        const Scalar gr = child.groupTarget.guideRate;
        totalGuideRate += gr;
        for (int c = 0; c < 3; ++c) {
            node.linearTerm[c] += gr * child.linearTerm[c];
        }
    }

    // Normalize
    if (totalGuideRate > Scalar(0)) {
        for (auto& lt : node.linearTerm) { lt /= totalGuideRate; }
    }

    // Group's own alpha: smallest among all GroupControlled children
    node.alphaToNextLimit = std::numeric_limits<Scalar>::max();
    node.nextLimitCtrl    = Well::ProducerCMode::CMODE_UNDEFINED;
    for (const auto& childName : node.children) {
        if (tree.count(childName) == 0) continue;
        const auto& child = tree.at(childName);
        if (child.ctrlStatus != ProdNodeCtrlStatus::GroupControlled) continue;
        if (child.alphaToNextLimit < node.alphaToNextLimit) {
            node.alphaToNextLimit = child.alphaToNextLimit;
            node.nextLimitCtrl    = child.nextLimitCtrl;
        }
    }
}

// ---------------------------------------------------------------------------

template<class Scalar>
std::pair<std::string, Scalar>
getNextLimitNode(const Tree<Scalar>& tree, const std::string& nodeName)
{
    if (tree.count(nodeName) == 0) {
        return {"", std::numeric_limits<Scalar>::max()};
    }
    const auto& node = tree.at(nodeName);

    if (node.type == ProdNodeType::Well) {
        return {nodeName, node.alphaToNextLimit};
    }

    std::string bestNode  = nodeName;
    Scalar      bestAlpha = node.alphaToNextLimit;

    for (const auto& childName : node.children) {
        if (tree.count(childName) == 0) continue;
        const auto& child = tree.at(childName);
        if (child.ctrlStatus != ProdNodeCtrlStatus::GroupControlled) continue;
        auto [cName, cAlpha] = getNextLimitNode(tree, childName);
        if (cAlpha < bestAlpha) {
            bestAlpha = cAlpha;
            bestNode  = cName;
        }
    }
    return {bestNode, bestAlpha};
}

// ---------------------------------------------------------------------------

template<class Scalar>
void stepAlpha(Tree<Scalar>& tree, const std::string& nodeName, Scalar alpha)
{
    if (tree.count(nodeName) == 0) return;
    auto& node = tree.at(nodeName);

    if (node.ctrlStatus == ProdNodeCtrlStatus::IndividualControlled) {
        return; // This node is fixed; do not change its rates
    }

    // Update rates: rates[c] -= alpha * linearTerm[c]  (production is negative)
    for (int c = 0; c < 3; ++c) {
        node.rates[c] -= alpha * node.linearTerm[c];
    }

    if (node.type == ProdNodeType::Well) {
        return;
    }

    // Recurse into GroupControlled children
    for (const auto& childName : node.children) {
        if (tree.count(childName) == 0) continue;
        auto& child = tree.at(childName);
        if (child.ctrlStatus == ProdNodeCtrlStatus::GroupControlled) {
            stepAlpha(tree, childName, alpha);
        }
    }

    // After updating children, recompute group rates as sum of children
    std::array<Scalar, 3> sumRates{};
    for (const auto& childName : node.children) {
        if (tree.count(childName) > 0) {
            const auto& child = tree.at(childName);
            for (int c = 0; c < 3; ++c) {
                sumRates[c] += child.rates[c];
            }
        }
    }
    node.rates = sumRates;
}

// ---------------------------------------------------------------------------

template<class Scalar>
void updateNode(Tree<Scalar>& tree, const std::string& nodeName)
{
    if (tree.count(nodeName) == 0) return;
    auto& node = tree.at(nodeName);

    // Switch to IndividualControlled
    node.ctrlStatus         = ProdNodeCtrlStatus::IndividualControlled;
    node.activeIndividualCtrl = node.nextLimitCtrl;
    node.alphaToNextLimit   = std::numeric_limits<Scalar>::max();

    // Walk up the tree and check if ancestors now have no GroupControlled children
    std::string parentName = node.parent;
    while (!parentName.empty() && tree.count(parentName) > 0) {
        auto& parent = tree.at(parentName);
        bool anyGroupControlled = false;
        for (const auto& childName : parent.children) {
            if (tree.count(childName) == 0) continue;
            if (tree.at(childName).ctrlStatus == ProdNodeCtrlStatus::GroupControlled) {
                anyGroupControlled = true;
                break;
            }
        }
        if (!anyGroupControlled) {
            if (parent.ctrlStatus != ProdNodeCtrlStatus::IndividualControlled) {
                parent.ctrlStatus = ProdNodeCtrlStatus::NoGroupChildren;
            }
        }
        parentName = parent.parent;
    }
}

// ---------------------------------------------------------------------------

template<class Scalar>
bool capIndividualAndSum(Tree<Scalar>& tree, const std::string& nodeName)
{
    if (tree.count(nodeName) == 0) return false;
    auto& node = tree.at(nodeName);

    bool anyChanged = false;

    if (node.type == ProdNodeType::Well) {
        // Cap each limit
        const std::array<Well::ProducerCMode, 5> modes = {
            Well::ProducerCMode::ORAT,
            Well::ProducerCMode::WRAT,
            Well::ProducerCMode::GRAT,
            Well::ProducerCMode::LRAT,
            Well::ProducerCMode::RESV
        };
        for (const auto& mode : modes) {
            const auto it = node.individualLimits.find(mode);
            if (it == node.individualLimits.end()) continue;
            const Scalar limit = it->second;
            const Scalar current = rateForMode(node.rates, mode, node.resvCoeff);
            if (current > limit * (Scalar(1) + kFeasibilityTolerance<Scalar>)) {
                // Scale rates down uniformly to satisfy the limit
                const Scalar scale = limit / current;
                for (int c = 0; c < 3; ++c) { node.rates[c] *= scale; }
                anyChanged = true;
                break; // One cap at a time is sufficient for feasibility
            }
        }
        // BHP total rate cap
        {
            const auto it = node.individualLimits.find(Well::ProducerCMode::BHP);
            if (it != node.individualLimits.end()) {
                const Scalar limit = it->second;
                Scalar currentTotal = Scalar(0);
                for (int c = 0; c < 3; ++c) { currentTotal += -node.rates[c]; }
                if (currentTotal > limit * (Scalar(1) + kFeasibilityTolerance<Scalar>)) {
                    const Scalar scale = limit / currentTotal;
                    for (int c = 0; c < 3; ++c) { node.rates[c] *= scale; }
                    anyChanged = true;
                }
            }
        }
        return anyChanged;
    }

    // Group: recurse into children, then recompute group rates
    for (const auto& childName : node.children) {
        if (tree.count(childName) > 0) {
            if (capIndividualAndSum(tree, childName)) {
                anyChanged = true;
            }
        }
    }

    // Recompute group rate as sum of children
    std::array<Scalar, 3> sumRates{};
    for (const auto& childName : node.children) {
        if (tree.count(childName) > 0) {
            for (int c = 0; c < 3; ++c) { sumRates[c] += tree.at(childName).rates[c]; }
        }
    }
    node.rates = sumRates;
    return anyChanged;
}

// ---------------------------------------------------------------------------

template<class Scalar>
void setAndUpdateTargets(Tree<Scalar>& tree, const std::string& nodeName)
{
    if (tree.count(nodeName) == 0) return;
    const auto& node = tree.at(nodeName);

    if (node.type == ProdNodeType::Well) return;

    // Total guide rate for GroupControlled children
    Scalar totalGuideRate = Scalar(0);
    for (const auto& childName : node.children) {
        if (tree.count(childName) == 0) continue;
        const auto& child = tree.at(childName);
        if (child.ctrlStatus == ProdNodeCtrlStatus::GroupControlled) {
            totalGuideRate += child.groupTarget.guideRate;
        }
    }
    if (totalGuideRate <= Scalar(0)) return;

    // Distribute the parent's total rate to GroupControlled children
    const Scalar parentTotal = std::accumulate(node.rates.begin(), node.rates.end(), Scalar(0),
                                                [](Scalar s, Scalar r){ return s + (-r); });

    for (const auto& childName : node.children) {
        if (tree.count(childName) == 0) continue;
        auto& child = tree.at(childName);
        if (child.ctrlStatus != ProdNodeCtrlStatus::GroupControlled) continue;

        const Scalar gr = child.groupTarget.guideRate;
        const Scalar childTarget = parentTotal * (gr / totalGuideRate);

        // Set the child's rates according to its fractions and the target
        for (int c = 0; c < 3; ++c) {
            child.rates[c] = -childTarget * child.wellRateFractions[c];
        }

        // Recurse
        setAndUpdateTargets(tree, childName);

        // After recursion, recompute group child rates from its children
        if (child.type == ProdNodeType::Group) {
            std::array<Scalar, 3> sumRates{};
            for (const auto& grandChildName : child.children) {
                if (tree.count(grandChildName) > 0) {
                    for (int c = 0; c < 3; ++c) {
                        sumRates[c] += tree.at(grandChildName).rates[c];
                    }
                }
            }
            child.rates = sumRates;
        }
    }
}

// ---------------------------------------------------------------------------

template<class Scalar>
bool makeFeasible(Tree<Scalar>& tree,
                  const std::string& topName,
                  Scalar tol,
                  int maxIter)
{
    if (tree.count(topName) == 0) return true;

    // If all rates are zero, seed with small initial rates
    auto& topNode = tree.at(topName);
    const Scalar initTotal = std::accumulate(topNode.rates.begin(), topNode.rates.end(), Scalar(0),
                                              [](Scalar s, Scalar r){ return s + (-r); });
    if (initTotal <= Scalar(0)) {
        // Set small rates in GroupControlled wells proportional to fractions
        std::function<void(const std::string&)> seedRates = [&](const std::string& n) {
            if (tree.count(n) == 0) return;
            auto& node = tree.at(n);
            if (node.type == ProdNodeType::Well) {
                for (int c = 0; c < 3; ++c) {
                    node.rates[c] = -kSmallRate<Scalar> * node.wellRateFractions[c];
                }
            } else {
                for (const auto& child : node.children) { seedRates(child); }
                // Recompute group rates as sum of children
                std::array<Scalar, 3> sumRates{};
                for (const auto& child : node.children) {
                    if (tree.count(child) > 0) {
                        for (int c = 0; c < 3; ++c) { sumRates[c] += tree.at(child).rates[c]; }
                    }
                }
                node.rates = sumRates;
            }
        };
        seedRates(topName);
    }

    for (int iter = 0; iter < maxIter; ++iter) {
        const bool changed = capIndividualAndSum(tree, topName);
        setAndUpdateTargets(tree, topName);
        if (!changed) return true;

        // Check convergence: all limits satisfied within tol
        bool converged = true;
        std::function<void(const std::string&)> checkTol = [&](const std::string& n) {
            if (!converged || tree.count(n) == 0) return;
            const auto& node = tree.at(n);
            for (const auto& [mode, limit] : node.individualLimits) {
                const Scalar current = (mode == Well::ProducerCMode::BHP)
                    ? std::accumulate(node.rates.begin(), node.rates.end(), Scalar(0),
                                       [](Scalar s, Scalar r){ return s + (-r); })
                    : rateForMode(node.rates, mode, node.resvCoeff);
                if (current > limit * (Scalar(1) + tol)) {
                    converged = false;
                    return;
                }
            }
            for (const auto& child : node.children) { checkTol(child); }
        };
        checkTol(topName);
        if (converged) return true;
    }
    return false; // Failed to converge
}

// ---------------------------------------------------------------------------

template<class Scalar>
void setFinalTargets(Tree<Scalar>& tree, const std::string& topName)
{
    if (tree.count(topName) == 0) return;
    auto& top = tree.at(topName);

    // Recursive target assignment.
    // The ctrlMode for each node's groupTarget is the parent GROUP's own control mode
    // (stored in ownCtrlMode), which was propagated from the schedule in buildTree.
    // The guide rate (groupTarget.guideRate) was pre-computed in buildTree for that mode.
    std::function<void(const std::string&, const std::string&, Scalar, Scalar)>
    assignTargets = [&](const std::string& parentGroupName,
                        const std::string& nodeName,
                        Scalar parentTarget,
                        Scalar parentGuideRateTotal)
    {
        if (tree.count(nodeName) == 0) return;
        auto& node = tree.at(nodeName);

        // Guide rate was set in buildTree for the parent group's control mode.
        const Scalar gr = node.groupTarget.guideRate;
        const Scalar myTarget = (parentGuideRateTotal > Scalar(0))
            ? parentTarget * (gr / parentGuideRateTotal) : Scalar(0);

        // Ctrl mode: use the parent group's own control mode.
        const Group::ProductionCMode parentCtrlMode =
            (tree.count(parentGroupName) > 0)
            ? tree.at(parentGroupName).ownCtrlMode
            : Group::ProductionCMode::ORAT; // fallback

        node.groupTarget.groupName      = parentGroupName;
        node.groupTarget.ctrlMode       = parentCtrlMode;
        node.groupTarget.value          = myTarget;
        node.groupTarget.guideRate      = gr; // unchanged
        node.groupTarget.guideRateRatio = (parentGuideRateTotal > Scalar(0)) ? (gr / parentGuideRateTotal) : Scalar(0);

        if (node.type == ProdNodeType::Group && node.ctrlStatus != ProdNodeCtrlStatus::IndividualControlled) {
            // Compute total guide rate of GroupControlled children
            Scalar childGuideTotal = Scalar(0);
            for (const auto& childName : node.children) {
                if (tree.count(childName) == 0) continue;
                if (tree.at(childName).ctrlStatus == ProdNodeCtrlStatus::GroupControlled) {
                    childGuideTotal += tree.at(childName).groupTarget.guideRate;
                }
            }
            for (const auto& childName : node.children) {
                if (tree.count(childName) == 0) continue;
                if (tree.at(childName).ctrlStatus == ProdNodeCtrlStatus::GroupControlled) {
                    assignTargets(nodeName, childName, myTarget, childGuideTotal);
                }
            }
        }
    };

    // Compute total guide rate of GroupControlled direct children of the top node
    Scalar topChildGuideTotal = Scalar(0);
    for (const auto& childName : top.children) {
        if (tree.count(childName) == 0) continue;
        if (tree.at(childName).ctrlStatus == ProdNodeCtrlStatus::GroupControlled) {
            topChildGuideTotal += tree.at(childName).groupTarget.guideRate;
        }
    }

    const Scalar topTotal = std::accumulate(top.rates.begin(), top.rates.end(), Scalar(0),
                                             [](Scalar s, Scalar r){ return s + (-r); });
    for (const auto& childName : top.children) {
        if (tree.count(childName) == 0) continue;
        if (tree.at(childName).ctrlStatus == ProdNodeCtrlStatus::GroupControlled) {
            assignTargets(topName, childName, topTotal, topChildGuideTotal);
        }
    }
}

// ---------------------------------------------------------------------------

template<class Scalar>
bool balanceGroupTree(Tree<Scalar>& tree, Scalar tol, int maxIter)
{
    bool allConverged = true;

    const auto ordering = getSubTreeOrdering(tree, "FIELD");
    for (const auto& topName : ordering) {
        if (tree.count(topName) == 0) continue;
        auto& topNode = tree.at(topName);

        // Make the sub-tree rates feasible first
        if (!makeFeasible(tree, topName, tol, maxIter)) {
            allConverged = false;
        }

        int iter = 0;
        // Iterate until no Undetermined nodes remain in this subtree
        while (topNode.ctrlStatus == ProdNodeCtrlStatus::Undetermined && iter < maxIter) {
            parametrizeTree(tree, topName);
            auto [limitNode, alpha] = getNextLimitNode(tree, topName);
            if (limitNode.empty() || alpha >= std::numeric_limits<Scalar>::max()) {
                break; // No limit found – stop
            }
            stepAlpha(tree, topName, alpha);
            updateNode(tree, limitNode);
            ++iter;
        }

        if (topNode.ctrlStatus == ProdNodeCtrlStatus::Undetermined) {
            allConverged = false;
        }

        setFinalTargets(tree, topName);
    }
    return allConverged;
}

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

    const bool converged = balanceGroupTree(tree, tol, maxIter);

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

template void parametrizeTree<double>(Tree<double>&, const std::string&);

template std::pair<std::string, double>
getNextLimitNode<double>(const Tree<double>&, const std::string&);

template void stepAlpha<double>(Tree<double>&, const std::string&, double);

template void updateNode<double>(Tree<double>&, const std::string&);

template bool capIndividualAndSum<double>(Tree<double>&, const std::string&);

template void setAndUpdateTargets<double>(Tree<double>&, const std::string&);

template bool makeFeasible<double>(Tree<double>&, const std::string&, double, int);

template void setFinalTargets<double>(Tree<double>&, const std::string&);

template bool balanceGroupTree<double>(Tree<double>&, double, int);

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

template void parametrizeTree<float>(Tree<float>&, const std::string&);

template std::pair<std::string, float>
getNextLimitNode<float>(const Tree<float>&, const std::string&);

template void stepAlpha<float>(Tree<float>&, const std::string&, float);

template void updateNode<float>(Tree<float>&, const std::string&);

template bool capIndividualAndSum<float>(Tree<float>&, const std::string&);

template void setAndUpdateTargets<float>(Tree<float>&, const std::string&);

template bool makeFeasible<float>(Tree<float>&, const std::string&, float, int);

template void setFinalTargets<float>(Tree<float>&, const std::string&);

template bool balanceGroupTree<float>(Tree<float>&, float, int);

template bool checkTreeValidity<float>(const Tree<float>&, const std::string&, float, DeferredLogger&);

template void applyTreeToState<float, BlackOilDefaultFluidSystemIndices>(
    const Tree<float>&,
    BlackoilWellModelGeneric<float, BlackOilDefaultFluidSystemIndices>&);

template bool runGroupTreeBalancer<float, BlackOilDefaultFluidSystemIndices>(
    BlackoilWellModelGeneric<float, BlackOilDefaultFluidSystemIndices>&,
    const SummaryState&, int, float, int, DeferredLogger&);

#endif

} // namespace Opm::ProdGroupTreeBalancer
