/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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
#include <opm/simulators/wells/BlackoilWellModelNetworkGeneric.hpp>

#include <opm/common/TimingMacros.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Network/Balance.hpp>

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/GroupStateHelper.hpp>
#include <opm/simulators/wells/BlackoilWellModelNetworkPressureComputation.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>

#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <cassert>
#include <map>
#include <string>
#include <utility>

namespace Opm {

namespace details {
    /// Helper to check if any network (production, gas injection, water injection) is active at a given time step.
    bool anyNetworkActive(const Schedule& schedule, const int timeStepIdx)
    {
        const auto& sstate = schedule[timeStepIdx];
        return sstate.network().active()
            || (sstate.injectionNetwork.get_ptr(Phase::GAS) != nullptr
                && sstate.injectionNetwork.get_ptr(Phase::GAS)->active())
            || (sstate.injectionNetwork.get_ptr(Phase::WATER) != nullptr
                && sstate.injectionNetwork.get_ptr(Phase::WATER)->active());
    }

    /// Helper to get all active networks (production, gas injection, water injection) at a given time step.
    std::vector<NetworkDescriptor>
    activeNetworks(const Schedule& schedule, const int timeStepIdx)
    {
        std::vector<NetworkDescriptor> active_networks;
        const auto& sstate = schedule[timeStepIdx];
        if (sstate.network().active()) {
            active_networks.push_back({NetworkDomain::Production, std::cref(sstate.network())});
        }
        if (sstate.injectionNetwork.get_ptr(Phase::GAS) != nullptr
            && sstate.injectionNetwork.get_ptr(Phase::GAS)->active()) {
            active_networks.push_back({NetworkDomain::InjectionGas, std::cref(*sstate.injectionNetwork.get_ptr(Phase::GAS))});
        }
        if (sstate.injectionNetwork.get_ptr(Phase::WATER) != nullptr
            && sstate.injectionNetwork.get_ptr(Phase::WATER)->active()) {
            active_networks.push_back({NetworkDomain::InjectionWater, std::cref(*sstate.injectionNetwork.get_ptr(Phase::WATER))});
        }
        return active_networks;
    }

    /// Helper to get all networks (production, gas injection, water injection) at a given time step,
    /// whether active or not.
    std::vector<NetworkDescriptor>
    networks(const Schedule& schedule, const int timeStepIdx)
    {
        std::vector<NetworkDescriptor> nw;
        const auto& sstate = schedule[timeStepIdx];
        nw.push_back({NetworkDomain::Production, std::cref(sstate.network())});
        if (sstate.injectionNetwork.get_ptr(Phase::GAS) != nullptr) {
            nw.push_back({NetworkDomain::InjectionGas, std::cref(*sstate.injectionNetwork.get_ptr(Phase::GAS))});
        }
        if (sstate.injectionNetwork.get_ptr(Phase::WATER) != nullptr) {
            nw.push_back({NetworkDomain::InjectionWater, std::cref(*sstate.injectionNetwork.get_ptr(Phase::WATER))});
        }
        return nw;
    }

} // namespace details

// The group-tree network solve: full steps, and from iteration 20 a line
// search, to break the cycles full steps can fall into.
template<class Scalar>
constexpr NetworkSolve::Parameters<Scalar> kGroupTreeSolveParams{1e-2, 50, 20};

// Above this max residual of its best iterate [bar or relative], an
// unconverged group-tree solve is reported as far from a network solution.
template<class Scalar>
constexpr Scalar kGroupTreeFarResidual = 5.0;


template<typename Scalar, typename IndexTraits>
BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
BlackoilWellModelNetworkGeneric(BlackoilWellModelGeneric<Scalar,IndexTraits>& well_model)
    : well_model_(well_model)
{
    this->setFromRestart(well_model_.eclState().getRestartNetworkPressures());
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
setFromRestart(const std::optional<std::map<std::string, double>>& node_pressures)
{
    if (node_pressures.has_value()) {
        if constexpr (std::is_same_v<Scalar,double>) {
            this->node_pressures_ = node_pressures.value();
        } else {
            for (const auto& it : node_pressures.value()) {
                this->node_pressures_[it.first] = it.second;
            }
        }
        this->syncProductionDomainState_();
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
updateActiveState(const int report_step)
{
    this->active_ = false;
    for (const auto& network : details::activeNetworks(well_model_.schedule(), report_step)) {
        updateActiveStateImpl(network.network.get());
    }
    this->active_ = well_model_.comm().max(active_);
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
updateActiveStateImpl(const Network::ExtNetwork& network)
{
    // Accumulates into active_ across the domains; an inactive network must not
    // clear what an earlier domain set.
    if (!network.active()) {
        return;
    }
    bool network_active = false;
    for (const auto& well : well_model_.genericWells()) {
        const bool is_partof_network = network.has_node(well->wellEcl().groupName());
        const bool prediction_mode = well->wellEcl().predictionMode();
        if (is_partof_network && prediction_mode) {
            network_active = true;
            break;
        }
    }
#ifdef RESERVOIR_COUPLING_ENABLED
    // A reservoir coupling master may have no local wells in the network
    // leaf groups, the leaf rates come from slave-reported
    // network_surface_rates instead.  Without this clause the master's
    // network solver short-circuits via active_=false and the leaf
    // pressures are never iterated.
    const auto& rescoup_proxy = well_model_.groupStateHelper().rescoup();
    if (!network_active && rescoup_proxy.isMaster()) {
        const auto& rescoup_master = rescoup_proxy.master();
        const auto num_slaves = rescoup_master.numSlaves();
        for (std::size_t s = 0; s < num_slaves && !network_active; ++s) {
            // Only a coupled slave supplies leaf rates this step; a master
            // group whose slave has not activated yet, or has already ended,
            // must not keep the network active (it would be solved against
            // missing/stale slave rates).
            if (!rescoup_master.slaveIsCoupled(s)) {
                continue;
            }
            for (const auto& master_group : rescoup_master.getMasterGroupNamesForSlave(s)) {
                if (network.has_node(master_group)) {
                    network_active = true;
                    break;
                }
            }
        }
    }
#endif
    this->active_ = this->active_ || network_active;
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
needPreStepRebalance(const int report_step) const
{
    const auto active_networks = details::activeNetworks(well_model_.schedule(), report_step);
    bool network_rebalance_necessary = false;
    for (const auto& well : well_model_.genericWells()) {
        const bool is_partof_network = std::any_of(active_networks.begin(),
                                                   active_networks.end(),
                                                   [&](const auto& network)
                                                   {
                                                       return network.network.get().has_node(well->wellEcl().groupName());
                                                   });
        // TODO: we might find more relevant events to be included here (including network change events?)
        const auto& events = well_model_.wellState().well(well->indexOfWell()).events;
        if (is_partof_network && events.hasEvent(ScheduleEvents::WELL_STATUS_CHANGE)) {
            network_rebalance_necessary = true;
            break;
        }
    }
    network_rebalance_necessary = well_model_.comm().max(network_rebalance_necessary);
    return network_rebalance_necessary;
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
shouldBalance(const int reportStepIdx) const
{
    // if network is not active, we do not need to balance the network
    if (!details::anyNetworkActive(well_model_.schedule(), reportStepIdx)) {
        return false;
    }

    const auto& balance = well_model_.schedule()[reportStepIdx].network_balance();
    const auto& iterCtx = well_model_.iterationContext();
    if (balance.mode() == Network::Balance::CalcMode::TimeStepStart) {
        return iterCtx.isFirstGlobalIteration();
    } else if (balance.mode() == Network::Balance::CalcMode::NUPCOL) {
        const int nupcol = well_model_.schedule()[reportStepIdx].nupcol();
        return iterCtx.withinNupcol(nupcol);
    } else {
        // We do not support any other rebalancing modes,
        // i.e. TimeInterval based rebalancing is not available.
        // This should be warned about elsewhere, so we choose to
        // avoid spamming with a warning here.
        return false;
    }
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
willBalanceOnNextIteration(const int reportStepIdx) const
{
    // if network is not active, we do not need to balance the network
    if (!details::anyNetworkActive(well_model_.schedule(), reportStepIdx)) {
        return false;
    }

    const auto& schedule_state = well_model_.schedule()[reportStepIdx];
    if (schedule_state.network_balance().mode() == Network::Balance::CalcMode::NUPCOL) {
        const int nupcol = schedule_state.nupcol();
        return well_model_.iterationContext().withinNupcol(nupcol - 1); // Note the -1 here!
    } else {
        // Any other rebalancing mode will only rebalance
        // at the start of the timestep.
        return false;
    }
}


template<typename Scalar, typename IndexTraits>
std::unordered_map<std::string, typename NetworkSolve::GroupTreeSystem<Scalar>::WellNetworkData>
BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
gatherWellNetworkDataForGroupTree(const int reportStepIdx) const
{
    using Sys = NetworkSolve::GroupTreeSystem<Scalar>;
    const auto& schedule = well_model_.schedule();
    const auto& summary_state = well_model_.summaryState();

    // Same reasoning as newtonProductionNodePressures(): every rank must end up
    // with the same answer, so the static data comes from the replicated
    // schedule and only what the well is currently doing is summed, contributed
    // once by the rank that owns it.
    std::map<std::string, const WellInterfaceGeneric<Scalar, IndexTraits>*> local;
    for (const auto& well : well_model_.genericWells()) {
        local.emplace(well->name(), well);
    }

    struct Candidate
    {
        std::string name;
        int vfp_table;
        Scalar efficiency;
    };
    std::vector<Candidate> candidates;
    for (const auto& name : schedule.wellNames(reportStepIdx)) {
        const auto& well = schedule.getWell(name, reportStepIdx);
        if (!well.isProducer() || !well.predictionMode()) {
            continue;
        }
        const auto controls = well.productionControls(summary_state);
        candidates.push_back({name, controls.vfp_table_number,
                              static_cast<Scalar>(well.getEfficiencyFactor(/*network=*/true))});
    }

    // Water, oil, gas positions in the well state's active-phase arrays --
    // this class's own [oil, water, gas] convention (Sys::kOil == 0), unlike
    // NetworkProductionSystem's [water, oil, gas] one.
    const auto& pu = well_model_.phaseUsage();
    std::array<int, Sys::NP> pos{};
    pos[Sys::kOil]   = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
    pos[Sys::kWater] = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
    pos[Sys::kGas]   = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
    if (std::any_of(pos.begin(), pos.end(), [](const int p) { return p < 0; })) {
        return {};   // needs all three phases active, same as the Newton builder
    }

    // Per candidate: present, usable ipr, three ipr_a, three ipr_b, efficiency
    // scaling, has a network-sourced (dynamic) thp limit, stopped, usable trial
    // ipr, three trial ipr_a, three trial ipr_b, tubing-table datum correction,
    // bhp the trial ipr was taken at, phase shares (set, three shares), three
    // reference rates, explicit WFR and GFR, current bhp, explicit lookup
    // (0 none, 1 below the first FLO value, 2 everywhere).
    constexpr int kEntries = 31;
    std::vector<Scalar> shared(candidates.size() * kEntries, Scalar{0});
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        const auto it = local.find(candidates[i].name);
        if (it == local.end() || !it->second->parallelWellInfo().isOwner()) {
            continue;
        }
        const auto& ws = well_model_.wellState()[it->second->indexOfWell()];
        // STOP is not SHUT: a stopped well is still part of the system (still
        // solved, still has a row) but pinned at zero rate, whereas a shut
        // well has no representation at all -- see the has_ipr handling below
        // and populateFromFlatNetwork()'s use of it. Only a genuinely shut
        // (or otherwise absent/non-owned) well is dropped from this gather
        // entirely.
        if (ws.status != WellStatus::OPEN && ws.status != WellStatus::STOP) {
            continue;
        }
        Scalar* e = &shared[i * kEntries];
        e[0] = Scalar{1};
        // A usable IPR: some phase's rate responds to bhp -- not necessarily
        // oil (a gas well has none).
        const auto responds = [&pos, &pu](const std::vector<Scalar>& b) {
            return static_cast<int>(b.size()) >= pu.numActivePhases()
                && std::any_of(pos.begin(), pos.end(), [&b](const int p) { return b[p] > Scalar{0}; });
        };
        if (responds(ws.implicit_ipr_b)) {
            e[1] = Scalar{1};
            for (int ph = 0; ph < Sys::NP; ++ph) {
                // ws holds q = b*bhp - a in opm's signed rates (production
                // negative); GroupTreeSystem wants production positive and
                // falling with bhp -- the same line negated.
                e[2 + ph] = ws.implicit_ipr_a[pos[ph]];
                e[5 + ph] = -ws.implicit_ipr_b[pos[ph]];
            }
        }
        e[8] = ws.efficiency_scaling_factor;
        // Where a flowing THP well starts in the network solve: its own
        // current solution, consistent with the IPR the solve uses.
        if (!it->second->wellIsStopped() && !it->second->awaitsFlow(ws)) {
            e[29] = ws.bhp;
        }
        // Individual + THP is the only case this ever matters for (see
        // extractFlatNetworkInput()'s own doc comment on networkThpWells);
        // harmless to compute unconditionally for every well.
        e[9] = it->second->getDynamicThpLimit().has_value() ? Scalar{1} : Scalar{0};
        // Only the well model's own dynamic stop (wellIsStopped() on an OPEN
        // well) makes a reopen candidate: a persistent STOP (deck, economic
        // limits) is not the network's to undo. So does an open well without
        // any rate (awaitsFlow()): the balancer leaves it out.
        // The tubing table's datum is not the well's reference depth: the
        // well's own solve applies this correction (see refreshWellNetworkData()).
        if (const auto dp = well_vfp_dp_.find(candidates[i].name); dp != well_vfp_dp_.end()) {
            e[18] = dp->second;
        }
        if (ws.status == WellStatus::OPEN && (it->second->wellIsStopped() || it->second->awaitsFlow(ws))) {
            e[10] = Scalar{1};
            if (responds(ws.stopped_ipr_b)) {
                e[11] = Scalar{1};
                for (int ph = 0; ph < Sys::NP; ++ph) {
                    e[12 + ph] = ws.stopped_ipr_a[pos[ph]];
                    e[15 + ph] = -ws.stopped_ipr_b[pos[ph]];
                }
                e[19] = ws.stopped_ipr_bhp;
            }
        }

        // Its phase shares (rates per unit of FLO) for the network solve
        // (GroupTreeSystem::WellNetworkData::phase_share): those of its IPR's
        // rates at its current bhp -- the linearisation the network uses, with
        // every phase as it flows there (not the well state's own rates, which
        // can be a single-phase guess set from a rate target before the well
        // was ever solved); for a stopped well with a trial IPR, the
        // composition its anchor judged its lift with (else the rates that IPR
        // gives at the bhp it was taken at); without any flow there, the
        // previous step's rates. Explicit lookups (below the table's first FLO
        // value, or everywhere for WVFPEXP or a well that switched to them)
        // change only the fractions its tubing is looked up with, not these
        // shares. Its reference rates
        // stand for it in a node's fallback fractions: those rates, else its
        // potentials. Its explicit fractions apply below the table's first FLO
        // value if it gets no shares.
        const int table_id = candidates[i].vfp_table;
        if (table_id > 0) {
            const auto& vfp = *it->second->vfpProperties();
            const auto& table = vfp.getProd()->getTable(table_id);
            const auto flo = [&table](const std::array<Scalar, Sys::NP>& r) {
                return detail::getFlo(table, r[Sys::kWater], r[Sys::kOil], r[Sys::kGas]);
            };
            std::array<Scalar, Sys::NP> q{};   // positive, [oil, water, gas]
            if (e[1] > Scalar{0}) {
                for (int ph = 0; ph < Sys::NP; ++ph) {
                    q[ph] = e[2 + ph] + e[5 + ph] * ws.bhp;   // ipr: a + b*bhp, b negative here
                }
            }
            if (e[11] > Scalar{0} && flo(q) < table.getFloAxis().front()) {
                // A stopped well with a trial IPR: the composition its anchor
                // judged its lift with, else the trial IPR's rates at its bhp.
                std::array<Scalar, Sys::NP> c{};
                if (static_cast<int>(ws.stopped_composition.size()) >= pu.numActivePhases()) {
                    for (int ph = 0; ph < Sys::NP; ++ph) {
                        c[ph] = -ws.stopped_composition[pos[ph]];
                    }
                }
                for (int ph = 0; ph < Sys::NP; ++ph) {
                    q[ph] = flo(c) > Scalar{0} ? c[ph]
                                               : e[12 + ph] + e[15 + ph] * e[19];   // trial ipr: a + b*bhp, b negative here
                }
            }
            // Explicit lookups (the deck's WVFPEXP, or the well's own switch to
            // them when it struggled to converge) at every FLO, else below the
            // table's first FLO value: they concern the tubing lookup only
            // (GroupTreeSystem::WellNetworkData::explicit_lookup), not the
            // well's phase rates, which stay those of its IPR.
            if (it->second->useVfpExplicit()) {
                e[30] = Scalar{2};
            } else if (flo(q) < table.getFloAxis().front()) {
                e[30] = Scalar{1};
            }
            // No flowing composition from the IPR: the previous step's rates.
            if (!(flo(q) > Scalar{0})) {
                for (int ph = 0; ph < Sys::NP; ++ph) {
                    q[ph] = -ws.prev_surface_rates[pos[ph]];
                }
            }
            if (flo(q) > Scalar{0}) {
                e[20] = Scalar{1};
                for (int ph = 0; ph < Sys::NP; ++ph) {
                    e[21 + ph] = std::max(q[ph], Scalar{0}) / flo(q);
                }
            } else {
                for (int ph = 0; ph < Sys::NP; ++ph) {
                    q[ph] = std::max(ws.well_potentials[pos[ph]], Scalar{0});
                }
            }
            for (int ph = 0; ph < Sys::NP; ++ph) {
                e[24 + ph] = std::max(q[ph], Scalar{0});
            }
            e[27] = vfp.getExplicitWFR(table_id, it->second->indexOfWell());
            e[28] = vfp.getExplicitGFR(table_id, it->second->indexOfWell());
        }
    }
    well_model_.comm().sum(shared.data(), shared.size());

    std::unordered_map<std::string, typename Sys::WellNetworkData> result;
    result.reserve(candidates.size());
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        const Scalar* e = &shared[i * kEntries];
        if (e[0] <= Scalar{0}) {
            continue;   // shut (or otherwise absent/not owned) on every rank --
                        // genuinely not part of the system, unlike a stopped well.
        }
        typename Sys::WellNetworkData data;
        // Open but no usable (linearised) inflow performance -- typically a
        // stopped well (updateIPRImplicit() itself warns its IPR is
        // "problematic" at zero rate) rather than a data-staleness gap, now
        // that the IPR refresh above runs for every predicting well. Kept in
        // the map anyway rather than dropped: populateFromFlatNetwork() reads
        // this flag to pin such a well at zero instead of tying it to a group
        // target or a THP row it has no real linearisation to solve.
        data.has_ipr = e[1] > Scalar{0};
        data.efficiency = candidates[i].efficiency * e[8];
        data.vfp_table = candidates[i].vfp_table;
        data.vfp_dp = e[18];
        // The alq a Thp well's tubing table sees is decided elsewhere (WLIFTOPT
        // gas lift, or a fixed deck value); left at its default (0) here, same
        // as every "unproblematic wells" test so far -- revisit once gas lift
        // needs to reach GroupTreeSystem's own solve.
        for (int ph = 0; ph < Sys::NP; ++ph) {
            data.ipr_a[ph] = e[2 + ph];
            data.ipr_b[ph] = e[5 + ph];
        }
        data.network_thp = e[9] > Scalar{0};
        data.stopped = e[10] > Scalar{0};
        data.has_trial_ipr = e[11] > Scalar{0};
        for (int ph = 0; ph < Sys::NP; ++ph) {
            data.trial_ipr_a[ph] = e[12 + ph];
            data.trial_ipr_b[ph] = e[15 + ph];
        }
        data.trial_bhp = e[19];
        if (e[20] > Scalar{0}) {
            data.phase_share = std::array<Scalar, Sys::NP>{e[21], e[22], e[23]};
        }
        for (int ph = 0; ph < Sys::NP; ++ph) {
            data.reference_q[ph] = e[24 + ph];
        }
        data.explicit_fractions = {e[27], e[28]};
        data.explicit_lookup = static_cast<int>(std::lround(e[30]));
        data.bhp = e[29];
        result.emplace(candidates[i].name, std::move(data));
    }
    return result;
}

template<typename Scalar, typename IndexTraits>
std::optional<typename BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::GroupTreeSolve>
BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
solveGroupTree(const Network::ExtNetwork& network,
              const int reportStepIdx,
              const Network::Node& root,
              const ProdGroupTreeBalancer::Tree<Scalar>& balancedTree,
              const bool offerReopenCandidates) const
{
    OPM_TIMEFUNCTION();
    using Sys = NetworkSolve::GroupTreeSystem<Scalar>;

    // Through the well model's deferred logger when there is one, so these
    // lines land in order with the balancer's own output (which is deferred).
    // Every rank reaches the same outcome; report it once.
    const auto log = [this](const std::string& msg) {
        if (well_model_.comm().rank() != 0) {
            return;
        }
        if (well_model_.groupStateHelper().hasDeferredLogger()) {
            well_model_.groupStateHelper().deferredLogger().debug(msg);
        } else {
            OpmLog::debug(msg);
        }
    };
    auto giveUp = [&](const std::string& why) {
        log(fmt::format("Network: solving the production network via the group-tree "
                        "balancer is not possible at report step {} ({}); using the "
                        "relaxed update.", reportStepIdx, why));
        // TEMPORARY: the relaxed fixed-point update is to go; make every use
        // of it visible.
        if (well_model_.comm().rank() == 0) {
            const auto msg = fmt::format("Network: the group-tree solve under {} gave up at report step {} "
                                         "({}); the relaxed fixed-point update is used.",
                                         root.name(), reportStepIdx, why);
            if (well_model_.groupStateHelper().hasDeferredLogger()) {
                well_model_.groupStateHelper().deferredLogger().warning("NETWORK_FIXED_POINT_FALLBACK", msg);
            } else {
                OpmLog::warning(msg);
            }
        }
        return std::optional<GroupTreeSolve>{};
    };

    if (!root.terminal_pressure().has_value()) {
        return giveUp(fmt::format("the tree under {} has no terminal pressure", root.name()));
    }
    // An empty tree: the balancer had no producing well (none open yet at
    // the start, or all closed). The network then has no balanced wells and
    // is solved for its nodes alone, plus any reopen candidates -- not
    // handed to the relaxed update.
    if (!balancedTree.empty() && balancedTree.count(root.name()) == 0) {
        return giveUp(fmt::format("{} is not a node in the balanced group tree", root.name()));
    }

    const auto& schedule = well_model_.schedule();
    Sys system(*well_model_.getVFPProperties().getProd());
    system.setTerminalPressure(*root.terminal_pressure());
    system.setAnalyticJacobian(analytic_jacobian_);

    // Nodes, parents before children -- same walk as
    // newtonProductionNodePressures(), minus autochoke (GroupTreeSystem has
    // no choke-target mechanism) and satellite-production network nodes
    // (GroupTreeSystem has no setNodeSource() equivalent; a satellite that
    // is not itself a network node is handled separately, as a group-tree
    // leaf -- see this function's own doc comment).
    std::map<std::string, int> index;
    std::vector<std::string> order{root.name()};
    // Node 0 is already the terminal -- GroupTreeSystem's own nodes_ starts
    // with a default-constructed placeholder there (see NetworkGroupTreeSystem.hpp's
    // nodes_ member), and setTerminalPressure() above already supplies its
    // pressure. Calling addNode() for it too, the way newtonProductionNodePressures()
    // does for NetworkProductionSystem (whose nodes_ starts empty), would push a
    // second entry ahead of every real node, shifting index/order's bookkeeping
    // one slot away from where each node actually lives in nodes_.
    index[order.front()] = 0;
    for (std::size_t at = 0; at < order.size(); ++at) {
        for (const auto& branch : network.downtree_branches(order[at])) {
            const auto& child = branch.downtree_node();
            if (index.count(child)) {
                continue;
            }
            index[child] = static_cast<int>(order.size());
            order.push_back(child);
            const auto& child_node = network.node(child);
            if (child_node.as_choke()) {
                return giveUp(fmt::format("{} is an autochoke node (not supported by the "
                                          "group-tree solver)", child));
            }
            // A network node need not be a group (BRANPROP/NODEPROP junctions);
            // only a group can carry satellite rates.
            if (schedule.hasGroup(child, reportStepIdx)
                && schedule.getGroup(child, reportStepIdx).hasSatelliteProduction()) {
                return giveUp(fmt::format("{} is a satellite-production network node (not "
                                          "supported by the group-tree solver)", child));
            }
            const auto& units = schedule.getUnits();
            Scalar alq = 0.0;
            if (branch.vfp_table().has_value()) {
                const auto& table = well_model_.getVFPProperties().getProd()->getTable(*branch.vfp_table());
                alq = branch.alq_value(VFPProdTable::ALQDimension(table.getALQType(), units)).value_or(0.0);
            }
            NetworkSolve::Node node{child, static_cast<int>(at), branch.vfp_table().value_or(NetworkSolve::NoTable)};
            node.efficiency = child_node.efficiency();
            const int node_idx = system.addNode(std::move(node), alq);
            // Another root of the network, below this one: its pressure is
            // fixed, its flow still counts upwards (see independentRoots()).
            if (child_node.terminal_pressure().has_value()) {
                system.setFixedPressure(node_idx, *child_node.terminal_pressure());
            }
        }
    }

    // The balancer's own tree does not need to coincide with the network's:
    // flatten it starting from this network root's own name, so a well
    // belonging to a different root's subtree cannot leak into this one's
    // well-data lookup below.
    std::unordered_set<std::string> networkThpWells;
    const auto wellNetworkData = gatherWellNetworkDataForGroupTree(reportStepIdx);
    for (const auto& [name, data] : wellNetworkData) {
        if (data.network_thp) {
            networkThpWells.insert(name);
        }
    }
    // The network's wells are those of its nodes that are groups (its leaf
    // nodes, in practice); its root need not be a group at all, nor its nodes
    // follow the group tree. So the whole balanced tree is flattened, and then
    // restricted to these wells: an Active group's equations are kept only if
    // all its wells are in this network.
    std::unordered_set<std::string> networkWells;
    for (const auto& name : order) {
        if (schedule.hasGroup(name, reportStepIdx)) {
            for (const auto& w : schedule.getGroup(name, reportStepIdx).wells()) {
                networkWells.insert(w);
            }
        }
    }
    const auto flat = ProdGroupTreeBalancer::restrictFlatNetworkInput(
        ProdGroupTreeBalancer::extractFlatNetworkInput(balancedTree, std::string("FIELD"),
                                                       well_model_.guideRate(), networkThpWells),
        balancedTree, networkWells);
    // No well of this network in the balanced tree (none producing: all
    // stopped, or none open yet): solved for its nodes and reopen
    // candidates alone, like an empty tree above.

    // Resolve each well's own network node (Well::groupName(), the same
    // convention newtonProductionNodePressures() uses) into this system's
    // own node index, so populateFromFlatNetwork() gets WellNetworkData with
    // a real node field rather than gatherWellNetworkDataForGroupTree()'s own
    // placeholder 0. Every well flat references (ownWells, plus every
    // top-level type == Well entry that is not a satellite) must resolve --
    // populateFromFlatNetwork() looks each one up unconditionally, and a
    // silently-excluded well would either throw there or, worse, just have
    // its rate go missing from whatever ancestor sum needed it.
    //
    // A genuinely absent well (shut, or simply not a predicting producer --
    // gatherWellNetworkDataForGroupTree() never puts one in wellNetworkData
    // at all, unlike a merely stopped one, which it keeps with has_ipr ==
    // false) does not enter the equations at all, the balancer's own
    // ownWells/GRUPTREE-derived membership has no awareness of well status to
    // have excluded it already. Rather than aborting the whole tree over it,
    // it gets the same zero-rate placeholder ensureWell() already builds for
    // a stopped well (WellNetworkData::has_ipr == false, forced here rather
    // than left at the struct's own true default -- see there for why a
    // stale positive target with no real ipr data must not be trusted). Node
    // 0 (the default) also keeps it out of every node's own physical flow
    // sum, the same as a satellite well's own "invisible" convention. Only a
    // well that IS present but whose own group can't be placed in this
    // system's topology is a real configuration gap worth giving up over.
    std::unordered_map<std::string, typename Sys::WellNetworkData> resolved;
    auto resolveWell = [&](const std::string& name) -> bool {
        if (resolved.count(name)) {
            return true;
        }
        const auto it = wellNetworkData.find(name);
        if (it == wellNetworkData.end() || !schedule.hasWell(name, reportStepIdx)) {
            typename Sys::WellNetworkData placeholder;
            placeholder.has_ipr = false;
            resolved.emplace(name, std::move(placeholder));
            return true;
        }
        const auto& node_name = schedule.getWell(name, reportStepIdx).groupName();
        const auto node_it = index.find(node_name);
        if (node_it == index.end()) {
            return false;   // a real topology gap -- give up, same as today
        }
        auto data = it->second;
        data.node = node_it->second;
        resolved.emplace(name, std::move(data));
        return true;
    };
    for (const auto& entry : flat) {
        for (const auto& w : entry.ownWells) {
            if (!resolveWell(w.name)) {
                return giveUp(fmt::format("{} has no usable network-side data", w.name));
            }
        }
        if (entry.type == ProdNodeType::Well && !entry.satelliteRates.has_value()
            && !resolveWell(entry.name)) {
            return giveUp(fmt::format("{} has no usable network-side data", entry.name));
        }
    }

    system.populateFromFlatNetwork(flat, resolved);

    // Reopen candidates: stopped wells under this root with a usable trial IPR.
    // A stopped well is never in the balanced tree, so it is added here, at its
    // own network node, and counted in the sum of its nearest Active ancestor
    // (with the efficiency accumulated up to it, as the balancer would).
    int num_candidates = 0;
    if (offerReopenCandidates) {
        for (const auto& [name, data] : wellNetworkData) {
            if (!data.stopped || !data.has_trial_ipr || resolved.count(name) > 0
                || !schedule.hasWell(name, reportStepIdx)) {
                continue;
            }
            const auto& well = schedule.getWell(name, reportStepIdx);
            const auto node_it = index.find(well.groupName());
            if (node_it == index.end()) {
                continue;   // not under this root
            }
            auto wd = data;
            wd.node = node_it->second;
            Scalar efficiency = well.getEfficiencyFactor()
                * well_model_.wellState().getGlobalEfficiencyScalingFactor(name);
            std::optional<int> active_node;
            std::string group = well.groupName();
            while (true) {
                active_node = system.activeNodeIndex(group);
                if (active_node.has_value() || group == "FIELD") {
                    break;
                }
                const auto& g = schedule.getGroup(group, reportStepIdx);
                efficiency *= g.getGroupEfficiencyFactor();
                group = g.parent();
            }
            system.addReopenCandidate(name, wd, active_node, efficiency);
            ++num_candidates;
        }
    }
    // DEBUG: once per network root and report step, how the wells under this
    // network's leaf groups end up (or do not end up) in the solve.
    if (group_tree_debug_logged_.insert({root.name(), reportStepIdx}).second) {
        const auto kindName = [](const auto kind) {
            switch (kind) {
            case Sys::WellKind::Pinned: return "pinned";
            case Sys::WellKind::Group: return "group";
            case Sys::WellKind::Thp: return "thp";
            }
            return "?";
        };
        const auto catName = [](const ProdNodeModeCategory c) {
            switch (c) {
            case ProdNodeModeCategory::Group: return "Group";
            case ProdNodeModeCategory::Individual: return "Individual";
            case ProdNodeModeCategory::None: return "None";
            case ProdNodeModeCategory::Transparent: return "Transparent";
            }
            return "?";
        };
        log(fmt::format("GroupTreeDebug root {} step {}: {} network nodes, {} flat entries, {} wells in the solve",
                        root.name(), reportStepIdx, order.size(), flat.size(), system.numWells()));
        for (std::size_t n = 0; n < order.size(); ++n) {
            std::string wells;
            const bool is_group = schedule.hasGroup(order[n], reportStepIdx);
            if (is_group) {
                for (const auto& w : schedule.getGroup(order[n], reportStepIdx).wells()) {
                    wells += " " + w;
                }
            }
            log(fmt::format("GroupTreeDebug   node {} {} ({}; leaf {}){}{}", n, order[n],
                            is_group ? "group" : "not a group",
                            network.downtree_branches(order[n]).empty() ? "yes" : "no",
                            wells.empty() ? "" : "; wells:", wells));
        }
        for (const auto& e : flat) {
            std::string own, children;
            for (const auto& w : e.ownWells) { own += " " + w.name; }
            for (const auto& c : e.activeChildren) { children += " " + c.name; }
            log(fmt::format("GroupTreeDebug   flat {} ({}, mode {}, target {:.5g}){}{}{}{}{}", e.name,
                            e.type == ProdNodeType::Well ? "well" : "group",
                            static_cast<int>(e.mode),
                            e.target * unit::day,
                            own.empty() ? "" : "; own:", own, children.empty() ? "" : "; children:", children,
                            e.satelliteRates.has_value() ? "; satellite" : ""));
        }
        std::set<std::string> in_solve;
        for (const auto& w : system.wells()) {
            in_solve.insert(w.name);
            log(fmt::format("GroupTreeDebug   solve well {}: {}, node {}{}", w.name, kindName(w.kind), w.node,
                            w.reopen_candidate ? ", candidate" : ""));
        }
        for (std::size_t n = 1; n < order.size(); ++n) {
            if (!schedule.hasGroup(order[n], reportStepIdx)) {
                continue;
            }
            for (const auto& w : schedule.getGroup(order[n], reportStepIdx).wells()) {
                if (in_solve.count(w) > 0) {
                    continue;
                }
                const auto t = balancedTree.find(w);
                const auto d = wellNetworkData.find(w);
                log(fmt::format("GroupTreeDebug   NOT in solve: {} at node {}: balanced tree {}, gather {}",
                                w, order[n],
                                t == balancedTree.end() ? std::string("absent")
                                    : fmt::format("{} (parent {}, rate O/W/G {:.4g}/{:.4g}/{:.4g})",
                                                  catName(t->second.modeCategory), t->second.parent,
                                                  -t->second.rates[0] * unit::day, -t->second.rates[1] * unit::day,
                                                  -t->second.rates[2] * unit::day),
                                d == wellNetworkData.end() ? std::string("absent")
                                    : fmt::format("has_ipr {}, stopped {}, trial {}", d->second.has_ipr,
                                                  d->second.stopped, d->second.has_trial_ipr)));
            }
        }
    }

    // Below a node table's first FLO value the fractions of its flow at the
    // previous converged solve apply; nodes without one get their wells'
    // reference fractions in finalize().
    for (std::size_t n = 1; n < order.size(); ++n) {
        if (const auto f = group_tree_node_fractions_.find(order[n]); f != group_tree_node_fractions_.end()) {
            system.setNodeExplicitFractions(static_cast<int>(n), f->second[0], f->second[1]);
        }
    }
    system.finalize();

    std::vector<Scalar> guess(system.numNodes());
    for (std::size_t n = 1; n < order.size(); ++n) {
        const auto it = this->nodePressures(details::NetworkDomain::Production).find(order[n]);
        guess[n - 1] = (it != this->nodePressures(details::NetworkDomain::Production).end())
            ? it->second : *root.terminal_pressure();
    }
    system.setBranchFlattening(branch_flattening_);

    // The step control (--group-tree-network-step): the full Newton step, or a
    // backtracking line search on the residual norm. The group-tree residual
    // has kinks (piecewise-linear tables, the slope-limited lookups' bridged
    // stretches), on which full steps can cycle although the solution is
    // unique (Norne, E-1H: a period-3 cycle with every step raising the
    // residual).
    const auto solveSystem = [this](Sys& sys, const std::vector<Scalar>& start) {
        return network_line_search_
            ? NetworkSolve::solve(sys, start, kGroupTreeSolveParams<Scalar>, NetworkSolve::LineSearch{})
            : NetworkSolve::solve(sys, start, kGroupTreeSolveParams<Scalar>, NetworkSolve::FullStep{});
    };
    auto result = solveSystem(system, guess);
    // DEBUG: cap switches (iterations on which some Thp well changed between
    // its cap and its tubing row) per solve.
    log(fmt::format("GroupTreeSwitches pass flat converged {} iterations {} switches {}",
                    result.converged ? 1 : 0, result.iterations, result.switches));
    if (!result.converged) {
        // TEMPORARY diagnostics of the last iterations (NetworkSolve::kDiagnosticIterations).
        log(fmt::format("Network diag: group-tree solve under {} at report step {} failed; "
                        "{} nodes, {} wells, {} unknowns{}",
                        root.name(), reportStepIdx, system.numNodes(), system.numWells(),
                        system.size(), system.usesAnalyticJacobian() ? ", analytic Jacobian" : ""));
        for (const auto& line : result.diagnostics) {
            log(line);
        }
        // A failure with reopen candidates: solve once more without them, so
        // that they stay stopped this round rather than the whole tree going
        // to the relaxed update.
        if (num_candidates > 0) {
            log(fmt::format("Network: the group-tree solve under {} at report step {} did not converge "
                            "with {} reopen candidates (max residual {:.3e}{}); solving again without them.",
                            root.name(), reportStepIdx, num_candidates, result.residual,
                            result.control_trace.empty()
                                ? std::string{}
                                : fmt::format(", controls {}", result.control_trace)));
            return solveGroupTree(network, reportStepIdx, root, balancedTree,
                                  /*offerReopenCandidates=*/false);
        }
        // Neither the flattened nor any other solve converged: continue with
        // the best iterate -- never worse than the start by its residual --
        // rather than the relaxed fixed-point update, and say so. No open-set
        // changes are made on it (updateGroupTreeOpenSet()).
        if (well_model_.comm().rank() == 0) {
            const bool far = result.best_residual > kGroupTreeFarResidual<Scalar>;
            const auto msg = fmt::format("Network: the group-tree solve under {} at report step {} did not "
                                         "converge in {} iterations; continuing with its best iterate (max "
                                         "residual {:.3g}{}, {} control switches{}).", root.name(),
                                         reportStepIdx, result.iterations - 1, result.best_residual,
                                         far ? ": FAR from a network solution" : "", result.switches,
                                         result.control_trace.empty()
                                             ? std::string{}
                                             : fmt::format(", controls {}", result.control_trace));
            if (well_model_.groupStateHelper().hasDeferredLogger()) {
                well_model_.groupStateHelper().deferredLogger().warning("NETWORK_NOT_CONVERGED", msg);
            } else {
                OpmLog::warning(msg);
            }
        }
        result.residual = result.best_residual;
        result.node_pressure = result.best_node_pressure;
        result.well_rate = result.best_well_rate;
        result.well_phase_rates = result.best_well_phase_rates;
        result.well_bhp = result.best_well_bhp;
    }
    // With flattened branch tables, a solution with some branch on its
    // flattened stretch is one of the flattened network only: solve the real
    // one from there (those branches unflattened, from this solution). If
    // that fails, test each such node (GroupTreeSystem::nodeMargin()): no
    // flowing state on its real table -- a cliff at network level -- stops
    // the Thp wells under it (updateGroupTreeOpenSet()); otherwise solve
    // again from the crossing. Still failing: warn and keep the flattened
    // solution.
    std::vector<std::string> network_cliff_wells;
    const auto flattened = result.converged ? system.flattenedNodes(result) : std::vector<int>{};
    if (!flattened.empty()) {
        std::string names;
        for (const int i : flattened) {
            names += (names.empty() ? "" : ", ") + system.node(i).name;
        }
        Sys refined = system;
        for (const int i : flattened) {
            refined.setNodeFlattening(i, false);
        }
        refined.setStartFromResult(result);
        std::vector<Scalar> from(result.node_pressure.begin() + 1, result.node_pressure.end());
        auto refined_result = solveSystem(refined, from);
        log(fmt::format("GroupTreeSwitches pass real converged {} iterations {} switches {}",
                        refined_result.converged ? 1 : 0, refined_result.iterations,
                        refined_result.switches));   // DEBUG
        bool cliff = false;
        if (!refined_result.converged) {
            bool retry = false;
            std::string margins;
            for (const int i : flattened) {
                const auto m = system.nodeMargin(i, result);
                if (!m.tested) {
                    continue;
                }
                margins += fmt::format(" {} margin {:.3f} bar{}", system.node(i).name, m.margin / unit::barsa,
                                       m.crossing.has_value()
                                           ? fmt::format(", crossing at {:.3f} bar", *m.crossing / unit::barsa)
                                           : std::string{});
                if (m.margin < Scalar{0}) {
                    cliff = true;
                    for (const int w : m.thp_wells) {
                        network_cliff_wells.push_back(system.wells()[w].name);
                    }
                } else if (m.crossing.has_value()) {
                    retry = true;
                    from[i - 1] = *m.crossing;
                    for (const int w : m.thp_wells) {
                        refined.setWellStartBhp(w, refined.thpWellAtNodePressure(w, *m.crossing).first);
                    }
                }
            }
            log(fmt::format("Network: real-table solve under {} at report step {} failed; node test:{}",
                            root.name(), reportStepIdx, margins));
            if (cliff) {
                if (well_model_.comm().rank() == 0) {
                    std::string wells;
                    for (const auto& w : network_cliff_wells) {
                        wells += (wells.empty() ? "" : ", ") + w;
                    }
                    const auto msg = fmt::format("Network: under {} at report step {} the wells {} cannot keep "
                                                 "their branch flowing on its real table (node test:{}); they "
                                                 "are stopped.", root.name(), reportStepIdx, wells, margins);
                    if (well_model_.groupStateHelper().hasDeferredLogger()) {
                        well_model_.groupStateHelper().deferredLogger().info(msg);
                    } else {
                        OpmLog::info(msg);
                    }
                }
            } else if (retry) {
                refined_result = solveSystem(refined, from);
                log(fmt::format("Network: real-table solve under {} at report step {} from the crossing: {} "
                                "in {} iterations.", root.name(), reportStepIdx,
                                refined_result.converged ? "converged" : "not converged",
                                refined_result.iterations));
            }
        }
        if (!cliff && !refined_result.converged) {
            // DEBUG: the failed real-table solve's diagnostics, as for a
            // failed first pass.
            log(fmt::format("Network diag: real-table solve under {} at report step {} failed (branches {})",
                            root.name(), reportStepIdx, names));
            for (const auto& line : refined_result.diagnostics) {
                log(line);
            }
            if (well_model_.comm().rank() == 0) {
                const auto msg = fmt::format("Network: branches {} under {} at report step {} ended on "
                                             "the flattened (falling, low-flow) part of their tables, and "
                                             "the network does not converge with their real tables from "
                                             "there (max residual {:.3e}). Keeping the flattened solution.",
                                             names, root.name(), reportStepIdx, refined_result.residual);
                if (well_model_.groupStateHelper().hasDeferredLogger()) {
                    well_model_.groupStateHelper().deferredLogger().warning("NETWORK_FLATTENED_BRANCH", msg);
                } else {
                    OpmLog::warning(msg);
                }
            }
        } else if (!cliff) {
            log(fmt::format("Network: branches {} under {} at report step {} ended on their flattened "
                            "table; solved again with their real tables in {} iterations.",
                            names, root.name(), reportStepIdx, refined_result.iterations));
            system = std::move(refined);
            result = std::move(refined_result);
        }
    }
    log(fmt::format("Network: {} the production network under {} via the group-tree "
                    "balancer at report step {} in {} iterations{}.",
                    result.converged ? "solved" : "did not solve (best iterate used)",
                    root.name(), reportStepIdx, result.iterations,
                    num_candidates > 0 ? fmt::format(" ({} reopen candidates)", num_candidates)
                                       : std::string{}));
    const auto fractions = system.nodeFractions(result.well_phase_rates);
    {   // DEBUG: the node flows and per-well rates of this converged solve
        // (wells under this root that are not in it are dropped, so they show
        // as absent rather than with an earlier solve's rates).
        const auto flows = system.nodeFlowsFromWellRates(result.well_phase_rates);
        for (std::size_t n = 1; n < order.size(); ++n) {
            if (schedule.hasGroup(order[n], reportStepIdx)) {
                for (const auto& wname : schedule.getGroup(order[n], reportStepIdx).wells()) {
                    group_tree_debug_well_rates_.erase(wname);
                }
            }
        }
        for (std::size_t n = 1; n < order.size(); ++n) {
            group_tree_debug_node_flows_[order[n]] = flows[n];
        }
        for (int w = 0; w < system.numWells(); ++w) {
            group_tree_debug_well_rates_[system.wells()[w].name] = result.well_phase_rates[w];
        }
        // DEBUG: the state this solve hands on -- node pressures (the wells'
        // THP limits at B4) with their liquid rate and GOR, and every well's
        // kind, oil rate and bhp.
        std::string nodes;
        for (std::size_t n = 1; n < order.size(); ++n) {
            const auto& q = flows[n];
            nodes += fmt::format(" {} {:.2f} bar (L {:.0f}, GOR {:.0f})", order[n],
                                 result.node_pressure[n] / unit::barsa,
                                 (q[Sys::kOil] + q[Sys::kWater]) * unit::day,
                                 q[Sys::kOil] > Scalar{0} ? q[Sys::kGas] / q[Sys::kOil] : Scalar{0});
        }
        std::string wells;
        for (int w = 0; w < system.numWells(); ++w) {
            const auto& well = system.wells()[w];
            const char* kind = well.reopen_candidate ? "cand"
                : (well.kind == Sys::WellKind::Thp ? "thp"
                   : (well.kind == Sys::WellKind::Group ? "grp" : "pin"));
            wells += fmt::format(" {} {} O/W/G {:.0f}/{:.0f}/{:.4g} bhp {:.2f};", well.name, kind,
                                 result.well_phase_rates[w][Sys::kOil] * unit::day,
                                 result.well_phase_rates[w][Sys::kWater] * unit::day,
                                 result.well_phase_rates[w][Sys::kGas] * unit::day,
                                 result.well_bhp[w] / unit::barsa);
        }
        log(fmt::format("GroupTreeState under {} at report step {}{}: nodes{}; wells{}", root.name(),
                        reportStepIdx, num_candidates > 0 ? " (with candidates)" : "", nodes, wells));
        // DEBUG: every Active (group-constrained) group of this solve: its
        // target against the sum its target row holds, and each member's
        // network-solve rates next to the well state's current ones.
        const auto& pu = well_model_.phaseUsage();
        const auto reports = system.activeNodeReports(result.well_phase_rates);
        // The same sums from the well state's rates (members not in the well
        // state count with their network-solve rates).
        const auto stateRates = [&](const std::string& name, const std::array<Scalar, 3>& fallback) {
            if (!well_model_.wellState().has(name)) {
                return fallback;
            }
            const auto& ws = well_model_.wellState().well(name);
            return std::array<Scalar, 3>{
                -ws.surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)],
                -ws.surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)],
                -ws.surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)]};
        };
        std::function<Scalar(const std::string&, const std::array<Scalar, 3>&)> stateTotal =
            [&](const std::string& name, const std::array<Scalar, 3>& weights) {
                Scalar sum = 0;
                for (const auto& r : reports) {
                    if (r.name != name) {
                        continue;
                    }
                    for (const auto& m : r.members) {
                        if (m.kind == 'c') {
                            sum += m.efficiency * stateTotal(m.name, weights);
                        } else {
                            const auto q = stateRates(m.name, m.q);
                            sum += m.efficiency * (weights[0] * q[0] + weights[1] * q[1] + weights[2] * q[2]);
                        }
                    }
                }
                return sum;
            };
        for (const auto& a : reports) {
            log(fmt::format("GroupTreeActive under {} at report step {}: {} {} target {:.6g} sum {:.6g} m3/d "
                            "(diff {:+.3g}): own wells {:.6g}, other members {:.6g}{}{}", root.name(),
                            reportStepIdx, a.name, WellProducerCMode2String(a.mode), a.target * unit::day,
                            a.total * unit::day, (a.total - a.target) * unit::day, a.own_total * unit::day,
                            (a.total - a.own_total) * unit::day, a.has_row ? "" : ", no row (no own wells)",
                            (a.has_row && a.own_total < Scalar{0}) ? ", NEGATIVE lambda" : ""));
            log(fmt::format("GroupTreeActive   {} sum from well-state rates {:.6g} m3/d (diff to target {:+.3g})",
                            a.name, stateTotal(a.name, a.weights) * unit::day,
                            (stateTotal(a.name, a.weights) - a.target) * unit::day));
            for (const auto& m : a.members) {
                if (m.kind == 'c') {
                    log(fmt::format("GroupTreeActive     group {} eff {:.4g}: sum {:.6g} m3/d", m.name,
                                    m.efficiency, m.mode_rate * unit::day));
                    continue;
                }
                std::string state = "not in well state";
                if (well_model_.wellState().has(m.name)) {
                    const auto& ws = well_model_.wellState().well(m.name);
                    state = fmt::format("well state O/W/G {:.5g}/{:.5g}/{:.5g} bhp {:.3f} thp {:.3f} {}",
                        -ws.surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)] * unit::day,
                        -ws.surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)] * unit::day,
                        -ws.surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)] * unit::day,
                        ws.bhp / unit::barsa, ws.thp / unit::barsa,
                        WellProducerCMode2String(ws.production_cmode));
                    // The network's bhp for it, and its IPR at the well
                    // state's bhp (= its rates there if the IPR is anchored
                    // at this state), with the IPR's slopes.
                    const auto& sw = system.wells()[m.well];
                    state += fmt::format("; net bhp {:.3f}; IPR at state bhp O/W/G {:.5g}/{:.5g}/{:.5g}, "
                                         "slope O/W/G {:.4g}/{:.4g}/{:.4g} m3/d/bar",
                                         result.well_bhp[m.well] / unit::barsa,
                                         (sw.ipr_a[Sys::kOil] + sw.ipr_b[Sys::kOil] * ws.bhp) * unit::day,
                                         (sw.ipr_a[Sys::kWater] + sw.ipr_b[Sys::kWater] * ws.bhp) * unit::day,
                                         (sw.ipr_a[Sys::kGas] + sw.ipr_b[Sys::kGas] * ws.bhp) * unit::day,
                                         sw.ipr_b[Sys::kOil] * unit::day * unit::barsa,
                                         sw.ipr_b[Sys::kWater] * unit::day * unit::barsa,
                                         sw.ipr_b[Sys::kGas] * unit::day * unit::barsa);
                    const auto& gt = ws.use_group_target_fallback ? ws.group_target_fallback : ws.group_target;
                    if (m.kind == 'g' || ws.production_cmode == Well::ProducerCMode::GRUP) {
                        state += gt.has_value()
                            ? fmt::format("; group target {} {} {:.6g} m3/d{}", gt->group_name,
                                          Group::ProductionCMode2String(gt->production_cmode),
                                          gt->target_value * unit::day,
                                          ws.use_group_target_fallback ? " (fallback)" : "")
                            : std::string("; no group target");
                    }
                }
                log(fmt::format("GroupTreeActive     well {} {} eff {:.4g}{}: O/W/G {:.5g}/{:.5g}/{:.5g}, "
                                "mode rate {:.6g} m3/d; {}", m.name,
                                m.kind == 'g' ? "grp" : (m.kind == 't' ? "thp" : "pin"), m.efficiency,
                                m.kind == 'g' ? fmt::format(" guide {:.4g}", m.guide_rate) : std::string{},
                                m.q[Sys::kOil] * unit::day, m.q[Sys::kWater] * unit::day,
                                m.q[Sys::kGas] * unit::day, m.mode_rate * unit::day, state));
            }
        }
    }
    for (std::size_t n = 1; n < order.size() && result.converged; ++n) {
        if (fractions[n].has_value()) {
            group_tree_node_fractions_[order[n]] = *fractions[n];
        }
    }
    GroupTreeSolve solved{std::move(system), std::move(order), std::move(result),
                          std::move(network_cliff_wells)};
    logCategorizationCheck(solved, reportStepIdx, balancedTree);
    return solved;
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
logCategorizationCheck(const GroupTreeSolve& solved,
                       const int reportStepIdx,
                       const ProdGroupTreeBalancer::Tree<Scalar>& balancedTree) const
{
    using Sys = NetworkSolve::GroupTreeSystem<Scalar>;
    using Issue = typename Sys::CategorizationIssue;
    using CMode = Well::ProducerCMode;
    // The solve is the same on every rank.
    if (well_model_.comm().rank() != 0) {
        return;
    }
    const auto& schedule = well_model_.schedule();
    const auto& summary_state = well_model_.summaryState();
    const auto& system = solved.system;
    const auto& result = solved.result;
    const typename Sys::CategorizationTolerance tol{};

    // The wells' own limits from the schedule (not their THP limit, which in
    // the network is the node pressure); RESV with the balancer's coefficients.
    std::unordered_map<std::string, typename Sys::IndividualLimits> limits;
    for (const auto& w : system.wells()) {
        if (!schedule.hasWell(w.name, reportStepIdx)) {
            continue;   // a satellite group
        }
        const auto& well = schedule.getWell(w.name, reportStepIdx);
        if (!well.isProducer()) {
            continue;
        }
        const auto controls = well.productionControls(summary_state);
        typename Sys::IndividualLimits lim;
        const auto add = [&controls, &lim](const CMode mode, const double value) {
            if (controls.hasControl(mode) && value > 0.0) {
                lim.rates.emplace_back(mode, static_cast<Scalar>(value));
            }
        };
        add(CMode::ORAT, controls.oil_rate);
        add(CMode::WRAT, controls.water_rate);
        add(CMode::GRAT, controls.gas_rate);
        add(CMode::LRAT, controls.liquid_rate);
        if (const auto t = balancedTree.find(w.name); t != balancedTree.end()) {
            lim.resv_coeff = t->second.resvCoeff;
            add(CMode::RESV, controls.resv_rate);
        }
        if (controls.hasControl(CMode::BHP)) {
            lim.bhp = static_cast<Scalar>(controls.bhp_limit);
        }
        limits.emplace(w.name, std::move(lim));
    }

    // A THP well's guide rate in an Active group's mode, as the balancer
    // looks it up for the group-controlled wells (extractFlatNetworkInput()).
    const auto& guide_rate = well_model_.guideRate();
    const auto guideRateOf = [&](const std::string& name, const CMode mode) -> std::optional<Scalar> {
        const auto t = balancedTree.find(name);
        if (t == balancedTree.end() || !t->second.availableForGroupControl
            || !(guide_rate.has(name) || guide_rate.hasPotentials(name))) {
            return std::nullopt;
        }
        GuideRateModel::Target target{GuideRateModel::Target::NONE};
        switch (mode) {
        case CMode::ORAT: target = GuideRateModel::Target::OIL; break;
        case CMode::WRAT: target = GuideRateModel::Target::WAT; break;
        case CMode::GRAT: target = GuideRateModel::Target::GAS; break;
        case CMode::LRAT: target = GuideRateModel::Target::LIQ; break;
        case CMode::RESV: target = GuideRateModel::Target::RES; break;
        default: return std::nullopt;
        }
        // RateVector is {oil, gas, water}; the tree's rates are negative for production.
        const auto& r = t->second.initialRates;
        const GuideRate::RateVector rv{-r[0], -r[2], -r[1]};
        const Scalar g = guide_rate.get(name, target, rv);
        return g > Scalar{0} ? std::optional<Scalar>(g) : std::nullopt;
    };

    auto issues = system.checkCategorization(result, limits, guideRateOf, tol);

    // Pass-through groups (in the balancer nothing of their own binds)
    // against their own limits: the network's rates for the wells in this
    // solve, the balancer's for the others. Groups without a well in this
    // solve are left out.
    std::unordered_map<std::string, std::array<Scalar, 3>> solved_rates;
    for (int w = 0; w < system.numWells(); ++w) {
        solved_rates[system.wells()[w].name] = result.well_phase_rates[w];
    }
    struct Total
    {
        std::array<Scalar, 3> q{};
        bool in_solve = false;
    };
    std::function<Total(const std::string&)> total = [&](const std::string& name) {
        Total t;
        const auto it = balancedTree.find(name);
        if (it == balancedTree.end()) {
            return t;
        }
        const auto& node = it->second;
        if (node.type == ProdNodeType::Well || node.isSatellite) {
            const auto s = solved_rates.find(name);
            if (node.type == ProdNodeType::Well && s != solved_rates.end()) {
                t.q = s->second;
                t.in_solve = true;
            } else {
                for (int c = 0; c < 3; ++c) {
                    t.q[c] = -node.rates[c];
                }
            }
            return t;
        }
        for (const auto& child : node.children) {
            const auto ct = total(child);
            const auto ch = balancedTree.find(child);
            const Scalar eff = ch != balancedTree.end() ? ch->second.efficiencyFactor : Scalar{1};
            for (int c = 0; c < 3; ++c) {
                t.q[c] += eff * ct.q[c];
            }
            t.in_solve = t.in_solve || ct.in_solve;
        }
        return t;
    };
    struct PassThroughIssue
    {
        std::string name;
        CMode mode;
        Scalar value;
        Scalar bound;
    };
    std::vector<PassThroughIssue> pass_through;
    for (const auto& [name, node] : balancedTree) {
        if (node.type != ProdNodeType::Group || node.isSatellite || node.Limits.empty()
            || node.modeCategory == ProdNodeModeCategory::Individual) {
            continue;
        }
        const auto t = total(name);
        if (!t.in_solve) {
            continue;
        }
        for (const auto& [mode, limit] : node.Limits) {
            std::array<Scalar, 3> weights{};
            switch (mode) {
            case CMode::ORAT: weights = {1, 0, 0}; break;
            case CMode::WRAT: weights = {0, 1, 0}; break;
            case CMode::GRAT: weights = {0, 0, 1}; break;
            case CMode::LRAT: weights = {1, 1, 0}; break;
            case CMode::RESV: weights = node.resvCoeff; break;
            default: continue;
            }
            const Scalar v = weights[0] * t.q[0] + weights[1] * t.q[1] + weights[2] * t.q[2];
            if (v > limit + tol.rel * std::abs(limit) + tol.abs) {
                pass_through.push_back({name, mode, v, limit});
            }
        }
    }

    // One line per solve (to count how often a solve ends with its
    // categorization contradicted), one per issue.
    const auto log = [this](const std::string& msg) {
        if (well_model_.groupStateHelper().hasDeferredLogger()) {
            well_model_.groupStateHelper().deferredLogger().debug(msg);
        } else {
            OpmLog::debug(msg);
        }
    };
    const auto& root = solved.order.front();
    const std::string state = result.converged ? "converged" : "best iterate";
    if (issues.empty() && pass_through.empty()) {
        log(fmt::format("GroupTreeCheck under {} at report step {} ({}): consistent", root, reportStepIdx,
                        state));
        return;
    }
    const std::array<std::pair<typename Issue::Kind, const char*>, 7> kinds{{
        {Issue::Kind::NegativeLambda, "negative lambda"},
        {Issue::Kind::GroupOverCapacity, "group-controlled over THP capacity"},
        {Issue::Kind::OverLimit, "over own limit"},
        {Issue::Kind::ThpOverShare, "THP over guide-rate share"},
        {Issue::Kind::PinnedOverCapacity, "fixed-rate over THP capacity"},
        {Issue::Kind::ActiveOverTarget, "Active group over target"},
        {Issue::Kind::ActiveCannotMeet, "Active group cannot meet target"},
    }};
    std::string counts;
    for (const auto& [kind, label] : kinds) {
        const auto n = std::count_if(issues.begin(), issues.end(),
                                     [kind = kind](const Issue& i) { return i.kind == kind; });
        if (n > 0) {
            counts += fmt::format("{}{} {}", counts.empty() ? "" : ", ", n, label);
        }
    }
    if (!pass_through.empty()) {
        counts += fmt::format("{}{} pass-through group over own limit", counts.empty() ? "" : ", ",
                              pass_through.size());
    }
    log(fmt::format("GroupTreeCheck under {} at report step {} ({}): inconsistent: {}", root, reportStepIdx,
                    state, counts));

    std::unordered_map<std::string, int> index;
    for (int w = 0; w < system.numWells(); ++w) {
        index[system.wells()[w].name] = w;
    }
    const auto wellKind = [&](const std::string& name) -> std::string {
        const auto it = index.find(name);
        if (it == index.end()) {
            return "";
        }
        const auto& w = system.wells()[it->second];
        const char* kind = w.reopen_candidate ? "candidate"
            : (w.kind == Sys::WellKind::Thp ? "THP" : (w.kind == Sys::WellKind::Group ? "group" : "fixed"));
        return fmt::format("{}, node {} at {:.3f} bar", kind, solved.order[w.node],
                           result.node_pressure[w.node] / unit::barsa);
    };
    const auto rate = [](const Scalar v) { return v * unit::day; };
    // WellProducerCMode2String() throws for CMODE_UNDEFINED (a comparison in FLO).
    const auto modeName = [](const CMode m) {
        return m == CMode::CMODE_UNDEFINED ? std::string("FLO") : WellProducerCMode2String(m);
    };
    for (const auto& i : issues) {
        const std::string mode = modeName(i.mode);
        switch (i.kind) {
        case Issue::Kind::NegativeLambda:
            log(fmt::format("GroupTreeCheck   negative lambda: {} ({}): its group-controlled wells total "
                            "{:.6g} m3/d", i.name, mode, rate(i.value)));
            break;
        case Issue::Kind::GroupOverCapacity:
            log(fmt::format("GroupTreeCheck   group-controlled over THP capacity: {} ({}): allocated {:.6g} "
                            "{} m3/d, its tubing lifts {:.6g}", i.name, wellKind(i.name), rate(i.value), mode,
                            rate(i.bound)));
            break;
        case Issue::Kind::OverLimit:
            if (i.mode == CMode::BHP) {
                log(fmt::format("GroupTreeCheck   over own limit: {} ({}): bhp {:.3f} bar below its BHP "
                                "limit {:.3f} bar", i.name, wellKind(i.name), i.value / unit::barsa,
                                i.bound / unit::barsa));
            } else {
                log(fmt::format("GroupTreeCheck   over own limit: {} ({}): {:.6g} {} m3/d above its limit "
                                "{:.6g}", i.name, wellKind(i.name), rate(i.value), mode, rate(i.bound)));
            }
            break;
        case Issue::Kind::ThpOverShare:
            log(fmt::format("GroupTreeCheck   THP over guide-rate share: {} ({}): {:.6g} {} m3/d above its "
                            "share {:.6g}", i.name, wellKind(i.name), rate(i.value), mode, rate(i.bound)));
            break;
        case Issue::Kind::PinnedOverCapacity:
            log(fmt::format("GroupTreeCheck   fixed-rate over THP capacity: {} ({}): FLO {:.6g} m3/d, its "
                            "tubing lifts {:.6g}", i.name, wellKind(i.name), rate(i.value), rate(i.bound)));
            break;
        case Issue::Kind::ActiveOverTarget:
            log(fmt::format("GroupTreeCheck   Active group over target: {} ({}, no group-controlled wells): "
                            "members total {:.6g} m3/d, target {:.6g}", i.name, mode, rate(i.value),
                            rate(i.bound)));
            break;
        case Issue::Kind::ActiveCannotMeet:
            log(fmt::format("GroupTreeCheck   Active group cannot meet target: {} ({}): members reach at "
                            "most {:.6g} m3/d, target {:.6g}", i.name, mode, rate(i.value), rate(i.bound)));
            break;
        }
    }
    for (const auto& p : pass_through) {
        log(fmt::format("GroupTreeCheck   pass-through group over own limit: {}: {:.6g} {} m3/d above its "
                        "limit {:.6g}", p.name, rate(p.value), modeName(p.mode),
                        rate(p.bound)));
    }
}

template<typename Scalar, typename IndexTraits>
std::optional<std::map<std::string, Scalar>>
BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
groupTreeProductionNodePressures(const Network::ExtNetwork& network,
                                 const int reportStepIdx,
                                 const Network::Node& root,
                                 const ProdGroupTreeBalancer::Tree<Scalar>& balancedTree) const
{
    const auto solved = solveGroupTree(network, reportStepIdx, root, balancedTree);
    if (!solved.has_value()) {
        return std::nullopt;
    }
    std::map<std::string, Scalar> pressures;
    for (std::size_t n = 0; n < solved->order.size(); ++n) {
        pressures[solved->order[n]] = solved->result.node_pressure[n];
    }
    return pressures;
}

template<typename Scalar, typename IndexTraits>
bool
BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
updateGroupTreeOpenSet(const Network::ExtNetwork& network,
                       const int reportStepIdx,
                       const Network::Node& root,
                       const ProdGroupTreeBalancer::Tree<Scalar>& balancedTree,
                       const bool offerReopenCandidates,
                       DeferredLogger& deferred_logger)
{
    using Sys = NetworkSolve::GroupTreeSystem<Scalar>;
    const auto solved = solveGroupTree(network, reportStepIdx, root, balancedTree,
                                       offerReopenCandidates);
    if (!solved.has_value()) {
        return false;   // the normal give-up/relaxed-update path handles this
    }
    const auto findWell = [this](const std::string& name) -> WellInterfaceGeneric<Scalar, IndexTraits>* {
        for (auto* well : well_model_.genericWells()) {
            if (well->name() == name) {
                return well;
            }
        }
        return nullptr;   // not on this rank
    };
    const bool report = well_model_.comm().rank() == 0;
    bool changed = false;
    // An unconverged solve (its best iterate, see solveGroupTree()) decides
    // nothing about the open set.
    if (!solved->result.converged) {
        return false;
    }

    // Reopen the candidates the solve lets flow on the real tubing curve.
    // The dynamic status (openWell()) is what the well model's own solves
    // read; the well state is seeded with the network's operating point and
    // the trial IPR the solve used, which is what the balancer and the next
    // network solve see until the well's own solve refreshes them.
    const auto& pu = well_model_.phaseUsage();
    std::array<int, Sys::NP> pos{};
    pos[Sys::kOil]   = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
    pos[Sys::kWater] = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
    pos[Sys::kGas]   = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
    const auto& cliff_wells = solved->network_cliff_wells;
    const auto atNetworkCliff = [&cliff_wells](const std::string& name) {
        return std::find(cliff_wells.begin(), cliff_wells.end(), name) != cliff_wells.end();
    };
    for (const auto& outcome : solved->system.reopenOutcomes(solved->result)) {
        if (!outcome.reopens || atNetworkCliff(outcome.name)) {
            // An open well without any rate (awaitsFlow()) that cannot flow
            // here is stopped, as if it had been stopped all along.
            if (auto* well = findWell(outcome.name);
                well != nullptr && well->awaitsFlow(well_model_.wellState()[well->indexOfWell()])) {
                well->stopWell(WellInterfaceGeneric<Scalar, IndexTraits>::StopReason::Network);
            }
            if (report) {
                deferred_logger.debug(fmt::format(
                    "Network: well {} stays stopped under the group-tree balancer at report "
                    "step {}: {} (bhp {:.3f} bar, thp {:.3f} bar).",
                    outcome.name, reportStepIdx,
                    outcome.at_shutin ? "no flow at this node pressure"
                                      : "it would flow only on the flattened tubing curve",
                    outcome.bhp / unit::barsa, outcome.thp / unit::barsa));
            }
            continue;
        }
        changed = true;
        if (auto* well = findWell(outcome.name)) {
            well->openWell();
            well->holdStopped(false);
            auto& ws = well_model_.wellState()[well->indexOfWell()];
            for (int ph = 0; ph < Sys::NP; ++ph) {
                ws.surface_rates[pos[ph]] = -outcome.q[ph];
            }
            ws.bhp = outcome.bhp;
            ws.thp = outcome.thp;
            ws.production_cmode = Well::ProducerCMode::THP;
            ws.implicit_ipr_a = ws.stopped_ipr_a;
            ws.implicit_ipr_b = ws.stopped_ipr_b;
        }
        if (report) {
            deferred_logger.debug(fmt::format(
                "Network: reopening well {} under the group-tree balancer at report step {} "
                "(oil {:.4g} m3/day at bhp {:.3f} bar, thp {:.3f} bar).",
                outcome.name, reportStepIdx, outcome.q[Sys::kOil] * unit::day,
                outcome.bhp / unit::barsa, outcome.thp / unit::barsa));
        }
    }

    // Stop the single worst open well whose converged point lies on the
    // flattened part of its tubing curve -- a point it cannot actually
    // operate at. The dynamic status is what the well model's own solves
    // re-evaluate; WellState::stopWell() (ws.status) would be persistent.
    // Wells that cannot keep their branch flowing on its real table (a cliff
    // at network level, solveGroupTree()): stopped, all of them at once.
    for (const auto& name : cliff_wells) {
        auto* well = findWell(name);
        if (well != nullptr && !well->wellIsStopped()) {
            well->stopWell(WellInterfaceGeneric<Scalar, IndexTraits>::StopReason::Network);
        }
        changed = true;
    }
    if (!cliff_wells.empty()) {
        return changed;
    }
    if (const auto worst = solved->system.worstCliffViolation(solved->result)) {
        changed = true;
        if (auto* well = findWell(*worst)) {
            well->stopWell(WellInterfaceGeneric<Scalar, IndexTraits>::StopReason::Network);
        }
        if (report) {
            deferred_logger.info(fmt::format(
                "Network: well {} converged past its own tubing-curve cliff under the group-tree "
                "balancer at report step {}; stopping it and rebalancing.", *worst, reportStepIdx));
        }
    }
    return changed;
}

template<typename Scalar, typename IndexTraits>
Scalar
BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
updatePressures(const int reportStepIdx,
                const Scalar damping_factor,
                const Scalar upper_update_bound,
                const details::DomainScope scope)
{
    OPM_TIMEFUNCTION();
    for (auto& invalid_nodes : this->domain_invalid_nodes_) {
        invalid_nodes.clear();
    }
    this->invalid_nodes_report_step_ = -1;

    if (!details::anyNetworkActive(well_model_.schedule(), reportStepIdx)) {
        return 0.0;
    }

    this->syncProductionDomainState_();
    const auto previous_node_pressures = this->domain_node_pressures_;

    // Production nodes placed by the group-tree solve: handed to the wells
    // undamped (see below).
    std::set<std::string> group_tree_nodes;
    for (const auto& network : details::activeNetworks(well_model_.schedule(), reportStepIdx)) {
        if (!details::inScope(scope, network.domain)) {
            continue;
        }
        NetworkPressures result;
        if (network.domain == details::NetworkDomain::Production) {
            result = this->computePressures(network.network.get(),
                                            *well_model_.getVFPProperties().getProd(),
                                            well_model_.schedule().getUnits(),
                                            reportStepIdx,
                                            well_model_.comm());
            if (this->group_tree_solver_ && this->balanced_group_tree_.has_value()) {
                // Per independent root; a tree the group-tree solve declines
                // keeps the fixed-point answer computed above.
                for (const auto& tree : independentRoots(network.network.get())) {
                    if (auto solved = this->groupTreeProductionNodePressures(
                            network.network.get(), reportStepIdx, tree.get(), *this->balanced_group_tree_)) {
                        for (const auto& [name, pressure] : *solved) {
                            result.node_pressures[name] = pressure;
                            result.invalid_nodes.erase(name);
                            group_tree_nodes.insert(name);
                        }
                    }
                }
            }
        } else {
            result = this->computePressures(network.network.get(),
                                            *well_model_.getVFPProperties().getInj(),
                                            well_model_.schedule().getUnits(),
                                            reportStepIdx,
                                            well_model_.comm());
        }
        this->nodePressures(network.domain) = std::move(result.node_pressures);
        this->branchData(network.domain) = std::move(result.branch_data);
        this->invalidNodes(network.domain) = std::move(result.invalid_nodes);
    }
    this->syncLegacyProductionState_();

    // here, the network imbalance is the difference between the previous nodal pressure and the new nodal pressure
    Scalar network_imbalance = 0.;
    if (!this->active()) {
        return network_imbalance;
    }

    for (const auto& network : details::activeNetworks(well_model_.schedule(), reportStepIdx)) {
        if (!details::inScope(scope, network.domain)) {
            continue;
        }
        auto& domain_pressures = this->nodePressures(network.domain);
        const auto& invalid = this->invalidNodes(network.domain);
        const auto& previous_domain_pressures = previous_node_pressures[details::domainIndex(network.domain)];

        if (!invalid.empty()) {
            // The VFP tables gave no pressure for these nodes (rate/pressure outside what
            // the tables can deliver). Keep the previous value and report the network as
            // unbalanced so that the wells get another chance to move into range.
            network_imbalance = std::max(network_imbalance, upper_update_bound);
            if (this->invalid_nodes_report_step_ != reportStepIdx) {
                this->invalid_nodes_report_step_ = reportStepIdx;
                OpmLog::warning(fmt::format("Network: no VFP solution for node(s) {} at report step {}; "
                                            "keeping the previous node pressure(s).",
                                            fmt::join(invalid, ", "), reportStepIdx + 1));
            }
        }

        if (!previous_domain_pressures.empty()) {
            for (auto& [name, new_pressure]: domain_pressures) {
                if (network.domain == details::NetworkDomain::Production
                    && this->node_pressure_overrides_.count(name) > 0) {
                    // Imposed by a reservoir coupling master: neither damped nor part of
                    // the imbalance, which only measures this network's own iteration.  It
                    // is already at the imposed value, see setNodePressureOverrides(), so
                    // what is computed from it before this update is not stale.
                    continue;
                }
                if (previous_domain_pressures.count(name) <= 0) {
                    if (std::abs(new_pressure) > network_imbalance) {
                        network_imbalance = std::abs(new_pressure);
                    }
                    continue;
                }

                const auto pressure = previous_domain_pressures.at(name);
                if (network.domain == details::NetworkDomain::Production
                    && group_tree_nodes.count(name) > 0) {
                    // The group-tree solve is the network's answer for the
                    // balancer's categorization: hand it to the wells as is. A
                    // well it puts past its tubing-curve cliff is the cliff
                    // diagnosis' business, not something to hide by creeping
                    // towards the answer.
                    network_imbalance = std::max(network_imbalance, std::abs(new_pressure - pressure));
                    continue;
                }
                if (invalid.count(name) > 0) {
                    new_pressure = pressure;
                    continue;
                }
                const Scalar change = (new_pressure - pressure);
                if (std::abs(change) > network_imbalance) {
                    network_imbalance = std::abs(change);
                }
                // We dampen the nodal pressure change during one iteration since our nodal pressure calculation
                // is somewhat explicit. There is a relative dampening factor applied to the update value, and also
                // the maximum update is limited (to 5 bar by default, can be changed with --network-max-pressure-update-in-bars).
                const Scalar damped_change = std::min(damping_factor * std::abs(change), upper_update_bound);
                const Scalar sign = change > 0 ? 1. : -1.;
                new_pressure = pressure + sign * damped_change;
            }
            continue;
        }

        for (const auto& [name, pressure]: domain_pressures) {
            if (std::abs(pressure) > network_imbalance) {
                network_imbalance = std::abs(pressure);
            }
        }
    }
    this->syncLegacyProductionState_();

    for (auto& well : well_model_.genericWells()) {
        if (!well->wellEcl().predictionMode()) {
            continue;
        }
        const auto domain = details::domainForWell(*well);
        if (!domain.has_value() || !details::inScope(scope, *domain)) {
            continue;
        }
        // The network cannot put a well without a VFP table under THP
        // control, so no THP limit is imposed on such a well, while its
        // rates still contribute to the network flows.
        if (well->wellEcl().vfp_table_number() > 0) {
            const auto it = this->nodePressures(*domain).find(well->wellEcl().groupName());
            if (it != this->nodePressures(*domain).end()) {
                if (this->invalidNodes(*domain).count(well->wellEcl().groupName()) > 0) {
                    // No valid leaf pressure this iteration; keep the well's current THP limit.
                    continue;
                }
                // The well belongs to a group that has a network pressure constraint;
                // set the dynamic THP constraint of the well accordingly.
                if (well->isProducer()) {
                    // For producers the leaf-node pressure is the group wellhead THP;
                    // apply it directly as a dynamic THP constraint.
                    this->imposeWellThpLimit(*well, it->second);
                    SingleWellState<Scalar, IndexTraits>& ws = well_model_.wellState()[well->indexOfWell()];
                    const bool thp_is_limit = ws.production_cmode == Well::ProducerCMode::THP;
                    // TODO: not sure why the thp is NOT updated properly elsewhere
                    if (thp_is_limit) {
                        ws.thp = well->getTHPConstraint(well_model_.summaryState());
                    }
                } else if (well->isInjector() && well->wellEcl().vfp_table_number() > 0) {
                    // For injectors, apply the network leaf-node pressure as a dynamic THP only
                    // if it falls within the individual well's VFPINJ table THP range.
                    // If P_leaf is outside the range, computeBhpAtThpLimitInj would extrapolate
                    // to invalid values and mark the well inoperable, causing a rate collapse.
                    const auto& inj_vfp = *well_model_.getVFPProperties().getInj();
                    const int table_id = well->wellEcl().injectionControls(
                        well_model_.summaryState()).vfp_table_number;
                    if (inj_vfp.hasTable(table_id)) {
                        const auto& thp_axis = inj_vfp.getTable(table_id).getTHPAxis();
                        const Scalar min_thp = static_cast<Scalar>(thp_axis.front());
                        const Scalar max_thp = static_cast<Scalar>(thp_axis.back());
                        const Scalar new_limit = it->second;
                        if (new_limit >= min_thp && new_limit <= max_thp) {
                            this->imposeWellThpLimit(*well, new_limit);
                            SingleWellState<Scalar, IndexTraits>& ws =
                                well_model_.wellState()[well->indexOfWell()];
                            const bool thp_is_limit =
                                ws.injection_cmode == Well::InjectorCMode::THP;
                            if (thp_is_limit) {
                                ws.thp = well->getTHPConstraint(well_model_.summaryState());
                            }
                        }
                    }
                }
            }
        }
    }
    return network_imbalance;
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
assignNodeAndBranchValues(std::map<std::string, data::NodeData>& nodevalues,
                          std::map<std::string, data::BranchData>& branchvalues,
                          std::map<std::string, data::BranchData>& converged_branchvalues,
                          const int reportStepIdx) const
{
    nodevalues.clear();
    branchvalues.clear();
    converged_branchvalues.clear();
    if (reportStepIdx < 0) return;
    for (const auto& [node, pressure] : node_pressures_) {
        nodevalues.emplace(node, data::NodeData{pressure});
        // Assign node values of well groups to GPR:WELLNAME
        const auto& sched = well_model_.schedule();
        if (!sched.hasGroup(node, reportStepIdx)) {
            continue;
        }
        const auto& group = sched.getGroup(node, reportStepIdx);
        for (const std::string& wellname : group.wells()) {
            nodevalues.emplace(wellname, data::NodeData{pressure});
        }
    }
    for (const auto& [branch, branch_data] : branch_data_) {
        branchvalues.emplace(branch, branch_data);
        // Skip wells (do not consider well->group a branch, at least not for now)
    }

    const auto& network = well_model_.schedule()[reportStepIdx].network();
    if (!network.active()) {
        return;
    }

    auto converged = this->computePressures(network,
                                            *well_model_.getVFPProperties().getProd(),
                                            well_model_.schedule().getUnits(),
                                            reportStepIdx,
                                            well_model_.comm());
    const auto& converged_pressures = converged.node_pressures;
    // DEBUG: where GPR (the last network update's pressures) and GNETPR (this
    // fixed-point recomputation from the final rates) differ, recompute the
    // node's branch lookup in variants to show which difference causes it.
    const auto& sched = well_model_.schedule();
    if (well_model_.comm().rank() == 0) {
        const auto& prod = *well_model_.getVFPProperties().getProd();
        for (const auto& [node, gnetpr] : converged_pressures) {
            const auto gpr_it = node_pressures_.find(node);
            const auto br_it = converged.branch_data.find(node);
            const auto up = network.uptree_branch(node);
            if (gpr_it == node_pressures_.end() || br_it == converged.branch_data.end() || !up
                || !up->vfp_table().has_value()) {
                continue;
            }
            const Scalar gpr = gpr_it->second;
            if (std::abs(gpr - gnetpr) < Scalar{0.5} * unit::barsa) {
                continue;
            }
            const int table_id = *up->vfp_table();
            const auto& table = prod.getTable(table_id);
            const auto& b = br_it->second;   // production positive
            const Scalar w = b.water_rate, o = b.oil_rate, g = b.gas_rate;
            const Scalar alq = up->alq_value(VFPProdTable::ALQDimension(table.getALQType(), sched.getUnits())).value_or(0.0);
            const std::string& up_node = up->uptree_node();
            const Scalar up_gnetpr = converged_pressures.count(up_node) ? converged_pressures.at(up_node) : Scalar{0};
            const Scalar up_gpr = node_pressures_.count(up_node) ? node_pressures_.at(up_node) : up_gnetpr;
            const auto fb = group_tree_node_fractions_.find(node);
            const Scalar fwfr = fb != group_tree_node_fractions_.end() ? fb->second[0] : Scalar{0};
            const Scalar fgfr = fb != group_tree_node_fractions_.end() ? fb->second[1] : Scalar{0};
            const Scalar flo = detail::getFlo(table, w, o, g);
            const Scalar p_noclamp = prod.bhp(table_id, -w, -o, -g, up_gnetpr, alq, 0.0, 0.0, false);
            const Scalar p_fallback = prod.bhp(table_id, -w, -o, -g, up_gnetpr, alq, fwfr, fgfr, false);
            const Scalar p_from_gpr = prod.bhp(table_id, -w, -o, -g, up_gpr, alq, fwfr, fgfr, false);
            std::string gt;
            if (const auto f = group_tree_debug_node_flows_.find(node); f != group_tree_debug_node_flows_.end()) {
                gt = fmt::format("; group-tree node flow O/W/G {:.5g}/{:.5g}/{:.5g} m3/d", f->second[0] * unit::day,
                                 f->second[1] * unit::day, f->second[2] * unit::day);
            }
            if (sched.hasGroup(node, reportStepIdx)) {
                for (const auto& wname : sched.getGroup(node, reportStepIdx).wells()) {
                    const auto r = group_tree_debug_well_rates_.find(wname);
                    gt += fmt::format("; {} group-tree O/W/G {}", wname, r == group_tree_debug_well_rates_.end()
                        ? std::string("absent")
                        : fmt::format("{:.5g}/{:.5g}/{:.5g}", r->second[0] * unit::day, r->second[1] * unit::day,
                                      r->second[2] * unit::day));
                    if (well_model_.wellState().has(wname)) {
                        const auto& ws = well_model_.wellState().well(wname);
                        const auto& pu = well_model_.phaseUsage();
                        gt += fmt::format(" final {:.5g}/{:.5g}/{:.5g}",
                            -ws.surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)] * unit::day,
                            -ws.surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)] * unit::day,
                            -ws.surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)] * unit::day);
                    }
                }
            }
            OpmLog::debug(gt.empty() ? std::string{} : "  " + gt.substr(2));
            OpmLog::debug(fmt::format(
                "Node pressure check {} (up {}): GPR {:.3f} vs GNETPR {:.3f} bar (diff {:+.3f}); "
                "unclamped {:.3f}, with fallback fractions {:.3f}, from upstream GPR {:.3f}; "
                "rates O/W/G {:.5g}/{:.5g}/{:.5g} m3/d, FLO {:.5g} (axis {:.5g}..{:.5g}); "
                "upstream GPR {:.3f} / GNETPR {:.3f} (THP axis {:.3f}..{:.3f}); fallback WFR/GFR {}; "
                "table {} ALQ {:.4g}; add gas lift gas {}",
                node, up_node, gpr / unit::barsa, gnetpr / unit::barsa, (gpr - gnetpr) / unit::barsa,
                p_noclamp / unit::barsa, p_fallback / unit::barsa, p_from_gpr / unit::barsa,
                o * unit::day, w * unit::day, g * unit::day, flo * unit::day,
                table.getFloAxis().front() * unit::day, table.getFloAxis().back() * unit::day,
                up_gpr / unit::barsa, up_gnetpr / unit::barsa,
                table.getTHPAxis().front() / unit::barsa, table.getTHPAxis().back() / unit::barsa,
                fb != group_tree_node_fractions_.end() ? fmt::format("{:.4g}/{:.4g}", fwfr, fgfr) : std::string("none"),
                table_id, alq, network.node(node).add_gas_lift_gas() ? "yes" : "no"));
        }
    }
    converged_branchvalues = std::move(converged.branch_data);
    for (const auto& [node, converged_pressure] : converged_pressures) {
        auto it = nodevalues.find(node);
        assert(it != nodevalues.end() );
        it->second.converged_pressure = converged_pressure;
        // Assign node values of group to GPR:WELLNAME
        const auto& sched = well_model_.schedule();
        if (!sched.hasGroup(node, reportStepIdx)) {
            continue;
        }
        const auto& group = sched.getGroup(node, reportStepIdx);
        for (const std::string& wellname : group.wells()) {
            auto it2 = nodevalues.find(wellname);
            assert(it2 != nodevalues.end());
            it2->second.converged_pressure = converged_pressure;
        }
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
initialize(const int report_step)
{
    const auto networks = details::networks(well_model_.schedule(), report_step);
    for (const auto& [domain, network] : networks) {
        // Discard pressures for nodes that are absent from the current network.
        // Retained per-well limits are kept because they can outlive the network.
        auto& node_pressures = this->nodePressures(domain);
        auto& branch_data = this->branchData(domain);
        auto& invalid_nodes = this->invalidNodes(domain);
        auto& last_valid_node_pressures = this->last_valid_domain_node_pressures_[details::domainIndex(domain)];
        auto& last_valid_branch_data = this->last_valid_domain_branch_data_[details::domainIndex(domain)];
        if (!network.get().active()) {
            node_pressures.clear();
            branch_data.clear();
            last_valid_node_pressures.clear();
            last_valid_branch_data.clear();
            invalid_nodes.clear();
        }
        else {
            const auto is_stale = [&network](const auto& node_pressure)
                { return !network.get().has_node(node_pressure.first); };

            std::erase_if(node_pressures, is_stale);
            std::erase_if(branch_data, is_stale);
            std::erase_if(last_valid_node_pressures, is_stale);
            std::erase_if(last_valid_branch_data, is_stale);
            std::erase_if(invalid_nodes, [&network](const auto& node)
                { return !network.get().has_node(node); });
        }
        if (domain == details::NetworkDomain::Production) {
            this->syncLegacyProductionState_();
        }
    }
    this->invalid_nodes_report_step_ = -1;

    // Retained THP limits can outlive network activity, so initialize every well.
    for (auto& well : well_model_.genericWells()) {
        this->initializeWell(*well);
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
initializeWell(WellInterfaceGeneric<Scalar,IndexTraits>& well)
{
    // The network cannot put a well without a VFP table under
    // THP control, so no THP limit is imposed on or retained for such a well,
    // while its rates still contribute to the network flows.
    if (well.wellEcl().vfp_table_number() <= 0) {
        return;
    }

    const auto domain = details::domainForWell(well);
    if (domain.has_value()) {
        const auto it = this->nodePressures(*domain).find(well.wellEcl().groupName());
        if (it != this->nodePressures(*domain).end() && well.isProducer()
            && this->invalidNodes(*domain).count(it->first) == 0) {
            // Carry the converged production-network pressure into the next
            // report step as the producer's starting THP constraint.
            this->imposeWellThpLimit(well, it->second);
        } else if (it == this->nodePressures(*domain).end()
                   || (well.isProducer() && this->invalidNodes(*domain).count(it->first) > 0)) {
            // Reapply a retained limit after network detachment or well reconstruction,
            // or when the current node pressure is only a no-solution placeholder.
            const auto& ws = well_model_.wellState().well(well.indexOfWell());
            if (ws.network_thp_limit.has_value()) {
                well.setDynamicThpLimit(*ws.network_thp_limit);
            }
        }
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
setNodePressureOverrides(std::map<std::string, Scalar> overrides)
{
    this->node_pressure_overrides_ = std::move(overrides);
    this->applyNodePressureOverrides_();
    this->syncLegacyProductionState_();
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
applyNodePressureOverrides_()
{
    auto& pressures = this->nodePressures(details::NetworkDomain::Production);
    for (const auto& [name, pressure] : this->node_pressure_overrides_) {
        if (const auto it = pressures.find(name); it != pressures.end()) {
            it->second = pressure;
        }
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
resetState()
{
    this->node_pressures_ = this->last_valid_node_pressures_;
    this->branch_data_ = this->last_valid_branch_data_;
    this->domain_node_pressures_ = this->last_valid_domain_node_pressures_;
    this->domain_branch_data_ = this->last_valid_domain_branch_data_;
    // The restored pressures predate the overrides in force now.
    this->applyNodePressureOverrides_();
    this->syncLegacyProductionState_();
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
imposeWellThpLimit(WellInterfaceGeneric<Scalar,IndexTraits>& well, const Scalar limit)
{
    well.setDynamicThpLimit(limit);
    well_model_.wellState().well(well.indexOfWell()).network_thp_limit = limit;
}

template <typename Scalar, typename IndexTraits>
typename BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::NetworkPressures
BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
computePressures(const Network::ExtNetwork& network,
                 const VFPProdProperties<Scalar>& vfp_prod_props,
                 const UnitSystem& unit_system,
                 const int reportStepIdx,
                 const Parallel::Communication& comm) const
{
    OPM_TIMEFUNCTION();
    if (!network.active()) {
        return {};
    }

    auto node_gas_consumption = std::map<std::string, Scalar>{};
    if (well_model_.phaseUsage().phaseIsActive(IndexTraits::gasPhaseIdx)) {
        node_gas_consumption = networkGasConsumptionRates<Scalar>(network,
                                                                  well_model_.schedule()[reportStepIdx],
                                                                  well_model_.summaryState());
    }

    NetworkPressureComputation<BlackoilWellModelGeneric<Scalar, IndexTraits>,
                               VFPProdProperties<Scalar>>
        network_pressure_computation(
            well_model_, network, vfp_prod_props, unit_system, reportStepIdx, comm,
            &this->node_pressure_overrides_,
            std::move(node_gas_consumption));

    auto [node_pressures, branch_data] = network_pressure_computation.run();
    return {std::move(node_pressures), std::move(branch_data), network_pressure_computation.invalidNodes()};
}

template <typename Scalar, typename IndexTraits>
typename BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::NetworkPressures
BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
computePressures(const Network::ExtNetwork& network,
                 const VFPInjProperties<Scalar>& vfp_inj_props,
                 const UnitSystem& unit_system,
                 const int reportStepIdx,
                 const Parallel::Communication& comm) const
{
    OPM_TIMEFUNCTION();
    if (!network.active()) {
        return {};
    }

    NetworkPressureComputation<BlackoilWellModelGeneric<Scalar, IndexTraits>,
                               VFPInjProperties<Scalar>>
        network_pressure_computation(
            well_model_, network, vfp_inj_props, unit_system, reportStepIdx, comm);

    auto [node_pressures, branch_data] = network_pressure_computation.run();
    return {std::move(node_pressures), std::move(branch_data), network_pressure_computation.invalidNodes()};
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
operator==(const BlackoilWellModelNetworkGeneric<Scalar,IndexTraits>& rhs) const
{
    return
           this->active_ == rhs.active_
        && this->node_pressures_ == rhs.node_pressures_
        && this->last_valid_node_pressures_ == rhs.last_valid_node_pressures_
        && this->branch_data_ == rhs.branch_data_
        && this->last_valid_branch_data_ == rhs.last_valid_branch_data_
        && this->domain_node_pressures_ == rhs.domain_node_pressures_
        && this->last_valid_domain_node_pressures_ == rhs.last_valid_domain_node_pressures_
        && this->domain_branch_data_ == rhs.domain_branch_data_
        && this->last_valid_domain_branch_data_ == rhs.last_valid_domain_branch_data_;
}

template class BlackoilWellModelNetworkGeneric<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class BlackoilWellModelNetworkGeneric<float, BlackOilDefaultFluidSystemIndices>;
#endif

}
