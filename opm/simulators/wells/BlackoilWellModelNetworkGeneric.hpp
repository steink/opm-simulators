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

#ifndef OPM_BLACKOILWELLMODEL_NETWORK_GENERIC_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_NETWORK_GENERIC_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/ScheduleTypes.hpp>

#include <opm/output/data/Groups.hpp>

#include <opm/simulators/flow/NewtonIterationContext.hpp>
#include <opm/simulators/wells/NetworkGroupTreeSystem.hpp>
#include <opm/simulators/wells/ProdGroupTreeBalancer.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <algorithm>
#include <array>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace Opm {
    class DeferredLogger;
    class Schedule;
    class UnitSystem;
    template<class Scalar, class IndexTraits> class BlackoilWellModelGeneric;
    template<typename Scalar, typename IndexTraits> class WellInterfaceGeneric;
    template<typename Scalar> class VFPInjProperties;
    template<typename Scalar> class VFPProdProperties;
}

namespace Opm {

namespace details {

    enum class NetworkDomain : std::size_t {
        Production = 0,
        InjectionGas,
        InjectionWater,
        Count
    };

    constexpr std::size_t domainIndex(const NetworkDomain domain)
    {
        return static_cast<std::size_t>(domain);
    }

    struct NetworkDescriptor {
        NetworkDomain domain;
        std::reference_wrapper<const Network::ExtNetwork> network;
    };

    /// The network domain corresponding to a well's type
    /// (producer, or injector and water/gas injection phase), or nullopt.
    template<class Well>
    std::optional<NetworkDomain> domainForWell(const Well& well)
    {
        if (well.isProducer()) {
            return NetworkDomain::Production;
        }
        if (well.isInjector()) {
            if (well.wellEcl().injectorType() == InjectorType::GAS) {
                return NetworkDomain::InjectionGas;
            }
            if (well.wellEcl().injectorType() == InjectorType::WATER) {
                return NetworkDomain::InjectionWater;
            }
        }
        return std::nullopt;
    }

    /// Helper to check if any network (production, gas injection, water injection) is active at a given time step.
    bool anyNetworkActive(const Schedule& schedule, const int timeStepIdx);

    /// Helper to get all active networks (production, gas injection, water injection) at a given time step.
    std::vector<NetworkDescriptor>
    activeNetworks(const Schedule& schedule, const int timeStepIdx);

    /// Helper to get all networks (production, gas injection, water injection) at a given time step,
    /// whether active or not.
    std::vector<NetworkDescriptor>
    networks(const Schedule& schedule, const int timeStepIdx);
    /// Which network domains a network update covers.
    enum class DomainScope { All, Production, Injection };

    /// True if \p domain is covered by \p scope.
    inline bool inScope(const DomainScope scope, const NetworkDomain domain)
    {
        switch (scope) {
        case DomainScope::Production: return domain == NetworkDomain::Production;
        case DomainScope::Injection:  return domain != NetworkDomain::Production;
        default:                      return true;
        }
    }

} // namespace details


/// Class for handling the blackoil well network model.
template<typename Scalar, typename IndexTraits>
class BlackoilWellModelNetworkGeneric
{
public:
    BlackoilWellModelNetworkGeneric(BlackoilWellModelGeneric<Scalar,IndexTraits>& well_model);

    virtual ~BlackoilWellModelNetworkGeneric() = default;

    /// return true if network is active (at least one network well in prediction mode)
    bool active() const
    { return active_; }

    const std::map<std::string, Scalar>&
    nodePressures() const { return node_pressures_; }

    // do not use, only needed for serialization testing
    void setNodePressures(const std::map<std::string, Scalar>& values)
    { node_pressures_ = values; }

    void setFromRestart(const std::optional<std::map<std::string, double>>& restart_pressures);

    //! \brief Pressures that a reservoir coupling master imposes on nodes of the
    //!        production network, by node (group) name.
    //!
    //! A reservoir coupling slave takes the pressure of a node that is a master
    //! group's leaf in the master's network from the master.  The pressure
    //! computation uses it as the fixed pressure of the node, instead of the one
    //! in the deck, and the nodes below it follow from it with the pressure drops
    //! in their branches.  So the THP limits of the slave's wells are the master's pressure
    //! plus those drops whenever they are set, and not only after the master's
    //! message is received.  The node itself is not damped by updatePressures().
    //! Only production network nodes with a fixed pressure in the deck are
    //! overridden.
    //! The cached pressures of the overridden nodes are set to the imposed values at once, so
    //! that what is computed from them before the next pressure update, such as the common THP
    //! of an auto-choke group below such a node, uses the current imposed pressure.
    void setNodePressureOverrides(std::map<std::string, Scalar> overrides);

    //! \brief Initialize wells according to network configuration.
    void initialize(const int report_step);

    //! \brief Initialize a single well according to network configuration.
    void initializeWell(WellInterfaceGeneric<Scalar,IndexTraits>& well);

    /// Checks if network is active (at least one network well on prediction).
    void updateActiveState(const int report_step);

    /// Checks if there are reasons to perform a pre-step network re-balance.
    /// (Currently, the only reasons are network well status changes.)
    /// (TODO: Consider if adding network change events would be helpful.)
    bool needPreStepRebalance(const int report_step) const;

    /// Checks if we shall perform a network re-balance.
    /// This is typically controlled by the NETBALAN keyword.
    bool shouldBalance(const int reportStepIndex) const;
    /// Checks if we will perform a network re-balance on the next Newton iteration.
    bool willBalanceOnNextIteration(const int reportStepIndex) const;

    /// Recompute the node pressures from the current leaf rates and move the applied
    /// pressures towards them, damped by damping_factor and capped at
    /// update_upper_bound. Returns the largest change over the nodes.
    ///
    /// \p scope restricts the update -- the evaluation, the pressure step and the
    /// THP-limit write-back -- to the production or the injection networks; the
    /// other domain's node pressures are left as they are.
    Scalar updatePressures(const int reportStepIdx,
                           const Scalar damping_factor,
                           const Scalar update_upper_bound,
                           const details::DomainScope scope = details::DomainScope::All);

    void assignNodeAndBranchValues(std::map<std::string, data::NodeData>& nodevalues,
                                   std::map<std::string, data::BranchData>& branchvalues,
                                   std::map<std::string, data::BranchData>& converged_branchvalues,
                                   const int reportStepIdx) const;

    void commitState()
    {
        this->last_valid_node_pressures_ = this->node_pressures_;
        this->last_valid_branch_data_ = this->branch_data_;
        this->last_valid_domain_node_pressures_ = this->domain_node_pressures_;
        this->last_valid_domain_branch_data_ = this->domain_branch_data_;
    }

    /// Restores the state saved by commitState().  The node pressure
    /// overrides are kept, as they are the current boundary condition, so
    /// they are applied to the restored pressures again.
    void resetState();

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(node_pressures_);
        serializer(last_valid_node_pressures_);
        serializer(branch_data_);
        serializer(last_valid_branch_data_);
        serializer(domain_node_pressures_);
        serializer(last_valid_domain_node_pressures_);
        serializer(domain_branch_data_);
        serializer(last_valid_domain_branch_data_);
        serializer(well_group_thp_calc_);
    }

    bool operator==(const BlackoilWellModelNetworkGeneric<Scalar,IndexTraits>& rhs) const;

    /// Solve the production network's group-tree-driven system (ProdGroupTreeBalancer
    /// + NetworkSolve::GroupTreeSystem) instead of relaxing the node pressures against
    /// the wells. Off by default; --network-solver=group-tree turns it on. The
    /// fixed-point answer is still computed first, and kept for a network the
    /// group-tree solve cannot model. Production only.
    void useGroupTreeSolver(const bool on) { group_tree_solver_ = on; }
    bool usesGroupTreeSolver() const { return group_tree_solver_; }

    /// The balanced group tree the next group-tree network solve is built for.
    /// Set by the caller when usesGroupTreeSolver() and a production network
    /// are both active; otherwise unused.
    void setBalancedGroupTree(ProdGroupTreeBalancer::Tree<Scalar> tree)
    {
        balanced_group_tree_ = std::move(tree);
    }

    /// The roots of \p network that head a tree of their own: each distinct
    /// fixed-pressure node without a branch above it. ExtNetwork::roots() also
    /// lists a fixed-pressure node that feeds a node further up (it is the
    /// uptree node of a branch with a fixed pressure), and lists a root once
    /// per branch below it; the group-tree solve treats the former as a
    /// fixed-pressure node inside its parent's tree, so it is solved once, with
    /// the group constraints spanning both.
    static std::vector<std::reference_wrapper<const Network::Node>>
    independentRoots(const Network::ExtNetwork& network)
    {
        std::vector<std::reference_wrapper<const Network::Node>> out;
        for (const auto& root : network.roots()) {
            const auto& name = root.get().name();
            if (network.uptree_branch(name).has_value()) {
                continue;
            }
            if (std::none_of(out.begin(), out.end(),
                             [&name](const auto& r) { return r.get().name() == name; })) {
                out.push_back(root);
            }
        }
        return out;
    }

    /// One B2/B3 step of the group-tree workflow: solve the network for
    /// \p balancedTree, then change the open set by what the solution says.
    ///  - With \p offerReopenCandidates (the first pass of B), stopped wells
    ///    with a usable trial IPR enter the solve as reopen candidates; those
    ///    that flow on the real tubing curve are reopened (openWell()) and their
    ///    well state seeded with the network's operating point.
    ///  - Wells that cannot keep their branch flowing on its real table are
    ///    stopped (stopWell(), the dynamic status).
    /// Returns true if the open set changed, so the caller rebalances and
    /// solves again; false if the solve did not converge or nothing changed.
    bool updateGroupTreeOpenSet(const Network::ExtNetwork& network,
                                const int reportStepIdx,
                                const Network::Node& root,
                                const ProdGroupTreeBalancer::Tree<Scalar>& balancedTree,
                                const bool offerReopenCandidates,
                                DeferredLogger& deferred_logger);

    /// Assemble the network Jacobian from the VFP table derivatives instead of
    /// differencing the residual.
    void useAnalyticJacobian(const bool on) { analytic_jacobian_ = on; }
    void useBranchFlattening(const bool on) { branch_flattening_ = on; }
    void useNetworkLineSearch(const bool on) { network_line_search_ = on; }
    /// Per local well, the hydrostatic correction its tubing table needs;
    /// computed on the typed side, where the well's density lives.
    void setWellVfpDp(const std::string& well, const Scalar dp) { well_vfp_dp_[well] = dp; }

protected:
    /// See setNodePressureOverrides().
    std::map<std::string, Scalar> node_pressure_overrides_;

    //! \brief Apply a network THP limit and retain it across well
    //! reconstruction or network detachment.
    void imposeWellThpLimit(WellInterfaceGeneric<Scalar,IndexTraits>& well,
                            const Scalar limit);

    /// Result of one network pressure evaluation for one network (domain).
    struct NetworkPressures
    {
        std::map<std::string, Scalar> node_pressures;
        std::map<std::string, data::BranchData> branch_data;
        // Nodes (and their descendants) whose VFP lookup has no solution; their
        // node_pressures entries are placeholders and must not be used.
        std::set<std::string> invalid_nodes;
    };

    NetworkPressures
    computePressures(const Network::ExtNetwork& network,
                     const VFPProdProperties<Scalar>& vfp_prod_props,
                     const UnitSystem& unit_system,
                     const int reportStepIdx,
                     const Parallel::Communication& comm) const;

    NetworkPressures
    computePressures(const Network::ExtNetwork& network,
                     const VFPInjProperties<Scalar>& vfp_inj_props,
                     const UnitSystem& unit_system,
                     const int reportStepIdx,
                     const Parallel::Communication& comm) const;

    void updateActiveStateImpl(const Network::ExtNetwork& network);

    static constexpr details::NetworkDomain productionNetworkDomain()
    {
        return details::NetworkDomain::Production;
    }

    const std::map<std::string, Scalar>& nodePressures(const details::NetworkDomain domain) const
    {
        return domain_node_pressures_[details::domainIndex(domain)];
    }

    std::map<std::string, Scalar>& nodePressures(const details::NetworkDomain domain)
    {
        return domain_node_pressures_[details::domainIndex(domain)];
    }

    const std::map<std::string, data::BranchData>& branchData(const details::NetworkDomain domain) const
    {
        return domain_branch_data_[details::domainIndex(domain)];
    }

    std::map<std::string, data::BranchData>& branchData(const details::NetworkDomain domain)
    {
        return domain_branch_data_[details::domainIndex(domain)];
    }

    const std::set<std::string>& invalidNodes(const details::NetworkDomain domain) const
    {
        return domain_invalid_nodes_[details::domainIndex(domain)];
    }

    std::set<std::string>& invalidNodes(const details::NetworkDomain domain)
    {
        return domain_invalid_nodes_[details::domainIndex(domain)];
    }

    /// Sets the cached production node pressures of the overridden nodes to
    /// the override values.  Does not sync the legacy state.
    void applyNodePressureOverrides_();

    void syncLegacyProductionState_()
    {
        this->node_pressures_ = this->nodePressures(productionNetworkDomain());
        this->branch_data_ = this->branchData(productionNetworkDomain());
    }

    void syncProductionDomainState_()
    {
        this->nodePressures(productionNetworkDomain()) = this->node_pressures_;
        this->branchData(productionNetworkDomain()) = this->branch_data_;
    }

    bool active_{false};
    BlackoilWellModelGeneric<Scalar,IndexTraits>& well_model_;

    // Network pressures for output and initialization
    std::map<std::string, Scalar> node_pressures_;
    // Network branch pressure drops and flow rates for output (outlet branch for production network, inlet branch for injection network)
    std::map<std::string, data::BranchData> branch_data_;
    // Domain-scoped pressure state to avoid collisions between production and injection networks.
    std::array<std::map<std::string, Scalar>, details::domainIndex(details::NetworkDomain::Count)> domain_node_pressures_;
    std::array<std::map<std::string, data::BranchData>, details::domainIndex(details::NetworkDomain::Count)> domain_branch_data_;
    /// Node pressures from the group-tree solve for one production network
    /// root, or nullopt if it cannot be done (logged at debug level) -- the
    /// caller then keeps the fixed-point result for that tree.
    std::optional<std::map<std::string, Scalar>>
    groupTreeProductionNodePressures(const Network::ExtNetwork& network,
                                     const int reportStepIdx,
                                     const Network::Node& root,
                                     const ProdGroupTreeBalancer::Tree<Scalar>& balancedTree) const;

    /// The per-well data NetworkSolve::GroupTreeSystem::WellNetworkData needs
    /// (vfp table, ipr_a/b, network-branch efficiency, ...), gathered MPI-safely:
    /// replicated static facts from the schedule, dynamic facts summed once
    /// from whichever rank owns the well, so every rank ends up with the same
    /// answer. Keyed by well name; the caller resolves each well's node once
    /// the topology is built. A well left out of the returned map is not an
    /// open (or stopped) producer in prediction mode.
    std::unordered_map<std::string, typename NetworkSolve::GroupTreeSystem<Scalar>::WellNetworkData>
    gatherWellNetworkDataForGroupTree(const int reportStepIdx) const;

    /// Everything through NetworkSolve::solve() for one production network
    /// root, keeping the system and the result (groupTreeProductionNodePressures()
    /// only needs the pressures; updateGroupTreeOpenSet() needs both). nullopt
    /// on anything the group-tree solve gives up on, logged at debug level.
    struct GroupTreeSolve {
        NetworkSolve::GroupTreeSystem<Scalar> system;
        std::vector<std::string> order;    // node names by index; order[0] == root.name()
        NetworkSolve::Result<Scalar> result;
        /// Thp wells under a node that cannot flow on its real branch table
        /// with them (GroupTreeSystem::nodeMargin() negative): stopped by
        /// updateGroupTreeOpenSet(), not reopened.
        std::vector<std::string> network_cliff_wells;
    };
    /// With \p offerReopenCandidates, stopped wells with a usable trial IPR
    /// enter the system as reopen candidates (see
    /// GroupTreeSystem::addReopenCandidate()).
    std::optional<GroupTreeSolve>
    solveGroupTree(const Network::ExtNetwork& network,
                   const int reportStepIdx,
                   const Network::Node& root,
                   const ProdGroupTreeBalancer::Tree<Scalar>& balancedTree,
                   const bool offerReopenCandidates = false) const;

    bool group_tree_solver_ = false;
    std::optional<ProdGroupTreeBalancer::Tree<Scalar>> balanced_group_tree_;
    bool analytic_jacobian_ = false;
    bool branch_flattening_ = true;   // group-tree solver: see GroupTreeSystem::setBranchFlattening()
    bool network_line_search_ = false;   // --group-tree-network-step=linesearch
    std::map<std::string, Scalar> well_vfp_dp_;
    /// Each production node's (WFR, GFR) at its last converged group-tree
    /// solve with a flow at or above its table's first FLO value: what its
    /// table is looked up with below that value at the next solve. Updated by
    /// solveGroupTree() (hence mutable).
    mutable std::map<std::string, std::array<Scalar, 2>> group_tree_node_fractions_;



    // Nodes without a valid VFP solution in the last evaluation (per domain); not serialized,
    // recomputed on every updatePressures().
    std::array<std::set<std::string>, details::domainIndex(details::NetworkDomain::Count)> domain_invalid_nodes_;
    int invalid_nodes_report_step_{-1};
    // Valid network pressures for output and initialization for safe restart after failed iterations
    std::map<std::string, Scalar> last_valid_node_pressures_;
    // Valid network branch pressure drops and flow rates for output (outlet branch for production network, inlet branch for injection network) for safe restart after failed iterations
    std::map<std::string, data::BranchData> last_valid_branch_data_;
    std::array<std::map<std::string, Scalar>, details::domainIndex(details::NetworkDomain::Count)> last_valid_domain_node_pressures_;
    std::array<std::map<std::string, data::BranchData>, details::domainIndex(details::NetworkDomain::Count)> last_valid_domain_branch_data_;
    // Cached auto-choke THP values used when balancing network groups.
    std::map<std::string, Scalar> well_group_thp_calc_;
};

} // namespace Opm

#endif
