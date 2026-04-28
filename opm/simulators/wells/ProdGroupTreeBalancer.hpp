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

#ifndef OPM_PROD_GROUP_TREE_BALANCER_HEADER_INCLUDED
#define OPM_PROD_GROUP_TREE_BALANCER_HEADER_INCLUDED

#include <opm/simulators/wells/ProdGroupTreeNode.hpp>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Opm {

class DeferredLogger;
class Group;
class GuideRate;
class Schedule;
class SummaryState;
template<class Scalar> class GroupState;
template<typename Scalar, typename IndexTraits> class BlackoilWellModelGeneric;
template<typename Scalar, typename IndexTraits> class WellState;

} // namespace Opm

namespace Opm::ProdGroupTreeBalancer {

/// Type alias for the tree map.
template<class Scalar>
using Tree = std::map<std::string, ProdGroupTreeNode<Scalar>>;

// ---------------------------------------------------------------------------
// Tree construction
// ---------------------------------------------------------------------------

/// Build the complete production group tree from FIELD down, populating all
/// well and group nodes with current rates, limits, guide rates, and initial
/// control status from existing production_cmode / group_state.production_control().
///
/// \param[in]  wellModel     Well model (for RESV coefficient computation)
/// \param[in]  summaryState  Summary state (for evaluating schedule quantities)
/// \param[in]  reportStep    Current report step index
/// \return     Map from node name to ProdGroupTreeNode
template<class Scalar, typename IndexTraits>
Tree<Scalar> buildTree(const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                       const SummaryState& summaryState,
                       int reportStep);

// ---------------------------------------------------------------------------
// Algorithm core
// ---------------------------------------------------------------------------

/// Compute a post-order (subtree top-nodes first, FIELD last) ordering of
/// "subtree roots" — the top-level nodes whose subtrees still need balancing.
/// This corresponds to Matlab get_sub_tree_ordering().
///
/// \param[in]  tree      The production group tree
/// \param[in]  rootName  Starting node (typically "FIELD")
/// \return     Ordered list of subtree-root names; FIELD last.
template<class Scalar>
std::vector<std::string> getSubTreeOrdering(const Tree<Scalar>& tree,
                                            const std::string& rootName);

/// Recursively compute linearTerm and alphaToNextLimit for the given subtree.
/// Corresponds to Matlab parametrize_tree().
///
/// \param[in,out]  tree      The production group tree
/// \param[in]      nodeName  Root of the subtree to parametrize
template<class Scalar>
void parametrizeTree(Tree<Scalar>& tree, const std::string& nodeName);

/// Find the node in the subtree rooted at nodeName with the smallest
/// alphaToNextLimit.  Corresponds to Matlab getNextLimitNode().
///
/// \param[in]  tree      The production group tree
/// \param[in]  nodeName  Root of the subtree to search
/// \return     Pair of (node_name, alpha_value) for the limiting node.
template<class Scalar>
std::pair<std::string, Scalar>
getNextLimitNode(const Tree<Scalar>& tree, const std::string& nodeName);

/// Step all rates in the GroupControlled subtree rooted at nodeName by alpha
/// along the linearTerm direction.  Corresponds to Matlab stepAlpha().
///
/// \param[in,out]  tree      The production group tree
/// \param[in]      nodeName  Root of the subtree
/// \param[in]      alpha     Step size
template<class Scalar>
void stepAlpha(Tree<Scalar>& tree, const std::string& nodeName, Scalar alpha);

/// Switch the given node to IndividualControlled, update linear terms upward,
/// and mark ancestor groups as NoGroupChildren if they have no remaining
/// GroupControlled children.  Corresponds to Matlab updateNode().
///
/// \param[in,out]  tree      The production group tree
/// \param[in]      nodeName  The node that has hit its individual limit
template<class Scalar>
void updateNode(Tree<Scalar>& tree, const std::string& nodeName);

/// Recursively cap individual rates at their limits, propagate sums upward,
/// and return whether any rate was capped (used inside makeFeasible).
///
/// \param[in,out]  tree      The production group tree
/// \param[in]      nodeName  Root of subtree
/// \return         true if any rate was modified
template<class Scalar>
bool capIndividualAndSum(Tree<Scalar>& tree, const std::string& nodeName);

/// Propagate a group target down through GroupControlled children using guide
/// rates.  Corresponds to part of Matlab make_feasible().
///
/// \param[in,out]  tree      The production group tree
/// \param[in]      nodeName  Root of the subtree
template<class Scalar>
void setAndUpdateTargets(Tree<Scalar>& tree, const std::string& nodeName);

/// Iterate capIndividualAndSum + setAndUpdateTargets until convergence or
/// \p maxIter is reached.  Sets initial tiny rates if current rates are zero.
/// Corresponds to Matlab make_feasible().
///
/// \param[in,out]  tree      The production group tree
/// \param[in]      topName   Name of the top-level subtree node
/// \param[in]      tol       Convergence tolerance (relative)
/// \param[in]      maxIter   Maximum number of iterations
/// \return         true if converged within maxIter
template<class Scalar>
bool makeFeasible(Tree<Scalar>& tree,
                  const std::string& topName,
                  Scalar tol,
                  int maxIter);

/// Assign final group targets recursively from topName downward.
/// Corresponds to Matlab set_targets().
///
/// \param[in,out]  tree      The production group tree
/// \param[in]      topName   Root of the subtree
template<class Scalar>
void setFinalTargets(Tree<Scalar>& tree, const std::string& topName);

/// Outer balancing loop.  Iterates over subtrees in getSubTreeOrdering() order,
/// and within each subtree calls parametrizeTree / getNextLimitNode / stepAlpha
/// / updateNode until no Undetermined nodes remain.
/// Corresponds to Matlab balance_group_tree_new().
///
/// \param[in,out]  tree      The production group tree
/// \param[in]      tol       Convergence tolerance
/// \param[in]      maxIter   Maximum iterations per subtree
/// \return         true if all subtrees converged
template<class Scalar>
bool balanceGroupTree(Tree<Scalar>& tree, Scalar tol, int maxIter);

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

/// Verify that every node is either:
///   (a) at its activeIndividualCtrl limit (within tol), or
///   (b) at its group target (within tol), or
///   (c) in NoGroupChildren / Undetermined state with no GroupControlled children.
/// Logs a warning for each violation and returns false if any are found.
///
/// \param[in]  tree      The production group tree
/// \param[in]  topName   Root of the subtree to check
/// \param[in]  tol       Relative tolerance for limit checks
/// \param[in]  logger    Deferred logger for warnings
/// \return     true if the tree is valid
template<class Scalar>
bool checkTreeValidity(const Tree<Scalar>& tree,
                       const std::string& topName,
                       Scalar tol,
                       DeferredLogger& logger);

// ---------------------------------------------------------------------------
// State write-back
// ---------------------------------------------------------------------------

/// Write the balanced tree state back to wellState and groupState:
///   - For each production well: set ws.surface_rates and ws.production_cmode
///   - For each group: set group_state.production_rates and group_state.production_control
/// Does NOT touch injection wells or groups.
///
/// \param[in]    tree        The balanced production group tree
/// \param[inout] wellModel   Well model (writes to its internal wellState and groupState)
template<class Scalar, typename IndexTraits>
void applyTreeToState(const Tree<Scalar>& tree,
                      BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel);

// ---------------------------------------------------------------------------
// Orchestrator
// ---------------------------------------------------------------------------

/// Top-level entry point: build tree, balance it, validate, and apply.
/// Logs timing and iteration counts; logs a warning on convergence failure.
///
/// \param[in]    wellModel     Well model (schedule, well/group state, guide rate, RESV coefficients)
/// \param[in]    summaryState  Summary state
/// \param[in]    reportStep    Current report step index
/// \param[in]    tol           Convergence tolerance (default: 1e-4)
/// \param[in]    maxIter       Maximum iterations per subtree (default: 100)
/// \param[in]    logger        Deferred logger
/// \return       true if the result passed checkTreeValidity
template<class Scalar, typename IndexTraits>
bool runGroupTreeBalancer(BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                          const SummaryState& summaryState,
                          int reportStep,
                          Scalar tol,
                          int maxIter,
                          DeferredLogger& logger);

} // namespace Opm::ProdGroupTreeBalancer

#endif // OPM_PROD_GROUP_TREE_BALANCER_HEADER_INCLUDED
