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

/// Get local tree descendants categorized by availability for group control.
/// Corresponds to Matlab getLocalTreeDescendants().
///
/// \param[in]  tree      The production group tree
/// \param[in]  nodeName  Root node
/// \return     Tuple of (c, c_fixed, c_trans) where:
///             c = nodes with guide-rate available for group control
///             c_fixed = nodes not available for group control
///             c_trans = transparent nodes (no guide-rate) between c and nodeName
template<class Scalar>
std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>>
getLocalTreeDescendants(const Tree<Scalar>& tree, const std::string& nodeName);

/// Reset rate_sums and guide_rate_sums for nodes from c up to origin.
/// Corresponds to Matlab reset_rates_and_guide_rate_sums().
///
/// \param[in,out]  tree    The production group tree
/// \param[in]      c       Nodes with guide rates
/// \param[in]      c_fixed Fixed nodes
/// \param[in]      c_trans Transparent nodes
/// \param[in]      origin  Origin node
/// \param[in]      mode    Target control mode
template<class Scalar>
void resetRatesAndGuideRateSums(Tree<Scalar>& tree,
                                 const std::vector<std::string>& c,
                                 const std::vector<std::string>& c_fixed,
                                 const std::vector<std::string>& c_trans,
                                 const std::string& origin,
                                 Well::ProducerCMode mode);

/// Compute limit-to-guide ratios for sorting children.
/// Corresponds to Matlab get_ratios_for_sorting().
///
/// \param[in]  tree  The production group tree
/// \param[in]  c     List of child node names
/// \return     Pair of (ratios, mode_indices) for each child
template<class Scalar>
std::pair<std::vector<Scalar>, std::vector<int>>
getRatiosForSorting(const Tree<Scalar>& tree, const std::vector<std::string>& c);

/// Check if there is a free path from node to node_control (all transparent).
/// Corresponds to Matlab hasFreePath().
///
/// \param[in]  tree          The production group tree
/// \param[in]  nodeName      Starting node
/// \param[in]  nodeControl   Controlling node
/// \return     true if path is free
template<class Scalar>
bool hasFreePath(const Tree<Scalar>& tree,
                 const std::string& nodeName,
                 const std::string& nodeControl);

/// Update transparent groups that reach their limits before the next child.
/// Corresponds to Matlab update_transparent_groups().
///
/// \param[in,out]  tree        The production group tree
/// \param[in]      nodeName    Origin node
/// \param[in,out]  c_trans     Transparent nodes (modified: updated nodes removed)
/// \param[in]      nextRatio   Next child's ratio
/// \param[in]      qm          Target rate for origin
/// \param[in]      mode        Control mode
/// \param[in]      tol         Tolerance
/// \return         List of updated transparent nodes
template<class Scalar>
std::vector<std::string>
updateTransparentGroups(Tree<Scalar>& tree,
                        const std::string& nodeName,
                        std::vector<std::string>& c_trans,
                        Scalar nextRatio,
                        Scalar qm,
                        Well::ProducerCMode mode,
                        Scalar tol);

/// Increment parent rate_sums from node up to origin.
/// Corresponds to Matlab increment_parent_rate_sums().
///
/// \param[in,out]  tree      The production group tree
/// \param[in]      nodeName  Starting node
/// \param[in]      origin    Origin node
/// \param[in]      rates     Rates to add (if not provided, uses node's rates)
template<class Scalar>
void incrementParentRateSums(Tree<Scalar>& tree,
                             const std::string& nodeName,
                             const std::string& origin,
                             const std::optional<std::array<Scalar, 3>>& rates = std::nullopt);

/// Decrement parent guide_rate_sums from node up to origin.
/// Corresponds to Matlab decrement_parent_guide_rate_sums().
///
/// \param[in,out]  tree            The production group tree
/// \param[in]      nodeName        Starting node
/// \param[in]      origin          Origin node
/// \param[in]      guideRateSums   Guide rate sums to subtract
template<class Scalar>
void decrementParentGuideRateSums(Tree<Scalar>& tree,
                                   const std::string& nodeName,
                                   const std::string& origin,
                                   const std::array<Scalar, 3>& guideRateSums);

/// Distribute rates to fallback children (tiny fractions in non-preferred mode).
/// Corresponds to Matlab distribute_fallback_rates().
///
/// \param[in,out]  tree      The production group tree
/// \param[in]      nodeName  Node distributing rates
/// \param[in]      c         All children
/// \param[in]      mode      Control mode
/// \param[in]      tol       Tolerance
template<class Scalar>
void distributeFallbackRates(Tree<Scalar>& tree,
                             const std::string& nodeName,
                             const std::vector<std::string>& c,
                             Well::ProducerCMode mode,
                             Scalar tol);

/// Recursively balance the group tree rooted at nodeName with target mode and rate.
/// This is the main balancing function that corresponds to Matlab balance_group_tree().
///
/// \param[in,out]  tree        The production group tree
/// \param[in]      nodeName    Root of the subtree to balance
/// \param[in]      targetMode  Target control mode
/// \param[in]      targetRate  Target rate for the control mode
/// \param[in]      tol         Convergence tolerance
template<class Scalar>
void balanceGroupTree(Tree<Scalar>& tree,
                      const std::string& nodeName,
                      Well::ProducerCMode targetMode,
                      Scalar targetRate,
                      Scalar tol);

/// Recursively cap individual rates at their limits, propagate sums upward,
/// and categorize nodes.  Corresponds to Matlab cap_individual_and_sum().
///
/// \param[in,out]  tree      The production group tree
/// \param[in]      nodeName  Root of subtree
/// \param[in]      tol       Tolerance for limit comparisons
/// \return         true if any category switched
template<class Scalar>
bool capIndividualAndSum(Tree<Scalar>& tree, const std::string& nodeName, Scalar tol);

/// Update node category based on rates and limits.
/// Corresponds to Matlab updateCategory().
///
/// \param[in,out]  tree                           The production group tree
/// \param[in]      nodeName                       Node to update
/// \param[in]      anyGroupControlledChildren     Whether node has any GRUP children
/// \param[in]      targetMode                     Target control mode (or CMODE_UNDEFINED)
/// \param[in]      tol                           Tolerance
template<class Scalar>
void updateCategory(Tree<Scalar>& tree,
                   const std::string& nodeName,
                   bool anyGroupControlledChildren,
                   Well::ProducerCMode targetMode,
                   Scalar tol);

/// Assign final group targets recursively from topName downward.
/// Corresponds to Matlab set_targets().
///
/// \param[in,out]  tree      The production group tree
/// \param[in]      topName   Root of the subtree
template<class Scalar>
void setTargets(Tree<Scalar>& tree, const std::string& topName);

/// Outer balancing loop for the entire tree.  Iterates over subtrees in
/// getSubTreeOrdering() order, balancing each one.
/// Corresponds to Matlab sortalgo2().
///
/// \param[in,out]  tree      The production group tree
/// \param[in]      tol       Convergence tolerance
/// \return         true if all subtrees balanced successfully
template<class Scalar>
bool runBalancingAlgorithm(Tree<Scalar>& tree, Scalar tol);

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
