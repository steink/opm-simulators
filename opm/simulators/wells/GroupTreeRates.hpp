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

#ifndef OPM_GROUP_TREE_RATES_HPP
#define OPM_GROUP_TREE_RATES_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace Opm {

/// Simplified node in the group/well hierarchy used for rate distribution.
///
/// Each node represents either a group or a well.  Wells are leaf nodes
/// (children is empty).  The root node typically represents the FIELD group.
///
/// The status field encodes the control state:
///   -1 : group-controlled (rate determined by parent allocation)
///    0 : fully determined by children (all children are fixed)
///    1 : individually rate-limited
template<class Scalar>
struct GroupTreeNode {
    int index{-1};                ///< Unique node index in the flat vector
    int parent{-1};               ///< Parent node index (-1 for root)
    Scalar rate{0};               ///< Current allocated rate
    Scalar limit{0};              ///< Rate limit for this node
    int status{-1};               ///< Control status (-1, 0, or 1)
    Scalar guideRate{0};          ///< Guide rate used for allocation
    std::string name;             ///< Group or well name
    std::vector<int> children;    ///< Indices of child nodes
};

/// Distribute rates through a group tree respecting individual node limits.
///
/// Starting from the root node (at @p rootIndex), rates are distributed
/// top-down to children according to their guide-rate fractions.  Whenever
/// a child's allocated rate exceeds its limit, that child is fixed at its
/// limit and the distribution is recomputed for the remaining
/// group-controlled children.  The process repeats until no limit
/// violations remain.
///
/// @param tree      Flat vector of tree nodes (modified in-place).
/// @param rootIndex Index of the subtree root in @p tree.
/// @param maxIter   Maximum number of outer iterations (safety limit).
/// @return          Number of outer iterations used.
template<class Scalar>
int distributeGroupTreeRates(std::vector<GroupTreeNode<Scalar>>& tree,
                             int rootIndex = 0,
                             int maxIter = 1000);

/// Distribute the rate of node @p nodeIndex to its children.
///
/// Children that are already individually limited (status != -1) keep
/// their current rate.  The remaining (available) rate is split among
/// group-controlled children proportionally to their guide rates.
/// If all children are fixed, the parent's status becomes 0 and its
/// rate is set to the sum of the children's rates.
///
/// @param tree      Flat vector of tree nodes (modified in-place).
/// @param nodeIndex Index of the node whose children are updated.
template<class Scalar>
void setSubRates(std::vector<GroupTreeNode<Scalar>>& tree,
                 int nodeIndex);

/// Find the child (in the subtree rooted at @p nodeIndex) whose rate
/// most exceeds its limit.
///
/// @param tree      Flat vector of tree nodes (read-only access).
/// @param nodeIndex Root of the subtree to search.
/// @param[out] worstExcess  The largest (rate - limit) value found.
/// @return          Index of the worst offending node, or -1 if no
///                  violation exists.
template<class Scalar>
int findWorstOffendingChild(const std::vector<GroupTreeNode<Scalar>>& tree,
                            int nodeIndex,
                            Scalar& worstExcess);

/// After fixing a child node, propagate upward: if all siblings are
/// also fixed, the parent becomes fully determined (status 0) with
/// rate equal to the sum of its children.
///
/// @param tree      Flat vector of tree nodes (modified in-place).
/// @param nodeIndex Index of the node whose parent is updated.
template<class Scalar>
void updateParentStatus(std::vector<GroupTreeNode<Scalar>>& tree,
                        int nodeIndex);

} // namespace Opm

#endif // OPM_GROUP_TREE_RATES_HPP
