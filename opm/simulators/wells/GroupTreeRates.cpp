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
#include <opm/simulators/wells/GroupTreeRates.hpp>

#include <cassert>
#include <numeric>

namespace Opm {

template<class Scalar>
int distributeGroupTreeRates(std::vector<GroupTreeNode<Scalar>>& tree,
                             int rootIndex,
                             int maxIter)
{
    if (rootIndex == 0) {
        // Initialise: all nodes group-controlled except root
        for (auto& node : tree) {
            node.status = -1;
        }
        tree[rootIndex].status = 1;
        tree[rootIndex].rate = tree[rootIndex].limit;
    }

    int iter = 0;
    for (; iter < maxIter; ++iter) {
        setSubRates(tree, rootIndex);

        Scalar worstExcess{0};
        const int ix = findWorstOffendingChild(tree, rootIndex, worstExcess);
        if (ix < 0) {
            break; // converged
        }

        // Fix the worst-offending node at its limit
        tree[ix].rate = tree[ix].limit;
        tree[ix].status = 1;

        // If it has children, recursively solve its subtree
        if (!tree[ix].children.empty()) {
            distributeGroupTreeRates(tree, ix, maxIter);
        }

        updateParentStatus(tree, ix);
    }
    return iter;
}

template<class Scalar>
void setSubRates(std::vector<GroupTreeNode<Scalar>>& tree,
                 int nodeIndex)
{
    const auto& children = tree[nodeIndex].children;
    if (children.empty()) {
        return; // leaf (well) node
    }

    // Separate children into fixed (status != -1) and group-controlled
    Scalar fixedRate{0};
    Scalar guideSum{0};
    bool allFixed = true;

    for (const int ci : children) {
        if (tree[ci].status != -1) {
            fixedRate += tree[ci].rate;
        } else {
            guideSum += tree[ci].guideRate;
            allFixed = false;
        }
    }

    if (allFixed) {
        tree[nodeIndex].status = 0;
        tree[nodeIndex].rate = fixedRate;
        return;
    }

    const Scalar availRate = tree[nodeIndex].rate - fixedRate;
    assert(availRate >= Scalar{0});
    assert(availRate > Scalar{0} || guideSum == Scalar{0});

    for (const int ci : children) {
        if (tree[ci].status == -1) {
            tree[ci].rate = (guideSum > Scalar{0})
                ? availRate * tree[ci].guideRate / guideSum
                : Scalar{0};
            setSubRates(tree, ci);
        }
    }
}

template<class Scalar>
int findWorstOffendingChild(const std::vector<GroupTreeNode<Scalar>>& tree,
                            int nodeIndex,
                            Scalar& worstExcess)
{
    const auto& children = tree[nodeIndex].children;
    int worstIndex = -1;

    for (const int ci : children) {
        const Scalar excess = tree[ci].rate - tree[ci].limit;
        if (excess > worstExcess) {
            worstExcess = excess;
            worstIndex = ci;
        }
    }

    // Recurse into children that have subtrees
    for (const int ci : children) {
        if (!tree[ci].children.empty()) {
            Scalar subtreeExcess = worstExcess;
            const int childWorst = findWorstOffendingChild(tree, ci, subtreeExcess);
            if (subtreeExcess > worstExcess) {
                worstExcess = subtreeExcess;
                worstIndex = childWorst;
            }
        }
    }

    if (worstExcess > Scalar{0} && worstIndex >= 0) {
        assert(tree[worstIndex].status == -1);
    } else {
        worstIndex = -1;
    }
    return worstIndex;
}

template<class Scalar>
void updateParentStatus(std::vector<GroupTreeNode<Scalar>>& tree,
                        int nodeIndex)
{
    const int pix = tree[nodeIndex].parent;
    if (pix < 0) {
        return;
    }

    const auto& siblings = tree[pix].children;
    bool allFixed = true;
    Scalar totalRate{0};

    for (const int ci : siblings) {
        if (tree[ci].status == -1) {
            allFixed = false;
            break;
        }
        totalRate += tree[ci].rate;
    }

    if (allFixed) {
        tree[pix].status = 0;
        tree[pix].rate = totalRate;
        updateParentStatus(tree, pix);
    }
}

// Explicit template instantiations
template struct GroupTreeNode<double>;
template struct GroupTreeNode<float>;

template int distributeGroupTreeRates<double>(
    std::vector<GroupTreeNode<double>>&, int, int);
template int distributeGroupTreeRates<float>(
    std::vector<GroupTreeNode<float>>&, int, int);

template void setSubRates<double>(
    std::vector<GroupTreeNode<double>>&, int);
template void setSubRates<float>(
    std::vector<GroupTreeNode<float>>&, int);

template int findWorstOffendingChild<double>(
    const std::vector<GroupTreeNode<double>>&, int, double&);
template int findWorstOffendingChild<float>(
    const std::vector<GroupTreeNode<float>>&, int, float&);

template void updateParentStatus<double>(
    std::vector<GroupTreeNode<double>>&, int);
template void updateParentStatus<float>(
    std::vector<GroupTreeNode<float>>&, int);

} // namespace Opm
