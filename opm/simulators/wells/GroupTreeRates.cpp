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

namespace Opm {

template<class Scalar>
int GroupTreeRates<Scalar>::
distribute(std::vector<GroupTreeNode<Scalar>>& tree,
           int root_index,
           int max_iter)
{
    if (root_index == 0) {
        // Initialise: all nodes group-controlled except root
        for (auto& node : tree) {
            node.status = -1;
        }
        tree[root_index].status = 1;
        tree[root_index].rate = tree[root_index].limit;
    }

    int iter = 0;
    for (; iter < max_iter; ++iter) {
        set_sub_rates(tree, root_index);

        Scalar worst_excess{0};
        const int ix = find_worst_offending_child(tree, root_index, worst_excess);
        if (ix < 0) {
            break; // converged
        }

        // Fix the worst-offending node at its limit
        tree[ix].rate = tree[ix].limit;
        tree[ix].status = 1;

        // If it has children, recursively solve its subtree
        if (!tree[ix].children.empty()) {
            distribute(tree, ix, max_iter);
        }

        update_parent_status(tree, ix);
    }
    return iter;
}

template<class Scalar>
void GroupTreeRates<Scalar>::
set_sub_rates(std::vector<GroupTreeNode<Scalar>>& tree,
              int node_index)
{
    const auto& children = tree[node_index].children;
    if (children.empty()) {
        return; // leaf (well) node
    }

    // Separate children into fixed (status != -1) and group-controlled
    Scalar fixed_rate{0};
    Scalar guide_sum{0};
    bool all_fixed = true;

    for (const int ci : children) {
        if (tree[ci].status != -1) {
            fixed_rate += tree[ci].rate;
        } else {
            guide_sum += tree[ci].guide_rate;
            all_fixed = false;
        }
    }

    if (all_fixed) {
        tree[node_index].status = 0;
        tree[node_index].rate = fixed_rate;
        return;
    }

    const Scalar avail_rate = tree[node_index].rate - fixed_rate;
    assert(avail_rate >= Scalar{0});
    assert(avail_rate > Scalar{0} || guide_sum == Scalar{0});

    for (const int ci : children) {
        if (tree[ci].status == -1) {
            tree[ci].rate = (guide_sum > Scalar{0})
                ? avail_rate * tree[ci].guide_rate / guide_sum
                : Scalar{0};
            set_sub_rates(tree, ci);
        }
    }
}

template<class Scalar>
int GroupTreeRates<Scalar>::
find_worst_offending_child(const std::vector<GroupTreeNode<Scalar>>& tree,
                           int node_index,
                           Scalar& worst_excess)
{
    const auto& children = tree[node_index].children;
    int worst_index = -1;

    for (const int ci : children) {
        const Scalar excess = tree[ci].rate - tree[ci].limit;
        if (excess > worst_excess) {
            worst_excess = excess;
            worst_index = ci;
        }
    }

    // Recurse into children that have subtrees
    for (const int ci : children) {
        if (!tree[ci].children.empty()) {
            Scalar subtree_excess = worst_excess;
            const int child_worst = find_worst_offending_child(tree, ci, subtree_excess);
            if (subtree_excess > worst_excess) {
                worst_excess = subtree_excess;
                worst_index = child_worst;
            }
        }
    }

    if (worst_excess > Scalar{0} && worst_index >= 0) {
        assert(tree[worst_index].status == -1);
    } else {
        worst_index = -1;
    }
    return worst_index;
}

template<class Scalar>
void GroupTreeRates<Scalar>::
update_parent_status(std::vector<GroupTreeNode<Scalar>>& tree,
                     int node_index)
{
    const int pix = tree[node_index].parent;
    if (pix < 0) {
        return;
    }

    const auto& siblings = tree[pix].children;
    bool all_fixed = true;
    Scalar total_rate{0};

    for (const int ci : siblings) {
        if (tree[ci].status == -1) {
            all_fixed = false;
            break;
        }
        total_rate += tree[ci].rate;
    }

    if (all_fixed) {
        tree[pix].status = 0;
        tree[pix].rate = total_rate;
        update_parent_status(tree, pix);
    }
}

// Explicit template instantiations
template struct GroupTreeNode<double>;
template struct GroupTreeNode<float>;

template class GroupTreeRates<double>;
template class GroupTreeRates<float>;

} // namespace Opm
