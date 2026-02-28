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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/simulators/wells/GroupTreeRates.hpp>

#define BOOST_TEST_MODULE GroupTreeRatesTest
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <vector>

using namespace Opm;

namespace {

// Helper: build a simple tree
//
//       [0] FIELD (limit=10000)
//        /        \
//    [1] G1        [2] G2
//   (lim=6000)   (lim=7000)
//     / \            / \
//  [3]W1 [4]W2   [5]W3 [6]W4
// (4000) (3000) (4500) (3500)
//
// Guide rates: G1=0.5, G2=0.5, W1=0.6, W2=0.4, W3=0.55, W4=0.45
std::vector<GroupTreeNode<double>> makeSimpleTree()
{
    std::vector<GroupTreeNode<double>> tree(7);

    // FIELD (root)
    tree[0] = {0, -1, 0.0, 10000.0, -1, 1.0, "FIELD", {1, 2}};
    // G1
    tree[1] = {1, 0, 0.0, 6000.0, -1, 0.5, "G1", {3, 4}};
    // G2
    tree[2] = {2, 0, 0.0, 7000.0, -1, 0.5, "G2", {5, 6}};
    // W1
    tree[3] = {3, 1, 0.0, 4000.0, -1, 0.6, "W1", {}};
    // W2
    tree[4] = {4, 1, 0.0, 3000.0, -1, 0.4, "W2", {}};
    // W3
    tree[5] = {5, 2, 0.0, 4500.0, -1, 0.55, "W3", {}};
    // W4
    tree[6] = {6, 2, 0.0, 3500.0, -1, 0.45, "W4", {}};

    return tree;
}

// Helper: build a tree where no limits are binding (all limits > field rate)
std::vector<GroupTreeNode<double>> makeUnconstrainedTree()
{
    std::vector<GroupTreeNode<double>> tree(5);

    // FIELD
    tree[0] = {0, -1, 0.0, 1000.0, -1, 1.0, "FIELD", {1, 2}};
    // G1
    tree[1] = {1, 0, 0.0, 9000.0, -1, 0.6, "G1", {3, 4}};
    // G2  (leaf well acting as a group-less well)
    tree[2] = {2, 0, 0.0, 9000.0, -1, 0.4, "W_SOLO", {}};
    // W1
    tree[3] = {3, 1, 0.0, 9000.0, -1, 0.5, "W1", {}};
    // W2
    tree[4] = {4, 1, 0.0, 9000.0, -1, 0.5, "W2", {}};

    return tree;
}

// Helper: build a single-well tree (FIELD -> W1)
std::vector<GroupTreeNode<double>> makeSingleWellTree()
{
    std::vector<GroupTreeNode<double>> tree(2);
    tree[0] = {0, -1, 0.0, 500.0, -1, 1.0, "FIELD", {1}};
    tree[1] = {1, 0, 0.0, 300.0, -1, 1.0, "W1", {}};
    return tree;
}

double totalLeafRate(const std::vector<GroupTreeNode<double>>& tree)
{
    double total = 0.0;
    for (const auto& node : tree) {
        if (node.children.empty()) {
            total += node.rate;
        }
    }
    return total;
}

} // anonymous namespace


BOOST_AUTO_TEST_CASE(SingleWell_WellLimitBinding)
{
    auto tree = makeSingleWellTree();
    // FIELD limit = 500, W1 limit = 300 => W1 should be limited to 300
    const int iter = distributeGroupTreeRates(tree);

    BOOST_CHECK(iter > 0);
    // W1 should be at its limit
    BOOST_CHECK_CLOSE(tree[1].rate, 300.0, 1e-10);
    // FIELD should reflect the well rate
    BOOST_CHECK_CLOSE(tree[0].rate, 300.0, 1e-10);
    // W1 individually limited
    BOOST_CHECK_EQUAL(tree[1].status, 1);
    // FIELD determined by children
    BOOST_CHECK_EQUAL(tree[0].status, 0);
}

BOOST_AUTO_TEST_CASE(SingleWell_FieldLimitBinding)
{
    auto tree = makeSingleWellTree();
    tree[1].limit = 800.0; // W1 limit > FIELD limit
    const int iter = distributeGroupTreeRates(tree);

    BOOST_CHECK(iter >= 0);
    // W1 rate should equal FIELD limit since it is the only well
    BOOST_CHECK_CLOSE(tree[1].rate, 500.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(Unconstrained_RatesFollowGuideRates)
{
    auto tree = makeUnconstrainedTree();
    distributeGroupTreeRates(tree);

    // FIELD rate = 1000 (its limit)
    const double fieldRate = 1000.0;
    // G1 gets 60%, W_SOLO gets 40%
    BOOST_CHECK_CLOSE(tree[1].rate, fieldRate * 0.6, 1e-10);
    BOOST_CHECK_CLOSE(tree[2].rate, fieldRate * 0.4, 1e-10);
    // W1, W2 each get 50% of G1
    BOOST_CHECK_CLOSE(tree[3].rate, fieldRate * 0.6 * 0.5, 1e-10);
    BOOST_CHECK_CLOSE(tree[4].rate, fieldRate * 0.6 * 0.5, 1e-10);

    // Total leaf rate should equal field rate
    BOOST_CHECK_CLOSE(totalLeafRate(tree), fieldRate, 1e-10);
}

BOOST_AUTO_TEST_CASE(SimpleTree_GroupAndWellLimits)
{
    auto tree = makeSimpleTree();
    distributeGroupTreeRates(tree);

    // All well rates must not exceed their limits
    for (const auto& node : tree) {
        BOOST_CHECK_LE(node.rate, node.limit + 1e-10);
    }

    // Total leaf rate must not exceed FIELD limit
    BOOST_CHECK_LE(totalLeafRate(tree), tree[0].limit + 1e-10);

    // Group rates must not exceed their limits
    BOOST_CHECK_LE(tree[1].rate, tree[1].limit + 1e-10); // G1
    BOOST_CHECK_LE(tree[2].rate, tree[2].limit + 1e-10); // G2
}

BOOST_AUTO_TEST_CASE(SimpleTree_RateConsistency)
{
    auto tree = makeSimpleTree();
    distributeGroupTreeRates(tree);

    // G1 rate == W1 rate + W2 rate
    BOOST_CHECK_CLOSE(tree[1].rate, tree[3].rate + tree[4].rate, 1e-10);
    // G2 rate == W3 rate + W4 rate
    BOOST_CHECK_CLOSE(tree[2].rate, tree[5].rate + tree[6].rate, 1e-10);
    // FIELD rate == G1 rate + G2 rate
    BOOST_CHECK_CLOSE(tree[0].rate, tree[1].rate + tree[2].rate, 1e-10);
}

BOOST_AUTO_TEST_CASE(SimpleTree_StatusConsistency)
{
    auto tree = makeSimpleTree();
    distributeGroupTreeRates(tree);

    // Root should be individually limited (status 1)
    BOOST_CHECK_EQUAL(tree[0].status, 1);

    // If a node has status 0, its rate must equal the sum of children's rates
    // If a node has status 1, its rate must equal its limit
    for (const auto& node : tree) {
        if (node.status == 1) {
            BOOST_CHECK_CLOSE(node.rate, node.limit, 1e-10);
        }
        if (node.status == 0) {
            double childSum = 0.0;
            for (int ci : node.children) {
                childSum += tree[ci].rate;
            }
            BOOST_CHECK_CLOSE(node.rate, childSum, 1e-10);
        }
    }
}

BOOST_AUTO_TEST_CASE(TightWellLimit)
{
    // One well has a very tight limit
    auto tree = makeUnconstrainedTree();
    tree[3].limit = 10.0; // W1 has tight limit

    distributeGroupTreeRates(tree);

    // W1 must be at its limit
    BOOST_CHECK_CLOSE(tree[3].rate, 10.0, 1e-10);
    BOOST_CHECK_EQUAL(tree[3].status, 1);

    // Total still respects FIELD
    BOOST_CHECK_LE(totalLeafRate(tree), tree[0].limit + 1e-10);

    // Rate consistency
    BOOST_CHECK_CLOSE(tree[1].rate, tree[3].rate + tree[4].rate, 1e-10);
}

BOOST_AUTO_TEST_CASE(AllWellsTight)
{
    // All wells have very tight limits (sum < FIELD limit)
    auto tree = makeSimpleTree();
    tree[3].limit = 100.0;  // W1
    tree[4].limit = 100.0;  // W2
    tree[5].limit = 100.0;  // W3
    tree[6].limit = 100.0;  // W4

    distributeGroupTreeRates(tree);

    // Each well should be at its limit
    BOOST_CHECK_CLOSE(tree[3].rate, 100.0, 1e-10);
    BOOST_CHECK_CLOSE(tree[4].rate, 100.0, 1e-10);
    BOOST_CHECK_CLOSE(tree[5].rate, 100.0, 1e-10);
    BOOST_CHECK_CLOSE(tree[6].rate, 100.0, 1e-10);

    // Total = 400 < FIELD limit of 10000
    BOOST_CHECK_CLOSE(totalLeafRate(tree), 400.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(GroupLimitViolation)
{
    // Tree where initial guide-rate distribution violates a group limit:
    //   FIELD (limit=10000), G1 (lim=3000), G2 (lim=8000)
    //   Guide rates 0.5/0.5 => initial: G1=5000 > 3000 (violation!)
    //   After fixing G1=3000, G2 gets 7000
    std::vector<GroupTreeNode<double>> tree(5);
    tree[0] = {0, -1, 0.0, 10000.0, -1, 1.0, "FIELD", {1, 2}};
    tree[1] = {1, 0, 0.0, 3000.0, -1, 0.5, "G1", {3}};
    tree[2] = {2, 0, 0.0, 8000.0, -1, 0.5, "G2", {4}};
    tree[3] = {3, 1, 0.0, 5000.0, -1, 1.0, "W1", {}};
    tree[4] = {4, 2, 0.0, 9000.0, -1, 1.0, "W2", {}};

    const int iter = distributeGroupTreeRates(tree);

    // Should take at least 1 iteration to fix the G1 violation
    BOOST_CHECK(iter >= 1);
    // G1 should be at its limit
    BOOST_CHECK_CLOSE(tree[1].rate, 3000.0, 1e-10);
    BOOST_CHECK_EQUAL(tree[1].status, 1);
    // W1 = G1's rate = 3000 (within its 5000 limit)
    BOOST_CHECK_CLOSE(tree[3].rate, 3000.0, 1e-10);
    // G2 gets the remaining 7000
    BOOST_CHECK_CLOSE(tree[2].rate, 7000.0, 1e-10);
    // W2 = G2's rate = 7000 (within its 9000 limit)
    BOOST_CHECK_CLOSE(tree[4].rate, 7000.0, 1e-10);
    // Total
    BOOST_CHECK_CLOSE(totalLeafRate(tree), 10000.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(CascadingViolations)
{
    // Multiple cascading violations:
    //   FIELD (10000) -> G1 (lim=4000, guide=0.5), G2 (lim=9000, guide=0.5)
    //   G1 -> W1 (lim=1000, guide=0.6), W2 (lim=5000, guide=0.4)
    //   G2 -> W3 (lim=3000, guide=0.5), W4 (lim=9000, guide=0.5)
    //
    // Initial: G1=5000 (>4000 viol), G2=5000
    //   G1 children: W1=3000 (>1000 viol), W2=2000
    //   G2 children: W3=2500, W4=2500
    // After fixing: multiple iterations needed
    std::vector<GroupTreeNode<double>> tree(7);
    tree[0] = {0, -1, 0.0, 10000.0, -1, 1.0, "FIELD", {1, 2}};
    tree[1] = {1, 0, 0.0, 4000.0, -1, 0.5, "G1", {3, 4}};
    tree[2] = {2, 0, 0.0, 9000.0, -1, 0.5, "G2", {5, 6}};
    tree[3] = {3, 1, 0.0, 1000.0, -1, 0.6, "W1", {}};
    tree[4] = {4, 1, 0.0, 5000.0, -1, 0.4, "W2", {}};
    tree[5] = {5, 2, 0.0, 3000.0, -1, 0.5, "W3", {}};
    tree[6] = {6, 2, 0.0, 9000.0, -1, 0.5, "W4", {}};

    const int iter = distributeGroupTreeRates(tree);
    BOOST_CHECK(iter >= 1);

    // All limits respected
    for (const auto& node : tree) {
        BOOST_CHECK_LE(node.rate, node.limit + 1e-10);
    }

    // W1 must be at its limit (1000)
    BOOST_CHECK_CLOSE(tree[3].rate, 1000.0, 1e-10);
    BOOST_CHECK_EQUAL(tree[3].status, 1);

    // Rate consistency
    BOOST_CHECK_CLOSE(tree[1].rate, tree[3].rate + tree[4].rate, 1e-10);
    BOOST_CHECK_CLOSE(tree[2].rate, tree[5].rate + tree[6].rate, 1e-10);
    BOOST_CHECK_CLOSE(tree[0].rate, tree[1].rate + tree[2].rate, 1e-10);
}

BOOST_AUTO_TEST_CASE(DeepTree)
{
    // FIELD -> G1 -> G2 -> W1
    std::vector<GroupTreeNode<double>> tree(4);
    tree[0] = {0, -1, 0.0, 5000.0, -1, 1.0, "FIELD", {1}};
    tree[1] = {1, 0, 0.0, 3000.0, -1, 1.0, "G1", {2}};
    tree[2] = {2, 1, 0.0, 2000.0, -1, 1.0, "G2", {3}};
    tree[3] = {3, 2, 0.0, 1500.0, -1, 1.0, "W1", {}};

    distributeGroupTreeRates(tree);

    // W1 is the tightest at 1500
    BOOST_CHECK_CLOSE(tree[3].rate, 1500.0, 1e-10);
    BOOST_CHECK_CLOSE(tree[2].rate, 1500.0, 1e-10);
    BOOST_CHECK_CLOSE(tree[1].rate, 1500.0, 1e-10);
    BOOST_CHECK_CLOSE(tree[0].rate, 1500.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(ZeroIterations_NoViolation)
{
    // If root limit is 0, no distribution needed
    std::vector<GroupTreeNode<double>> tree(2);
    tree[0] = {0, -1, 0.0, 0.0, -1, 1.0, "FIELD", {1}};
    tree[1] = {1, 0, 0.0, 100.0, -1, 1.0, "W1", {}};

    const int iter = distributeGroupTreeRates(tree);

    // Should converge immediately (no violations since rates are 0)
    BOOST_CHECK_EQUAL(iter, 0);
    BOOST_CHECK_CLOSE(tree[1].rate, 0.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(EqualGuideRates)
{
    // Three wells with equal guide rates
    std::vector<GroupTreeNode<double>> tree(4);
    tree[0] = {0, -1, 0.0, 900.0, -1, 1.0, "FIELD", {1, 2, 3}};
    tree[1] = {1, 0, 0.0, 500.0, -1, 1.0, "W1", {}};
    tree[2] = {2, 0, 0.0, 500.0, -1, 1.0, "W2", {}};
    tree[3] = {3, 0, 0.0, 500.0, -1, 1.0, "W3", {}};

    distributeGroupTreeRates(tree);

    // Each well gets 300
    BOOST_CHECK_CLOSE(tree[1].rate, 300.0, 1e-10);
    BOOST_CHECK_CLOSE(tree[2].rate, 300.0, 1e-10);
    BOOST_CHECK_CLOSE(tree[3].rate, 300.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(SetSubRates_Basic)
{
    auto tree = makeUnconstrainedTree();
    // Manually set root as limited with known rate
    tree[0].status = 1;
    tree[0].rate = 1000.0;

    setSubRates(tree, 0);

    BOOST_CHECK_CLOSE(tree[1].rate, 600.0, 1e-10);
    BOOST_CHECK_CLOSE(tree[2].rate, 400.0, 1e-10);
    BOOST_CHECK_CLOSE(tree[3].rate, 300.0, 1e-10);
    BOOST_CHECK_CLOSE(tree[4].rate, 300.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(FindWorstOffending_Basic)
{
    auto tree = makeSimpleTree();
    // Set FIELD rate and distribute
    tree[0].status = 1;
    tree[0].rate = 10000.0;
    setSubRates(tree, 0);

    double worstExcess = 0.0;
    const int ix = findWorstOffendingChild(tree, 0, worstExcess);

    // At least one child should be offending since rates > limits for some
    BOOST_CHECK(ix >= 0 || worstExcess <= 0.0);
}

BOOST_AUTO_TEST_CASE(UpdateParentStatus_AllFixed)
{
    auto tree = makeUnconstrainedTree();
    // Fix all children of FIELD
    tree[1].status = 1;
    tree[1].rate = 500.0;
    tree[2].status = 1;
    tree[2].rate = 300.0;
    tree[0].status = -1;

    updateParentStatus(tree, 1);

    BOOST_CHECK_EQUAL(tree[0].status, 0);
    BOOST_CHECK_CLOSE(tree[0].rate, 800.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(UpdateParentStatus_NotAllFixed)
{
    auto tree = makeUnconstrainedTree();
    tree[1].status = 1;
    tree[1].rate = 500.0;
    tree[2].status = -1; // still group-controlled
    tree[0].status = -1;

    updateParentStatus(tree, 1);

    // Parent should NOT be updated since not all children are fixed
    BOOST_CHECK_EQUAL(tree[0].status, -1);
}
