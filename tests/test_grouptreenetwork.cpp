/*
  Copyright 2026 SINTEF Digital

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

#define BOOST_TEST_MODULE GroupTreeNetworkTests

#include <opm/simulators/wells/ProdGroupTreeBalancer.hpp>
#include <opm/simulators/wells/ProdGroupTreeNode.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>

#include <fmt/format.h>

#include <array>
#include <map>
#include <optional>
#include <string>
#include <unordered_set>
#include <vector>

using namespace Opm;

namespace {

// A deck with nothing in it, just enough for GuideRate's constructor. Deck,
// EclipseState and Schedule all need to stay alive together in the caller's
// own scope -- Schedule holds references into the other two, so returning
// one of these by value from a helper leaves it holding dangling references
// the moment the helper returns.
// Same minimal-grid skeleton as tests/test_wellprodindexcalculator.cpp's
// createWell(), which is known to build a Deck/EclipseState/Schedule without
// tripping EclipseState's region-statistics setup -- a hand-rolled DX/DY/DZ
// grid with no PROPS/REGIONS section left EclipseState's destructor
// segfaulting on a half-initialized std::optional<FIPRegionStatistics>.
constexpr const char* kEmptyDeck = R"(RUNSPEC
DIMENS
  10 10 3 /
START
 8 OCT 2020 /
GRID
DXV
  10*100.0 /
DYV
  10*100.0 /
DZV
  3*10.0 /
DEPTHZ
  121*2000.0 /
PERMX
  300*100.0 /
PERMY
  300*100.0 /
PERMZ
  300*10.0 /
PORO
  300*0.3 /
SCHEDULE
TSTEP
  10
/
END
)";

// The same three-level tree used to validate the balancer against the
// equations elsewhere: FIELD -> PLAT -> {GP1 -> {W1, W2}, GP2 -> {W3}}.
// Spike for Part 1 of groups_and_network_clean.md: what does
// categorizeBalancedNode() actually leave on each node once the balancer
// has converged, so the flat-extraction pass can be written against real
// behaviour instead of a guess from reading the source.
ProdGroupTreeBalancer::Tree<double>
buildTree(double target, const std::array<double, 3>& cap, const std::array<double, 3>& guide)
{
    const std::map<std::string, std::string> parent{
        {"PLAT", "FIELD"}, {"GP1", "PLAT"}, {"GP2", "PLAT"},
        {"W1", "GP1"}, {"W2", "GP1"}, {"W3", "GP2"}};
    const std::vector<std::string> wnames{"W1", "W2", "W3"};

    ProdGroupTreeBalancer::Tree<double> tree;
    auto addGroup = [&](const std::string& name, const std::vector<std::string>& kids, const bool has_target) {
        ProdGroupTreeNode<double> n;
        n.name = name;
        n.type = ProdNodeType::Group;
        n.parent = parent.count(name) ? parent.at(name) : std::string{};
        n.children = kids;
        n.availableForGroupControl = true;
        // Left false, matching the balancer-vs-equations test this tree is
        // copied from: setting it true makes propagateGuideRatesAndMode() look
        // the group up in `guideRate` too, which then needs guide_rate.compute()
        // called for every group name as well as every well -- worth remembering
        // for Part 4's real integration, not needed for this spike's purpose.
        n.modeCategory = has_target ? ProdNodeModeCategory::Individual : ProdNodeModeCategory::Group;
        if (has_target) {
            n.mode = Well::ProducerCMode::ORAT;
            n.preferredMode = Group::ProductionCMode::ORAT;
            n.Limits[Well::ProducerCMode::ORAT] = target / 86400.0;
        }
        tree.emplace(name, std::move(n));
    };
    addGroup("FIELD", {"PLAT"}, false);
    addGroup("PLAT", {"GP1", "GP2"}, true);
    addGroup("GP1", {"W1", "W2"}, false);
    addGroup("GP2", {"W3"}, false);
    for (std::size_t i = 0; i < wnames.size(); ++i) {
        ProdGroupTreeNode<double> n;
        n.name = wnames[i];
        n.type = ProdNodeType::Well;
        n.parent = parent.at(wnames[i]);
        n.availableForGroupControl = true;
        n.mode = Well::ProducerCMode::GRUP;
        n.hasGuideRate = true;
        n.Limits[Well::ProducerCMode::ORAT] = cap[i] / 86400.0;
        n.rates = {-cap[i] / 86400.0, 0.0, 0.0};
        n.initialRates = n.rates;
        tree.emplace(wnames[i], std::move(n));
    }
    (void)guide;   // guide rates are supplied to the GuideRate object, not the tree
    return tree;
}

void logNode(const ProdGroupTreeBalancer::Tree<double>& tree, const std::string& name)
{
    const auto& n = tree.at(name);
    const char* cat = [&] {
        switch (n.modeCategory) {
        case ProdNodeModeCategory::Group: return "Group";
        case ProdNodeModeCategory::Individual: return "Individual";
        case ProdNodeModeCategory::None: return "None";
        case ProdNodeModeCategory::Transparent: return "Transparent";
        }
        return "?";
    }();
    const double oil = -n.rates[0] * 86400.0;
    BOOST_TEST_MESSAGE(fmt::format("  {:6s} category={:11s} mode={:3d} oil={:.1f}",
                                    name, cat, static_cast<int>(n.mode), oil));
}

void logFlat(const ProdGroupTreeBalancer::FlatNetworkInput<double>& flat)
{
    for (const auto& n : flat) {
        std::string wells;
        for (const auto& w : n.ownWells) {
            wells += fmt::format(" {}(w={:.6g})", w.name, w.weight);
        }
        std::string thp;
        for (const auto& w : n.thpWells) {
            thp += fmt::format(" {}(eff={:.3f})", w.name, w.weight);
        }
        std::string children;
        for (const auto& c : n.activeChildren) {
            children += fmt::format(" {}(eff={:.3f})", c.name, c.efficiency);
        }
        BOOST_TEST_MESSAGE(fmt::format("  {:6s} target={:.1f}  ownWells:{}  thpWells:{}  activeChildren:{}",
                                        n.name, n.target * 86400.0,
                                        wells.empty() ? " (none)" : wells,
                                        thp.empty() ? " (none)" : thp,
                                        children.empty() ? " (none)" : children));
    }
}

// A tree in the *converged* state extractFlatNetworkInput() should see once
// two limits bind at different levels: PLAT's own ORAT target (3000) and,
// nested inside it, GP1's own tighter ORAT limit (1800) -- W1/W2 split that
// 900/900, GP2/W3 (pass-through) absorbs PLAT's remaining budget (1200).
// This is the case the plan's first cut at Part 1 got wrong: PLAT's own sum
// must still see GP1's (and so W1/W2's) production even though W1/W2 answer
// to GP1's lambda, not PLAT's -- hence activeChildren rather than a
// further-flattened well list.
//
// Built directly in this final state rather than by running the balancer:
// extractFlatNetworkInput() is a pure function of whatever categorization
// and rates the tree holds, and coercing the real (sorted, guide-rate-ratio
// driven) distribution algorithm into landing on one specific hand-picked
// split is its own exercise, independent of what this is testing -- see the
// modecategory_semantics/flat_extraction_on_a_pass_through_tree tests above
// for that (they run the real balancer and check its output directly).
ProdGroupTreeBalancer::Tree<double> buildNestedActiveTree()
{
    const std::map<std::string, std::string> parent{
        {"PLAT", "FIELD"}, {"GP1", "PLAT"}, {"GP2", "PLAT"},
        {"W1", "GP1"}, {"W2", "GP1"}, {"W3", "GP2"}};

    ProdGroupTreeBalancer::Tree<double> tree;
    auto addGroup = [&](const std::string& name, const std::vector<std::string>& kids,
                        ProdNodeModeCategory category, double oil,
                        const std::optional<double>& target) {
        ProdGroupTreeNode<double> n;
        n.name = name;
        n.type = ProdNodeType::Group;
        n.parent = parent.count(name) ? parent.at(name) : std::string{};
        n.children = kids;
        n.availableForGroupControl = true;
        n.modeCategory = category;
        n.rates = {-oil / 86400.0, 0.0, 0.0};
        n.initialRates = n.rates;
        if (target) {
            n.mode = Well::ProducerCMode::ORAT;
            n.preferredMode = Group::ProductionCMode::ORAT;
            n.Limits[Well::ProducerCMode::ORAT] = *target / 86400.0;
        }
        tree.emplace(name, std::move(n));
    };
    addGroup("FIELD", {"PLAT"}, ProdNodeModeCategory::None, 3000.0, std::nullopt);
    addGroup("PLAT", {"GP1", "GP2"}, ProdNodeModeCategory::Individual, 3000.0, 3000.0);
    addGroup("GP1", {"W1", "W2"}, ProdNodeModeCategory::Individual, 1800.0, 1800.0);
    addGroup("GP2", {"W3"}, ProdNodeModeCategory::Transparent, 1200.0, std::nullopt);

    auto addWell = [&](const std::string& name, double oil, double guide) {
        ProdGroupTreeNode<double> n;
        n.name = name;
        n.type = ProdNodeType::Well;
        n.parent = parent.at(name);
        n.availableForGroupControl = true;
        n.mode = Well::ProducerCMode::GRUP;
        n.modeCategory = ProdNodeModeCategory::Group;
        n.hasGuideRate = true;
        n.rates = {-oil / 86400.0, 0.0, 0.0};
        n.initialRates = n.rates;
        n.groupTarget.guideRate = guide / 86400.0;   // unused by extraction; belt and suspenders
        tree.emplace(name, std::move(n));
    };
    addWell("W1", 900.0, 1.0);
    addWell("W2", 900.0, 1.0);
    addWell("W3", 1200.0, 1.0);
    return tree;
}

} // namespace

BOOST_AUTO_TEST_CASE(modecategory_semantics_on_a_three_level_tree)
{
    // "target binds, one well capped": PLAT's own ORAT target (3000) is what
    // actually binds; W3 is capped below its share (400 < 1000), so GP1
    // (W1+W2) carries the remaining 2600, split evenly (1300/1300).
    Opm::Parser parser;
    const auto deck = parser.parseString(kEmptyDeck);
    const EclipseState es{deck};
    const Schedule schedule{deck, es};
    GuideRate guide_rate{schedule};
    DeferredLogger logger;
    const std::vector<std::string> wnames{"W1", "W2", "W3"};
    const std::array<double, 3> guide{1.0, 1.0, 1.0};
    for (std::size_t i = 0; i < wnames.size(); ++i) {
        guide_rate.compute(wnames[i], /*report_step=*/0, /*sim_time=*/0.0, guide[i] / 86400.0, 0.0, 0.0);
    }

    auto tree = buildTree(/*target=*/3000.0, /*cap=*/{2000.0, 2000.0, 400.0}, guide);
    const bool ok = ProdGroupTreeBalancer::balanceTreeForTesting(tree, guide_rate, 1e-8, logger);
    BOOST_REQUIRE(ok);

    BOOST_TEST_MESSAGE("Tree after balancing (target 3000, caps 2000/2000/400):");
    for (const auto& name : {std::string("FIELD"), std::string("PLAT"), std::string("GP1"),
                             std::string("GP2"), std::string("W1"), std::string("W2"), std::string("W3")}) {
        logNode(tree, name);
    }

    // What Part 1 of groups_and_network_clean.md needs to know: PLAT's own
    // limit is what's actually binding, so PLAT must be the node the flat
    // extraction treats as "Active" (its own target T, a new lambda).
    BOOST_CHECK(tree.at("PLAT").modeCategory == ProdNodeModeCategory::Individual);
    BOOST_CHECK_CLOSE(-tree.at("PLAT").rates[0] * 86400.0, 3000.0, 1e-6);

    // W3 is capped by its own limit -- pinned, excluded from any group sum.
    BOOST_CHECK(tree.at("W3").modeCategory == ProdNodeModeCategory::Individual);
    BOOST_CHECK_CLOSE(-tree.at("W3").rates[0] * 86400.0, 400.0, 1e-6);

    // W1, W2 are fully under PLAT's (flattened, through GP1) control.
    BOOST_CHECK(tree.at("W1").modeCategory == ProdNodeModeCategory::Group);
    BOOST_CHECK(tree.at("W2").modeCategory == ProdNodeModeCategory::Group);
    BOOST_CHECK_CLOSE(-tree.at("W1").rates[0] * 86400.0, 1300.0, 1e-6);
    BOOST_CHECK_CLOSE(-tree.at("W2").rates[0] * 86400.0, 1300.0, 1e-6);
}

// Part 1 proper: the same tree flattened. PLAT is the top-level Active node
// (GP1/GP2 are pure pass-through); it owns W1 and W2 directly, and references
// W3 -- individually limited, its own separate entry -- as an active child,
// since W3's production still counts toward PLAT's own 3000 target even
// though W3's own row is what actually determines its rate.
BOOST_AUTO_TEST_CASE(flat_extraction_on_a_pass_through_tree)
{
    Opm::Parser parser;
    const auto deck = parser.parseString(kEmptyDeck);
    const EclipseState es{deck};
    const Schedule schedule{deck, es};
    GuideRate guide_rate{schedule};
    DeferredLogger logger;
    const std::vector<std::string> wnames{"W1", "W2", "W3"};
    const std::array<double, 3> guide{3.0, 2.0, 1.0};   // unequal, so a mismatched weight would show up
    for (std::size_t i = 0; i < wnames.size(); ++i) {
        guide_rate.compute(wnames[i], 0, 0.0, guide[i] / 86400.0, 0.0, 0.0);
    }

    auto tree = buildTree(/*target=*/3000.0, /*cap=*/{2000.0, 2000.0, 400.0}, guide);
    BOOST_REQUIRE(ProdGroupTreeBalancer::balanceTreeForTesting(tree, guide_rate, 1e-8, logger));

    const auto flat = ProdGroupTreeBalancer::extractFlatNetworkInput(tree, std::string("FIELD"), guide_rate);
    BOOST_TEST_MESSAGE("Flattened (pass-through) tree:");
    logFlat(flat);

    // Exactly two Active nodes: PLAT (a top-level root) and W3 (referenced
    // from PLAT's activeChildren, not itself a root). FIELD/GP1/GP2 never
    // appear at all.
    BOOST_REQUIRE_EQUAL(flat.size(), 2U);
    const auto& plat = (flat[0].name == "PLAT") ? flat[0] : flat[1];
    const auto& w3 = (flat[0].name == "W3") ? flat[0] : flat[1];

    BOOST_CHECK_EQUAL(plat.name, "PLAT");
    BOOST_CHECK(plat.type == ProdNodeType::Group);
    BOOST_CHECK_CLOSE(plat.target * 86400.0, 3000.0, 1e-6);
    BOOST_REQUIRE_EQUAL(plat.ownWells.size(), 2U);
    // Both W1 and W2 sit directly under GP1, one flattened hop from PLAT, so
    // their weight is exactly their own guide rate (efficiency 1 throughout).
    for (const auto& w : plat.ownWells) {
        const double expected_guide = (w.name == "W1") ? guide[0] : (w.name == "W2") ? guide[1] : -1.0;
        BOOST_CHECK(expected_guide > 0.0);
        BOOST_CHECK_CLOSE(w.weight, expected_guide / 86400.0, 1e-6);
    }
    BOOST_REQUIRE_EQUAL(plat.activeChildren.size(), 1U);
    BOOST_CHECK_EQUAL(plat.activeChildren.front().name, "W3");
    BOOST_CHECK_CLOSE(plat.activeChildren.front().efficiency, 1.0, 1e-6);

    BOOST_CHECK_EQUAL(w3.name, "W3");
    BOOST_CHECK(w3.type == ProdNodeType::Well);
    BOOST_CHECK_CLOSE(w3.target * 86400.0, 400.0, 1e-6);
    BOOST_CHECK(w3.ownWells.empty());
    BOOST_CHECK(w3.activeChildren.empty());
}

// The case Part 1's first draft got wrong: GP1 has its own tighter limit, so
// it is Individual too, nested inside PLAT's. PLAT's own sum must still
// account for GP1's production even though W1/W2 now answer to GP1's lambda,
// not PLAT's -- so PLAT should list GP1 as an active child, NOT list W1/W2 in
// its own wells, and GP1 should appear as its own separate entry owning them.
BOOST_AUTO_TEST_CASE(flat_extraction_with_a_nested_active_group)
{
    Opm::Parser parser;
    const auto deck = parser.parseString(kEmptyDeck);
    const EclipseState es{deck};
    const Schedule schedule{deck, es};
    GuideRate guide_rate{schedule};
    for (const auto& name : {std::string("W1"), std::string("W2"), std::string("W3")}) {
        guide_rate.compute(name, 0, 0.0, 1.0 / 86400.0, 0.0, 0.0);
    }

    auto tree = buildNestedActiveTree();
    const auto flat = ProdGroupTreeBalancer::extractFlatNetworkInput(tree, std::string("FIELD"), guide_rate);
    BOOST_TEST_MESSAGE("Flattened (nested-active) tree:");
    logFlat(flat);

    // Two Active nodes: PLAT and GP1 (W3 is Group-category, tied to PLAT
    // through the pass-through GP2, same as W1/W2 were in the first case).
    BOOST_REQUIRE_EQUAL(flat.size(), 2U);
    const auto& plat = (flat[0].name == "PLAT") ? flat[0] : flat[1];
    const auto& gp1 = (flat[0].name == "GP1") ? flat[0] : flat[1];
    BOOST_CHECK_EQUAL(plat.name, "PLAT");
    BOOST_CHECK_EQUAL(gp1.name, "GP1");

    // PLAT: GP1 is an active child, not folded into ownWells; W3 (through
    // GP2) is PLAT's only direct well.
    BOOST_REQUIRE_EQUAL(plat.activeChildren.size(), 1U);
    BOOST_CHECK_EQUAL(plat.activeChildren.front().name, "GP1");
    BOOST_REQUIRE_EQUAL(plat.ownWells.size(), 1U);
    BOOST_CHECK_EQUAL(plat.ownWells.front().name, "W3");

    // GP1: owns W1 and W2 directly, no active children of its own.
    BOOST_CHECK(gp1.activeChildren.empty());
    BOOST_REQUIRE_EQUAL(gp1.ownWells.size(), 2U);
    for (const auto& w : gp1.ownWells) {
        BOOST_CHECK(w.name == "W1" || w.name == "W2");
    }
    BOOST_CHECK_CLOSE(gp1.target * 86400.0, 1800.0, 1e-6);
}

// The balancer has no notion of THP control -- it categorized W2 as Group,
// same as W1 -- but the caller knows (from WellState, in the real
// integration) that W2 is actually on the network's THP control. Extraction
// must move it from ownWells (tied to PLAT's lambda) to thpWells (still
// counted in PLAT's sum, but via efficiency alone, since its rate answers to
// its own bhp/thp/IPR row, not to a guide-rate share of PLAT's target).
BOOST_AUTO_TEST_CASE(thp_controlled_well_is_not_tied_to_the_lambda)
{
    Opm::Parser parser;
    const auto deck = parser.parseString(kEmptyDeck);
    const EclipseState es{deck};
    const Schedule schedule{deck, es};
    GuideRate guide_rate{schedule};
    DeferredLogger logger;
    const std::vector<std::string> wnames{"W1", "W2", "W3"};
    const std::array<double, 3> guide{1.0, 1.0, 1.0};
    for (std::size_t i = 0; i < wnames.size(); ++i) {
        guide_rate.compute(wnames[i], 0, 0.0, guide[i] / 86400.0, 0.0, 0.0);
    }

    auto tree = buildTree(/*target=*/3000.0, /*cap=*/{2000.0, 2000.0, 400.0}, guide);
    BOOST_REQUIRE(ProdGroupTreeBalancer::balanceTreeForTesting(tree, guide_rate, 1e-8, logger));
    BOOST_REQUIRE(tree.at("W2").modeCategory == ProdNodeModeCategory::Group);   // the balancer's own view

    const std::unordered_set<std::string> thpControlled{"W2"};
    const auto flat = ProdGroupTreeBalancer::extractFlatNetworkInput(tree, std::string("FIELD"),
                                                                      guide_rate, thpControlled);
    BOOST_TEST_MESSAGE("Flattened (W2 reclassified as THP-controlled):");
    logFlat(flat);

    BOOST_REQUIRE_EQUAL(flat.size(), 2U);
    const auto& plat = (flat[0].name == "PLAT") ? flat[0] : flat[1];
    BOOST_CHECK_EQUAL(plat.name, "PLAT");

    // W1 stays a normal ownWell; W2 moves to thpWells with weight ==
    // cumulative efficiency alone (1.0 here), not guide-rate*efficiency.
    BOOST_REQUIRE_EQUAL(plat.ownWells.size(), 1U);
    BOOST_CHECK_EQUAL(plat.ownWells.front().name, "W1");
    BOOST_REQUIRE_EQUAL(plat.thpWells.size(), 1U);
    BOOST_CHECK_EQUAL(plat.thpWells.front().name, "W2");
    BOOST_CHECK_CLOSE(plat.thpWells.front().weight, 1.0, 1e-6);
}
