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

#define BOOST_TEST_MODULE NetworkGroupTreeSystemTests

#include <opm/simulators/wells/NetworkGroupTreeSystem.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Schedule/VFPProdTable.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <fmt/format.h>

#include <array>
#include <cmath>
#include <string>

using namespace Opm;
using namespace Opm::NetworkSolve;

namespace {

// The same small, hand-checkable VFPPROD table used in test_networksolve.cpp
// (originally from test_networkpressure.cpp): table 3, LIQ/WCT/GOR/THP/GRAT,
// METRIC. FLO axis 20/100/1000/2000, THP axis 10/30 bar, WFR axis 0/0.5/1,
// GFR axis 100 (one point), ALQ axis 0 (one point).
const std::string kVfpProd = R"(
VFPPROD
     3     250.00      LIQ        WCT         GOR         THP        GRAT      METRIC   BHP      /
       20.0       100.0    1000.0     2000.0 /
      10.00      30.00 /
      0.000      0.5      1.0 /
       100.0 /
        0.0 /
  1  1  1  1    12.0   15.0   20.0   30.0 /
  1  2  1  1    13.0   16.0   21.0   31.0 /
  1  3  1  1    14.0   17.0   22.0   32.0 /
  2  1  1  1    32.0   35.0   40.0   50.0 /
  2  2  1  1    33.0   36.0   41.0   51.0 /
  2  3  1  1    34.0   37.0   42.0   52.0 /
)";

const NetworkSolve::Parameters<double> kParams{1e-7, 50};

struct Fixture
{
    Fixture()
        : deck(Opm::Parser{}.parseString(kVfpProd))
        , table(deck["VFPPROD"].front(), /*gaslift_opt_active=*/false, UnitSystem{})
    {
        props.addTable(table);
    }
    const Opm::Deck deck;
    const VFPProdTable table;
    VFPProdProperties<double> props;
};

} // namespace

// The simplest possible case: one node, one Pinned well (a fixed, already-
// known rate -- no group, no THP). Validates the node-pressure row and the
// VFP table wiring (tableBhp's water/oil/gas reordering in particular)
// against a direct call to the same table, with nothing else in the system
// to get in the way.
BOOST_FIXTURE_TEST_CASE(single_pinned_well, Fixture)
{
    GroupTreeSystem<double> sys(props);
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/3, /*efficiency=*/1.0});
    const double terminal = 10.0 * unit::barsa;
    sys.setTerminalPressure(terminal);

    const double oil_rate = 500.0 * unit::cubic(unit::meter) / unit::day;
    GroupTreeSystem<double>::Well w;
    w.name = "W1";
    w.node = 1;
    w.kind = GroupTreeSystem<double>::WellKind::Pinned;
    w.fixed_q = {oil_rate, 0.0, 0.0};   // oil, water, gas
    sys.addWell(w);

    const std::vector<double> guess{15.0 * unit::barsa};
    sys.finalize();
    const auto result = NetworkSolve::solve(sys, guess, kParams, FullStep{});
    BOOST_REQUIRE(result.converged);

    // A direct call to the same table, at the same rate and terminal
    // pressure, must give exactly the node pressure the system converged to.
    const double expected = props.bhp(3, /*aqua=*/0.0, /*liquid=*/-oil_rate, /*vapour=*/0.0,
                                      terminal, /*alq=*/0.0, 0.0, 0.0, false);
    BOOST_REQUIRE_EQUAL(result.node_pressure.size(), 2U);   // terminal + N1
    BOOST_CHECK_CLOSE(result.node_pressure[0], terminal, 1e-9);
    BOOST_CHECK_CLOSE(result.node_pressure[1], expected, 1e-6);
    BOOST_REQUIRE_EQUAL(result.well_rate.size(), 1U);
    BOOST_CHECK_CLOSE(result.well_rate[0], oil_rate, 1e-9);
}

// Same shape as single_pinned_well, but the node's own branch carries a
// nonzero ALQ (BRANPROP/GRUPNET) -- addNode()'s second argument, stored in
// node_alq_ and threaded into the node row's own tableBhp() call instead of
// a hardcoded zero. The shared Fixture's table has only one ALQ axis point
// (a query at any alq would silently ignore it and this test would pass by
// accident), so this one builds its own two-point-ALQ table specifically to
// tell "alq threaded through" apart from "alq silently dropped".
BOOST_AUTO_TEST_CASE(single_pinned_well_with_node_alq)
{
    const std::string vfp_two_alq_points = R"(
VFPPROD
     7     250.00      LIQ        WCT         GOR         THP        GRAT      METRIC   BHP      /
       500.0 /
       10.00 /
        0.000 /
       100.0 /
        0.0       10.0 /
  1  1  1  1    20.0 /
  1  1  1  2    30.0 /
)";
    const auto deck = Opm::Parser{}.parseString(vfp_two_alq_points);
    const VFPProdTable table(deck["VFPPROD"].front(), /*gaslift_opt_active=*/false, UnitSystem{});
    VFPProdProperties<double> local_props;
    local_props.addTable(table);

    GroupTreeSystem<double> sys(local_props);
    const double alq = 7.5;
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/7, /*efficiency=*/1.0}, alq);
    const double terminal = 10.0 * unit::barsa;
    sys.setTerminalPressure(terminal);

    const double oil_rate = 500.0 * unit::cubic(unit::meter) / unit::day;
    GroupTreeSystem<double>::Well w;
    w.name = "W1";
    w.node = 1;
    w.kind = GroupTreeSystem<double>::WellKind::Pinned;
    w.fixed_q = {oil_rate, 0.0, 0.0};
    sys.addWell(w);

    const std::vector<double> guess{15.0 * unit::barsa};
    sys.finalize();
    const auto result = NetworkSolve::solve(sys, guess, kParams, FullStep{});
    BOOST_REQUIRE(result.converged);

    const double expected_with_alq = local_props.bhp(7, 0.0, -oil_rate, 0.0, terminal, alq, 0.0, 0.0, false);
    const double expected_without_alq = local_props.bhp(7, 0.0, -oil_rate, 0.0, terminal, 0.0, 0.0, 0.0, false);
    BOOST_REQUIRE_EQUAL(result.node_pressure.size(), 2U);
    BOOST_CHECK_CLOSE(result.node_pressure[1], expected_with_alq, 1e-6);
    BOOST_CHECK(std::abs(result.node_pressure[1] - expected_without_alq) > 1e-3 * unit::barsa);
}

// One well, Group-controlled, tied to one Active node's own lambda -- the
// single-well case makes lambda closed-form (no THP member under this node to
// couple it to pressure) regardless of the well's own guide rate: whatever it
// is, a lone well gets the whole target. Exercises the ownWells guideRate/
// efficiency split (Part 1's fix) and the affine-IPR-inversion (bhpFromTarget)
// for the first time, not just the pass-through node row Test 1 covered.
BOOST_FIXTURE_TEST_CASE(single_group_controlled_well, Fixture)
{
    GroupTreeSystem<double> sys(props);
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/3, /*efficiency=*/1.0});
    const double terminal = 10.0 * unit::barsa;
    sys.setTerminalPressure(terminal);

    const double target = 300.0 * unit::cubic(unit::meter) / unit::day;
    GroupTreeSystem<double>::ActiveNode a;
    a.mode = Opm::Well::ProducerCMode::ORAT;
    a.target = target;
    const int active_idx = sys.addActiveNode(a);

    // ipr_a/ipr_b picked so the well can deliver the target at all (the
    // pivot point doesn't need to be anywhere in particular -- the test
    // checks internal consistency against a direct table call, not a
    // hand-computed bhp).
    GroupTreeSystem<double>::Well w;
    w.name = "W1";
    w.node = 1;
    w.kind = GroupTreeSystem<double>::WellKind::Group;
    w.active_node = active_idx;
    w.guide_rate = 150.0 * unit::cubic(unit::meter) / unit::day;   // different from target: lambda != 1
    w.ipr_b[0] = -target / (50.0 * unit::barsa);
    w.ipr_a[0] = target - w.ipr_b[0] * (20.0 * unit::barsa);       // q(oil) = target at bhp = 20 bar
    const int well_idx = sys.addWell(w);
    sys.activeNode(active_idx).own_wells.push_back(well_idx);

    const std::vector<double> guess{15.0 * unit::barsa};
    sys.finalize();
    const auto result = NetworkSolve::solve(sys, guess, kParams, FullStep{});
    BOOST_REQUIRE(result.converged);

    // The lone well under this Active node gets the whole target, whatever
    // its own guide rate is.
    BOOST_REQUIRE_EQUAL(result.well_rate.size(), 1U);
    BOOST_CHECK_CLOSE(result.well_rate[0], target, 1e-6);

    // Internal consistency: the converged node pressure must be exactly what
    // the same table gives directly for the rates the system converged to.
    BOOST_REQUIRE_EQUAL(result.well_phase_rates.size(), 1U);
    const auto& q = result.well_phase_rates[0];
    const double expected = props.bhp(3, /*aqua=*/-q[1], /*liquid=*/-q[0], /*vapour=*/-q[2],
                                      terminal, /*alq=*/0.0, 0.0, 0.0, false);
    BOOST_REQUIRE_EQUAL(result.node_pressure.size(), 2U);
    BOOST_CHECK_CLOSE(result.node_pressure[1], expected, 1e-6);
}

// One well, on the network's THP control -- the one kind not covered above,
// and the only one with a genuine unknown of its own (bhp) and real two-way
// pressure coupling: the well's own row needs its node's pressure (its thp)
// to get its bhp, and the node's own row needs the well's rate (from that
// same bhp) to get its pressure. Reusing the same VFP table for both the
// well's tubing curve and the node's branch curve is physically odd but
// numerically exactly what is needed to exercise that coupling -- nothing
// about the class cares what a table number "means".
BOOST_FIXTURE_TEST_CASE(single_thp_controlled_well, Fixture)
{
    GroupTreeSystem<double> sys(props);
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/3, /*efficiency=*/1.0});
    const double terminal = 10.0 * unit::barsa;
    sys.setTerminalPressure(terminal);

    GroupTreeSystem<double>::Well w;
    w.name = "W1";
    w.node = 1;
    w.kind = GroupTreeSystem<double>::WellKind::Thp;
    w.vfp_table = 3;
    const double q_pivot = 500.0 * unit::cubic(unit::meter) / unit::day;
    w.ipr_b[0] = -q_pivot / (50.0 * unit::barsa);
    w.ipr_a[0] = q_pivot - w.ipr_b[0] * (15.0 * unit::barsa);   // q(oil) = q_pivot at bhp = 15 bar
    sys.addWell(w);

    const std::vector<double> guess{15.0 * unit::barsa};
    sys.finalize();
    const auto result = NetworkSolve::solve(sys, guess, kParams, FullStep{});
    BOOST_REQUIRE(result.converged);
    BOOST_TEST_MESSAGE(fmt::format("converged in {} iterations, node_pressure={:.4f} bar, well_bhp={:.4f} bar",
                                    result.iterations, result.node_pressure[1] / unit::barsa,
                                    result.well_bhp[0] / unit::barsa));
    // The starting guess (15 bar for both node pressure and bhp) is not
    // already a solution -- at that guess the well's own row alone wants
    // bhp = tableBhp(3, thp=15bar, q(bhp=15bar)) != 15 bar, so a real Newton
    // step was needed, not just a lucky initial guess.
    BOOST_CHECK(std::abs(result.iterations) >= 1);
    BOOST_CHECK(std::abs(result.node_pressure[1] - guess[0]) > 1e-3 * unit::barsa
                || std::abs(result.well_bhp[0] - guess[0]) > 1e-3 * unit::barsa);

    // Both of the coupled equations must hold at the converged state: the
    // well's own bhp against the table at its node's (converged) pressure,
    // and that node's pressure against the table at the terminal, for the
    // rate the converged bhp implies.
    BOOST_REQUIRE_EQUAL(result.well_bhp.size(), 1U);
    BOOST_REQUIRE_EQUAL(result.well_phase_rates.size(), 1U);
    BOOST_REQUIRE_EQUAL(result.node_pressure.size(), 2U);
    const double bhp = result.well_bhp[0];
    const auto& q = result.well_phase_rates[0];
    BOOST_CHECK_CLOSE(q[0], w.ipr_a[0] + w.ipr_b[0] * bhp, 1e-9);   // ipr itself, sanity check

    const double expected_bhp = props.bhp(3, /*aqua=*/-q[1], /*liquid=*/-q[0], /*vapour=*/-q[2],
                                          result.node_pressure[1], /*alq=*/0.0, 0.0, 0.0, false);
    BOOST_CHECK_CLOSE(bhp, expected_bhp, 1e-6);

    const double expected_node_pressure = props.bhp(3, /*aqua=*/-q[1], /*liquid=*/-q[0], /*vapour=*/-q[2],
                                                    terminal, /*alq=*/0.0, 0.0, 0.0, false);
    BOOST_CHECK_CLOSE(result.node_pressure[1], expected_node_pressure, 1e-6);
}

// The safety net a Thp well needs before ~Lambda even exists: this well's ipr
// is so depleted (bhp_shutin = 5 bar) that no bhp in the table's own range
// (its lowest corner, thp=10/flo=20, is already 12 bar) can ever match the
// tubing curve -- there is no valid, positive-flow operating point at all.
// Left unguarded, Newton would chase a crossing that does not exist and push
// bhp past bhp_shutin (negative phase rates, meaningless for production).
// updateControls()/limitStep() must instead settle it exactly at bhp_shutin,
// zero rate, decoupled from thp -- the node has no table of its own here
// specifically so its row cannot mask the well's row settling correctly.
BOOST_FIXTURE_TEST_CASE(thp_well_capped_at_its_own_shut_in_bhp, Fixture)
{
    GroupTreeSystem<double> sys(props);
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/NoTable, /*efficiency=*/1.0});
    const double terminal = 10.0 * unit::barsa;
    sys.setTerminalPressure(terminal);

    GroupTreeSystem<double>::Well w;
    w.name = "W1";
    w.node = 1;
    w.kind = GroupTreeSystem<double>::WellKind::Thp;
    w.vfp_table = 3;
    const double bhp_shutin = 5.0 * unit::barsa;
    const double q_pivot = 500.0 * unit::cubic(unit::meter) / unit::day;
    w.ipr_b[0] = -q_pivot / (50.0 * unit::barsa);
    w.ipr_a[0] = -w.ipr_b[0] * bhp_shutin;   // q(bhp_shutin) = 0 exactly
    sys.addWell(w);

    const std::vector<double> guess{15.0 * unit::barsa};
    sys.finalize();
    const auto result = NetworkSolve::solve(sys, guess, kParams, FullStep{});
    BOOST_REQUIRE(result.converged);
    BOOST_TEST_MESSAGE(fmt::format("converged in {} iterations: bhp={:.4f} bar, rate={:.4f} m3/d",
                                    result.iterations, result.well_bhp[0] / unit::barsa,
                                    result.well_rate[0] / (unit::cubic(unit::meter) / unit::day)));

    BOOST_REQUIRE_EQUAL(result.well_bhp.size(), 1U);
    BOOST_REQUIRE_EQUAL(result.well_rate.size(), 1U);
    BOOST_CHECK_CLOSE(result.well_bhp[0], bhp_shutin, 1e-6);
    BOOST_CHECK_SMALL(result.well_rate[0], 1e-6 * q_pivot);
}

// The case the whole "keep lambda explicit" design point is about: one
// Active node with both a group-controlled well (tied to its lambda) and a
// THP-controlled well (a genuine unknown, pressure-dependent) under it. Their
// rates are coupled *through* that shared lambda -- checked directly via the
// node's own target equation, sum(both wells' oil) == target, which only
// holds if lambda genuinely adjusted for whatever the THP well's converged
// rate turned out to be.
BOOST_FIXTURE_TEST_CASE(mixed_group_and_thp_wells_share_one_lambda, Fixture)
{
    GroupTreeSystem<double> sys(props);
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/3, /*efficiency=*/1.0});
    const double terminal = 10.0 * unit::barsa;
    sys.setTerminalPressure(terminal);

    const double target = 800.0 * unit::cubic(unit::meter) / unit::day;
    GroupTreeSystem<double>::ActiveNode a;
    a.mode = Opm::Well::ProducerCMode::ORAT;
    a.target = target;
    const int active_idx = sys.addActiveNode(a);

    GroupTreeSystem<double>::Well g;
    g.name = "G1";
    g.node = 1;
    g.kind = GroupTreeSystem<double>::WellKind::Group;
    g.active_node = active_idx;
    g.guide_rate = 300.0 * unit::cubic(unit::meter) / unit::day;
    const double g_pivot = 400.0 * unit::cubic(unit::meter) / unit::day;
    g.ipr_b[0] = -g_pivot / (50.0 * unit::barsa);
    g.ipr_a[0] = g_pivot - g.ipr_b[0] * (20.0 * unit::barsa);
    const int g_idx = sys.addWell(g);

    GroupTreeSystem<double>::Well t;
    t.name = "T1";
    t.node = 1;
    t.kind = GroupTreeSystem<double>::WellKind::Thp;
    t.vfp_table = 3;
    t.efficiency = 0.9;   // exercised here too: efficiency scales this well's
                          // contribution to both the node and the active-node sum
    const double t_pivot = 300.0 * unit::cubic(unit::meter) / unit::day;
    t.ipr_b[0] = -t_pivot / (50.0 * unit::barsa);
    t.ipr_a[0] = t_pivot - t.ipr_b[0] * (18.0 * unit::barsa);
    const int t_idx = sys.addWell(t);

    sys.activeNode(active_idx).own_wells.push_back(g_idx);
    sys.activeNode(active_idx).member_wells.emplace_back(t_idx, t.efficiency);

    const std::vector<double> guess{15.0 * unit::barsa};
    sys.finalize();
    const auto result = NetworkSolve::solve(sys, guess, kParams, FullStep{});
    BOOST_REQUIRE(result.converged);
    BOOST_TEST_MESSAGE(fmt::format("converged in {} iterations: G1={:.2f} T1={:.2f} m3/d, node_pressure={:.4f} bar",
                                    result.iterations, result.well_rate[0] / (unit::cubic(unit::meter) / unit::day),
                                    result.well_rate[1] / (unit::cubic(unit::meter) / unit::day),
                                    result.node_pressure[1] / unit::barsa));

    // The active node's own target equation: this only holds if lambda
    // genuinely responded to the THP well's (pressure-dependent) rate.
    BOOST_REQUIRE_EQUAL(result.well_rate.size(), 2U);
    const double g1_oil = result.well_rate[0];
    const double t1_oil = result.well_rate[1];
    BOOST_CHECK_CLOSE(g1_oil + t.efficiency * t1_oil, target, 1e-6);

    // Node-pressure self-consistency, both wells' full phase triples,
    // efficiency-scaled (G1's efficiency is the default 1.0).
    BOOST_REQUIRE_EQUAL(result.well_phase_rates.size(), 2U);
    std::array<double, 3> q_node{};
    for (int p = 0; p < 3; ++p) {
        q_node[p] = result.well_phase_rates[0][p] + t.efficiency * result.well_phase_rates[1][p];
    }
    const double expected_node_pressure = props.bhp(3, -q_node[1], -q_node[0], -q_node[2],
                                                    terminal, 0.0, 0.0, 0.0, false);
    BOOST_CHECK_CLOSE(result.node_pressure[1], expected_node_pressure, 1e-6);

    // T1's own row, same self-consistency check as the single-THP-well test.
    const double expected_t1_bhp = props.bhp(3, -result.well_phase_rates[1][1], -result.well_phase_rates[1][0],
                                             -result.well_phase_rates[1][2], result.node_pressure[1],
                                             0.0, 0.0, 0.0, false);
    BOOST_CHECK_CLOSE(result.well_bhp[1], expected_t1_bhp, 1e-6);
}

// A nested Active node: PLAT (parent, target Tp) has its own well *and* a
// nested child GP1 (its own target Tc) beneath it. PLAT's row must pick up
// GP1's *target* -- a constant, efficiency-scaled -- not any of GP1's own
// wells or lambda: lambda_PLAT = (Tp - eff*Tc) / (PLAT's own guide rates)
// should hold exactly regardless of anything happening inside GP1, which is
// the whole point of activeChildren referencing by name/target rather than
// flattening a second binding constraint away.
BOOST_FIXTURE_TEST_CASE(nested_active_node_uses_the_childs_target_not_its_wells, Fixture)
{
    GroupTreeSystem<double> sys(props);
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/3, /*efficiency=*/1.0});
    const double terminal = 10.0 * unit::barsa;
    sys.setTerminalPressure(terminal);

    const double Tc = 300.0 * unit::cubic(unit::meter) / unit::day;
    const double Tp = 1000.0 * unit::cubic(unit::meter) / unit::day;
    const double eff_child = 0.9;

    GroupTreeSystem<double>::ActiveNode child;
    child.mode = Opm::Well::ProducerCMode::ORAT;
    child.target = Tc;
    const int child_idx = sys.addActiveNode(child);

    GroupTreeSystem<double>::ActiveNode parent;
    parent.mode = Opm::Well::ProducerCMode::ORAT;
    parent.target = Tp;
    const int parent_idx = sys.addActiveNode(parent);

    GroupTreeSystem<double>::Well wc;
    wc.name = "WC";
    wc.node = 1;
    wc.kind = GroupTreeSystem<double>::WellKind::Group;
    wc.active_node = child_idx;
    wc.guide_rate = 150.0 * unit::cubic(unit::meter) / unit::day;   // irrelevant: lone well gets the whole target
    const double wc_pivot = 400.0 * unit::cubic(unit::meter) / unit::day;
    wc.ipr_b[0] = -wc_pivot / (50.0 * unit::barsa);
    wc.ipr_a[0] = wc_pivot - wc.ipr_b[0] * (20.0 * unit::barsa);
    const int wc_idx = sys.addWell(wc);

    GroupTreeSystem<double>::Well wp;
    wp.name = "WP";
    wp.node = 1;
    wp.kind = GroupTreeSystem<double>::WellKind::Group;
    wp.active_node = parent_idx;
    wp.guide_rate = 400.0 * unit::cubic(unit::meter) / unit::day;
    const double wp_pivot = 900.0 * unit::cubic(unit::meter) / unit::day;
    wp.ipr_b[0] = -wp_pivot / (50.0 * unit::barsa);
    wp.ipr_a[0] = wp_pivot - wp.ipr_b[0] * (20.0 * unit::barsa);
    const int wp_idx = sys.addWell(wp);

    sys.activeNode(child_idx).own_wells.push_back(wc_idx);
    sys.activeNode(parent_idx).own_wells.push_back(wp_idx);
    sys.activeNode(parent_idx).active_children.emplace_back(child_idx, eff_child);

    const std::vector<double> guess{15.0 * unit::barsa};
    sys.finalize();
    const auto result = NetworkSolve::solve(sys, guess, kParams, FullStep{});
    BOOST_REQUIRE(result.converged);

    // GP1 (child) is a lone well: it gets its whole target, exactly, no
    // matter what PLAT (parent) is doing.
    BOOST_REQUIRE_EQUAL(result.well_rate.size(), 2U);
    const double wc_oil = result.well_rate[0];
    const double wp_oil = result.well_rate[1];
    BOOST_CHECK_CLOSE(wc_oil, Tc, 1e-6);

    // PLAT's own target equation: its own well plus eff_child * GP1's actual
    // oil total. Both modes are ORAT here, so that total is numerically the
    // same as GP1's target Tc -- this test alone cannot tell "uses GP1's
    // target" and "uses GP1's actual oil total" apart; see
    // nested_active_node_with_different_modes_must_project_not_copy_target
    // for the case where the two genuinely differ.
    BOOST_CHECK_CLOSE(wp_oil + eff_child * Tc, Tp, 1e-6);

    // Node-pressure self-consistency: both wells feed the same physical node,
    // each at its own (default 1.0) network efficiency -- eff_child never
    // enters the physical flow balance, only PLAT's own target bookkeeping.
    BOOST_REQUIRE_EQUAL(result.well_phase_rates.size(), 2U);
    std::array<double, 3> q_node{};
    for (int p = 0; p < 3; ++p) {
        q_node[p] = result.well_phase_rates[0][p] + result.well_phase_rates[1][p];
    }
    const double expected_node_pressure = props.bhp(3, -q_node[1], -q_node[0], -q_node[2],
                                                    terminal, 0.0, 0.0, 0.0, false);
    BOOST_CHECK_CLOSE(result.node_pressure[1], expected_node_pressure, 1e-6);
}

// The bug the previous test could not have caught: PLAT is ORAT, GP1 (nested
// beneath it) is WRAT. GP1's own row only pins its *water* total to its
// target -- it says nothing about GP1's oil total, which is a different
// number entirely (GP1's well has separate oil and water IPR coefficients,
// both evaluated at the one bhp GP1's water equation settles on). PLAT's own
// row must sum that actual oil total, not GP1's (water-valued) target --
// using the target directly would silently add a water-rate number into an
// oil-rate sum.
BOOST_FIXTURE_TEST_CASE(nested_active_node_with_different_modes_must_project_not_copy_target, Fixture)
{
    GroupTreeSystem<double> sys(props);
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/3, /*efficiency=*/1.0});
    const double terminal = 10.0 * unit::barsa;
    sys.setTerminalPressure(terminal);

    const double m3d = unit::cubic(unit::meter) / unit::day;
    const double Tc_water = 200.0 * m3d;    // GP1's own target: WATER
    const double Tp_oil = 1000.0 * m3d;     // PLAT's own target: OIL
    const double eff_child = 0.8;

    GroupTreeSystem<double>::ActiveNode child;
    child.mode = Opm::Well::ProducerCMode::WRAT;
    child.target = Tc_water;
    const int child_idx = sys.addActiveNode(child);

    GroupTreeSystem<double>::ActiveNode parent;
    parent.mode = Opm::Well::ProducerCMode::ORAT;
    parent.target = Tp_oil;
    const int parent_idx = sys.addActiveNode(parent);

    // GP1's one well: its own row ties WATER to lambda_child (bhp comes out
    // at 20 bar, the pivot, since a lone well gets the whole target exactly);
    // its OIL rate at that same bhp is a completely separate number, fixed by
    // separate ipr_a/ipr_b[oil] coefficients unrelated to the water ones.
    GroupTreeSystem<double>::Well wc;
    wc.name = "GP1_WELL";
    wc.node = 1;
    wc.kind = GroupTreeSystem<double>::WellKind::Group;
    wc.active_node = child_idx;
    wc.guide_rate = 150.0 * m3d;   // irrelevant, lone well under GP1
    wc.ipr_b[GroupTreeSystem<double>::kWater] = -Tc_water / (50.0 * unit::barsa);
    wc.ipr_a[GroupTreeSystem<double>::kWater] =
        Tc_water - wc.ipr_b[GroupTreeSystem<double>::kWater] * (20.0 * unit::barsa);
    const double wc_oil_at_pivot = 150.0 * m3d;   // deliberately unrelated to Tc_water
    wc.ipr_b[GroupTreeSystem<double>::kOil] = -wc_oil_at_pivot / (80.0 * unit::barsa);
    wc.ipr_a[GroupTreeSystem<double>::kOil] =
        wc_oil_at_pivot - wc.ipr_b[GroupTreeSystem<double>::kOil] * (20.0 * unit::barsa);
    const int wc_idx = sys.addWell(wc);

    // PLAT's own well: whatever is left of Tp_oil after GP1's (efficiency-
    // scaled) oil contribution.
    GroupTreeSystem<double>::Well wp;
    wp.name = "PLAT_WELL";
    wp.node = 1;
    wp.kind = GroupTreeSystem<double>::WellKind::Group;
    wp.active_node = parent_idx;
    wp.guide_rate = 400.0 * m3d;   // irrelevant, lone well under PLAT
    const double expected_wp_oil = Tp_oil - eff_child * wc_oil_at_pivot;   // 1000 - 0.8*150 = 880
    wp.ipr_b[GroupTreeSystem<double>::kOil] = -expected_wp_oil / (50.0 * unit::barsa);
    wp.ipr_a[GroupTreeSystem<double>::kOil] =
        expected_wp_oil - wp.ipr_b[GroupTreeSystem<double>::kOil] * (20.0 * unit::barsa);
    const int wp_idx = sys.addWell(wp);

    sys.activeNode(child_idx).own_wells.push_back(wc_idx);
    sys.activeNode(parent_idx).own_wells.push_back(wp_idx);
    sys.activeNode(parent_idx).active_children.emplace_back(child_idx, eff_child);

    const std::vector<double> guess{15.0 * unit::barsa};
    sys.finalize();
    const auto result = NetworkSolve::solve(sys, guess, kParams, FullStep{});
    BOOST_REQUIRE(result.converged);
    BOOST_REQUIRE_EQUAL(result.well_phase_rates.size(), 2U);
    const auto& q_wc = result.well_phase_rates[0];
    const auto& q_wp = result.well_phase_rates[1];
    BOOST_TEST_MESSAGE(fmt::format("GP1: water={:.2f} oil={:.2f}  PLAT: oil={:.2f} m3/d",
                                    q_wc[GroupTreeSystem<double>::kWater] / m3d,
                                    q_wc[GroupTreeSystem<double>::kOil] / m3d,
                                    q_wp[GroupTreeSystem<double>::kOil] / m3d));

    // GP1's own row: water, exactly its target (unaffected by the bug either way).
    BOOST_CHECK_CLOSE(q_wc[GroupTreeSystem<double>::kWater], Tc_water, 1e-6);
    // GP1's oil is a different number entirely -- the bug would never check
    // this value at all, but it is what PLAT's row needs.
    BOOST_CHECK_CLOSE(q_wc[GroupTreeSystem<double>::kOil], wc_oil_at_pivot, 1e-6);
    // The assertion the bug fails: with the old code, PLAT's row would have
    // solved wp_oil = Tp_oil - eff_child * Tc_water = 1000 - 0.8*200 = 840,
    // not the correct 1000 - 0.8*150 = 880.
    BOOST_CHECK_CLOSE(q_wp[GroupTreeSystem<double>::kOil], expected_wp_oil, 1e-6);
    BOOST_CHECK_CLOSE(q_wp[GroupTreeSystem<double>::kOil] + eff_child * q_wc[GroupTreeSystem<double>::kOil],
                      Tp_oil, 1e-6);
}

// Coverage for the case discussed after the fix above: PLAT (parent) and GP1
// (child) are on the *same* mode (both ORAT), but GP1's only member is a THP
// well, not a group-controlled one. A THP well's rate has a genuine nonzero
// pressure derivative; PLAT's row needs to see that derivative through GP1
// (activeNodeTotal recurses into member_wells exactly the same as own_wells, so
// it does, uniformly, regardless of mode -- this is what makes the fix from
// the previous test also correct here, not just for the cross-mode case).
// Unlike that test, a "read GP1's target directly" version would reach the
// *same* converged numbers here (GP1's own row still forces its actual oil
// total to its target either way), so this cannot demonstrate the old
// formulation was wrong by checking a value -- what it could have been wrong
// about is the quality of the local linearization Newton uses to get there,
// which is not something a converged-value assertion observes. This test is
// coverage for the configuration (it does converge, to the physically
// expected numbers), not a regression test for this specific point.
BOOST_FIXTURE_TEST_CASE(nested_active_node_with_a_thp_well_under_the_child, Fixture)
{
    GroupTreeSystem<double> sys(props);
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/3, /*efficiency=*/1.0});
    const double terminal = 10.0 * unit::barsa;
    sys.setTerminalPressure(terminal);

    const double m3d = unit::cubic(unit::meter) / unit::day;
    const double Tc = 300.0 * m3d;    // GP1's own target, ORAT -- same mode as PLAT
    const double Tp = 1000.0 * m3d;   // PLAT's own target, ORAT
    const double eff_child = 0.85;

    GroupTreeSystem<double>::ActiveNode child;
    child.mode = Opm::Well::ProducerCMode::ORAT;
    child.target = Tc;
    const int child_idx = sys.addActiveNode(child);

    GroupTreeSystem<double>::ActiveNode parent;
    parent.mode = Opm::Well::ProducerCMode::ORAT;
    parent.target = Tp;
    const int parent_idx = sys.addActiveNode(parent);

    // GP1 needs a group-controlled member too, not just the THP well: with
    // nothing tied to lambda_child, GP1's own row would have zero Jacobian
    // column for lambda_child (a real edge case this test tripped over, see
    // the discussion) -- a singular system, not something this test is about.
    GroupTreeSystem<double>::Well wt;
    wt.name = "GP1_THP_WELL";
    wt.node = 1;
    wt.kind = GroupTreeSystem<double>::WellKind::Thp;
    wt.vfp_table = 3;
    const double t_pivot = 120.0 * m3d;
    wt.ipr_b[GroupTreeSystem<double>::kOil] = -t_pivot / (50.0 * unit::barsa);
    wt.ipr_a[GroupTreeSystem<double>::kOil] =
        t_pivot - wt.ipr_b[GroupTreeSystem<double>::kOil] * (18.0 * unit::barsa);
    const int wt_idx = sys.addWell(wt);

    GroupTreeSystem<double>::Well wg;
    wg.name = "GP1_GROUP_WELL";
    wg.node = 1;
    wg.kind = GroupTreeSystem<double>::WellKind::Group;
    wg.active_node = child_idx;
    wg.guide_rate = 200.0 * m3d;
    const double wg_pivot = 250.0 * m3d;   // plenty of headroom over Tc - t_pivot
    wg.ipr_b[GroupTreeSystem<double>::kOil] = -wg_pivot / (50.0 * unit::barsa);
    wg.ipr_a[GroupTreeSystem<double>::kOil] =
        wg_pivot - wg.ipr_b[GroupTreeSystem<double>::kOil] * (20.0 * unit::barsa);
    const int wg_idx = sys.addWell(wg);

    GroupTreeSystem<double>::Well wp;
    wp.name = "PLAT_WELL";
    wp.node = 1;
    wp.kind = GroupTreeSystem<double>::WellKind::Group;
    wp.active_node = parent_idx;
    wp.guide_rate = 400.0 * m3d;   // irrelevant, lone well under PLAT
    const double expected_wp_oil = Tp - eff_child * Tc;   // 1000 - 0.85*300 = 745
    wp.ipr_b[GroupTreeSystem<double>::kOil] = -expected_wp_oil / (50.0 * unit::barsa);
    wp.ipr_a[GroupTreeSystem<double>::kOil] =
        expected_wp_oil - wp.ipr_b[GroupTreeSystem<double>::kOil] * (20.0 * unit::barsa);
    const int wp_idx = sys.addWell(wp);

    sys.activeNode(child_idx).member_wells.emplace_back(wt_idx, wt.efficiency);
    sys.activeNode(child_idx).own_wells.push_back(wg_idx);
    sys.activeNode(parent_idx).own_wells.push_back(wp_idx);
    sys.activeNode(parent_idx).active_children.emplace_back(child_idx, eff_child);

    const std::vector<double> guess{15.0 * unit::barsa};
    sys.finalize();
    const auto result = NetworkSolve::solve(sys, guess, kParams, FullStep{});
    BOOST_REQUIRE(result.converged);
    BOOST_REQUIRE_EQUAL(result.well_rate.size(), 3U);
    const double wt_oil = result.well_rate[0];
    const double wg_oil = result.well_rate[1];
    const double wp_oil = result.well_rate[2];
    BOOST_TEST_MESSAGE(fmt::format("GP1: thp={:.2f} group={:.2f}  PLAT={:.2f} m3/d",
                                    wt_oil / m3d, wg_oil / m3d, wp_oil / m3d));

    // GP1's own row: its group well fills in whatever the THP well doesn't
    // deliver, so the two together hit GP1's target exactly -- lambda_child
    // genuinely responds to the THP well's (pressure-dependent) rate here.
    BOOST_CHECK_CLOSE(wt_oil + wg_oil, Tc, 1e-6);
    // PLAT's row: its own well plus eff_child * GP1's actual total (== Tc).
    BOOST_CHECK_CLOSE(wp_oil, expected_wp_oil, 1e-6);
    BOOST_CHECK_CLOSE(wp_oil + eff_child * (wt_oil + wg_oil), Tp, 1e-6);
}

// Case B from the design discussion: GP1's *only* member is a THP well --
// own_wells is empty, so finalize() must drop GP1's row/lambda entirely (a
// naive implementation gives lambda_GP1 a zero Jacobian column: singular,
// result.converged == false). GP1's target is then unenforceable -- nothing
// left, network-side, can move its total toward it -- so unlike every other
// test here, GP1's THP well is *not* expected to land anywhere near Tc; it
// settles wherever its own pressure-coupled physics puts it. PLAT still
// needs GP1's actual (whatever that turns out to be) oil total, via
// activeNodeTotal(), which works the same whether or not GP1 has a row of
// its own.
BOOST_FIXTURE_TEST_CASE(active_node_with_empty_own_wells_has_its_row_dropped, Fixture)
{
    GroupTreeSystem<double> sys(props);
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/3, /*efficiency=*/1.0});
    const double terminal = 10.0 * unit::barsa;
    sys.setTerminalPressure(terminal);

    const double m3d = unit::cubic(unit::meter) / unit::day;
    const double Tc = 300.0 * m3d;    // GP1's target: unenforceable, kept only to show it is ignored
    const double Tp = 1000.0 * m3d;
    const double eff_child = 0.8;

    GroupTreeSystem<double>::ActiveNode child;
    child.mode = Opm::Well::ProducerCMode::ORAT;
    child.target = Tc;
    const int child_idx = sys.addActiveNode(child);

    GroupTreeSystem<double>::ActiveNode parent;
    parent.mode = Opm::Well::ProducerCMode::ORAT;
    parent.target = Tp;
    const int parent_idx = sys.addActiveNode(parent);

    GroupTreeSystem<double>::Well wt;
    wt.name = "GP1_ONLY_WELL";
    wt.node = 1;
    wt.kind = GroupTreeSystem<double>::WellKind::Thp;
    wt.vfp_table = 3;
    const double t_pivot = 150.0 * m3d;   // deliberately far from Tc=300
    wt.ipr_b[GroupTreeSystem<double>::kOil] = -t_pivot / (50.0 * unit::barsa);
    wt.ipr_a[GroupTreeSystem<double>::kOil] =
        t_pivot - wt.ipr_b[GroupTreeSystem<double>::kOil] * (18.0 * unit::barsa);
    const int wt_idx = sys.addWell(wt);

    GroupTreeSystem<double>::Well wp;
    wp.name = "PLAT_WELL";
    wp.node = 1;
    wp.kind = GroupTreeSystem<double>::WellKind::Group;
    wp.active_node = parent_idx;
    wp.guide_rate = 400.0 * m3d;
    // PLAT's own well doesn't know what GP1 will actually deliver either --
    // give it enough headroom to cover whatever's left of Tp.
    const double wp_pivot = 900.0 * m3d;
    wp.ipr_b[GroupTreeSystem<double>::kOil] = -wp_pivot / (50.0 * unit::barsa);
    wp.ipr_a[GroupTreeSystem<double>::kOil] =
        wp_pivot - wp.ipr_b[GroupTreeSystem<double>::kOil] * (20.0 * unit::barsa);
    const int wp_idx = sys.addWell(wp);

    sys.activeNode(child_idx).member_wells.emplace_back(wt_idx, wt.efficiency);
    sys.activeNode(parent_idx).own_wells.push_back(wp_idx);
    sys.activeNode(parent_idx).active_children.emplace_back(child_idx, eff_child);

    sys.finalize();
    // The point of finalize(): GP1 (empty own_wells) gets no lambda; only
    // PLAT does. One fewer unknown than a naive "one lambda per active node"
    // count would give.
    BOOST_CHECK_EQUAL(sys.numActiveNodes(), 2);
    BOOST_CHECK_EQUAL(sys.numActiveLambdas(), 1);
    BOOST_CHECK_EQUAL(sys.size(), sys.numNodes() + sys.numActiveLambdas() + 1 /*one thp well*/);

    const std::vector<double> guess{15.0 * unit::barsa};
    const auto result = NetworkSolve::solve(sys, guess, kParams, FullStep{});
    BOOST_REQUIRE(result.converged);
    BOOST_REQUIRE_EQUAL(result.well_rate.size(), 2U);
    const double gp1_oil = result.well_rate[0];
    const double plat_oil = result.well_rate[1];
    BOOST_TEST_MESSAGE(fmt::format("GP1(thp, unenforced target {:.0f})={:.2f}  PLAT={:.2f} m3/d",
                                    Tc / m3d, gp1_oil / m3d, plat_oil / m3d));

    // GP1's target is not met -- nothing here enforces it any more.
    BOOST_CHECK(std::abs(gp1_oil - Tc) > 1.0 * m3d);
    // PLAT's row still holds, using GP1's actual (unenforced) total.
    BOOST_CHECK_CLOSE(plat_oil + eff_child * gp1_oil, Tp, 1e-6);

    // Node-pressure self-consistency, as in every other test.
    BOOST_REQUIRE_EQUAL(result.well_phase_rates.size(), 2U);
    std::array<double, 3> q_node{};
    for (int p = 0; p < 3; ++p) {
        q_node[p] = result.well_phase_rates[0][p] + result.well_phase_rates[1][p];
    }
    const double expected_node_pressure = props.bhp(3, -q_node[1], -q_node[0], -q_node[2],
                                                    terminal, 0.0, 0.0, 0.0, false);
    BOOST_CHECK_CLOSE(result.node_pressure[1], expected_node_pressure, 1e-6);
}

// Coverage for member_wells holding a Pinned well, not just a Thp one -- every
// other test above only ever exercised this list with WellKind::Thp. A Pinned
// well is an ordinary Individual well referenced by its parent (Part 1's
// activeChildren, for a type == Well entry with no networkThp flag): its
// contribution to the parent's sum is its fixed rate, efficiency-scaled, with
// no bhp unknown and no row of its own -- wellPhaseRatesOwn() already returns
// fixed_q directly for WellKind::Pinned, so activeNodeTotal()'s member_wells
// loop needs no kind-specific handling to get this right.
BOOST_FIXTURE_TEST_CASE(pinned_well_counts_toward_its_parents_sum_via_member_wells, Fixture)
{
    GroupTreeSystem<double> sys(props);
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/3, /*efficiency=*/1.0});
    const double terminal = 10.0 * unit::barsa;
    sys.setTerminalPressure(terminal);

    const double m3d = unit::cubic(unit::meter) / unit::day;
    const double Tp = 1000.0 * m3d;
    const double fixed_oil = 200.0 * m3d;
    const double eff_pinned = 0.75;

    GroupTreeSystem<double>::ActiveNode plat;
    plat.mode = Opm::Well::ProducerCMode::ORAT;
    plat.target = Tp;
    const int plat_idx = sys.addActiveNode(plat);

    GroupTreeSystem<double>::Well w_fixed;
    w_fixed.name = "W_FIXED";
    w_fixed.node = 1;
    w_fixed.kind = GroupTreeSystem<double>::WellKind::Pinned;
    w_fixed.fixed_q = {fixed_oil, 0.0, 0.0};
    const int w_fixed_idx = sys.addWell(w_fixed);

    GroupTreeSystem<double>::Well wp;
    wp.name = "PLAT_WELL";
    wp.node = 1;
    wp.kind = GroupTreeSystem<double>::WellKind::Group;
    wp.active_node = plat_idx;
    wp.guide_rate = 400.0 * m3d;
    const double expected_wp_oil = Tp - eff_pinned * fixed_oil;   // 1000 - 0.75*200 = 850
    wp.ipr_b[0] = -expected_wp_oil / (50.0 * unit::barsa);
    wp.ipr_a[0] = expected_wp_oil - wp.ipr_b[0] * (20.0 * unit::barsa);
    const int wp_idx = sys.addWell(wp);

    sys.activeNode(plat_idx).own_wells.push_back(wp_idx);
    sys.activeNode(plat_idx).member_wells.emplace_back(w_fixed_idx, eff_pinned);

    const std::vector<double> guess{15.0 * unit::barsa};
    sys.finalize();
    const auto result = NetworkSolve::solve(sys, guess, kParams, FullStep{});
    BOOST_REQUIRE(result.converged);

    BOOST_REQUIRE_EQUAL(result.well_rate.size(), 2U);
    const double w_fixed_oil = result.well_rate[0];
    const double wp_oil = result.well_rate[1];
    BOOST_CHECK_CLOSE(w_fixed_oil, fixed_oil, 1e-9);   // pinned: never moves
    BOOST_CHECK_CLOSE(wp_oil, expected_wp_oil, 1e-6);
    BOOST_CHECK_CLOSE(wp_oil + eff_pinned * w_fixed_oil, Tp, 1e-6);

    BOOST_REQUIRE_EQUAL(result.well_phase_rates.size(), 2U);
    std::array<double, 3> q_node{};
    for (int p = 0; p < 3; ++p) {
        q_node[p] = result.well_phase_rates[0][p] + result.well_phase_rates[1][p];
    }
    const double expected_node_pressure = props.bhp(3, -q_node[1], -q_node[0], -q_node[2],
                                                    terminal, 0.0, 0.0, 0.0, false);
    BOOST_CHECK_CLOSE(result.node_pressure[1], expected_node_pressure, 1e-6);
}

namespace {

// A single-thp-row table with a hand-designed low-flow "liquid loading" hump
// (the same one test_flattenedtubingcurve.cpp derives its own expectations
// against): bhp drops from 30 to 20 to 15 bar as flo rises 10->20->30 m3/d,
// then climbs normally to 25 and 40 bar at 40/50 m3/d.
const std::string kHumpedVfpProd = R"(
VFPPROD
     7     250.00      OIL        WCT         GOR         THP        GRAT      METRIC   BHP      /
       10.0  20.0  30.0  40.0  50.0 /
      20.00 /
      0.000 /
       100.0 /
        0.0 /
  1  1  1  1    30.0   20.0   15.0   25.0   40.0 /
)";

} // namespace

// Proof that Well::ipr_slope_limit actually changes what Newton converges to,
// not just that it compiles: the node has no table of its own and the
// terminal pressure sits exactly at the humped table's one thp row (20 bar),
// so the well's own bhp equation is solved against, respectively, the raw
// table and its slope-limited form, with nothing else in the system to blur
// the comparison.
//
// The table's one thp row, in bar against flo in m3/d, is
//   (10,30) (20,20) (30,15) (40,25) (50,40)
// so its segment slopes are -1, -0.5, +1, +1.5 bar/(m3/d): a liquid-loading
// branch falling steeply out of flo=10, a minimum near flo=30, then the
// friction-dominated rise.
//
// PI = 0.96 (m3/d)/bar puts the IPR's own slope at -1/0.96 = -1.0417
// bar/(m3/d), so the slope limit is 0.95 * that = -0.9896: the [10,20]
// segment (-1) is steeper than that and gets flattened, [20,30] (-0.5) and
// everything above it does not. A query on [10,20] therefore extrapolates
// [20,30] backwards, which is the line bhp = 30 - 0.5*flo -- and that line
// *continues* [20,30] itself, so the slope-limited curve is one straight
// segment over the whole of [10,30].
//
// bhp_shutin = 40.6 bar is then picked so that each curve has exactly *one*
// crossing with IPR(flo) = 40.6 - flo/0.96, which keeps the comparison from
// depending on which root Newton happens to fall into:
//   - against the raw table, only the falling [10,20] branch is crossed, at
//     (flo=14.4, bhp=25.6) -- an unstable operating point, the very thing the
//     slope limit exists to keep a solve away from;
//   - against the slope-limited curve, only bhp = 30 - 0.5*flo is crossed, at
//     flo = 254.4/13, bhp = 262.8/13.
// Both roots sit inside [10,20], i.e. inside the flattened stretch, so the
// limit is genuinely in play at the solution and not merely on the way there.
BOOST_AUTO_TEST_CASE(slope_limit_moves_the_thp_wells_converged_point_off_the_cliff)
{
    const auto deck = Opm::Parser{}.parseString(kHumpedVfpProd);
    const VFPProdTable table(deck["VFPPROD"].front(), /*gaslift_opt_active=*/false, UnitSystem{});
    VFPProdProperties<double> humped_props;
    humped_props.addTable(table);

    const double m3d = unit::cubic(unit::meter) / unit::day;
    const double PI = 0.96 * m3d / unit::barsa;
    const double bhp_shutin = 40.6 * unit::barsa;

    GroupTreeSystem<double>::Well w;
    w.name = "W1";
    w.node = 1;
    w.kind = GroupTreeSystem<double>::WellKind::Thp;
    w.vfp_table = table.getTableNum();
    w.ipr_b[GroupTreeSystem<double>::kOil] = -PI;
    w.ipr_a[GroupTreeSystem<double>::kOil] = PI * bhp_shutin;

    auto buildAndSolve = [&](const bool use_slope_limit) {
        GroupTreeSystem<double> sys(humped_props);
        sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/NoTable, /*efficiency=*/1.0});
        sys.setTerminalPressure(20.0 * unit::barsa);   // exactly the table's one thp row

        auto well = w;
        if (use_slope_limit) {
            // The same limit populateFromFlatNetwork() would have given it.
            well.ipr_slope_limit = sys.iprSlopeLimit(well);
            BOOST_REQUIRE(well.ipr_slope_limit.has_value());
            BOOST_CHECK_CLOSE(*well.ipr_slope_limit, (-0.95 / 0.96) * unit::barsa / m3d, 1e-8);
        }
        sys.addWell(well);

        const std::vector<double> guess{20.0 * unit::barsa};
        sys.finalize();
        return NetworkSolve::solve(sys, guess, kParams, FullStep{});
    };

    const auto real = buildAndSolve(false);
    const auto flat = buildAndSolve(true);
    BOOST_REQUIRE(real.converged);
    BOOST_REQUIRE(flat.converged);

    BOOST_REQUIRE_EQUAL(real.well_bhp.size(), 1U);
    BOOST_REQUIRE_EQUAL(flat.well_bhp.size(), 1U);
    BOOST_TEST_MESSAGE(fmt::format("real: bhp={:.4f} bar oil={:.4f} m3/d   limited: bhp={:.4f} bar oil={:.4f} m3/d",
                                    real.well_bhp[0] / unit::barsa, real.well_rate[0] / m3d,
                                    flat.well_bhp[0] / unit::barsa, flat.well_rate[0] / m3d));

    // Against the raw table: the unstable crossing on the falling branch.
    BOOST_CHECK_CLOSE(real.well_bhp[0], 25.6 * unit::barsa, 1e-4);
    BOOST_CHECK_CLOSE(real.well_rate[0], 14.4 * m3d, 1e-4);

    // Against the slope-limited curve: a genuinely different point, off that
    // branch -- proof the limit is what actually got solved against, not the
    // raw table silently reused.
    BOOST_CHECK_CLOSE(flat.well_bhp[0], (262.8 / 13.0) * unit::barsa, 1e-4);
    BOOST_CHECK_CLOSE(flat.well_rate[0], (254.4 / 13.0) * m3d, 1e-4);
}
