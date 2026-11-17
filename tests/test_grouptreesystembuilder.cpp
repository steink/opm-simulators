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

// Part 4's builder: GroupTreeSystem::populateFromFlatNetwork() translates
// ProdGroupTreeBalancer::FlatNetworkInput (Part 1's output, group-tree-only,
// no network topology) into GroupTreeSystem's own wells_/activeNodes_ (Part
// 2's Newton system). These tests build FlatNetworkInput by hand -- Part 1's
// own tests already cover extractFlatNetworkInput() itself -- and check that
// the translator resolves names into the right kind of reference (Pinned vs.
// Thp well, own_wells vs. member_wells vs. active_children) and, for the
// nested-active-group case, reproduces the same physics already validated
// directly against GroupTreeSystem's own API in test_networkgrouptreesystem.cpp.

#include <config.h>

#define BOOST_TEST_MODULE GroupTreeSystemBuilderTests

#include <opm/simulators/wells/NetworkGroupTreeSystem.hpp>
#include <opm/simulators/wells/ProdGroupTreeBalancer.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Schedule/VFPProdTable.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <array>
#include <string>
#include <unordered_map>

using namespace Opm;
using namespace Opm::NetworkSolve;

namespace {

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

using Sys = GroupTreeSystem<double>;
using Flat = ProdGroupTreeBalancer::FlatNetworkInput<double>;
using FlatNode = ProdGroupTreeBalancer::FlatActiveNode<double>;

} // namespace

// A single Individual well, pinned at its own ORAT target: the builder must
// invert that target into fixed_q via the well's own ipr, exactly the way
// bhpFromTarget already does for WellKind::Group -- checked here against a
// direct call to the same inversion, not a hand-derived number, since the
// inversion itself is already covered elsewhere.
BOOST_FIXTURE_TEST_CASE(builder_pins_an_individual_well_at_its_target, Fixture)
{
    Sys sys(props);
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/3, /*efficiency=*/1.0});
    sys.setTerminalPressure(10.0 * unit::barsa);

    const double m3d = unit::cubic(unit::meter) / unit::day;
    const double target = 400.0 * m3d;

    FlatNode w;
    w.name = "W1";
    w.type = ProdNodeType::Well;
    w.mode = Well::ProducerCMode::ORAT;
    w.target = target;
    const Flat flat{w};

    Sys::WellNetworkData wd;
    wd.node = 1;
    wd.efficiency = 0.95;
    wd.vfp_table = 3;
    wd.ipr_b[0] = -target / (60.0 * unit::barsa);
    wd.ipr_a[0] = target - wd.ipr_b[0] * (25.0 * unit::barsa);
    const std::unordered_map<std::string, Sys::WellNetworkData> wellData{{"W1", wd}};

    sys.populateFromFlatNetwork(flat, wellData);
    BOOST_REQUIRE_EQUAL(sys.numWells(), 1);
    BOOST_CHECK_EQUAL(sys.wells()[0].name, "W1");
    BOOST_CHECK(sys.wells()[0].kind == Sys::WellKind::Pinned);
    BOOST_CHECK_EQUAL(sys.wells()[0].node, 1);
    BOOST_CHECK_CLOSE(sys.wells()[0].efficiency, 0.95, 1e-9);
    BOOST_CHECK_CLOSE(sys.wells()[0].fixed_q[0], target, 1e-6);   // by construction of wd's ipr

    sys.finalize();
    const std::vector<double> guess{15.0 * unit::barsa};
    const auto result = NetworkSolve::solve(sys, guess, NetworkSolve::Parameters<double>{1e-7, 50}, FullStep{});
    BOOST_REQUIRE(result.converged);
    BOOST_REQUIRE_EQUAL(result.well_rate.size(), 1U);
    BOOST_CHECK_CLOSE(result.well_rate[0], target, 1e-6);
}

// target <= 0 (stopped, or not yet meaningful before Part 1b) must give a
// genuine full stop -- fixed_q all zero -- not bhpFromTarget's single-phase
// zero-crossing, which would leave the other phases flowing at whatever the
// pivot bhp happens to give them.
BOOST_FIXTURE_TEST_CASE(builder_fully_stops_a_zero_target_well, Fixture)
{
    Sys sys(props);
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/3, /*efficiency=*/1.0});
    sys.setTerminalPressure(10.0 * unit::barsa);

    FlatNode w;
    w.name = "W1";
    w.type = ProdNodeType::Well;
    w.mode = Well::ProducerCMode::ORAT;
    w.target = 0.0;
    const Flat flat{w};

    Sys::WellNetworkData wd;
    wd.node = 1;
    wd.vfp_table = 3;
    // Deliberately nonzero ipr: if the builder mishandled target <= 0 by
    // still calling bhpFromTarget, water/gas could come out nonzero even
    // though oil is forced to zero.
    const double m3d = unit::cubic(unit::meter) / unit::day;
    wd.ipr_b = {-5.0 * m3d / unit::barsa, -3.0 * m3d / unit::barsa, -1.0 * m3d / unit::barsa};
    wd.ipr_a = {100.0 * m3d, 80.0 * m3d, 20.0 * m3d};
    const std::unordered_map<std::string, Sys::WellNetworkData> wellData{{"W1", wd}};

    sys.populateFromFlatNetwork(flat, wellData);
    BOOST_REQUIRE_EQUAL(sys.numWells(), 1);
    for (int p = 0; p < 3; ++p) {
        BOOST_CHECK_EQUAL(sys.wells()[0].fixed_q[p], 0.0);
    }
}

// An Individual well flagged networkThp must become a live Thp well, not a
// pin: no fixed_q, a genuine bhp unknown tied to its own node's pressure.
BOOST_FIXTURE_TEST_CASE(builder_gives_a_network_thp_well_a_live_row, Fixture)
{
    Sys sys(props);
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/3, /*efficiency=*/1.0});
    const double terminal = 10.0 * unit::barsa;
    sys.setTerminalPressure(terminal);

    FlatNode w;
    w.name = "W1";
    w.type = ProdNodeType::Well;
    w.mode = Well::ProducerCMode::THP;
    w.networkThp = true;
    const Flat flat{w};

    const double m3d = unit::cubic(unit::meter) / unit::day;
    Sys::WellNetworkData wd;
    wd.node = 1;
    wd.vfp_table = 3;
    const double pivot = 300.0 * m3d;
    wd.ipr_b[0] = -pivot / (50.0 * unit::barsa);
    wd.ipr_a[0] = pivot - wd.ipr_b[0] * (18.0 * unit::barsa);
    const std::unordered_map<std::string, Sys::WellNetworkData> wellData{{"W1", wd}};

    sys.populateFromFlatNetwork(flat, wellData);
    BOOST_REQUIRE_EQUAL(sys.numWells(), 1);
    BOOST_CHECK(sys.wells()[0].kind == Sys::WellKind::Thp);

    sys.finalize();
    const std::vector<double> guess{15.0 * unit::barsa};
    const auto result = NetworkSolve::solve(sys, guess, NetworkSolve::Parameters<double>{1e-7, 50}, FullStep{});
    BOOST_REQUIRE(result.converged);

    // Self-consistency: converged bhp against a direct table call at the
    // converged node pressure (thp) and rate, same check style as
    // test_networkgrouptreesystem.cpp's own THP test.
    BOOST_REQUIRE_EQUAL(result.well_phase_rates.size(), 1U);
    const auto& q = result.well_phase_rates[0];
    const double expected_bhp = props.bhp(3, -q[1], -q[0], -q[2], result.node_pressure[1], 0.0, 0.0, 0.0, false);
    BOOST_CHECK_CLOSE(result.well_bhp[0], expected_bhp, 1e-6);
}

// The cross-mode nested-active-group case from
// test_networkgrouptreesystem.cpp's nested_active_node_with_different_modes_
// must_project_not_copy_target, this time built through the translator from a
// hand-built FlatNetworkInput (PLAT: ORAT, GP1: WRAT, referenced via
// activeChildren) instead of directly through GroupTreeSystem's own API --
// exercising ensureActive()'s activeChildren resolution and its distinction
// between own_wells (built fresh from FlatWellShare) and active_children
// (looked up by name, recursively built).
BOOST_FIXTURE_TEST_CASE(builder_resolves_a_nested_active_group_across_modes, Fixture)
{
    Sys sys(props);
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/3, /*efficiency=*/1.0});
    const double terminal = 10.0 * unit::barsa;
    sys.setTerminalPressure(terminal);

    const double m3d = unit::cubic(unit::meter) / unit::day;
    const double Tc_water = 200.0 * m3d;
    const double Tp_oil = 1000.0 * m3d;
    const double eff_child = 0.8;

    FlatNode gp1;
    gp1.name = "GP1";
    gp1.type = ProdNodeType::Group;
    gp1.mode = Well::ProducerCMode::WRAT;
    gp1.target = Tc_water;
    gp1.ownWells.push_back({"GP1_WELL", /*guideRate=*/150.0 * m3d, /*efficiency=*/1.0});

    FlatNode plat;
    plat.name = "PLAT";
    plat.type = ProdNodeType::Group;
    plat.mode = Well::ProducerCMode::ORAT;
    plat.target = Tp_oil;
    plat.ownWells.push_back({"PLAT_WELL", /*guideRate=*/400.0 * m3d, /*efficiency=*/1.0});
    plat.activeChildren.push_back({"GP1", eff_child});

    const Flat flat{gp1, plat};

    std::unordered_map<std::string, Sys::WellNetworkData> wellData;
    {
        Sys::WellNetworkData wd;
        wd.node = 1;
        wd.vfp_table = 3;
        wd.ipr_b[Sys::kWater] = -Tc_water / (50.0 * unit::barsa);
        wd.ipr_a[Sys::kWater] = Tc_water - wd.ipr_b[Sys::kWater] * (20.0 * unit::barsa);
        const double gp1_oil_at_pivot = 150.0 * m3d;   // deliberately unrelated to Tc_water
        wd.ipr_b[Sys::kOil] = -gp1_oil_at_pivot / (80.0 * unit::barsa);
        wd.ipr_a[Sys::kOil] = gp1_oil_at_pivot - wd.ipr_b[Sys::kOil] * (20.0 * unit::barsa);
        wellData["GP1_WELL"] = wd;
    }
    const double expected_plat_oil = Tp_oil - eff_child * 150.0 * m3d;   // 1000 - 0.8*150 = 880
    {
        Sys::WellNetworkData wd;
        wd.node = 1;
        wd.vfp_table = 3;
        wd.ipr_b[Sys::kOil] = -expected_plat_oil / (50.0 * unit::barsa);
        wd.ipr_a[Sys::kOil] = expected_plat_oil - wd.ipr_b[Sys::kOil] * (20.0 * unit::barsa);
        wellData["PLAT_WELL"] = wd;
    }

    sys.populateFromFlatNetwork(flat, wellData);
    BOOST_REQUIRE_EQUAL(sys.numActiveNodes(), 2);
    BOOST_REQUIRE_EQUAL(sys.numWells(), 2);

    sys.finalize();
    const std::vector<double> guess{15.0 * unit::barsa};
    const auto result = NetworkSolve::solve(sys, guess, NetworkSolve::Parameters<double>{1e-7, 50}, FullStep{});
    BOOST_REQUIRE(result.converged);

    BOOST_REQUIRE_EQUAL(result.well_phase_rates.size(), 2U);
    const auto& q_gp1 = (sys.wells()[0].name == "GP1_WELL") ? result.well_phase_rates[0] : result.well_phase_rates[1];
    const auto& q_plat = (sys.wells()[0].name == "PLAT_WELL") ? result.well_phase_rates[0] : result.well_phase_rates[1];

    BOOST_CHECK_CLOSE(q_gp1[Sys::kWater], Tc_water, 1e-6);
    // The bug the original (pre-translator) test caught: using GP1's target
    // (a water number) directly in PLAT's oil sum would give 840, not 880.
    BOOST_CHECK_CLOSE(q_plat[Sys::kOil] + eff_child * q_gp1[Sys::kOil], Tp_oil, 1e-6);
    BOOST_CHECK_CLOSE(q_plat[Sys::kOil], expected_plat_oil, 1e-6);
}

// A GSATPROD satellite group referenced via activeChildren: no ipr/vfp data
// at all (it is not in wellNetworkData -- the builder must not even look it
// up there), just a fixed three-phase rate that still has to reach PLAT's
// own oil-target row through activeNodeTotal(), exactly like a pinned well's
// target would, but carrying all three phases rather than one projected
// scalar (see FlatActiveNode::satelliteRates's own doc comment).
BOOST_FIXTURE_TEST_CASE(builder_resolves_a_satellite_group_as_a_fixed_source, Fixture)
{
    Sys sys(props);
    sys.addNode(Node{"N1", /*parent=*/0, /*vfp_table=*/3, /*efficiency=*/1.0});
    sys.setTerminalPressure(10.0 * unit::barsa);

    const double m3d = unit::cubic(unit::meter) / unit::day;
    const double Tp_oil = 1000.0 * m3d;
    const double eff_sat = 0.9;
    const double sat_oil = 150.0 * m3d, sat_water = 20.0 * m3d, sat_gas = 5.0 * m3d;

    FlatNode sat;
    sat.name = "SAT";
    sat.type = ProdNodeType::Group;
    sat.satelliteRates = std::array<double, 3>{sat_oil, sat_water, sat_gas};

    FlatNode plat;
    plat.name = "PLAT";
    plat.type = ProdNodeType::Group;
    plat.mode = Well::ProducerCMode::ORAT;
    plat.target = Tp_oil;
    plat.ownWells.push_back({"PLAT_WELL", /*guideRate=*/400.0 * m3d, /*efficiency=*/1.0});
    plat.activeChildren.push_back({"SAT", eff_sat});

    const Flat flat{sat, plat};

    std::unordered_map<std::string, Sys::WellNetworkData> wellData;
    const double expected_plat_oil = Tp_oil - eff_sat * sat_oil;   // 1000 - 0.9*150 = 865
    Sys::WellNetworkData wd;
    wd.node = 1;
    wd.vfp_table = 3;
    wd.ipr_b[0] = -expected_plat_oil / (50.0 * unit::barsa);
    wd.ipr_a[0] = expected_plat_oil - wd.ipr_b[0] * (20.0 * unit::barsa);
    wellData["PLAT_WELL"] = wd;
    // Deliberately no entry for "SAT" -- if the builder ever tried
    // wellNetworkData.at("SAT") this test would throw before reaching solve().

    sys.populateFromFlatNetwork(flat, wellData);
    BOOST_REQUIRE_EQUAL(sys.numActiveNodes(), 1);
    BOOST_REQUIRE_EQUAL(sys.numWells(), 2);

    sys.finalize();
    const std::vector<double> guess{15.0 * unit::barsa};
    const auto result = NetworkSolve::solve(sys, guess, NetworkSolve::Parameters<double>{1e-7, 50}, FullStep{});
    BOOST_REQUIRE(result.converged);

    BOOST_REQUIRE_EQUAL(result.well_phase_rates.size(), 2U);
    const auto& q_sat = (sys.wells()[0].name == "SAT") ? result.well_phase_rates[0] : result.well_phase_rates[1];
    const auto& q_plat = (sys.wells()[0].name == "PLAT_WELL") ? result.well_phase_rates[0] : result.well_phase_rates[1];

    // The satellite's own rate is exactly its fixed source, on every phase.
    BOOST_CHECK_CLOSE(q_sat[Sys::kOil], sat_oil, 1e-9);
    BOOST_CHECK_CLOSE(q_sat[Sys::kWater], sat_water, 1e-9);
    BOOST_CHECK_CLOSE(q_sat[Sys::kGas], sat_gas, 1e-9);
    // PLAT's own row correctly used it, efficiency-scaled.
    BOOST_CHECK_CLOSE(q_plat[Sys::kOil], expected_plat_oil, 1e-6);
    BOOST_CHECK_CLOSE(q_plat[Sys::kOil] + eff_sat * q_sat[Sys::kOil], Tp_oil, 1e-6);

    // The satellite has no network-topology role: node 0 (the terminal),
    // never matched by any node's own flow sum.
    BOOST_CHECK_EQUAL(sys.wells()[(sys.wells()[0].name == "SAT") ? 0 : 1].node, 0);
}
