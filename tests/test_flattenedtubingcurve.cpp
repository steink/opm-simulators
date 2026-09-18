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

#define BOOST_TEST_MODULE FlattenedTubingCurveTests

#include <opm/simulators/wells/FlattenedTubingCurve.hpp>

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

using namespace Opm;

namespace {

// A single-thp-row table with a hand-designed low-flow "liquid loading" hump:
// bhp drops from 30 to 15 bar as flo rises from 10 to 30 m3/d (the cliff),
// then climbs normally (friction-dominated) to 25 and 40 bar at 40/50 m3/d.
// FLO_TYPE = OIL, so getFlo() just returns the oil argument directly --
// chosen so the well's own ipr (oil-only below) sweeps flo with no unit
// surprises.
const std::string kHumpedVfpProd = R"(
VFPPROD
     5     250.00      OIL        WCT         GOR         THP        GRAT      METRIC   BHP      /
       10.0  20.0  30.0  40.0  50.0 /
      20.00 /
      0.000 /
       100.0 /
        0.0 /
  1  1  1  1    30.0   20.0   15.0   25.0   40.0 /
)";

struct Fixture
{
    Fixture()
        : deck(Opm::Parser{}.parseString(kHumpedVfpProd))
        , table(deck["VFPPROD"].front(), /*gaslift_opt_active=*/false, UnitSystem{})
    {
        props.addTable(table);
    }
    const Opm::Deck deck;
    const VFPProdTable table;
    VFPProdProperties<double> props;
};

const double bar = unit::barsa;
const double m3d = unit::cubic(unit::meter) / unit::day;

} // namespace

// The hand-worked example from groups_and_network_clean.md's derivation:
// PI = 1 (m3/d)/bar, bhp_shutin = 50 bar, eps = 0.05 bar/(m3/d), so
// IPR'(FLO) = -1 bar/(m3/d) and the required threshold is -0.95 bar/(m3/d).
// Walking the table's own flo grid {10,20,30,40,50} against bhp
// {30,20,15,25,40}: the chord (10->20) has slope -1.0 (fails, discarded);
// (10->30) has slope -0.75 (clears -0.95, accepted); (30->40) is 1.0 and
// (40->50) is 1.5 (both clear it easily). So flo=20 is the only discarded
// point, bridged by the chord from (10,30) to (30,15), and flo=30 is where
// the real curve becomes trustworthy again.
BOOST_FIXTURE_TEST_CASE(low_flow_hump_gets_bridged_by_a_single_chord, Fixture)
{
    std::array<double, 3> ipr_b{};
    const double PI = 1.0 * m3d / bar;
    ipr_b[FlattenedTubingCurve<double>::kOil] = -PI;
    const double eps = 0.05 * bar / m3d;

    const FlattenedTubingCurve<double> tilde(table, props, table.getTableNum(), ipr_b, /*alq=*/0.0, eps);
    BOOST_REQUIRE(!tilde.degenerate());

    // The bridged interval's own endpoint: real data trusted again from here.
    BOOST_CHECK_CLOSE(tilde.firstUndiscardedFlo(20.0 * bar), 30.0 * m3d, 1e-6);

    // bhp() takes q in NetworkGroupTreeSystem's own [oil, water, gas]
    // convention; this table's FLO_TYPE is OIL, so oil alone is flo.
    auto q_oil = [](const double flo) { return std::array<double, 3>{flo, 0.0, 0.0}; };

    // Grid points either side of the discard are untouched.
    BOOST_CHECK_CLOSE(tilde.bhp(20.0 * bar, q_oil(10.0 * m3d)), 30.0 * bar, 1e-6);
    BOOST_CHECK_CLOSE(tilde.bhp(20.0 * bar, q_oil(30.0 * m3d)), 15.0 * bar, 1e-6);
    BOOST_CHECK_CLOSE(tilde.bhp(20.0 * bar, q_oil(40.0 * m3d)), 25.0 * bar, 1e-6);
    BOOST_CHECK_CLOSE(tilde.bhp(20.0 * bar, q_oil(50.0 * m3d)), 40.0 * bar, 1e-6);

    // The discarded point itself: chord value, not the real (and lower) 20 bar.
    BOOST_CHECK_CLOSE(tilde.bhp(20.0 * bar, q_oil(20.0 * m3d)), 22.5 * bar, 1e-6);

    // Interior queries: linear between whatever data_ actually holds at the
    // bracketing grid points (bridged-to-real for the first, both-real for
    // the second) -- confirms bhp() itself does not special-case which.
    BOOST_CHECK_CLOSE(tilde.bhp(20.0 * bar, q_oil(25.0 * m3d)), 18.75 * bar, 1e-6);
    BOOST_CHECK_CLOSE(tilde.bhp(20.0 * bar, q_oil(45.0 * m3d)), 32.5 * bar, 1e-6);
}

// A well whose own phase mix the table's FLO_TYPE cannot see at all (this
// table is FLO_OIL; an all-water ipr has zero oil component) -- a table/well
// mismatch, not an ordinary "cannot flow at this thp". Distinct from the
// per-row case (an ordinary, expected outcome handled by bhp()'s own flat
// tail plus NetworkGroupTreeSystem's bhp_shutin cap, no flag needed).
BOOST_FIXTURE_TEST_CASE(phase_mix_the_table_cannot_see_is_degenerate, Fixture)
{
    std::array<double, 3> ipr_b{};
    const double PI = 1.0 * m3d / bar;
    ipr_b[FlattenedTubingCurve<double>::kWater] = -PI;

    const FlattenedTubingCurve<double> tilde(table, props, table.getTableNum(), ipr_b,
                                            /*alq=*/0.0, /*eps=*/0.05 * bar / m3d);
    BOOST_CHECK(tilde.degenerate());
}
