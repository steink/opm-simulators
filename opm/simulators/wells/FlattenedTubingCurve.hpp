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

#ifndef OPM_FLATTENED_TUBING_CURVE_HEADER_INCLUDED
#define OPM_FLATTENED_TUBING_CURVE_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/VFPProdTable.hpp>

#include <opm/simulators/wells/VFPProdProperties.hpp>

#include <algorithm>
#include <array>
#include <vector>

namespace Opm {

/// $\widetilde\Lambda$ from groups_and_network_clean.md's "$\widetilde\Lambda$
/// and the cliff diagnosis" section: a well's own tubing curve (real $\Lambda$,
/// from its VFP table), flattened against its own constant-fraction IPR so
/// that the two curves have at most one crossing -- the mechanism that lets
/// NetworkGroupTreeSystem solve a Thp well's row to full Newton convergence
/// with no discrete branching, even for a well whose real tubing curve has a
/// low-flow liquid-loading hump.
///
/// Built once per well whenever its ipr_a/ipr_b (and so its own IPR line) are
/// (re)computed -- the same cadence Part 1b's stopped-well trial IPR already
/// runs at -- not once per Newton iterate. thp is a genuine table axis here
/// (the network node pressure driving it moves during the solve); wfr/gfr/alq
/// are fixed at construction, matching the "constant well phase fractions"
/// assumption this whole design already relies on for the well's own
/// guide-rate/VFP bookkeeping elsewhere.
template<class Scalar>
class FlattenedTubingCurve
{
public:
    static constexpr int kOil = 0, kWater = 1, kGas = 2;

    /// \param table     The well's own VFP table (its axes are reused as-is:
    ///                  discarded points get their *data* overwritten, never
    ///                  the axes themselves).
    /// \param props     For the actual table lookups during construction;
    ///                  never touched again afterwards -- bhp() below is a
    ///                  self-contained 2-D lookup over the flattened data.
    /// \param table_id  table's own number, as props.bhp() needs it.
    /// \param ipr_b     The well's own affine IPR's slope, q_p = ipr_a[p] +
    ///                  ipr_b[p]*bhp, [oil, water, gas] order, positive =
    ///                  production (matching NetworkGroupTreeSystem's own
    ///                  convention). Only the slope is needed here: it fixes
    ///                  both the well's phase-fraction direction and (via
    ///                  flo_per_unit_drop below) the IPR's own slope in
    ///                  FLO-space -- the flattening condition compares
    ///                  slopes only, never an absolute bhp, so ipr_a/
    ///                  bhp_shutin play no part in this class at all (the
    ///                  bhp_shutin cap itself lives in NetworkGroupTreeSystem).
    /// \param alq       Fixed artificial-lift value for every lookup.
    /// \param eps       Minimum required margin of Lambda' over IPR' (see
    ///                  the doc's derivation of the flattening condition).
    FlattenedTubingCurve(const VFPProdTable& table,
                        const VFPProdProperties<Scalar>& props,
                        const int table_id,
                        const std::array<Scalar, 3>& ipr_b,
                        const Scalar alq,
                        const Scalar eps)
        : flo_type_(table.getFloType())
        , thp_axis_(table.getTHPAxis())
        , flo_axis_(table.getFloAxis())
    {
        const std::size_t nt = thp_axis_.size();
        const std::size_t nf = flo_axis_.size();
        data_.assign(nt * nf, Scalar{0});
        first_undiscarded_flo_.assign(nt, Scalar{0});

        // IPR expressed as a function of FLO rather than bhp: at one
        // (SI-)unit of drawdown below bhp_shutin every phase rate is exactly
        // -ipr_b[p] (since q_p(bhp) = ipr_b[p]*(bhp - bhp_shutin) under the
        // constant-fraction assumption), and FLO is linear and homogeneous in
        // the phase rates (a sum/selection of them, whatever the table's own
        // FLO_TYPE is) -- so this one lookup gives dFLO/dbhp directly,
        // without this class ever needing to know FLO_TYPE's own weighting,
        // or bhp_shutin's actual value (it is a derivative, so where exactly
        // "one unit of drawdown" is taken from cancels out below).
        const std::array<Scalar, 3> unit_drop_q{-ipr_b[kOil], -ipr_b[kWater], -ipr_b[kGas]};
        const Scalar flo_per_unit_drop = floOf(flo_type_, unit_drop_q[kWater], unit_drop_q[kOil], unit_drop_q[kGas]);

        // flo_per_unit_drop <= 0 means this well's own phase mix does not
        // correspond to any positive FLO at all under this table's FLO_TYPE
        // (e.g. a phase combination the table's own rate definition simply
        // cannot see) -- a table/well mismatch, not an ordinary "this well
        // cannot flow at the current thp" (that case is handled per-thp-row
        // below, and by NetworkGroupTreeSystem's bhp_shutin cap regardless of
        // anything this class does).
        if (!(flo_per_unit_drop > Scalar{0})) {
            degenerate_ = true;
            return;
        }
        // ipr' in FLO-space: IPR(FLO) = bhp_shutin - FLO/flo_per_unit_drop.
        const Scalar ipr_prime = -Scalar{1} / flo_per_unit_drop;
        const Scalar threshold = ipr_prime + eps;

        for (std::size_t t = 0; t < nt; ++t) {
            buildRow(t, props, table_id, unit_drop_q, flo_per_unit_drop, alq, threshold);
        }
    }

    /// $\widetilde\Lambda(thp, q)$: same calling convention as
    /// NetworkGroupTreeSystem::tableBhp() -- \p q is [oil, water, gas],
    /// positive = production -- so a Thp well's row can switch between the
    /// two with a one-line change and nothing else. FLO is derived from \p q
    /// exactly as the real table's own FLO_TYPE would (floOf() duplicates
    /// VFPHelpers::getFlo()'s three-case switch rather than keeping the
    /// table itself alive just for this, the same tradeoff
    /// NetworkGroupTreeSystem::phaseWeights() already makes for its own
    /// small, table-independent duplicate of ProdGroupTreeBalancer's
    /// projectOnMode()). What follows is then a plain bilinear lookup over
    /// the flattened data; queries outside the table's own axis range are
    /// clamped to the boundary rather than extrapolated -- deliberately
    /// conservative for a first version; revisit if a real deck needs
    /// otherwise.
    Scalar bhp(const Scalar thp, const std::array<Scalar, 3>& q) const
    {
        const Scalar flo = floOf(flo_type_, q[kWater], q[kOil], q[kGas]);
        const auto [it0, iw] = bracket(thp_axis_, thp);
        const auto [if0, fw] = bracket(flo_axis_, flo);
        const Scalar v00 = at(it0, if0);
        const Scalar v01 = at(it0, if0 + 1);
        const Scalar v10 = at(it0 + 1, if0);
        const Scalar v11 = at(it0 + 1, if0 + 1);
        const Scalar v0 = v00 + fw * (v01 - v00);
        const Scalar v1 = v10 + fw * (v11 - v10);
        return v0 + iw * (v1 - v0);
    }

    /// The FLO value up to which this thp row's curve was bridged by a
    /// chord -- see the class-level doc comment's "at most one crossing"
    /// point. Beyond it, bhp() returns exactly what the real table would
    /// (no discarding happened there), so the post-convergence diagnosis can
    /// skip evaluating the real VFP table entirely whenever the converged
    /// FLO already clears this threshold. Linearly interpolated between the
    /// table's own thp rows, clamped at the ends.
    ///
    /// This assumes the low-flow liquid-loading hump is the *only* discarded
    /// region on a row -- the physically expected case this whole mechanism
    /// targets. bhp() itself has no such assumption (its own construction
    /// below discards as many separate intervals as the slope bound
    /// actually requires, anywhere on the curve); it is only this shortcut
    /// that would miss a second, later discarded interval, if a real deck's
    /// table ever had one.
    Scalar firstUndiscardedFlo(const Scalar thp) const
    {
        const auto [it0, iw] = bracket(thp_axis_, thp);
        const Scalar v0 = first_undiscarded_flo_[it0];
        const Scalar v1 = first_undiscarded_flo_[std::min(it0 + 1, thp_axis_.size() - 1)];
        return v0 + iw * (v1 - v0);
    }

    /// True if this well's own phase mix cannot be expressed as any positive
    /// FLO on this table at all -- see the constructor's own doc comment.
    /// Distinct from an ordinary per-thp-row "cannot flow at this thp",
    /// which bhp()/firstUndiscardedFlo() already handle gracefully (a flat
    /// row, always inside the "discarded" region) with no need for the
    /// caller to treat it specially.
    bool degenerate() const { return degenerate_; }

private:
    std::size_t idx(const std::size_t t, const std::size_t f) const { return t * flo_axis_.size() + f; }
    Scalar at(const std::size_t t, const std::size_t f) const { return data_[idx(t, f)]; }

    /// Duplicates VFPHelpers::getFlo()'s own three-case switch (a pure,
    /// sign-agnostic selection/sum of its arguments -- positive input gives
    /// positive FLO, matching this class's own convention throughout, with
    /// no separate sign handling needed here) so this class never has to
    /// keep the VFPProdTable itself alive past construction.
    static Scalar floOf(const VFPProdTable::FLO_TYPE type, const Scalar aqua, const Scalar liquid,
                        const Scalar vapour)
    {
        switch (type) {
        case VFPProdTable::FLO_TYPE::FLO_OIL: return liquid;
        case VFPProdTable::FLO_TYPE::FLO_LIQ: return aqua + liquid;
        case VFPProdTable::FLO_TYPE::FLO_GAS: return vapour;
        }
        return Scalar{0};
    }

    /// The grid interval containing \p v (clamped to the axis's own range)
    /// and the fractional weight within it, for a plain linear interpolation.
    static std::pair<std::size_t, Scalar> bracket(const std::vector<double>& axis, const Scalar v)
    {
        if (axis.size() == 1) {
            return {0, Scalar{0}};
        }
        if (v <= axis.front()) {
            return {0, Scalar{0}};
        }
        if (v >= axis.back()) {
            return {axis.size() - 2, Scalar{1}};
        }
        const auto it = std::upper_bound(axis.begin(), axis.end(), v);
        const std::size_t i1 = static_cast<std::size_t>(it - axis.begin());
        const std::size_t i0 = i1 - 1;
        const Scalar w = (v - axis[i0]) / (axis[i1] - axis[i0]);
        return {i0, w};
    }

    /// One thp row's worth of $\Lambda$ samples at the table's own FLO grid,
    /// then the greedy discard/bridge pass described in the class-level doc
    /// comment: a monotone stack of "kept" indices, extended left to right,
    /// popping back whenever even the widest chord so far still fails the
    /// slope bound. A point that never clears the bound (checked all the way
    /// back to the row's very first kept index) is simply left out and
    /// retried against the next candidate -- never removed from
    /// consideration by shrinking the search, since only kept.back() (not
    /// kept itself) ever gets compared against.
    void buildRow(const std::size_t t, const VFPProdProperties<Scalar>& props,
                 const int table_id, const std::array<Scalar, 3>& unit_drop_q, const Scalar flo_per_unit_drop,
                 const Scalar alq, const Scalar threshold)
    {
        const std::size_t nf = flo_axis_.size();
        std::vector<Scalar> lam(nf);
        for (std::size_t f = 0; f < nf; ++f) {
            const Scalar s = flo_axis_[f] / flo_per_unit_drop;
            // props.bhp() wants negative-for-production aqua/liquid/vapour
            // (the same convention GroupTreeSystem::tableBhp() already
            // negates for) -- unit_drop_q, like this whole class's own
            // ipr_b, is positive-for-production, so it is negated here, at
            // the one place the two conventions meet.
            lam[f] = props.bhp(table_id, -s * unit_drop_q[kWater], -s * unit_drop_q[kOil], -s * unit_drop_q[kGas],
                               thp_axis_[t], alq, Scalar{0}, Scalar{0}, false);
        }

        std::vector<std::size_t> kept{0};
        for (std::size_t f = 1; f < nf; ++f) {
            bool accepted = false;
            while (true) {
                const std::size_t i = kept.back();
                const Scalar slope = (lam[f] - lam[i]) / (flo_axis_[f] - flo_axis_[i]);
                if (slope > threshold) {
                    accepted = true;
                    break;
                }
                if (kept.size() == 1) {
                    break;   // f cannot be validly reached yet; leave it out, try f+1 next
                }
                kept.pop_back();
            }
            if (accepted) {
                kept.push_back(f);
            }
        }

        for (std::size_t s = 0; s + 1 < kept.size(); ++s) {
            const std::size_t i0 = kept[s];
            const std::size_t i1 = kept[s + 1];
            for (std::size_t f = i0; f <= i1; ++f) {
                const Scalar w = (flo_axis_[f] - flo_axis_[i0]) / (flo_axis_[i1] - flo_axis_[i0]);
                data_[idx(t, f)] = lam[i0] + w * (lam[i1] - lam[i0]);
            }
        }
        // Tail beyond the last kept point (including the whole row, if kept
        // never grew past its first entry): held flat, never guessed at --
        // see the class-level doc comment on why this is still safe (f' = 0
        // there, still clearing the bound against IPR's own negative slope).
        for (std::size_t f = kept.back(); f < nf; ++f) {
            data_[idx(t, f)] = lam[kept.back()];
        }

        first_undiscarded_flo_[t] = (kept.size() >= 2) ? flo_axis_[kept[1]] : flo_axis_.back();
    }

    VFPProdTable::FLO_TYPE flo_type_;
    std::vector<double> thp_axis_;
    std::vector<double> flo_axis_;
    std::vector<Scalar> data_;                     // [t][f], row-major
    std::vector<Scalar> first_undiscarded_flo_;    // per thp row
    bool degenerate_ = false;
};

} // namespace Opm

#endif // OPM_FLATTENED_TUBING_CURVE_HEADER_INCLUDED
