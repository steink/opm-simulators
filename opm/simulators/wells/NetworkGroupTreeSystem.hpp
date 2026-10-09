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
#ifndef OPM_NETWORK_GROUP_TREE_SYSTEM_HEADER_INCLUDED
#define OPM_NETWORK_GROUP_TREE_SYSTEM_HEADER_INCLUDED

#include <opm/simulators/wells/NetworkSolve.hpp>
#include <opm/simulators/wells/ProdGroupTreeBalancer.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace Opm::NetworkSolve {

/// The compact, group-tree-driven production network system from
/// groups_and_network_clean.md (Part 2 of the implementation plan): the
/// balancer (ProdGroupTreeBalancer) decides which wells are on group control,
/// individually limited, or on THP, once per outer iterate; this system has
/// no active-set state of its own at all -- updateControls() never moves
/// anything -- it just assembles and solves the resulting smooth system for
/// node pressures (and, for THP wells, bhp).
///
/// First version: "unproblematic" wells only, i.e. no stopped-well trial IPR
/// (Part 1b) yet. Lift-cliff handling (Part 3) is opt-in per Thp well: one
/// with a Well::ipr_slope_limit set solves its row against the slope-limited
/// tubing curve (VFPProdProperties::bhp_with_slope_limit()) instead of the
/// real table -- see thpWellResidualRow(). No autochoke either. Every well
/// has one of three kinds:
///  - Pinned:  a fixed, already-known rate (an Individual well, its own limit
///             binds) -- no unknown of its own at all.
///  - Group:   tied to an Active node's own lambda, q_target_mode = guideRate
///             * lambda (see ProdGroupTreeBalancer::FlatWellShare for why
///             guideRate and efficiency have to stay separate here).
///  - Thp:     its own bhp is a genuine unknown, tied to the network node's
///             pressure (its thp) via the real tubing curve.
/// One lambda unknown per Active node (see the plan's fill-in note: this is
/// deliberate, not something to eliminate -- a node with both Group and Thp
/// members needs it to stay sparse).
template<class Scalar>
class GroupTreeSystem : public SystemBase<Scalar>
{
public:
    using State = std::vector<Scalar>;
    using ScalarType = Scalar;
    static constexpr int NP = 3;                    // oil, water, gas
    static constexpr int kOil = 0, kWater = 1, kGas = 2;   // ProdGroupTreeBalancer's order

    enum class WellKind { Pinned, Group, Thp };

    struct Well
    {
        std::string name;
        int node = 0;                     // index into nodes_: this well's own network node
        WellKind kind = WellKind::Pinned;
        Scalar efficiency = 1;            // WEFAC as it applies to the network branch

        std::array<Scalar, NP> fixed_q{}; // Pinned only: the well's own (pre-efficiency) rate

        int active_node = -1;             // Group only: index into activeNodes_
        Scalar guide_rate = 0;            // Group only: q_target_mode = guide_rate * lambda

        // Group and Thp both need the tubing/IPR data to fill in bhp and the
        // other two phases; Pinned does not (its q is already fully known).
        int vfp_table = NoTable;
        Scalar alq = 0;
        // The tubing table's datum depth vs. the well's reference depth: the
        // well's bhp is the table's value minus this (as in the well's own
        // solve and NetworkProductionSystem). Thp wells only.
        Scalar vfp_dp = 0;
        std::array<Scalar, NP> ipr_a{};
        std::array<Scalar, NP> ipr_b{};   // q_p = ipr_a[p] + ipr_b[p] * bhp

        // Thp only: the bhp at which every phase's IPR rate is exactly zero
        // (the reservoir's own shut-in pressure) -- set by finalize(), not the
        // caller. A Newton step is never allowed to push this well's bhp past
        // it (limitStep()): beyond this point every phase rate is negative,
        // which is not a real operating point for a production well and is
        // not something the tubing table is meaningful for either. Once the
        // *current* iterate sits at or past it, updateControls() switches
        // this well's own row to bhp - bhp_shutin = 0 (decoupled from thp)
        // instead of the usual tubing-curve equation -- see residual().
        Scalar bhp_shutin = 0;

        // Thp only, and optional even then: Part 3's cliff protection, as the
        // steepest d(bhp)/d(FLO) this well's own tubing curve is allowed to
        // have before it is flattened -- the well's own IPR slope in that
        // same (table-oriented, positive-FLO) sense, plus a margin. Passed
        // straight to VFPProdProperties::bhp_with_slope_limit(); see
        // detail::SlopeLimit for what the flattening does and why.
        //
        // Unset means "use the real table directly, unflattened" -- the only
        // behaviour before this field existed, still exercised by every
        // earlier Thp test in test_networkgrouptreesystem.cpp, and still a
        // legitimate choice for a well known not to have a lift cliff. See
        // thpWellResidualRow().
        std::optional<Scalar> ipr_slope_limit;

        // Thp only: a stopped well offered for reopening (addReopenCandidate()).
        // Its ipr_a/b are its trial IPR; it starts at bhp_shutin (rate 0), so
        // the solve itself decides whether its node pressure lets it flow --
        // see reopenOutcomes(). Never ranked by worstCliffViolation(): a
        // candidate on the flattened curve is simply not reopened.
        bool reopen_candidate = false;
        // Thp wells: where it starts in the solve (0 if unknown) -- a flowing
        // well's current bhp, a reopen candidate's trial-IPR bhp. A candidate
        // starts there rather than at bhp_shutin: for a tubing curve with a
        // minimum, the curve at zero rate needs more than the shut-in bhp, so
        // a start at the cap would never leave it.
        Scalar start_bhp = 0;

        // Set by applyPhaseShare(): the water and gas fractions (the table's
        // WFR and GFR) of the well's fixed phase proportions, which its tubing
        // lookups use at every FLO. Unset: the fractions of the iterate's own
        // rates, and explicit_fractions below the table's first FLO value (as
        // in the well's own solve).
        std::optional<std::array<Scalar, 2>> fractions;
        std::array<Scalar, 2> explicit_fractions{};

        // Rates (positive, [oil, water, gas]) that stand for this well when a
        // node's fallback fractions are mixed from its wells -- see finalize().
        std::array<Scalar, NP> reference_q{};
    };

    /// One Active node from ProdGroupTreeBalancer::extractFlatNetworkInput():
    /// its own lambda and target, and the (already name-to-index-resolved)
    /// members contributing to its sum.
    struct ActiveNode
    {
        std::string name;
        ::Opm::Well::ProducerCMode mode{::Opm::Well::ProducerCMode::CMODE_UNDEFINED};
        Scalar target = 0;
        std::array<Scalar, NP> resv_coeff{};             // RESV mode only
        std::vector<int> own_wells;                      // indices into wells_, WellKind::Group

        // Standalone wells (Pinned or Thp) whose own row, if they have one,
        // lives elsewhere -- keyed globally by kind in residual(), never per
        // ActiveNode -- but that still count toward this node's sum, via
        // efficiency alone rather than a guide-rate allocation of this node's
        // own lambda. A Pinned well here is an Individual well whose own
        // limit already fixes it (or a confirmed-stopped one); a Thp well is
        // on the network's live THP control. Either way activeNodeTotal()
        // treats them identically: wellPhaseRatesOwn(w, x) is already correct
        // for both kinds.
        std::vector<std::pair<int, Scalar>> member_wells;       // (index into wells_, efficiency)
        std::vector<std::pair<int, Scalar>> active_children;   // (index into activeNodes_, efficiency)
    };

    explicit GroupTreeSystem(const VFPProdProperties<Scalar>& props) : props_(&props) {}

    /// Nodes must be added parent-before-child (node 0 is the terminal,
    /// Node::parent == -1; every other node's parent must already have a
    /// smaller index) -- residual() sums flow bottom-up in one pass relying
    /// on that order, the same convention NetworkSolve::Node itself documents.
    ///
    /// Building order: addActiveNode() before the wells/nested active nodes
    /// that reference it via Well::active_node (they need its index), but its
    /// own own_wells/member_wells/active_children can only be filled in *after*
    /// those exist (they need well/active-node indices in turn) -- hence
    /// activeNode(i), a mutable accessor, rather than trying to hand the
    /// whole ActiveNode over complete in one call.
    // \p alq is quoted in the units of the branch's own table's alq type
    // (same convention as NetworkProductionSystem::addNode) -- a genuine
    // node-level unknown-free constant, since it comes straight off the
    // branch (BRANPROP/GRUPNET), not from anything this system solves for.
    int addNode(Node n, const Scalar alq = Scalar{0})
    {
        nodes_.push_back(std::move(n));
        node_alq_.push_back(alq);
        fixed_pressure_.push_back(std::nullopt);
        node_explicit_fractions_.push_back(std::nullopt);
        node_flattened_.push_back(false);
        return static_cast<int>(nodes_.size()) - 1;
    }

    /// Pin node \p node's pressure at \p pressure instead of computing it from
    /// the branch above: a fixed-pressure node below the root (another root of
    /// the network that also feeds a node further up). Its subtree is solved
    /// from that pressure, and its flow still counts in every node above it --
    /// the same treatment the fixed-point network computation gives it.
    void setFixedPressure(const int node, const Scalar pressure)
    {
        fixed_pressure_[node] = pressure;
    }

    /// The water and gas fractions (WFR, GFR) node \p node's table is looked
    /// up with when its flow is below the table's first FLO value -- e.g. the
    /// node's fractions at the previous converged solve. Above it the
    /// fractions of the node's own flow apply. A node left unset gets the
    /// fractions of its wells' reference rates (Well::reference_q) in
    /// finalize(), if they have any.
    void setNodeExplicitFractions(const int node, const Scalar wfr, const Scalar gfr)
    {
        node_explicit_fractions_[node] = std::array<Scalar, 2>{wfr, gfr};
    }

    /// Each node's (WFR, GFR) for the given well rates (e.g. a converged
    /// result's well_phase_rates), for nodes with a table whose flow is at
    /// least the table's first FLO value; nullopt for the others and the
    /// terminal. What setNodeExplicitFractions() takes at the next solve.
    std::vector<std::optional<std::array<Scalar, 2>>>
    nodeFractions(const std::vector<std::array<Scalar, NP>>& well_q) const
    {
        const auto flows = nodeFlowsFromWellRates(well_q);
        std::vector<std::optional<std::array<Scalar, 2>>> out(nodes_.size());
        for (int i = 1; i <= numNodes(); ++i) {
            if (!hasTable(nodes_[i])) {
                continue;
            }
            const auto& table = props_->getTable(nodes_[i].vfp_table);
            const auto& q = flows[i];
            if (detail::getFlo(table, q[kWater], q[kOil], q[kGas]) >= table.getFloAxis().front()) {
                // getWFR()/getGFR() take opm's production-negative rates.
                out[i] = std::array<Scalar, 2>{detail::getWFR(table, -q[kWater], -q[kOil], -q[kGas]),
                                               detail::getGFR(table, -q[kWater], -q[kOil], -q[kGas])};
            }
        }
        return out;
    }

    /// Every node's flow (efficiency-scaled, [oil, water, gas], index 0
    /// unused) for the given well rates, e.g. a converged result's
    /// well_phase_rates.
    std::vector<std::array<Scalar, NP>>
    nodeFlowsFromWellRates(const std::vector<std::array<Scalar, NP>>& well_q) const
    {
        std::vector<std::array<Scalar, NP>> q(nodes_.size());
        for (int i = numNodes(); i >= 1; --i) {
            std::array<Scalar, NP> qi{};
            for (const int w : wells_at_[i]) {
                for (int p = 0; p < NP; ++p) { qi[p] += wells_[w].efficiency * well_q[w][p]; }
            }
            for (const int c : children_[i]) {
                for (int p = 0; p < NP; ++p) { qi[p] += nodes_[c].efficiency * q[c][p]; }
            }
            q[i] = qi;
        }
        return q;
    }
    int addWell(Well w) { wells_.push_back(std::move(w)); return static_cast<int>(wells_.size()) - 1; }
    int addActiveNode(ActiveNode a) { activeNodes_.push_back(std::move(a)); return static_cast<int>(activeNodes_.size()) - 1; }
    void setTerminalPressure(const Scalar p) { terminal_pressure_ = p; }

    ActiveNode& activeNode(const int i) { return activeNodes_[i]; }

    /// Must be called once, after every node/well/active-node and its
    /// membership lists (own_wells/member_wells/active_children) are fully
    /// populated, before the system is used at all.
    ///
    /// Decides which active nodes actually get a lambda unknown and row: one
    /// with empty own_wells has nothing of its own to adjust -- every one of
    /// its immediate members turned out to be individually constrained
    /// itself, either by its own limit (a rare coincidence: several siblings'
    /// own limits summing exactly to the parent's) or, far more commonly, by
    /// the network's THP control (an ordinary small group whose one or two
    /// wells currently sit on THP) -- so its target is unenforceable here:
    /// there is no lever left, network-side, that could move its total
    /// toward it. Its row is dropped entirely rather than left as a phantom,
    /// unconstrained unknown (a zero Jacobian column, which is a singular
    /// system, not a harmless one). Ancestors can still read its actual
    /// subtree total via activeNodeTotal() regardless of whether it has a
    /// row of its own.
    /// The physical/network-side facts about a well that ProdGroupTreeBalancer's
    /// FlatNetworkInput has no way to know, since it only ever sees the group
    /// tree, not the network topology: which network node the well feeds into
    /// (must already be in this system via addNode()), its tubing curve, and
    /// its own network-branch efficiency (WEFAC as the *network* topology sees
    /// it -- a separate scaling from a FlatWellShare/FlatChildRef's efficiency,
    /// which is cumulative up the *group* tree; the two trees need not
    /// coincide, so in general these differ). One entry per well named
    /// anywhere in the FlatNetworkInput passed to populateFromFlatNetwork().
    struct WellNetworkData
    {
        int node = 0;
        Scalar efficiency = 1;
        int vfp_table = NoTable;
        Scalar alq = 0;
        Scalar vfp_dp = 0;   // see Well::vfp_dp
        std::array<Scalar, NP> ipr_a{};
        std::array<Scalar, NP> ipr_b{};

        // False for an open-but-not-flowing well -- most commonly one the well
        // model has temporarily stopped (Well::Status::STOP, distinct from
        // shut: still part of the system, still solved, just pinned at zero
        // rate for now) -- with no usable linearised IPR to tie a group target
        // or a THP row to. populateFromFlatNetwork() pins such a well at zero
        // instead of building a Group/Thp row it has no real data for. A
        // genuinely shut (or otherwise absent) well never reaches this struct
        // at all -- see gatherWellNetworkDataForGroupTree()'s own doc comment.
        bool has_ipr = true;

        // Not used by populateFromFlatNetwork() itself -- this rides along
        // because it comes from the same MPI-safe gather as everything else
        // above (well->getDynamicThpLimit().has_value()), for the caller's
        // own use building extractFlatNetworkInput()'s networkThpWells set
        // before this struct even comes into play.
        bool network_thp = false;

        // Stopped (persistently or dynamically), and -- when has_trial_ipr --
        // the flowing trial IPR the well model computed for it (stopped_ipr_a/b,
        // in this struct's own sign convention), for addReopenCandidate(). Not
        // used by populateFromFlatNetwork(): a stopped well is never in the
        // balanced tree.
        bool stopped = false;
        bool has_trial_ipr = false;
        std::array<Scalar, NP> trial_ipr_a{};
        std::array<Scalar, NP> trial_ipr_b{};
        Scalar trial_bhp = 0;   // the bhp the trial IPR was taken at; 0 if unknown

        // The well's phase rates per unit of the table's FLO ([oil, water,
        // gas], production positive), fixed for the solve: see
        // applyPhaseShare(). Unset: per-phase IPRs as given.
        std::optional<std::array<Scalar, NP>> phase_share;
        // See Well::explicit_fractions, Well::reference_q.
        std::array<Scalar, 2> explicit_fractions{};
        std::array<Scalar, NP> reference_q{};
        // The well's current bhp, if it flows (0 otherwise): where a Thp
        // well starts (Well::start_bhp).
        Scalar bhp = 0;
    };

    /// Add a stopped well as a reopen candidate (Well::reopen_candidate): a Thp
    /// well on its trial IPR at its own network node, counted in the sum of
    /// \p active_node (with \p efficiency) when it has an Active ancestor. Call
    /// after populateFromFlatNetwork() and before finalize().
    void addReopenCandidate(const std::string& name,
                            const WellNetworkData& wd,
                            const std::optional<int> active_node,
                            const Scalar efficiency)
    {
        Well well;
        well.name = name;
        well.node = wd.node;
        well.efficiency = wd.efficiency;
        well.vfp_table = wd.vfp_table;
        well.alq = wd.alq;
        well.vfp_dp = wd.vfp_dp;
        well.ipr_a = wd.trial_ipr_a;
        well.ipr_b = wd.trial_ipr_b;
        well.explicit_fractions = wd.explicit_fractions;
        well.reference_q = wd.reference_q;
        if (wd.phase_share.has_value()) {
            applyPhaseShare(well, *wd.phase_share);
        }
        well.kind = WellKind::Thp;
        well.ipr_slope_limit = iprSlopeLimit(well);
        well.reopen_candidate = true;
        well.start_bhp = wd.trial_bhp;
        const int idx = addWell(std::move(well));
        if (active_node.has_value()) {
            activeNodes_[*active_node].member_wells.emplace_back(idx, efficiency);
        }
    }

    /// Index of the Active node named \p name, if there is one.
    std::optional<int> activeNodeIndex(const std::string& name) const
    {
        for (std::size_t k = 0; k < activeNodes_.size(); ++k) {
            if (activeNodes_[k].name == name) {
                return static_cast<int>(k);
            }
        }
        return std::nullopt;
    }

    /// What the converged solve says about each reopen candidate: it reopens
    /// if it flows (off its shut-in cap) on the real tubing curve -- one that
    /// only crosses on the flattened curve would be stopped again straight
    /// away, so it stays stopped. \p q is its own (pre-efficiency) rate, oil /
    /// water / gas, production positive.
    struct ReopenOutcome
    {
        std::string name;
        bool reopens = false;
        bool at_shutin = false;       // settled on its cap: no flow at this node pressure
        bool on_real_curve = true;    // false: crosses only on the flattened curve
        std::array<Scalar, NP> q{};
        Scalar bhp = 0;
        Scalar thp = 0;
    };
    std::vector<ReopenOutcome> reopenOutcomes(const Result<Scalar>& result) const
    {
        std::vector<ReopenOutcome> out;
        for (std::size_t w = 0; w < wells_.size(); ++w) {
            const auto& well = wells_[w];
            if (!well.reopen_candidate) {
                continue;
            }
            ReopenOutcome o;
            o.name = well.name;
            o.q = result.well_phase_rates[w];
            o.bhp = result.well_bhp[w];
            o.thp = result.node_pressure[well.node];
            const bool flows = std::any_of(o.q.begin(), o.q.end(),
                                           [](const Scalar r) { return r > Scalar{0}; });
            o.at_shutin = !flows || o.bhp >= well.bhp_shutin - Scalar{1.0e-3} * unit::barsa;
            o.on_real_curve = !well.ipr_slope_limit.has_value()
                || slopeLimitedBhp(well, o.thp, o.q).limit == detail::SlopeLimit::Unflattened;
            o.reopens = !o.at_shutin && o.on_real_curve;
            out.push_back(std::move(o));
        }
        return out;
    }

    /// Populate wells_/activeNodes_ from an already-balanced, already-flattened
    /// group tree (ProdGroupTreeBalancer::extractFlatNetworkInput()'s output),
    /// resolving every name into the indices this class actually needs. Must
    /// be called after the network topology (addNode()/setTerminalPressure())
    /// is already in place -- \p wellNetworkData's node indices refer to it --
    /// and before finalize().
    ///
    /// Each FlatActiveNode becomes exactly one of:
    ///  - type == Well, !networkThp: a Pinned Well, fixed_q inverted from its
    ///    own target via its own mode/resvCoeff (bhpFromTarget) -- except a
    ///    target <= 0 (stopped, or not yet meaningful before Part 1b), which
    ///    is fixed_q = 0 directly: bhpFromTarget would only zero the single
    ///    projected mode, not actually stop the well on every phase.
    ///  - type == Well, networkThp: a Thp Well -- no fixed_q, its bhp is a
    ///    genuine unknown (residual() adds its row automatically, keyed by
    ///    kind, not by anything this function does).
    ///  - type == Group: an ActiveNode, own_wells built fresh from ownWells
    ///    (these never have their own FlatActiveNode entry -- see
    ///    FlatWellShare's own doc comment) via wellNetworkData, and
    ///    active_children/member_wells resolved from activeChildren by
    ///    looking up the referenced entry's own type.
    /// Memoized by name so a FlatNetworkInput's actual ordering (children can
    /// precede or follow the parent that references them) never matters.
    void populateFromFlatNetwork(const ProdGroupTreeBalancer::FlatNetworkInput<Scalar>& flat,
                                 const std::unordered_map<std::string, WellNetworkData>& wellNetworkData)
    {
        using FlatNode = ProdGroupTreeBalancer::FlatActiveNode<Scalar>;
        std::unordered_map<std::string, const FlatNode*> byName;
        for (const auto& e : flat) {
            byName[e.name] = &e;
        }
        std::unordered_map<std::string, int> wellIdx, activeIdx;

        // The fields every well takes from wellNetworkData regardless of kind.
        auto baseWell = [&](const std::string& name) {
            const auto& wd = wellNetworkData.at(name);
            Well well;
            well.name = name;
            well.node = wd.node;
            well.efficiency = wd.efficiency;
            well.vfp_table = wd.vfp_table;
            well.alq = wd.alq;
            well.vfp_dp = wd.vfp_dp;
            well.ipr_a = wd.ipr_a;
            well.ipr_b = wd.ipr_b;
            well.explicit_fractions = wd.explicit_fractions;
            well.reference_q = wd.reference_q;
            well.start_bhp = wd.bhp;
            if (wd.phase_share.has_value()) {
                applyPhaseShare(well, *wd.phase_share);
            }
            return well;
        };

        // A GSATPROD satellite group (FlatActiveNode::satelliteRates set):
        // no IPR/VFP/network-node data at all -- it is not a real well, so
        // it is never in wellNetworkData -- just a fixed three-phase rate,
        // pinned exactly like an ordinary Pinned well's fixed_q, but with
        // node left at 0 (the terminal): a satellite has no network-topology
        // role, only a group-tree bookkeeping one, and node 0 is never
        // matched by residual()'s node-flow loop (i <- 1..numNodes()), so
        // this well is correctly invisible to every physical node-pressure
        // equation while still counting toward whichever ancestor's
        // activeNodeTotal() references it via member_wells.
        auto ensureSatellite = [&](const std::string& name) -> int {
            if (const auto it = wellIdx.find(name); it != wellIdx.end()) {
                return it->second;
            }
            const auto& src = *byName.at(name);
            Well well;
            well.name = name;
            well.kind = WellKind::Pinned;
            well.fixed_q = *src.satelliteRates;
            const int idx = addWell(std::move(well));
            wellIdx[name] = idx;
            return idx;
        };

        std::function<int(const std::string&)> ensureWell;
        std::function<int(const std::string&)> ensureActive;

        ensureWell = [&](const std::string& name) -> int {
            if (const auto it = wellIdx.find(name); it != wellIdx.end()) {
                return it->second;
            }
            const auto& src = *byName.at(name);
            Well well = baseWell(name);
            // A well with no usable IPR right now (WellNetworkData::has_ipr) --
            // typically stopped, not shut: still part of the system, but with
            // nothing to tie a THP row or a positive target to -- is pinned at
            // zero regardless of what the balancer's own mode/target say,
            // exactly like a genuine full stop below.
            const bool has_ipr = wellNetworkData.at(name).has_ipr;
            if (src.networkThp && has_ipr) {
                well.kind = WellKind::Thp;
                well.ipr_slope_limit = iprSlopeLimit(well);
            } else {
                well.kind = WellKind::Pinned;
                if (has_ipr && src.target > Scalar{0}) {
                    const auto weights = phaseWeights(src.mode, src.resvCoeff);
                    const Scalar bhp = bhpFromTarget(well, weights, src.target);
                    for (int p = 0; p < NP; ++p) {
                        well.fixed_q[p] = well.ipr_a[p] + well.ipr_b[p] * bhp;
                    }
                }
                // !has_ipr, or target <= 0: stopped (or not yet meaningful) --
                // well.fixed_q stays value-initialized zero, a genuine full stop.
            }
            const int idx = addWell(std::move(well));
            wellIdx[name] = idx;
            return idx;
        };

        ensureActive = [&](const std::string& name) -> int {
            if (const auto it = activeIdx.find(name); it != activeIdx.end()) {
                return it->second;
            }
            const auto& src = *byName.at(name);
            ActiveNode a;
            a.name = src.name;
            a.mode = src.mode;
            a.target = src.target;
            a.resv_coeff = src.resvCoeff;
            const int idx = addActiveNode(std::move(a));
            activeIdx[name] = idx;   // register before recursing

            for (const auto& w : src.ownWells) {
                Well well = baseWell(w.name);
                if (wellNetworkData.at(w.name).has_ipr) {
                    well.kind = WellKind::Group;
                    well.active_node = idx;
                    well.guide_rate = w.guideRate;
                    const int wIdx = addWell(std::move(well));
                    activeNodes_[idx].own_wells.push_back(wIdx);
                } else {
                    // No usable IPR (typically stopped, not shut -- see
                    // WellNetworkData::has_ipr): nothing to tie this node's
                    // lambda to for this well, so it is pinned at zero and
                    // counted as a fixed (zero) contribution instead, the same
                    // as any other member_wells entry -- not left in own_wells,
                    // where finalize() would otherwise wrongly see a real
                    // lambda-tied well and keep a row with nothing to adjust
                    // if this turns out to be the node's only "own" well.
                    well.kind = WellKind::Pinned;
                    const Scalar efficiency = well.efficiency;
                    const int wIdx = addWell(std::move(well));
                    activeNodes_[idx].member_wells.emplace_back(wIdx, efficiency);
                }
            }
            for (const auto& child : src.activeChildren) {
                const auto& childSrc = *byName.at(child.name);
                if (childSrc.satelliteRates.has_value()) {
                    const int wIdx = ensureSatellite(child.name);
                    activeNodes_[idx].member_wells.emplace_back(wIdx, child.efficiency);
                } else if (childSrc.type == ProdNodeType::Well) {
                    const int wIdx = ensureWell(child.name);
                    activeNodes_[idx].member_wells.emplace_back(wIdx, child.efficiency);
                } else {
                    const int aIdx = ensureActive(child.name);
                    activeNodes_[idx].active_children.emplace_back(aIdx, child.efficiency);
                }
            }
            return idx;
        };

        for (const auto& e : flat) {
            if (e.satelliteRates.has_value()) {
                ensureSatellite(e.name);
            } else if (e.type == ProdNodeType::Well) {
                ensureWell(e.name);
            } else {
                ensureActive(e.name);
            }
        }
    }

    void finalize()
    {
        lambda_slot_.assign(activeNodes_.size(), -1);
        num_lambdas_ = 0;
        for (std::size_t k = 0; k < activeNodes_.size(); ++k) {
            if (!activeNodes_[k].own_wells.empty()) {
                lambda_slot_[k] = num_lambdas_++;
            }
        }
        thp_capped_.assign(wells_.size(), false);
        for (auto& well : wells_) {
            if (well.kind == WellKind::Thp) {
                well.bhp_shutin = shutInBhp(well);
            }
        }

        // Topology and column lookups, so residual() and jacobian() need not
        // search for them on every call.
        wells_at_.assign(nodes_.size(), {});
        children_.assign(nodes_.size(), {});
        thp_col_.assign(wells_.size(), -1);
        num_thp_ = 0;
        for (int w = 0; w < numWells(); ++w) {
            if (wells_[w].node > 0) {
                wells_at_[wells_[w].node].push_back(w);
            }
            if (wells_[w].kind == WellKind::Thp) {
                thp_col_[w] = numNodes() + numActiveLambdas() + num_thp_++;
            }
        }
        for (std::size_t c = 1; c < nodes_.size(); ++c) {
            if (nodes_[c].parent > 0) {
                children_[nodes_[c].parent].push_back(static_cast<int>(c));
            }
        }

        // Nodes without fallback fractions: those of their wells' reference
        // rates, mixed as the node would carry them.
        std::vector<std::array<Scalar, NP>> reference(wells_.size());
        for (int w = 0; w < numWells(); ++w) {
            reference[w] = wells_[w].reference_q;
        }
        const auto mixed = nodeFlowsFromWellRates(reference);
        for (int i = 1; i <= numNodes(); ++i) {
            if (node_explicit_fractions_[i].has_value() || !hasTable(nodes_[i])) {
                continue;
            }
            const auto& table = props_->getTable(nodes_[i].vfp_table);
            const auto& q = mixed[i];
            if (detail::getFlo(table, q[kWater], q[kOil], q[kGas]) > Scalar{0}) {
                // getWFR()/getGFR() take opm's production-negative rates.
                node_explicit_fractions_[i] = std::array<Scalar, 2>{
                    detail::getWFR(table, -q[kWater], -q[kOil], -q[kGas]),
                    detail::getGFR(table, -q[kWater], -q[kOil], -q[kGas])};
            }
        }
        buildRateDerivatives();
    }

    /// Assemble jacobian() from the table derivatives instead of letting
    /// solve() difference residual(). See jacobian().
    void setAnalyticJacobian(const bool on) { analytic_jacobian_ = on; }

    /// Look up node \p node's branch table flattened so that its pressure
    /// never falls with flow (slope limit 0), or the table as it is.
    void setNodeFlattening(const int node, const bool on) { node_flattened_[node] = on; }

    /// setNodeFlattening() for every node with a table and no fixed pressure.
    /// A branch table on its falling low-flow branch makes the network
    /// non-monotone (more than one solution); flattened, it is not.
    void setBranchFlattening(const bool on)
    {
        for (int i = 1; i <= numNodes(); ++i) {
            node_flattened_[i] = on && hasTable(nodes_[i]) && !fixed_pressure_[i].has_value();
        }
    }

    /// The flattened nodes whose lookup at \p result's converged flows and
    /// pressures is not on the real table: there the solution is one of the
    /// flattened network, not of the real one.
    std::vector<int> flattenedNodes(const Result<Scalar>& result) const
    {
        std::vector<int> out;
        const auto flows = nodeFlowsFromWellRates(result.well_phase_rates);
        for (int i = 1; i <= numNodes(); ++i) {
            if (!node_flattened_[i]) {
                continue;
            }
            const auto& q = flows[i];
            const auto f = nodeFractions(i);
            const auto limited = props_->bhp_with_slope_limit(nodes_[i].vfp_table, -q[kWater], -q[kOil],
                                                              -q[kGas], result.node_pressure[nodes_[i].parent],
                                                              node_alq_[i], f.wfr, f.gfr, f.fixed, Scalar{0});
            if (limited.limit != detail::SlopeLimit::Unflattened) {
                out.push_back(i);
            }
        }
        return out;
    }

    /// Start every Thp well at its bhp in \p result (a converged solve of
    /// this system).
    void setStartFromResult(const Result<Scalar>& result)
    {
        for (std::size_t w = 0; w < wells_.size(); ++w) {
            if (wells_[w].kind == WellKind::Thp) {
                wells_[w].start_bhp = result.well_bhp[w];
            }
        }
    }

    const Node& node(const int i) const { return nodes_[i]; }

    int numNodes() const { return static_cast<int>(nodes_.size()) - 1; }   // excludes the terminal
    int numActiveNodes() const { return static_cast<int>(activeNodes_.size()); }
    /// Active nodes that actually got a lambda unknown -- see finalize().
    int numActiveLambdas() const { return num_lambdas_; }
    int numWells() const override { return static_cast<int>(wells_.size()); }
    const std::vector<Well>& wells() const { return wells_; }
    const std::vector<ActiveNode>& activeNodes() const { return activeNodes_; }
    Scalar terminalPressure() const { return terminal_pressure_; }

    /// The steepest d(bhp)/d(FLO) \p well's tubing curve may have before
    /// bhp_with_slope_limit() flattens it -- the well's own IPR slope in that
    /// same sense, with a margin. Along the IPR, FLO = flo(ipr_a + ipr_b*bhp),
    /// and flo() is a selection/sum of the phase rates, so d(FLO)/d(bhp) is
    /// just flo(ipr_b) -- negative for a producer, hence a negative limit.
    ///
    /// The margin is relative rather than an absolute pressure-per-rate
    /// epsilon: a curve exactly parallel to the IPR admits no unique crossing
    /// either, and "within 5% of the IPR's own slope" scales with the well
    /// instead of needing to be tuned per deck -- and needs no unit
    /// conversion to get wrong. nullopt for a well whose own phase mix
    /// registers as no FLO at all on this table: nothing to limit against, so
    /// it keeps the unflattened curve.
    ///
    /// populateFromFlatNetwork() sets Well::ipr_slope_limit from this for
    /// every Thp well it builds; it is public so that a hand-built well (the
    /// tests) can be given the same limit the real path would compute.
    /// Give \p well one IPR in its table's FLO with fixed phase proportions:
    /// FLO = a_FLO + b_FLO * bhp (the per-phase coefficients combined as the
    /// table's FLO combines rates) and q_p = share[p] * FLO, with share the
    /// phase rates per unit of FLO. Its rates, its node's flow and its tubing
    /// lookups (Well::fractions) then all have one composition at every bhp,
    /// so the FLO-based slope limit keeps its row monotone, and all phases
    /// reach zero at one shut-in bhp. The price is that the composition does
    /// not change with the rate, as the well's own inflow can -- the outer
    /// rounds correct it from a new start point.
    void applyPhaseShare(Well& well, const std::array<Scalar, NP>& share) const
    {
        const auto& table = props_->getTable(well.vfp_table);
        const Scalar a_flo = detail::getFlo(table, well.ipr_a[kWater], well.ipr_a[kOil], well.ipr_a[kGas]);
        const Scalar b_flo = detail::getFlo(table, well.ipr_b[kWater], well.ipr_b[kOil], well.ipr_b[kGas]);
        for (int p = 0; p < NP; ++p) {
            well.ipr_a[p] = share[p] * a_flo;
            well.ipr_b[p] = share[p] * b_flo;
        }
        // getWFR()/getGFR() take opm's production-negative rates.
        well.fractions = std::array<Scalar, 2>{
            detail::getWFR(table, -share[kWater], -share[kOil], -share[kGas]),
            detail::getGFR(table, -share[kWater], -share[kOil], -share[kGas])};
    }

    std::optional<Scalar> iprSlopeLimit(const Well& well) const
    {
        const auto& table = props_->getTable(well.vfp_table);
        const Scalar dflo_dbhp = detail::getFlo(table, well.ipr_b[kWater],
                                                well.ipr_b[kOil], well.ipr_b[kGas]);
        if (dflo_dbhp == Scalar{0}) {
            return std::nullopt;
        }
        constexpr Scalar margin = Scalar{0.95};
        return margin / dflo_dbhp;
    }

    /// After a converged solve, the name of the Thp well whose converged
    /// operating point most needs correcting -- see Well::ipr_slope_limit for
    /// why: a well solved against the slope-limited curve converges cleanly
    /// even where its real tubing curve falls faster than its own IPR, a
    /// region it cannot actually operate in. Re-evaluating the lookup at the
    /// converged point says whether any flattening was needed there at all
    /// (detail::SlopeLimit::Unflattened means the real table was used as is,
    /// so the answer stands); anything else means this well's row needs to be
    /// pinned at zero and the whole tree rebalanced -- the caller's job, not
    /// this class's, which knows nothing of the balancer or WellState.
    ///
    /// Among the wells that did need flattening, the one whose flattened bhp
    /// departs furthest from what the real table says at the same point is
    /// returned (nullopt if there are none), matching the "stop the worst
    /// offender, one at a time" anti-oscillation choice. A Thp well with no
    /// slope limit set is never flagged -- there is nothing to diagnose
    /// against.
    std::optional<std::string> worstCliffViolation(const Result<Scalar>& result) const
    {
        std::optional<std::string> worst;
        Scalar worst_gap = Scalar{0};
        for (std::size_t w = 0; w < wells_.size(); ++w) {
            const auto& well = wells_[w];
            if (well.kind != WellKind::Thp || !well.ipr_slope_limit.has_value()
                || well.reopen_candidate) {
                continue;
            }
            const Scalar thp = result.node_pressure[well.node];
            const auto& q = result.well_phase_rates[w];
            const auto limited = slopeLimitedBhp(well, thp, q);
            if (limited.limit == detail::SlopeLimit::Unflattened) {
                continue;
            }
            const Scalar gap = std::abs(limited.evaluation.value
                                        - tableBhp(well.vfp_table, thp, q, well.alq,
                                                   wellFractions(well, q)));
            if (!worst.has_value() || gap > worst_gap) {
                worst_gap = gap;
                worst = well.name;
            }
        }
        return worst;
    }

    int size() const override { return numNodes() + numActiveLambdas() + numThpWells(); }

    Scalar columnScale(const int i) const override
    {
        // Node pressures and thp-well bhps (everything but the middle
        // numActiveLambdas() entries) are bar-scale. A lambda is dimensionless
        // -- guideRate is rate-valued (GuideRate::get() falls back to raw
        // potentials, in the same units as a rate), so target/sum(guideRate)
        // is an O(1) ratio -- but this is the piece most likely to need
        // revisiting once this runs against a real deck's guide rates.
        const bool is_lambda = i >= numNodes() && i < numNodes() + numActiveLambdas();
        return is_lambda ? Scalar{1} : unit::barsa;
    }

    State start(const State& node_pressure_guess) const override
    {
        State x(size(), Scalar{0});
        const int n = numNodes();
        for (int i = 0; i < n; ++i) {
            x[i] = fixed_pressure_[i + 1].value_or(node_pressure_guess[i]);
        }
        for (int k = 0; k < numActiveNodes(); ++k) {
            if (lambda_slot_[k] < 0) {
                continue;   // dropped: see finalize()
            }
            const auto& a = activeNodes_[k];
            Scalar guide_sum = 0;
            for (const int w : a.own_wells) {
                guide_sum += wells_[w].guide_rate;
            }
            x[lambdaIdx(k)] = (guide_sum > Scalar{0}) ? a.target / guide_sum : Scalar{0};
        }
        for (int w = 0; w < numWells(); ++w) {
            if (wells_[w].kind == WellKind::Thp) {
                // A flowing well starts at its current bhp, a reopen candidate
                // where its trial IPR was taken (on the flowing side of a
                // tubing curve with a minimum). Without either: a candidate at
                // shut-in (rate 0); a flowing well at the node pressure guess,
                // as before (a low bhp, a high rate).
                const auto& well = wells_[w];
                x[thpBhpIdx(w)] = well.start_bhp > Scalar{0} ? std::min(well.start_bhp, well.bhp_shutin)
                    : (well.reopen_candidate ? well.bhp_shutin : node_pressure_guess[well.node - 1]);
            }
        }
        return x;
    }

    /// Every phase of well w's own (pre-efficiency) rate, at state x.
    std::array<Scalar, NP> wellPhaseRatesOwn(const int w, const State& x) const
    {
        const auto& well = wells_[w];
        Scalar bhp{};
        switch (well.kind) {
        case WellKind::Pinned:
            return well.fixed_q;
        case WellKind::Group: {
            const Scalar q_target_mode = well.guide_rate * x[lambdaIdx(well.active_node)];
            const auto weights = phaseWeights(activeNodes_[well.active_node].mode,
                                              activeNodes_[well.active_node].resv_coeff);
            bhp = bhpFromTarget(well, weights, q_target_mode);
            break;
        }
        case WellKind::Thp:
            bhp = x[thpBhpIdx(w)];
            break;
        }
        std::array<Scalar, NP> q{};
        for (int p = 0; p < NP; ++p) {
            q[p] = well.ipr_a[p] + well.ipr_b[p] * bhp;
        }
        return q;
    }

    /// Active node k's whole subtree (its own wells, plus every nested active
    /// child, recursively), projected onto `weights` -- which is *not*
    /// necessarily node k's own mode. It has to be given by the caller: a
    /// parent whose mode differs from a nested child's cannot use the child's
    /// target as a stand-in for "the child's contribution on the parent's
    /// mode" (a water target is not an oil rate, whatever its value), so it
    /// asks this same subtree for its own total projected onto the parent's
    /// weights instead. Only ever the same as node k's own target when the
    /// caller happens to ask with node k's own weights (i.e. this is node k's
    /// own row calling with its own mode).
    Scalar activeNodeTotal(const int k, const State& x, const std::array<Scalar, NP>& weights) const
    {
        const auto& a = activeNodes_[k];
        Scalar sum = 0;
        for (const int w : a.own_wells) {
            sum += wells_[w].efficiency * projectOnMode(wellPhaseRatesOwn(w, x), weights);
        }
        for (const auto& [w, eff] : a.member_wells) {
            sum += eff * projectOnMode(wellPhaseRatesOwn(w, x), weights);
        }
        for (const auto& [child, eff] : a.active_children) {
            sum += eff * activeNodeTotal(child, x, weights);   // same weights all the way down
        }
        return sum;
    }

    State residual(const State& x) const override
    {
        State r(size(), Scalar{0});
        const int n = numNodes();

        // Node rows, bottom-up: q at node i is its own wells' (efficiency-
        // scaled) rates plus its children's (efficiency-scaled) q, and the
        // node's pressure is what the tubing table says that q needs at the
        // upstream pressure -- or just the upstream pressure, node-to-node,
        // where there is no table.
        const auto q = nodeFlows(x);
        for (int i = 1; i <= n; ++i) {
            if (fixed_pressure_[i].has_value()) {
                r[pIdx(i)] = (x[pIdx(i)] - *fixed_pressure_[i]) / unit::barsa;
                continue;
            }
            const Scalar upstream = (nodes_[i].parent == 0) ? terminal_pressure_ : x[pIdx(nodes_[i].parent)];
            const Scalar computed = hasTable(nodes_[i]) ? nodeTableBhp(i, upstream, q[i]) : upstream;
            r[pIdx(i)] = (x[pIdx(i)] - computed) / unit::barsa;
        }

        // Active-node rows: own_wells (tied to this node's own lambda) +
        // member_wells + activeChildren's subtrees, all projected onto this
        // node's own mode -- see activeNodeTotal() for why a nested child's
        // contribution has to be computed that way rather than read off its
        // own target when the two modes differ.
        for (int k = 0; k < numActiveNodes(); ++k) {
            if (lambda_slot_[k] < 0) {
                continue;   // dropped: empty own_wells, nothing to adjust -- see finalize()
            }
            const auto& a = activeNodes_[k];
            const auto weights = phaseWeights(a.mode, a.resv_coeff);
            const Scalar sum = activeNodeTotal(k, x, weights);
            const Scalar scale = (a.target > Scalar{0}) ? a.target : Scalar{1};
            r[lambdaIdx(k)] = (sum - a.target) / scale;
        }

        // Thp-well rows -- see thpWellResidualRow() for what each one is.
        for (int w = 0; w < numWells(); ++w) {
            if (wells_[w].kind != WellKind::Thp) { continue; }
            r[thpBhpIdx(w)] = thpWellResidualRow(w, x);
        }
        return r;
    }

    /// This well's own bhp against what the tubing table says it needs at
    /// its node's pressure -- its thp -- for the rate its own IPR gives it
    /// at that same bhp. The network's only role here is exactly this:
    /// supplying the thp a well's own row is solved against. Except when
    /// updateControls() found the *current* iterate's bhp at or past
    /// bhp_shutin: every phase rate would be negative past that point, not a
    /// real operating point, so the row switches to a plain pin at
    /// bhp_shutin, decoupled from thp entirely -- see WellKind::Thp members'
    /// own doc comment on bhp_shutin for why the tubing-curve equation
    /// itself, left as is, would not converge there (its own root and
    /// bhp_shutin are generally different numbers).
    Scalar thpWellResidualRow(const int w, const State& x) const
    {
        const auto& well = wells_[w];
        const Scalar bhp = x[thpBhpIdx(w)];
        if (thp_capped_[w]) {
            return (bhp - well.bhp_shutin) / unit::barsa;
        }
        const Scalar thp = x[pIdx(well.node)];
        const auto qw = wellPhaseRatesOwn(w, x);
        const Scalar computed = slopeLimitedBhp(well, thp, qw).evaluation.value - well.vfp_dp;
        return (bhp - computed) / unit::barsa;
    }

    /// This well's own tubing curve at (thp, q), flattened against its own
    /// IPR where Well::ipr_slope_limit says to -- see that field. Without a
    /// limit set this is exactly tableBhp(), reported as Unflattened.
    detail::SlopeLimitedEvaluation<Scalar>
    slopeLimitedBhp(const Well& well, const Scalar thp, const std::array<Scalar, NP>& q) const
    {
        // props_->bhp*() want water, oil, gas (aqua, liquid, vapour) and
        // negative-for-production; q here is oil, water, gas and positive --
        // the same two conversions tableBhp() makes, at the same one place.
        if (!well.ipr_slope_limit.has_value()) {
            detail::SlopeLimitedEvaluation<Scalar> plain;
            plain.evaluation.value = tableBhp(well.vfp_table, thp, q, well.alq, wellFractions(well, q));
            return plain;
        }
        const auto f = wellFractions(well, q);
        return props_->bhp_with_slope_limit(well.vfp_table, -q[kWater], -q[kOil], -q[kGas],
                                            thp, well.alq, f.wfr, f.gfr, f.fixed,
                                            *well.ipr_slope_limit);
    }

    /// The Jacobian of residual(), from the table derivatives. Every well rate
    /// is affine in the unknowns (Pinned: constant; Group: in its node's
    /// lambda; Thp: in its own bhp), so node flows and target rows have
    /// constant derivatives, precomputed by finalize(); only the tubing tables
    /// need derivatives at x (tableLookup()). Rows are scaled exactly as
    /// residual() scales them.
    DenseMatrix<Scalar> jacobian(const State& x) const override
    {
        DenseMatrix<Scalar> J(size());
        const Scalar bar = unit::barsa;
        const int n = numNodes();
        const auto q = nodeFlows(x);

        // Node rows: (p_i - table(p_parent, q_i)) / bar.
        for (int i = 1; i <= n; ++i) {
            const int row = pIdx(i);
            const int parent = nodes_[i].parent;
            J(row, row) += Scalar{1} / bar;
            if (fixed_pressure_[i].has_value()) {
                continue;
            }
            if (hasTable(nodes_[i])) {
                const Scalar upstream = (parent == 0) ? terminal_pressure_ : x[pIdx(parent)];
                const auto t = tableLookup(nodes_[i].vfp_table, upstream, q[i], node_alq_[i],
                                           nodeSlopeLimit(i), nodeFractions(i));
                if (parent != 0) {
                    J(row, pIdx(parent)) -= t.dthp / bar;
                }
                for (const auto& [col, dq] : node_flow_jac_[i]) {
                    J(row, col) -= (t.dq[kOil] * dq[kOil] + t.dq[kWater] * dq[kWater]
                                    + t.dq[kGas] * dq[kGas]) / bar;
                }
            } else if (parent != 0) {
                J(row, pIdx(parent)) -= Scalar{1} / bar;
            }
        }

        // Target rows: linear, see buildRateDerivatives().
        for (int k = 0; k < numActiveNodes(); ++k) {
            if (lambda_slot_[k] < 0) {
                continue;
            }
            for (const auto& [col, c] : target_jac_[k]) {
                J(lambdaIdx(k), col) += c;
            }
        }

        // Thp-well rows: (bhp - table~(p_node, q_w(bhp))) / bar, or the
        // shut-in pin (bhp - bhp_shutin) / bar while capped.
        for (int w = 0; w < numWells(); ++w) {
            const auto& well = wells_[w];
            if (well.kind != WellKind::Thp) {
                continue;
            }
            const int row = thpBhpIdx(w);
            if (thp_capped_[w]) {
                J(row, row) += Scalar{1} / bar;
                continue;
            }
            const Scalar thp = x[pIdx(well.node)];
            const auto qw = wellPhaseRatesOwn(w, x);
            const auto t = tableLookup(well.vfp_table, thp, qw, well.alq, well.ipr_slope_limit,
                                       wellFractions(well, qw));
            const Scalar dtable_dbhp = t.dq[kOil] * well.ipr_b[kOil]
                + t.dq[kWater] * well.ipr_b[kWater] + t.dq[kGas] * well.ipr_b[kGas];
            J(row, row) += (Scalar{1} - dtable_dbhp) / bar;
            J(row, pIdx(well.node)) -= t.dthp / bar;
        }
        return J;
    }
    bool usesAnalyticJacobian() const override { return analytic_jacobian_; }

    /// The one piece of per-iterate state this class has: whether each Thp
    /// well's *current* bhp has reached its own bhp_shutin (see that field's
    /// doc comment). Deliberately not the kind of combinatorial, multi-well,
    /// hysteresis-carrying switching this design otherwise avoids inside
    /// Newton -- it is a single scalar bound on one well's own unknown,
    /// decided fresh from the current iterate alone every time, the same way
    /// a standard bound-constrained ("projected") Newton method treats an
    /// active box constraint.
    ///
    /// At the cap the rate is zero, and the well stays capped only while its
    /// tubing curve still needs at least bhp_shutin there; if it needs less,
    /// the free row would pull bhp below the cap -- the well can flow -- so the
    /// cap is released (the projected-Newton rule for an active bound). Without
    /// this a well that once reached its cap, or a reopen candidate starting
    /// at it, could never flow again within the solve.
    bool updateControls(const State& x) override
    {
        bool moved = false;
        for (int w = 0; w < numWells(); ++w) {
            if (wells_[w].kind != WellKind::Thp) { continue; }
            const auto& well = wells_[w];
            bool capped = x[thpBhpIdx(w)] >= well.bhp_shutin;
            if (capped) {
                // Just off the cap, so the phase fractions are the well's own
                // rather than 0/0 at exactly zero rate.
                const Scalar bhp = well.bhp_shutin - Scalar{1.0e-3} * unit::barsa;
                std::array<Scalar, NP> q{};
                for (int p = 0; p < NP; ++p) {
                    q[p] = well.ipr_a[p] + well.ipr_b[p] * bhp;
                }
                const Scalar thp = x[pIdx(well.node)];
                capped = slopeLimitedBhp(well, thp, q).evaluation.value - well.vfp_dp
                    >= well.bhp_shutin;
            }
            moved = moved || (capped != thp_capped_[w]);
            thp_capped_[w] = capped;
        }
        return moved;
    }
    char controlLetter(const int w) const override
    {
        switch (wells_[w].kind) {
        case WellKind::Pinned: return 'P';
        case WellKind::Group:  return 'G';
        case WellKind::Thp:    return thp_capped_[w] ? 'C' : 'T';
        }
        return '?';
    }

    /// No step limiting except the one hard physical bound: a Thp well's bhp
    /// must never be pushed past its own bhp_shutin, where every phase's IPR
    /// rate turns negative -- not a real operating point, and not one the
    /// tubing table (queried with that negative rate) has any reason to
    /// behave sensibly at either.
    State limitStep(const State& x, const State& dx) const override
    {
        State limited = dx;
        for (int w = 0; w < numWells(); ++w) {
            const auto& well = wells_[w];
            if (well.kind != WellKind::Thp) { continue; }
            const int idx = thpBhpIdx(w);
            if (x[idx] + limited[idx] > well.bhp_shutin) {
                limited[idx] = well.bhp_shutin - x[idx];
            }
        }
        return limited;
    }


    // ------------------------------------------------------------------
    // TEMPORARY diagnostics for solves that fail to converge (see
    // NetworkSolve::SystemBase::describeIteration()).
    // ------------------------------------------------------------------

    /// What unknown (column) or equation (row) \p i is: the layout is node
    /// pressures, then one lambda per active node that has one, then the Thp
    /// wells' bhp.
    std::string describeIndex(const int i) const
    {
        if (i < numNodes()) {
            return fmt::format("node {}", nodes_[i + 1].name);
        }
        if (i < numNodes() + numActiveLambdas()) {
            for (int k = 0; k < numActiveNodes(); ++k) {
                if (lambda_slot_[k] >= 0 && lambdaIdx(k) == i) {
                    return fmt::format("target {} ({})", activeNodes_[k].name,
                                       ::Opm::WellProducerCMode2String(activeNodes_[k].mode));
                }
            }
            return "lambda ?";
        }
        for (int w = 0; w < numWells(); ++w) {
            if (wells_[w].kind == WellKind::Thp && thpBhpIdx(w) == i) {
                return fmt::format("well {}", wells_[w].name);
            }
        }
        return "?";
    }

    static std::string slopeLimitName(const detail::SlopeLimit limit)
    {
        switch (limit) {
        case detail::SlopeLimit::Unflattened: return "unflattened";
        case detail::SlopeLimit::Bridged: return "bridged";
        case detail::SlopeLimit::NonMonotone: return "NON-MONOTONE";
        case detail::SlopeLimit::Clamped: return "CLAMPED";
        }
        return "?";
    }

    std::string describeWell(const int w, const State& x) const
    {
        const auto& well = wells_[w];
        const auto q = wellPhaseRatesOwn(w, x);
        const auto rates = fmt::format("O/W/G {:.1f}/{:.1f}/{:.1f} m3/d", q[kOil] * unit::day,
                                       q[kWater] * unit::day, q[kGas] * unit::day);
        if (well.kind != WellKind::Thp) {
            return fmt::format("  well {} [{}]: {}", well.name, controlLetter(w), rates);
        }
        const Scalar thp = x[pIdx(well.node)];
        const Scalar bhp = x[thpBhpIdx(w)];
        const auto tubing = slopeLimitedBhp(well, thp, q);
        const auto f = wellFractions(well, q);
        return fmt::format("  well {} [{}{}]: thp {:.3f}, bhp {:.3f}, shut-in {:.3f}, tubing {:.3f} ({}), "
                           "row {:.3e}, {}; fractions WFR/GFR fixed {}, explicit {:.4g}/{:.4g}, used {:.4g}/{:.4g}{}",
                           well.name, controlLetter(w), well.reopen_candidate ? ", candidate" : "",
                           thp / unit::barsa, bhp / unit::barsa, well.bhp_shutin / unit::barsa,
                           (tubing.evaluation.value - well.vfp_dp) / unit::barsa,
                           slopeLimitName(tubing.limit), thpWellResidualRow(w, x), rates,
                           well.fractions.has_value()
                               ? fmt::format("{:.4g}/{:.4g}", (*well.fractions)[0], (*well.fractions)[1])
                               : std::string("none"),
                           well.explicit_fractions[0], well.explicit_fractions[1], f.wfr, f.gfr,
                           f.fixed ? "" : " (below first FLO only)");
    }

    std::string describeRow(const int i, const State& x, const State& r,
                            const std::vector<std::array<Scalar, NP>>& flows) const
    {
        if (i < numNodes()) {
            const int node = i + 1;
            const Scalar upstream = (nodes_[node].parent == 0) ? terminal_pressure_
                                                               : x[pIdx(nodes_[node].parent)];
            const Scalar computed = fixed_pressure_[node].has_value() ? *fixed_pressure_[node]
                : (hasTable(nodes_[node]) ? nodeTableBhp(node, upstream, flows[node]) : upstream);
            return fmt::format("  node {}: r {:.3e}, p {:.3f}, needed {:.3f} (upstream {:.3f}{}), "
                               "O/W/G {:.1f}/{:.1f}/{:.1f} m3/d",
                               nodes_[node].name, r[i], x[i] / unit::barsa, computed / unit::barsa,
                               upstream / unit::barsa,
                               fixed_pressure_[node].has_value() ? ", fixed"
                                   : (hasTable(nodes_[node]) ? "" : ", no table"),
                               flows[node][kOil] * unit::day, flows[node][kWater] * unit::day,
                               flows[node][kGas] * unit::day);
        }
        if (i < numNodes() + numActiveLambdas()) {
            for (int k = 0; k < numActiveNodes(); ++k) {
                if (lambda_slot_[k] < 0 || lambdaIdx(k) != i) {
                    continue;
                }
                const auto& a = activeNodes_[k];
                const auto weights = phaseWeights(a.mode, a.resv_coeff);
                int capped = 0;
                for (const auto& [w, eff] : a.member_wells) {
                    capped += (wells_[w].kind == WellKind::Thp && thp_capped_[w]) ? 1 : 0;
                }
                return fmt::format("  target {} ({}): r {:.3e}, lambda {:.4e}, target {:.3f}, sum {:.3f} "
                                   "(own {}, members {} of which {} capped, children {})",
                                   a.name, ::Opm::WellProducerCMode2String(a.mode), r[i], x[i],
                                   a.target * unit::day, activeNodeTotal(k, x, weights) * unit::day,
                                   a.own_wells.size(), a.member_wells.size(), capped,
                                   a.active_children.size());
            }
        }
        for (int w = 0; w < numWells(); ++w) {
            if (wells_[w].kind == WellKind::Thp && thpBhpIdx(w) == i) {
                return describeWell(w, x);
            }
        }
        return fmt::format("  row {}: r {:.3e}", i, r[i]);
    }

    std::vector<std::string> describeIteration(const int it, const State& x, const State& r,
                                               const State& dx_raw, const State& dx_limited,
                                               const bool controls_moved) const override
    {
        std::vector<std::string> out;
        const int n = size();
        auto norm2 = [](const State& v) {
            Scalar s = 0;
            for (const auto e : v) { s += e * e; }
            return std::sqrt(s);
        };
        Scalar worst = 0;
        for (const auto e : r) { worst = std::max(worst, std::abs(e)); }

        std::string set;
        for (int w = 0; w < numWells(); ++w) { set += controlLetter(w); }

        // The step, in units of each column's scale, and what it does to the
        // residual if taken in full.
        std::string step = "no step (singular Jacobian)";
        if (static_cast<int>(dx_limited.size()) == n) {
            int i_max = 0;
            Scalar s_max = -1;
            for (int i = 0; i < n; ++i) {
                const Scalar s = std::abs(dx_raw[i]) / columnScale(i);
                if (s > s_max) { s_max = s; i_max = i; }
            }
            State next = x;
            for (int i = 0; i < n; ++i) { next[i] += dx_limited[i]; }
            const auto r_next = residual(next);
            Scalar worst_next = 0;
            for (const auto e : r_next) { worst_next = std::max(worst_next, std::abs(e)); }
            step = fmt::format("max |dx|/scale {:.3e} at {} (dx {:.4e}, limited to {:.4e}); "
                               "full step gives max|r| {:.3e}, |r|2 {:.3e}",
                               s_max, describeIndex(i_max), dx_raw[i_max], dx_limited[i_max],
                               worst_next, norm2(r_next));
        }
        out.push_back(fmt::format("Network diag it {}: max|r| {:.3e}, |r|2 {:.3e}, controls {}, set {}; {}",
                                  it, worst, norm2(r), controls_moved ? "moved" : "fixed", set, step));

        // The three worst rows.
        std::vector<int> order(n);
        for (int i = 0; i < n; ++i) { order[i] = i; }
        const int top = std::min(n, 3);
        std::partial_sort(order.begin(), order.begin() + top, order.end(),
                          [&r](const int a, const int b) { return std::abs(r[a]) > std::abs(r[b]); });
        const auto flows = nodeFlows(x);
        for (int j = 0; j < top; ++j) {
            out.push_back(describeRow(order[j], x, r, flows));
        }

        // Wells whose control letter changed since the last described iteration.
        if (diag_letters_.size() == set.size()) {
            for (int w = 0; w < numWells(); ++w) {
                if (diag_letters_[w] != set[w]) {
                    out.push_back(fmt::format("  changed {} -> {}:", diag_letters_[w], set[w])
                                  + describeWell(w, x));
                }
            }
        }
        diag_letters_ = set;

        // The analytic Jacobian against a forward difference of the residual,
        // with the same steps NetworkSolve::solve() takes when it differences.
        if (usesAnalyticJacobian()) {
            const auto J = jacobian(x);
            Scalar worst_mismatch = 0;
            int wi = 0;
            int wj = 0;
            Scalar a_val = 0;
            Scalar fd_val = 0;
            for (int j = 0; j < n; ++j) {
                const Scalar h_nominal = Scalar{1e-2} * columnScale(j);
                State unit_dx(n, Scalar{0});
                unit_dx[j] = h_nominal;
                auto limited = limitStep(x, unit_dx);
                if (std::abs(limited[j]) < Scalar{1e-6} * std::abs(h_nominal)) {
                    unit_dx[j] = -h_nominal;
                    limited = limitStep(x, unit_dx);
                }
                const Scalar h = limited[j];
                if (h == Scalar{0}) { continue; }
                State shifted = x;
                for (int k = 0; k < n; ++k) { shifted[k] += limited[k]; }
                const auto rj = residual(shifted);
                for (int i = 0; i < n; ++i) {
                    const Scalar fd = (rj[i] - r[i]) / h;
                    // Mismatch relative to the column's scale, so columns in
                    // Pa and dimensionless lambdas compare alike.
                    const Scalar mismatch = std::abs(J(i, j) - fd) * columnScale(j);
                    if (mismatch > worst_mismatch) {
                        worst_mismatch = mismatch;
                        wi = i; wj = j; a_val = J(i, j); fd_val = fd;
                    }
                }
            }
            // fd is a forward difference with a step of 1% of the column
            // scale: a relative error of about 1e-3 or less is its own
            // truncation error, a large one is a wrong (or discontinuous)
            // derivative.
            const Scalar rel = std::abs(a_val - fd_val)
                / std::max({std::abs(a_val), std::abs(fd_val), std::numeric_limits<Scalar>::min()});
            out.push_back(fmt::format("  jacobian check: worst |analytic - fd| * scale {:.3e} (relative {:.2e}) "
                                      "at row {}, column {} (analytic {:.4e}, fd {:.4e})",
                                      worst_mismatch, rel, describeIndex(wi), describeIndex(wj),
                                      a_val, fd_val));
        }
        return out;
    }

    std::vector<std::string> describeState(const State& x) const override
    {
        std::vector<std::string> out;
        out.push_back(fmt::format("Network diag final state: terminal {:.3f} bar",
                                  terminal_pressure_ / unit::barsa));
        const auto flows = nodeFlows(x);
        const auto r = residual(x);
        for (int i = 0; i < size(); ++i) {
            if (i < numNodes() + numActiveLambdas()) {
                out.push_back(describeRow(i, x, r, flows));
            }
        }
        for (int w = 0; w < numWells(); ++w) {
            out.push_back(describeWell(w, x));
        }
        return out;
    }

    State pressures(const State& x) const override
    {
        State p(nodes_.size());
        p[0] = terminal_pressure_;
        for (int i = 1; i <= numNodes(); ++i) { p[i] = x[pIdx(i)]; }
        return p;
    }
    State wellRates(const State& x) const override
    {
        State q(numWells());
        for (int w = 0; w < numWells(); ++w) { q[w] = wellPhaseRatesOwn(w, x)[kOil]; }
        return q;
    }
    std::vector<std::array<Scalar, 3>> wellPhaseRates(const State& x) const override
    {
        std::vector<std::array<Scalar, 3>> q(numWells());
        for (int w = 0; w < numWells(); ++w) { q[w] = wellPhaseRatesOwn(w, x); }
        return q;
    }
    State wellBhps(const State& x) const override
    {
        State bhp(numWells(), Scalar{0});
        for (int w = 0; w < numWells(); ++w) {
            const auto& well = wells_[w];
            if (well.kind == WellKind::Thp) {
                bhp[w] = x[thpBhpIdx(w)];
            } else {
                // From the phase whose rate responds most to bhp (not oil
                // as such: a gas well has none).
                const auto qw = wellPhaseRatesOwn(w, x);
                int p = 0;
                for (int k = 1; k < NP; ++k) {
                    if (std::abs(well.ipr_b[k]) > std::abs(well.ipr_b[p])) { p = k; }
                }
                bhp[w] = (well.ipr_b[p] != Scalar{0}) ? (qw[p] - well.ipr_a[p]) / well.ipr_b[p] : Scalar{0};
            }
        }
        return bhp;
    }

private:
    mutable std::string diag_letters_;   // TEMPORARY: see describeIteration()

    int numThpWells() const { return num_thp_; }
    int pIdx(const int node) const { return node - 1; }
    // Only ever called for a node with lambda_slot_[active_node] >= 0: the
    // only way to reach it is via a well in that node's own own_wells (the
    // Group case in wellPhaseRatesOwn()) or the node's own row in residual(),
    // both of which are exactly the two places that only exist because
    // own_wells is non-empty.
    int lambdaIdx(const int active_node) const { return numNodes() + lambda_slot_[active_node]; }
    int thpBhpIdx(const int w) const { return thp_col_[w]; }

    /// Every node's (efficiency-scaled) inflow at x, bottom-up: its own wells'
    /// rates plus its children's flows. Index 0 (the terminal) is unused.
    std::vector<std::array<Scalar, NP>> nodeFlows(const State& x) const
    {
        std::vector<std::array<Scalar, NP>> q(nodes_.size());
        for (int i = numNodes(); i >= 1; --i) {
            std::array<Scalar, NP> qi{};
            for (const int w : wells_at_[i]) {
                const auto qw = wellPhaseRatesOwn(w, x);
                for (int p = 0; p < NP; ++p) { qi[p] += wells_[w].efficiency * qw[p]; }
            }
            for (const int c : children_[i]) {
                for (int p = 0; p < NP; ++p) { qi[p] += nodes_[c].efficiency * q[c][p]; }
            }
            q[i] = qi;
        }
        return q;
    }

    /// A tubing-table lookup with its derivatives: the value and d/dthp from
    /// bhp_with_slope_limit() (the partials of the very line the value comes
    /// from), d/dq by automatic differentiation of the same lookup in the three
    /// phase rates, which also carries how the water and gas fractions move
    /// with them. With no limit the lookup is the plain table (the limit can
    /// never trigger), but with the unclipped FLO derivative -- the templated
    /// bhp() clips it at zero, which would describe a different function than
    /// the value on a downward-sloping part of the table.
    /// Explicit WFR/GFR handed to the VFP lookup, and whether they apply
    /// everywhere (fixed) or only below the table's first FLO value.
    struct LookupFractions
    {
        Scalar wfr = 0;
        Scalar gfr = 0;
        bool fixed = false;
    };

    struct TableDerivatives
    {
        Scalar value = 0;
        Scalar dthp = 0;
        std::array<Scalar, NP> dq{};   // d/d(own positive rate), oil, water, gas
    };
    TableDerivatives tableLookup(const int table, const Scalar thp,
                                 const std::array<Scalar, NP>& q, const Scalar alq,
                                 const std::optional<Scalar> max_slope,
                                 const LookupFractions& f) const
    {
        const Scalar limit = max_slope.value_or(std::numeric_limits<Scalar>::lowest());
        // props_->bhp*() want water, oil, gas and negative-for-production.
        const auto plain = props_->bhp_with_slope_limit(table, -q[kWater], -q[kOil], -q[kGas],
                                                        thp, alq, f.wfr, f.gfr, f.fixed,
                                                        limit);
        using Eval = DenseAd::Evaluation<Scalar, NP>;
        const Eval aqua = Eval::createVariable(-q[kWater], kWater);
        const Eval liquid = Eval::createVariable(-q[kOil], kOil);
        const Eval vapour = Eval::createVariable(-q[kGas], kGas);
        const Eval ad = props_->bhp_with_slope_limit(table, aqua, liquid, vapour, thp, alq,
                                                     f.wfr, f.gfr, f.fixed, limit);
        TableDerivatives out;
        out.value = plain.evaluation.value;
        out.dthp = plain.evaluation.dthp;
        for (int p = 0; p < NP; ++p) {
            out.dq[p] = -ad.derivative(p);   // d/dq = -d/d(-q)
        }
        return out;
    }

    /// Precompute the constant derivatives jacobian() needs (called by
    /// finalize()): each well's rates with respect to the unknowns, each
    /// node's inflow from those, bottom-up with efficiencies, and each target
    /// row -- already divided by the row's scale, as residual() does.
    void buildRateDerivatives()
    {
        using Sparse = std::map<int, std::array<Scalar, NP>>;
        std::vector<Sparse> rate(wells_.size());
        for (int w = 0; w < numWells(); ++w) {
            const auto& well = wells_[w];
            if (well.kind == WellKind::Group) {
                // q_p = a_p + b_p * (g * lambda - A) / B, A and B the IPR
                // projected on the node's mode (bhpFromTarget()).
                const auto& a = activeNodes_[well.active_node];
                const Scalar B = projectOnMode(well.ipr_b, phaseWeights(a.mode, a.resv_coeff));
                if (B != Scalar{0}) {
                    auto& d = rate[w][lambdaIdx(well.active_node)];
                    for (int p = 0; p < NP; ++p) { d[p] = well.ipr_b[p] * well.guide_rate / B; }
                }
            } else if (well.kind == WellKind::Thp) {
                rate[w][thpBhpIdx(w)] = well.ipr_b;
            }
        }

        std::vector<Sparse> flow(nodes_.size());
        for (int i = numNodes(); i >= 1; --i) {
            for (const int w : wells_at_[i]) {
                for (const auto& [col, d] : rate[w]) {
                    auto& f = flow[i][col];
                    for (int p = 0; p < NP; ++p) { f[p] += wells_[w].efficiency * d[p]; }
                }
            }
            for (const int c : children_[i]) {
                for (const auto& [col, d] : flow[c]) {
                    auto& f = flow[i][col];
                    for (int p = 0; p < NP; ++p) { f[p] += nodes_[c].efficiency * d[p]; }
                }
            }
        }
        node_flow_jac_.assign(nodes_.size(), {});
        for (std::size_t i = 1; i < nodes_.size(); ++i) {
            node_flow_jac_[i].assign(flow[i].begin(), flow[i].end());
        }

        // Target rows: activeNodeTotal() with the node's own weights, divided
        // by the same scale residual() uses.
        target_jac_.assign(activeNodes_.size(), {});
        for (int k = 0; k < numActiveNodes(); ++k) {
            if (lambda_slot_[k] < 0) {
                continue;
            }
            const auto& a = activeNodes_[k];
            const auto weights = phaseWeights(a.mode, a.resv_coeff);
            const Scalar scale = (a.target > Scalar{0}) ? a.target : Scalar{1};
            std::map<int, Scalar> row;
            std::function<void(int, Scalar)> add = [&](const int node, const Scalar factor) {
                const auto& an = activeNodes_[node];
                auto addWell = [&](const int w, const Scalar eff) {
                    for (const auto& [col, d] : rate[w]) {
                        row[col] += factor * eff * projectOnMode(d, weights) / scale;
                    }
                };
                for (const int w : an.own_wells) { addWell(w, wells_[w].efficiency); }
                for (const auto& [w, eff] : an.member_wells) { addWell(w, eff); }
                for (const auto& [child, eff] : an.active_children) { add(child, factor * eff); }
            };
            add(k, Scalar{1});
            target_jac_[k].assign(row.begin(), row.end());
        }
    }
    bool hasTable(const Node& n) const { return n.vfp_table != NoTable; }

    /// Node \p i's table at (upstream, q), flattened if setNodeFlattening().
    Scalar nodeTableBhp(const int i, const Scalar upstream, const std::array<Scalar, NP>& q) const
    {
        if (!node_flattened_[i]) {
            return tableBhp(nodes_[i].vfp_table, upstream, q, node_alq_[i], nodeFractions(i));
        }
        const auto f = nodeFractions(i);
        return props_->bhp_with_slope_limit(nodes_[i].vfp_table, -q[kWater], -q[kOil], -q[kGas],
                                            upstream, node_alq_[i], f.wfr, f.gfr, f.fixed, Scalar{0})
            .evaluation.value;
    }

    std::optional<Scalar> nodeSlopeLimit(const int i) const
    {
        return node_flattened_[i] ? std::optional<Scalar>{Scalar{0}} : std::nullopt;
    }

    Scalar tableBhp(const int table, const Scalar thp, const std::array<Scalar, NP>& q, const Scalar alq,
                    const LookupFractions& f) const
    {
        // props_->bhp() wants water, oil, gas (aqua, liquid, vapour); q here is
        // oil, water, gas (ProdGroupTreeBalancer's order) -- reordered right here,
        // the one place the two conventions meet.
        return props_->bhp(table, -q[kWater], -q[kOil], -q[kGas], thp, alq, f.wfr, f.gfr, f.fixed);
    }

    /// What a well's tubing lookup is evaluated with: the fractions of its
    /// fixed phase proportions (Well::fractions) at every FLO -- the same as
    /// its rates' own wherever they are defined; without fixed proportions,
    /// the iterate's own, and its explicit ones below the first FLO value.
    LookupFractions wellFractions(const Well& well, const std::array<Scalar, NP>& /*q*/) const
    {
        if (!well.fractions.has_value()) {
            return {well.explicit_fractions[0], well.explicit_fractions[1], false};
        }
        return {(*well.fractions)[0], (*well.fractions)[1], true};
    }

    /// What node \p i's table is evaluated with: the fractions of its own
    /// flow, and below the table's first FLO value its fallback fractions.
    LookupFractions nodeFractions(const int i) const
    {
        const auto& f = node_explicit_fractions_[i];
        return f.has_value() ? LookupFractions{(*f)[0], (*f)[1], false} : LookupFractions{};
    }

    /// The phase-rate weights a given target mode measures, in this class's
    /// [oil, water, gas] order -- mirrors ProdGroupTreeBalancer.cpp's
    /// projectOnMode()/modeWeights(), duplicated rather than shared since it
    /// is a two-line function and the two files intentionally don't depend on
    /// each other beyond ProdGroupTreeBalancer's own public header.
    static std::array<Scalar, NP> phaseWeights(const ::Opm::Well::ProducerCMode mode,
                                               const std::array<Scalar, NP>& resv_coeff)
    {
        switch (mode) {
        case ::Opm::Well::ProducerCMode::ORAT: return {1, 0, 0};
        case ::Opm::Well::ProducerCMode::WRAT: return {0, 1, 0};
        case ::Opm::Well::ProducerCMode::GRAT: return {0, 0, 1};
        case ::Opm::Well::ProducerCMode::LRAT: return {1, 1, 0};
        case ::Opm::Well::ProducerCMode::RESV: return resv_coeff;
        default: return {0, 0, 0};
        }
    }
    static Scalar projectOnMode(const std::array<Scalar, NP>& q, const std::array<Scalar, NP>& weights)
    {
        return weights[kOil] * q[kOil] + weights[kWater] * q[kWater] + weights[kGas] * q[kGas];
    }

    /// Invert q_p = ipr_a[p] + ipr_b[p]*bhp, projected onto `weights`, for the
    /// bhp that gives target on that projection -- closed form since the
    /// projection of an affine function is affine. Used for every one of the
    /// five modes phaseWeights() can return, not just a single phase.
    static Scalar bhpFromTarget(const Well& well, const std::array<Scalar, NP>& weights, const Scalar target)
    {
        const Scalar a = projectOnMode(well.ipr_a, weights);
        const Scalar b = projectOnMode(well.ipr_b, weights);
        return (b != Scalar{0}) ? (target - a) / b : Scalar{0};
    }

    /// The bhp at which every phase's IPR rate is exactly zero -- see
    /// Well::bhp_shutin's own doc comment. Picks the phase with the largest
    /// |ipr_b| to invert, for robustness against a phase that happens to
    /// carry a near-zero (but not exactly zero) coefficient; the constant-
    /// fraction IPR assumption this whole design already relies on (see
    /// groups_and_network_clean.md) means every active phase agrees on this
    /// bhp anyway, so which one is used is only a numerical-conditioning
    /// choice, not a modeling one.
    static Scalar shutInBhp(const Well& well)
    {
        int best = 0;
        for (int p = 1; p < NP; ++p) {
            if (std::abs(well.ipr_b[p]) > std::abs(well.ipr_b[best])) {
                best = p;
            }
        }
        return (well.ipr_b[best] != Scalar{0}) ? -well.ipr_a[best] / well.ipr_b[best] : Scalar{0};
    }

    const VFPProdProperties<Scalar>* props_;
    std::vector<Node> nodes_{Node{}};   // index 0: the terminal (parent == -1)
    std::vector<Scalar> node_alq_{Scalar{0}};   // parallel to nodes_; see addNode()
    std::vector<bool> node_flattened_{false};   // parallel to nodes_; see setNodeFlattening()
    std::vector<std::optional<std::array<Scalar, 2>>> node_explicit_fractions_{std::nullopt};   // see setNodeExplicitFractions()
    std::vector<std::optional<Scalar>> fixed_pressure_{std::nullopt};   // parallel to nodes_; see setFixedPressure()
    std::vector<Well> wells_;
    std::vector<ActiveNode> activeNodes_;
    Scalar terminal_pressure_ = 0;

    // Set by finalize(): lambda_slot_[k] is k's position in the lambda block
    // of the state vector, or -1 if k has no lambda unknown at all (empty
    // own_wells -- see finalize()'s own comment).
    std::vector<int> lambda_slot_;
    int num_lambdas_ = 0;

    // Set by finalize() (sized, all false) and updateControls() (updated
    // every iterate): thp_capped_[w] is only meaningful for a Thp-kind well,
    // true once its current bhp has reached bhp_shutin -- see that field's
    // own doc comment and residual()'s use of it.
    std::vector<bool> thp_capped_;

    // Set by finalize(): wells feeding each node, each node's children, each
    // well's Thp column (-1 if not Thp), and the constant derivatives
    // jacobian() uses (see buildRateDerivatives()).
    std::vector<std::vector<int>> wells_at_;
    std::vector<std::vector<int>> children_;
    std::vector<int> thp_col_;
    int num_thp_ = 0;
    std::vector<std::vector<std::pair<int, std::array<Scalar, NP>>>> node_flow_jac_;
    std::vector<std::vector<std::pair<int, Scalar>>> target_jac_;
    bool analytic_jacobian_ = false;
};

} // namespace Opm::NetworkSolve

#endif // OPM_NETWORK_GROUP_TREE_SYSTEM_HEADER_INCLUDED
