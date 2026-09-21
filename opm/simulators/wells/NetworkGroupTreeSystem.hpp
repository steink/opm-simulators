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

#include <opm/simulators/wells/FlattenedTubingCurve.hpp>
#include <opm/simulators/wells/NetworkSolve.hpp>
#include <opm/simulators/wells/ProdGroupTreeBalancer.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <array>
#include <cmath>
#include <functional>
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
/// with a Well::tilde_lambda set solves its row against that flattened curve
/// instead of the real table -- see thpWellResidualRow(). No autochoke
/// either. Every well has one of three kinds:
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

        // Thp only, and optional even then: Part 3's flattened tubing curve,
        // built by the caller from this well's own table/ipr whenever those
        // are (re)computed (not by this class). Unset means "use the real
        // table directly, unflattened" -- the only behaviour before this
        // field existed, still exercised by every earlier Thp test in
        // test_networkgrouptreesystem.cpp, and still a legitimate choice for
        // a well known not to have a lift cliff. See thpWellResidualRow().
        std::optional<FlattenedTubingCurve<Scalar>> tilde_lambda;
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
        return static_cast<int>(nodes_.size()) - 1;
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
    };

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
            well.ipr_a = wd.ipr_a;
            well.ipr_b = wd.ipr_b;
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
    }

    int numNodes() const { return static_cast<int>(nodes_.size()) - 1; }   // excludes the terminal
    int numActiveNodes() const { return static_cast<int>(activeNodes_.size()); }
    /// Active nodes that actually got a lambda unknown -- see finalize().
    int numActiveLambdas() const { return num_lambdas_; }
    int numWells() const override { return static_cast<int>(wells_.size()); }
    const std::vector<Well>& wells() const { return wells_; }
    const std::vector<ActiveNode>& activeNodes() const { return activeNodes_; }
    Scalar terminalPressure() const { return terminal_pressure_; }

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
            x[i] = node_pressure_guess[i];
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
                x[thpBhpIdx(w)] = node_pressure_guess[wells_[w].node - 1];
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
        std::vector<std::array<Scalar, NP>> q(nodes_.size());
        for (int i = n; i >= 1; --i) {
            std::array<Scalar, NP> qi{};
            for (int w = 0; w < numWells(); ++w) {
                if (wells_[w].node != i) { continue; }
                const auto qw = wellPhaseRatesOwn(w, x);
                for (int p = 0; p < NP; ++p) { qi[p] += wells_[w].efficiency * qw[p]; }
            }
            for (std::size_t c = 1; c < nodes_.size(); ++c) {
                if (nodes_[c].parent != i) { continue; }
                for (int p = 0; p < NP; ++p) { qi[p] += nodes_[c].efficiency * q[c][p]; }
            }
            q[i] = qi;
        }
        for (int i = 1; i <= n; ++i) {
            const Scalar upstream = (nodes_[i].parent == 0) ? terminal_pressure_ : x[pIdx(nodes_[i].parent)];
            const Scalar computed = hasTable(nodes_[i])
                ? tableBhp(nodes_[i].vfp_table, upstream, q[i], node_alq_[i]) : upstream;
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
        const Scalar computed = well.tilde_lambda ? well.tilde_lambda->bhp(thp, qw)
                                                  : tableBhp(well.vfp_table, thp, qw, well.alq);
        return (bhp - computed) / unit::barsa;
    }

    DenseMatrix<Scalar> jacobian(const State&) const override
    {
        // No analytic Jacobian yet -- solve() differences residual() instead
        // (see usesAnalyticJacobian()). Never called; the trivial empty
        // matrix here is just to satisfy the pure virtual.
        return DenseMatrix<Scalar>(size());
    }
    bool usesAnalyticJacobian() const override { return false; }

    /// The one piece of per-iterate state this class has: whether each Thp
    /// well's *current* bhp has reached its own bhp_shutin (see that field's
    /// doc comment). Deliberately not the kind of combinatorial, multi-well,
    /// hysteresis-carrying switching this design otherwise avoids inside
    /// Newton -- it is a single scalar bound on one well's own unknown,
    /// decided fresh from the current iterate alone every time, the same way
    /// a standard bound-constrained ("projected") Newton method treats an
    /// active box constraint.
    bool updateControls(const State& x) override
    {
        bool moved = false;
        for (int w = 0; w < numWells(); ++w) {
            if (wells_[w].kind != WellKind::Thp) { continue; }
            const bool capped = x[thpBhpIdx(w)] >= wells_[w].bhp_shutin;
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
                const auto qw = wellPhaseRatesOwn(w, x);
                bhp[w] = (well.ipr_b[kOil] != Scalar{0}) ? (qw[kOil] - well.ipr_a[kOil]) / well.ipr_b[kOil] : Scalar{0};
            }
        }
        return bhp;
    }

private:
    int numThpWells() const
    {
        int n = 0;
        for (const auto& w : wells_) { n += (w.kind == WellKind::Thp) ? 1 : 0; }
        return n;
    }
    int pIdx(const int node) const { return node - 1; }
    // Only ever called for a node with lambda_slot_[active_node] >= 0: the
    // only way to reach it is via a well in that node's own own_wells (the
    // Group case in wellPhaseRatesOwn()) or the node's own row in residual(),
    // both of which are exactly the two places that only exist because
    // own_wells is non-empty.
    int lambdaIdx(const int active_node) const { return numNodes() + lambda_slot_[active_node]; }
    int thpBhpIdx(const int w) const
    {
        int idx = numNodes() + numActiveLambdas();
        for (int i = 0; i < w; ++i) { idx += (wells_[i].kind == WellKind::Thp) ? 1 : 0; }
        return idx;
    }
    bool hasTable(const Node& n) const { return n.vfp_table != NoTable; }

    Scalar tableBhp(const int table, const Scalar thp, const std::array<Scalar, NP>& q, const Scalar alq) const
    {
        // props_->bhp() wants water, oil, gas (aqua, liquid, vapour); q here is
        // oil, water, gas (ProdGroupTreeBalancer's order) -- reordered right here,
        // the one place the two conventions meet.
        return props_->bhp(table, -q[kWater], -q[kOil], -q[kGas], thp, alq, Scalar{0}, Scalar{0}, false);
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
};

} // namespace Opm::NetworkSolve

#endif // OPM_NETWORK_GROUP_TREE_SYSTEM_HEADER_INCLUDED
