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

#ifndef OPM_PROD_GROUP_TREE_BALANCER_HEADER_INCLUDED
#define OPM_PROD_GROUP_TREE_BALANCER_HEADER_INCLUDED

#include <opm/simulators/wells/ProdGroupTreeNode.hpp>

#include <array>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace Opm {

class DeferredLogger;
class GuideRate;
class SummaryState;
template<typename Scalar, typename IndexTraits> class BlackoilWellModelGeneric;

} // namespace Opm

namespace Opm::ProdGroupTreeBalancer {

/// Type alias for the tree map.
template<class Scalar>
using Tree = std::map<std::string, ProdGroupTreeNode<Scalar>>;

/// A well tied directly to an Active node's own lambda: the *node's* own
/// equation sees efficiency * guideRate * lambda (efficiency-scaled, like any
/// other contribution to its sum), but the *well's* own rate -- what its own
/// IPR row actually has to invert to a bhp -- is guideRate * lambda alone:
/// efficiency describes how much of the well's own production reaches the
/// Active node, not the well's own physical rate, so it must not be baked
/// into a single product the way it safely can be for thpWells/activeChildren
/// (their own rate comes from somewhere else entirely -- their own bhp row,
/// or their own separate target -- so nothing there ever needs unscaling).
template<class Scalar>
struct FlatWellShare
{
    std::string name;
    Scalar guideRate;
    Scalar efficiency;
};

/// A reference from one Active node's own sum to another Active node found
/// elsewhere in the same FlatNetworkInput (by name) -- either a well whose
/// own limit binds (or, Part 1b, is stopped, or is on the network's live THP
/// control -- see FlatActiveNode::networkThp) or a nested group with its own
/// separate target. Either way its own total is still physically part of
/// this node's subtree and must count toward this node's sum, just scaled by
/// the cumulative efficiency factor between the two, and NOT via this node's
/// own lambda (the referenced node has its own row -- pinned, a live THP row,
/// or its own group equation -- that already determines it).
template<class Scalar>
struct FlatChildRef
{
    std::string name;
    Scalar efficiency;
};

/// One "Active" node -- a well or group whose own limit is what is actually
/// binding (ProdNodeModeCategory::Individual) -- in the tree extractFlatNetworkInput()
/// produces from an already-balanced Tree. Active nodes are the only ones with an
/// equation of their own; everything strictly between two of them (a group whose
/// own modeCategory is Group, None or Transparent -- nothing of its own binds) is
/// pass-through and never appears here at all: its wells are folded straight into
/// the nearest Active ancestor's ownWells.
///
/// A group's own target equation sums two kinds of contribution: ownWells
/// (tied to this node's own lambda -- always ProdNodeModeCategory::Group wells;
/// GRUP and THP are mutually exclusive on a well's control mode, so a well
/// counted here is never also a THP well) and activeChildren -- every other
/// Active node (well or nested group) found within this node's original
/// subtree, referenced by name into the same FlatNetworkInput rather than
/// further flattened, since each one's own total is already determined by its
/// own separate entry. The flattening only ever collapses pass-through
/// layers; a genuine second binding constraint nested inside the first (of
/// either kind) keeps its own entry.
///
/// A well-type entry (type == Well) has empty ownWells/activeChildren, and is
/// exactly one of:
///  - genuinely limit-bound (target > 0, pinned at that limit);
///  - currently stopped (target == 0, from a zero Limits entry -- Part 1b);
///  - on the network's live THP control (mode == THP && networkThp == true):
///    target is not meaningful (there is nothing to pin -- its rate comes
///    from its own bhp/thp/IPR row, genuinely responding to pressure). The
///    balancer itself has no notion of this distinction (mode == THP alone
///    covers both a fixed deck THP limit and a live network one); the caller
///    supplies it via extractFlatNetworkInput()'s networkThpWells set.
/// Which of the three applies is a fact the caller reads off target/mode/
/// networkThp, since this function has no access to well state or IPR data
/// to decide it itself.
///
/// A satellite-group entry (type == Group, satelliteRates has a value) is a
/// fourth, well-*like* case: a GSATPROD group has no wells of its own and no
/// network node either (a satellite rate stands in for wells never modelled
/// at all) -- its own subtree contributes a fixed, known rate to whichever
/// ancestor references it, exactly like a pinned well's target does, just
/// carried as the full three-phase rate rather than a single mode-projected
/// scalar (a satellite's ancestor can be on any mode, not necessarily the
/// satellite's own). ownWells/activeChildren are empty for this case too.
template<class Scalar>
struct FlatActiveNode
{
    std::string name;
    ProdNodeType type{ProdNodeType::Well};
    Well::ProducerCMode mode{Well::ProducerCMode::CMODE_UNDEFINED};
    Scalar target{0};
    bool networkThp{false};   // only meaningful when type == Well && mode == THP
    std::array<Scalar, 3> resvCoeff{};   // [oil, water, gas]; RESV mode only, own or (if type == Well) ancestor's
    std::optional<std::array<Scalar, 3>> satelliteRates;   // [oil, water, gas], positive = production
    std::vector<FlatWellShare<Scalar>> ownWells;
    std::vector<FlatChildRef<Scalar>> activeChildren;
};

/// One entry per Active node in the tree, roots (nodes whose nearest Active
/// ancestor is none, i.e. nothing above them binds either) first.
template<class Scalar>
using FlatNetworkInput = std::vector<FlatActiveNode<Scalar>>;

/// Flatten an already-balanced tree (as produced by balanceTreeForTesting() or
/// runGroupTreeBalancer()) into the reduced tree of Active nodes described above,
/// starting the search from \p rootName. A read-only pass over the tree; the only
/// reason it also takes \p guideRate is that a well's own guide rate -- needed for
/// ownWells' allocation weight -- is not reliably left on the tree node itself by
/// balancing and has to be looked up the same way the balancer looks it up.
///
/// \p networkThpWells names every Individual, mode == THP well the caller knows
/// has a *network*-sourced (as opposed to deck-sourced) THP limit -- in the real
/// integration, well->getDynamicThpLimit().has_value(). The balancer itself has
/// no way to know this: it only ever populates mode == THP from the well's
/// current strictest limit, with no notion of where that limit came from. A
/// well named here gets FlatActiveNode::networkThp == true instead of an
/// ordinary pin. Deliberately a plain name set rather than a WellInterface
/// reference, so this header stays free of any well-model-specific type; the
/// caller does that one lookup itself. Never consulted for a Group-category
/// well: GRUP and THP are mutually exclusive on a well's control mode, so that
/// case cannot arise.
template<class Scalar>
FlatNetworkInput<Scalar> extractFlatNetworkInput(const Tree<Scalar>& tree,
                                                 const std::string& rootName,
                                                 const GuideRate& guideRate,
                                                 const std::unordered_set<std::string>& networkThpWells = {});

/// Top-level entry point: build tree, balance it, validate, and apply.
/// All internal functions are implementation details not exposed through this interface.
///
/// \param[in]    wellModel     Well model (schedule, well/group state, guide rate, RESV coefficients)
/// \param[in]    summaryState  Summary state
/// \param[in]    reportStep    Current report step index
/// \param[in]    tol           Convergence tolerance
/// \param[in]    limits        Globally gathered well limits (from prepareWellsForBalancing_*)
/// \param[in]    logger        Deferred logger
/// \return       true if the result passed checkTreeValidity
/// Balance a tree that was built by hand rather than from a well model, so the
/// algorithm can be used as an oracle in a standalone test. Everything else in
/// this file is untouched; this is the only addition.
template<class Scalar>
bool balanceTreeForTesting(Tree<Scalar>& tree,
                           const GuideRate& guideRate,
                           Scalar tol,
                           DeferredLogger& logger);

/// The balanced tree runGroupTreeBalancer() builds internally, exposed so a
/// caller that needs the tree itself (e.g. to flatten it via
/// extractFlatNetworkInput() for a real network-pressure solve) can get at it
/// without also going through runGroupTreeBalancer()'s own applyTreeToState()
/// write-back, which is only appropriate for the no-network case.
template<class Scalar>
struct BalancedTree
{
    Tree<Scalar> tree;
    bool success{true};   ///< runBalancingAlgorithm() converged (true if there was nothing to balance)
    bool valid{true};     ///< checkTreeValidity() passed (true if there was nothing to balance)
};

/// buildTree() + runBalancingAlgorithm() + checkTreeValidity(), with the same
/// rank-0 logging runGroupTreeBalancer() itself does -- everything
/// runGroupTreeBalancer() does short of the final applyTreeToState() write-back.
/// An empty \p limits (no active wells / no wells with positive potentials)
/// short-circuits to an empty, trivially-valid result, same as
/// runGroupTreeBalancer()'s own early return.
template<class Scalar, typename IndexTraits>
BalancedTree<Scalar> balanceGroupTree(BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                                      const SummaryState& summaryState,
                                      int reportStep,
                                      Scalar tol,
                                      const std::unordered_map<std::string, std::pair<int, Scalar>>& limits,
                                      DeferredLogger& logger);

/// Commit a balanced tree to WellState/GroupState: each producer's
/// production_cmode (GRUP for a Group-category well, its own limiting mode for
/// an Individual one), its group target and fallback target, and each group's
/// production control mode.
///
/// \p commitRates additionally writes the tree's well and group rates. That is
/// right when nothing else decides the rates (no production network); with a
/// network the rates come from the network and local well solves, and the
/// tree's rates for a THP well are only the balancer's estimate of its limit,
/// so the caller passes false.
///
/// \return true if any local well's or any group's production control mode
/// changed. The well part is rank-local, so the caller must reduce it over
/// ranks.
template<class Scalar, typename IndexTraits>
bool applyTreeToState(const Tree<Scalar>& tree,
                      BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                      DeferredLogger& logger,
                      bool commitRates = true);

template<class Scalar, typename IndexTraits>
bool runGroupTreeBalancer(BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                          const SummaryState& summaryState,
                          int reportStep,
                          Scalar tol,
                          const std::unordered_map<std::string, std::pair<int, Scalar>>& limits,
                          DeferredLogger& logger);

} // namespace Opm::ProdGroupTreeBalancer

#endif // OPM_PROD_GROUP_TREE_BALANCER_HEADER_INCLUDED
