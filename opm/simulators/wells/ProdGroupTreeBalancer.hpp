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

#include <map>
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

/// A well tied directly to an Active node's own lambda: q_w = weight * lambda,
/// weight = the well's guide rate scaled by its cumulative efficiency factor
/// up to the Active node.
template<class Scalar>
struct FlatWellShare
{
    std::string name;
    Scalar weight;
};

/// A reference from one Active node's own sum to another Active node found
/// elsewhere in the same FlatNetworkInput (by name) -- either a well whose
/// own limit binds (or, Part 1b, is stopped) or a nested group with its own
/// separate target. Either way its own total is still physically part of
/// this node's subtree and must count toward this node's sum, just scaled by
/// the cumulative efficiency factor between the two, and NOT via this node's
/// own lambda (the referenced node has its own row -- pinned, a THP row on a
/// trial IPR, or its own group equation -- that already determines it).
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
/// A group's own target equation sums three kinds of contribution: ownWells
/// (tied to this node's own lambda), thpWells (see below), and activeChildren
/// -- every other Active node (well or nested group) found within this node's
/// original subtree, referenced by name into the same FlatNetworkInput rather
/// than further flattened, since each one's own total is already determined
/// by its own separate entry. The flattening only ever collapses pass-through
/// layers; a genuine second binding constraint nested inside the first (of
/// either kind) keeps its own entry.
///
/// thpWells are wells this node's own guide-rate allocation would otherwise
/// claim (ProdNodeModeCategory::Group) but that are actually on the network's
/// THP control -- their rate is not g_w*lambda at all, it comes from their own
/// bhp/thp/IPR row, genuinely responding to the network's pressures. They
/// still count toward this node's sum (their production is still physically
/// part of this node's subtree), just via weight = cumulative efficiency alone
/// (no guide-rate allocation applies to a rate this node doesn't get to set).
/// The balancer itself has no notion of THP control -- see extractFlatNetworkInput().
///
/// A well-type entry (type == Well) has empty ownWells/thpWells/activeChildren:
/// it is either genuinely limit-bound (target > 0, pinned) or currently stopped
/// (target == 0, from a zero Limits entry) -- which one is a fact the caller
/// reads off target, since this function has no access to well state or IPR
/// data to decide it itself.
template<class Scalar>
struct FlatActiveNode
{
    std::string name;
    ProdNodeType type{ProdNodeType::Well};
    Well::ProducerCMode mode{Well::ProducerCMode::CMODE_UNDEFINED};
    Scalar target{0};
    std::vector<FlatWellShare<Scalar>> ownWells;
    std::vector<FlatWellShare<Scalar>> thpWells;
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
/// \p thpControlledWells names every well the caller knows is on the network's
/// THP control (e.g. from WellState's production_cmode) -- a fact the balancer
/// itself has no way to know, since it was written without the network in mind
/// and only ever distinguishes Individual from Group. Deliberately a plain name
/// set rather than a WellState reference, so this header stays free of any
/// network-specific type; the caller does that one lookup itself.
template<class Scalar>
FlatNetworkInput<Scalar> extractFlatNetworkInput(const Tree<Scalar>& tree,
                                                 const std::string& rootName,
                                                 const GuideRate& guideRate,
                                                 const std::unordered_set<std::string>& thpControlledWells = {});

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

template<class Scalar, typename IndexTraits>
bool runGroupTreeBalancer(BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                          const SummaryState& summaryState,
                          int reportStep,
                          Scalar tol,
                          const std::unordered_map<std::string, std::pair<int, Scalar>>& limits,
                          DeferredLogger& logger);

} // namespace Opm::ProdGroupTreeBalancer

#endif // OPM_PROD_GROUP_TREE_BALANCER_HEADER_INCLUDED
