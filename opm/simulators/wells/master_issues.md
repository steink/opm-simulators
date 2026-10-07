# Problems found in master

Issues in existing `opm-simulators` code, found while working on the
group-tree network branch, that are not specific to that branch. Each entry
states how it was found, what goes wrong, the cause, and a proposed fix.
Nothing here has been fixed in master unless stated.

---

## 1. Fixed-point network computation leaves part of a multi-root network at 0 bar

**Where:** `NetworkPressureComputation::run()` / `computeNodePressures()`
([BlackoilWellModelNetworkPressureComputation.hpp](BlackoilWellModelNetworkPressureComputation.hpp)).

**Found with:** `opm-tests/network/NETWORK-01-MULTIROOT.DATA`, default
network solver (fixed-point). The deck is not among the regression tests in
`regressionTests.cmake` (only NETWORK-01, NETWORK-01_STANDARD,
NETWORK-01-REROUTE, NETWORK-01-REROUTE_STD and NETWORK-01-WTEST are), so the
topology may not be regarded as fully supported.

**Topology:** FIELD (fixed at 80 bar) ← GRPA ← PRODA (PROD1, PROD3), and
GRPA ← GRPB (fixed at 82 bar) ← PRODB (PROD2). GRPB is a fixed-pressure node
that also feeds GRPA. The branches PRODA→GRPA and GRPA→FIELD have table 9999
(no pressure loss).

**What goes wrong:** PRODA and GRPA get a pressure of **0 bar** throughout the
run (GPR:PRODA, GPR:GRPA, and the "Network node … pressure" debug lines).
PROD1 therefore has no network THP limit and produces its full 5000 m3/day at
a THP of about 72 bar, below the 80 bar that its node, connected to FIELD
without pressure loss, must have. The saved run in `opm-tests/network/tmp_master`
shows the same zeros, so this is master behaviour, not something introduced
on the branch.

**Cause:**
1. `ExtNetwork::roots()` returns the uptree node of every branch whose uptree
   node has a fixed pressure, in branch order. Here that is GRPB before FIELD
   (and a root is listed once per branch below it).
2. `run()` processes the roots in that order. For GRPB, `computeNodePressures()`
   sets GRPB to its terminal pressure and, because GRPB has an uptree branch,
   computes that branch's data from `node_pressures_[(*upbranch).uptree_node()]`.
   `operator[]` inserts GRPA with the value 0.
3. When FIELD's tree is processed next, the "do not traverse subtree more than
   once" check finds GRPA already in `node_pressures_` and skips it. PRODA is
   then computed from GRPA's 0 bar through its no-loss branch.

**Proposed fix:** process the roots that have no branch above them first
(a fixed-pressure node with an uptree branch is then reached, and computed,
from its parent's tree, and its own pass finds everything done). Also read the
upstream pressure with `find()` rather than `operator[]`, so a missing upstream
value cannot be created as 0. The branch uses the same root selection for the
group-tree solve (`BlackoilWellModelNetworkGeneric::independentRoots()`).

**Effect of fixing:** changes fixed-point results for any deck with a
fixed-pressure node below the root. With the group-tree solver, which handles
this topology correctly, NETWORK-01-MULTIROOT has PRODA = GRPA = 80 bar and
PROD1 on THP control at about 4300 m3/day; FOPT is 6.5% lower than the
fixed-point run.

## Network node pressures below the first FLO knot use WFR = GFR = 0

`NetworkVfpPressureCalculator::compute()`
(`BlackoilWellModelNetworkPressureComputation.hpp`) calls
`VFPProdProperties::bhp()` with explicit WFR = GFR = 0 and
`use_expvfp = false` ("we dont support explicit lookup"). Below a table's
first FLO knot `bhp()` switches to the explicit fractions, so a branch
carrying little or no flow is evaluated as dry, gas-free oil: a much heavier
column and a jump in the curve at the knot. Affects nodes whose wells are
shut in or nearly so. A fix needs fractions for the node, e.g. its wells'
explicit fractions weighted by their rates or potentials.

