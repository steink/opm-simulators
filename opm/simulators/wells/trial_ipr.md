# Trial IPRs for stopped wells

How a stopped producer gets a *trial IPR* in the group-tree workflow, and how
the network solve uses it to decide whether the well reopens. State as of the
working tree after "flowing trial IPR" (on top of 200dcf962).

Code: `WellInterface::updateStoppedWellTrialIpr()`, `computeAnchor()`,
`trialSolveAtThp()` (`WellInterface_impl.hpp`);
`WellBhpThpCalculator::maxFlowingThp()`, `VFPHelpers::maxFlowingThp()` /
`liftMargin()`; `BlackoilWellModelNetwork::refreshWellNetworkData()`,
`BlackoilWellModelNetworkGeneric::lowestNodePressure()`,
`gatherWellNetworkDataForGroupTree()`, `solveGroupTree()`,
`updateGroupTreeOpenSet()`; `GroupTreeSystem::addReopenCandidate()`,
`reopenOutcomes()` (`NetworkGroupTreeSystem.hpp`).

---

## 1. Why

A stopped well's own implicit IPR is linearised at zero rate, a degenerate
state, and is no use for deciding whether it could flow again. The network
solve needs an IPR the well would have *if it flowed*, at a point near where
it would run, to treat it as a candidate: a THP well that either flows at the
node pressure the solve finds or stays at its shut-in cap.

## 2. When it is computed

In `refreshWellNetworkData()`, at the start of every group-tree workflow round
(before B), for every predicting well, right after its implicit IPR refresh.
`updateStoppedWellTrialIpr()` first clears the previous trial IPR
(`SingleWellState::stopped_ipr_a/b`, `stopped_ipr_bhp`), then returns without
one unless the well is:

- a producer (injectors never get one);
- not SHUT, and stopped -- persistently (`ws.status == STOP`) or, far more
  often, dynamically (`wellIsStopped()`);
- not held stopped for the whole timestep (`isHeldStoppedForTimestep()`: its
  initial solve at the start of the step failed).

The well object and the well state are never changed by the trial: all solves
run on scratch copies, and the well's status, stop reason and operability are
restored afterwards. Only `stopped_ipr_a/b` and `stopped_ipr_bhp` are written.

The caller also passes the **lowest pressure the well's network node can
have** (`lowestNodePressure()`): the nearest fixed (terminal) pressure at or
above the well's group in the production network. Nullopt if the well's group
is not a network node (no network, or not in it).

## 3. Without a tubing table

`estimateOperableBhp()` on scratch copies, with explicit fractions: the
minimum bhp on the tubing curve at the THP limit, a BHP-controlled solve there,
and `estimateStableBhp()` from the IPR of that solve. If it finds an operable
bhp, the IPR of that solve is the trial IPR and that bhp its point; otherwise
no trial IPR. (Without a table there is no THP to speak of; this path is
mainly for completeness.)

## 4. With a tubing table: the anchor

`computeAnchor(start_thp = THP limit)` finds the well's **maximum flowing
THP** -- the highest THP at which its IPR still meets its tubing curve -- by
alternating well solves and table-only lift-margin computations.

**Start.** A stopped well is not flowing, so the first solve is at the bhp on
the tubing curve at the THP limit: `bhp0 = max(calculateMinimumBhpAtThp(THP
limit), bhp limit)` (the curve's minimum at that THP). The well is opened on the
scratch copy. It has no rates yet, so this lookup uses its explicit fractions
(`use_vfpexplicit = true`).

**Composition.** The lift margin looks up the tubing curve with one water and
gas fraction for the whole FLO range. The anchor takes them from its **first
flowing solve** and keeps them for all its rounds:
- Before that solve the well has no rates, so the explicit fractions (rates at
  the start of the timestep, or when it last flowed) are all there is: the
  start bhp uses them.
- After it, explicit fractions would only add a lag (C-1H in model5: water cut
  0.025 explicit, 0.195 actual).
- Each round's own rates would not do either: rounds approach the touching
  point, near shut-in for a well without a cliff, where the composition is that
  of its last producing layer (C-1H: water cut towards 1).

Below the table's first FLO value the lookup still uses the explicit fractions.
The first flowing solve is the anchor's lowest bhp, so its highest rate: the
most realistic composition it sees, whatever the case. It is stored with the
trial IPR (`SingleWellState::stopped_composition`), and the network gives the
candidate the same composition (section 6).

**Rounds** (at most 12):
1. `solveWellWithBhp(bhp)`. Track a bracket: `flows_at` (highest bhp seen to
   flow) and `stops_at` (lowest bhp seen to stop).
   - Not converged before anything flowed: status *SolveFailed*, stop.
   - Stopped, or not converged after a flow: treated as a stop; if nothing has
     flowed yet, status *NoFlow*, stop; otherwise bisect the bracket and
     continue (stop when the bracket is narrower than 0.1 bar).
2. The solve flowed: refresh the implicit IPR and compute the maximum flowing
   THP for that IPR (`maxFlowingThp()`, below). The THP, touching bhp and
   touching FLO become the anchor's current values; the first round's IPR is
   kept separately (`start_flows`, `start_bhp`, `start_ipr_a/b`).
3. Settled if the table says *capped* (flows even at the table's highest THP),
   or the touching bhp is within 0.1 bar of this round's bhp, or -- when it
   touches at zero rate -- the maximum flowing THP moved less than 0.1 bar.
   Otherwise the next bhp is the touching bhp (bisected inside the bracket if it
   is past a bhp already seen to stop).

**A final solve** at the touching bhp (or, failing that, at `flows_at`) gives
the anchor's IPR, rates and the THP the table gives for them.

**Status:** *Flows* (found), *Capped* (flows at the table's highest THP; the
THP is a lower bound), *NoFlow* (cannot flow anywhere in the table's THP range),
*SolveFailed* (numerical).

### The lift margin and maximum flowing THP (tables only, no solves)

For a THP, `liftMargin()` is the minimum over FLO of (bhp the tubing curve
requires − bhp the IPR has available), over FLO = 0, the table's FLO knots below
the FLO at the bhp limit, and that FLO itself; required bhps are corrected for
the datum depth and WVFPDP. A margin ≤ 0 means the well can flow at that THP.
`maxFlowingThp()` scans the THP knots from the top for the highest one where it
can flow, then bisects up to the next knot. The FLO where the minimum is taken
is the *touching FLO*: positive for a well with a **lift cliff** (curve with a
minimum the IPR touches from below), zero when the highest THP is reached at
zero rate (no cliff: e.g. the left branch of a U-shaped curve, typical for wet
gas wells).

## 5. Choosing the trial IPR

`margin` is `--group-tree-reopen-thp-margin` (bar, default 0), `limit` the
well's current THP limit (the node pressure in a network, otherwise its deck
THP limit), `lowest` the node's lowest possible pressure (`limit` without a
network).

| Case | Trial IPR |
|---|---|
| anchor *NoFlow* or *SolveFailed* | none |
| **can flow at the current pressure**: max flowing THP ≥ limit + margin | |
|   touches at FLO > 0 (lift cliff) | the anchor's IPR, at its touching bhp |
|   touches at zero rate, and the first solve (at the THP limit) flowed | the first solve's IPR, at its bhp |
|   touches at zero rate, first solve did not flow | falls through to the next part |
| **cannot flow at the current pressure** | |
|   max flowing THP < lowest + margin | none: it cannot flow at any pressure its node can have |
|   touches at FLO > 0 | the anchor's IPR at its maximum flowing THP (below the limit) |
|   touches at zero rate | `trialSolveAtThp()` at THPs stepping down from top = min(max flowing THP, limit) towards lowest: top − f·(top − lowest), f = 0.25, 0.5, 1; the first that flows. None flows: none |

Without a network `lowest` equals `limit`, so a well that cannot flow at its
fixed THP limit gets no trial IPR -- as before this change.

`trialSolveAtThp(thp)`: on scratch copies, a BHP-controlled solve at
`max(calculateMinimumBhpAtThp(thp), bhp limit)` (that lookup with explicit
fractions, as for the anchor's start); if it converges with the well flowing,
that bhp and the IPR there.

### Why these choices

- **Current pressure first.** In most cases the node pressure of the last
  network solve is a good guess of where the well would run.
- **Lower pressures only if needed.** The current node pressure may include
  flow that will not be there when the well reopens -- in particular the well's
  own flow before it stopped (e.g. a well put on a group target the network
  back-pressure would not let it deliver, which then stopped). Judging it
  against that pressure would keep it stopped for good; the network solve
  decides, with the candidate's own THP row, at the pressure the node actually
  gets.
- **The anchor's IPR at a cliff, not at zero rate.** With a cliff, a reopened
  well runs near its lift limit, where the anchor's IPR is linearised. Touching
  at zero rate, the touching point is the shut-in point: never an operating
  point, often far off the well's productivity where it would run, and a
  candidate there starts on its own cap in the network solve (T/C chatter).
- **`lowest` as the gate.** A node's pressure in a production network is never
  below the fixed pressure upstream of it; below that, offering the well only
  costs solves.

## 6. Into the network solve

**Gather** (`gatherWellNetworkDataForGroupTree()`): a stopped well
(`ws.status == OPEN && wellIsStopped()`) has a trial IPR if some phase's
`stopped_ipr_b` is positive (any phase, not just oil). It carries the trial
IPR (sign-converted to production positive), `stopped_ipr_bhp`, and phase
shares (rates per unit of FLO) from `stopped_composition`: the composition the
anchor judged its lift with (its first flowing solve, section 4). So the
anchor's gate and the network solve see the same fluid, whichever point the
trial IPR comes from. Fallbacks: without a stored composition (no tubing
table), the rates the trial IPR gives at `stopped_ipr_bhp`; below the table's
first FLO value, the rates behind its explicit fractions.

**Candidates** (`solveGroupTree()`, only when B offers them -- the first pass
of its settle loop): every stopped well under this network with a trial IPR,
not already in the solve, is added with `addReopenCandidate()`:
- a THP well at its own network node, with the trial IPR made proportional in
  FLO with its phase shares (`applyPhaseShare()`), and the slope limit of that
  IPR;
- counted in the sum of its nearest Active ancestor group (efficiency
  accumulated up to it), if it has one;
- starting at `min(stopped_ipr_bhp, shut-in bhp)` (at the shut-in bhp if
  `stopped_ipr_bhp` is unknown).

**Outcome** (`reopenOutcomes()`): after a converged solve a candidate reopens
if it flows -- some phase rate positive and its bhp more than 1e-3 bar below
its shut-in bhp -- *and* its crossing is on the real tubing curve, not only on
the flattened one. Otherwise it stays stopped ("no flow at this node pressure"
or "it would flow only on the flattened tubing curve").

**Reopening** (`updateGroupTreeOpenSet()`): the well is opened (dynamic
status), its stop hold released, and its well state seeded with the network's
operating point (rates, bhp, THP = node pressure), THP control, and the trial
IPR as its implicit IPR -- what the balancer and the next solve see until the
well's own solve refreshes them. The A3 check then compares the well's own
solve with that IPR; a reopened well that stops again counts as a mismatch of
1 (another round).

**Failure:** a solve with candidates that does not converge is solved again
without them; the candidates then stay stopped this round.

## 7. Caveats and open points

- **Explicit fractions** remain only where a stopped well has no rates (the
  start bhp of the anchor and of `trialSolveAtThp()`) and below a table's first
  FLO value. Whether they should be the start-of-step rates or those of the
  most recent converged flowing solve is open (`next_steps.md`, Step 5c).
- **Accuracy of the IPR.** A trial IPR is only as good as the well solve it is
  taken from; strongly nonlinear inflow (A11: gas rate per bar of drawdown
  changing fourfold between 74 and 110 bar) makes any single linearisation
  local, and a loose well-solve tolerance makes it worse.
- **The descending solves** (zero-rate anchors that cannot flow at the current
  pressure) were not exercised by the test decks so far.
- **Cost.** Every stopped predicting producer gets an anchor (up to 12 solves
  plus a final one) every round, and possibly up to three more solves.
- **Reopening outside the network.** The wells' own operability checks can
  still reopen a stopped well at zero rate with explicit fractions, at a node
  pressure computed without it; that path is part of the open/stop
  oscillation work, not of the trial IPR.
