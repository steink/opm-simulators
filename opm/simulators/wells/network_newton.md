# Making the group-tree network Newton converge

Everything the group-tree network solve (`--network-solver=group-tree`) does
so that its Newton iteration converges, and why. State of the branch after
`7ba3e4bcc` (2026-10-15). Code: `NetworkGroupTreeSystem.hpp`
(`GroupTreeSystem`), `NetworkSolve.hpp` (`solve()`),
`BlackoilWellModelNetworkGeneric.cpp` (`solveGroupTree()`,
`updateGroupTreeOpenSet()`).

## 1. The system

For the configuration the balancer hands over (which wells are on group
control, which groups hold a limit), one Newton solve for:

| unknown | row |
|---|---|
| node pressure p_i (below the fixed-pressure root) | p_i - T_branch(p_parent, q_i) = 0, in bar |
| one lambda per active group node with group wells | (sum of its subtree on its mode - target) / target = 0 |
| bhp of each THP well | bhp - (T_well(p_node, q(bhp)) - datum correction) = 0, in bar; or bhp - shut-in = 0 on its cap |

Group wells have no unknown of their own: their rate in the node's mode is
guide rate x lambda. Wells at a fixed rate are pinned. Converged when the
max-norm of the scaled residual is below 1e-2 (0.01 bar, 1 %) and no cap
changed in the last iteration; at most 50 iterations.

The solve runs in **two passes**: first on flattened branch tables, then,
where needed, on the real ones (sections 2 and 3). Section 4 lists what both
passes share; section 5 what happens when the second pass fails; section 6
what the workflow around the solve does.

## 2. Pass 1: the flattened network

**Branch tables flattened** (`--group-tree-branch-flattening`, default on).
Every branch table is looked up with slope limit 0: where the required
pressure falls with flow (the low-flow, liquid-loading part), the table is
replaced by a flat continuation. Together with the slope-limited well tubing
(4.2), every row is monotone in its own unknown and the network has a unique
solution. This is the robust pass: since the cap test (4.3) it almost never
fails (0-11 failures per Norne run against about 1500 solves).

## 3. Pass 2: the real branch tables

Only where the flattened solution has some branch on its flattened part
(`flattenedNodes()`): those branches are put back on their real tables and
the network solved again, starting from the flattened solution (node
pressures and the THP wells' bhps, `setStartFromResult()`). Usually 2-3
iterations. Most solves need it at some point on Norne (real operating points
often sit on the falling part of a flowline curve).

`--group-tree-real-table-pass=false` skips it and accepts the flattened
solution. Not recommended: on model5 with iteration hold it breaks down
(STDW 59 -> 130 timesteps, MSW 11 -> 128), and on Norne it produces about
4 % more oil from flattened-only states.

## 4. Used in both passes

### 4.1 One FLO-IPR per well, fixed phase proportions
Each well's per-phase IPR is replaced by one IPR in its tubing table's FLO,
with its phases in fixed proportions (`applyPhaseShare()`). Its rates, its
contribution to its node and its tubing lookups then have one composition at
every bhp, all phases reach zero at one shut-in bhp, and the FLO-based slope
limit (4.2) keeps the row monotone. The proportions: the well's own rates
from its IPR at its current bhp; for a stopped well with a trial IPR, the
composition its anchor judged its lift with (`trial_ipr.md`). Without it the
lookups and the rates disagreed (0/0 fractions, inconsistent shares).

### 4.2 Well tubing curves slope-limited against their own IPR
A THP well's tubing curve is looked up with a slope limit of 0.95 x its
FLO-IPR slope (`iprSlopeLimit()`): where it falls faster than that, it is
replaced by the line of the first interval flat enough. The row then rises
with bhp (d/dbhp >= 0.05), so it has at most one root: no cliffs, no
multiple crossings inside the solve. Whether the converged point is on the
real curve is checked afterwards (6.3). No limit for an IPR whose FLO does
not fall with bhp (`floIprValid()`, `trial_ipr.md`).

### 4.3 The shut-in cap, decided from the node pressure
A THP well's bhp is never stepped past its shut-in bhp (`limitStep()`). It is
on its cap (row bhp = shut-in, rate zero) when its slope-limited tubing needs
at least the shut-in bhp at (next to) zero rate at the current node pressure
-- decided at every iterate (`--group-tree-cap-test=node-pressure`, default).
On the slope-limited curve the row rises with bhp, so the well has a root
below shut-in exactly when this test fails: the cap follows the node
pressure, not the bhp iterate. The old rule (`bhp`: only once the bhp
iterate reaches shut-in) made Newton push a well that cannot lift across
tubing-table corners to reach its cap (Norne B-3H), where full steps cycle
and a line search stalls. Cap switches are rare: at most 0.2 per solve on
average on Norne.

### 4.4 Node fractions below a table's first FLO value
A branch table below its first FLO value is evaluated with fixed water and
gas fractions: those of the node's flow at the previous converged solve,
else its wells' reference rates (`setNodeExplicitFractions()`,
`finalize()`). Master's node calculation uses WFR = GFR = 0 there (dry,
gas-free oil; `master_issues.md`).

### 4.5 Analytic Jacobian
From the tables' one-sided derivatives at the point itself
(`--network-analytic-jacobian`, default true since `120966a32`). The
finite-difference Jacobian steps a bhp by 0.01 bar; for a steep-IPR well that
moves its rate across a table knot (B-3H: 100 m3/d), and the derivative comes
out as an average of two segments (9.4 against 19.6), so the Newton direction
is no descent direction. The finite-difference fallback still clamps its
exploratory step at the shut-in bound.

### 4.6 Full Newton steps
`--group-tree-network-step=full` (default). A step on which some cap changed
is always taken in full. The residual has kinks (table knots, the
slope-limited curves' bridged stretches), and a line search
(`linesearch`: halve the step until the 2-norm of the scaled residual drops,
at most 12 times) can creep up to a kink and stall when the solution lies
across it or behind a cap switch (B-3H before the cap test 4.3).

- Before the cap test (4.3) and the node test (5) the line search was much
  worse than full steps (Norne real-table failures 4 -> 51, 9 -> 66).
- With both: identical to full steps in 8 of 12 runs (Norne, model5,
  FLOW-CGC). Where it differs it fails pass 2 more often (Norne 5 -> 14,
  13 -> 21; model5 STDW 3 -> 9), but the node test catches every one of
  those (cliffs, or retries from the crossing that converge), so no solve
  reaches the fallback; slightly more chops in two Norne runs (5 -> 11,
  14 -> 15). No gain anywhere.
- Possible later: a line search that kicks in only when full steps run into
  trouble -- after a number of iterations, or when the residual has risen
  several times in a row (the cycles of E-1H and E-2H) -- so that it breaks
  cycles without stalling the solves full steps converge.

### 4.7 Scaling and start
Pressure and bhp rows and columns in bar, target rows relative to the
target. Start: node pressures from the previous network solution (terminal
pressure where there is none); THP wells at their current bhp (capped at
shut-in), reopen candidates at their trial bhp (at most shut-in); each
lambda at target / sum of its group wells' guide rates.

### 4.8 No phantom unknowns
An active group node whose own wells are all individually constrained has
nothing for its lambda to move: its lambda and row are dropped
(`finalize()`), else the Jacobian is singular.

## 5. When pass 2 fails: the node test

`nodeMargin()` for each node whose branch ended on its flattened part:

- q(p): the flow the wells under the node deliver at node pressure p -- every
  THP well on its own row at p (`thpWellAtNodePressure()`, the unique root or
  its cap), the other wells and the child nodes at their solved rates.
- g(p) = p - T_branch(upstream, q(p)) on the real table, sampled from the
  upstream pressure to where the last THP well under the node is capped.
- **margin < 0** (g negative wherever a THP well flows): no flowing state
  exists -- a cliff at network level. The wells cannot keep the branch
  flowing; the only real state has them shut in. The THP wells under the node
  are stopped (stop reason Network, not reopened that round) and the tree
  rebalanced. Norne: E-TPL on E-2H alone, g peaks at -2.1 bar, the only root
  is E-2H shut in at 61.5 bar.
- **margin >= 0**: a flowing state exists and Newton missed it; the
  real-table solve is retried from the crossing (node pressure and the wells'
  bhps there).
- Still failing: `NETWORK_FLATTENED_BRANCH` warning and the relaxed
  fixed-point update for that network (or the flattened solution, with
  `--group-tree-flattened-branch-fallback=flattened`). None left on the Norne
  runs.

Limits: the upstream pressure is held at its flattened value (exact for
Norne's templates under a fixed node), and child nodes at their flows.

## 6. Around the solve

1. **Well cliffs:** a THP well whose converged point is on the flattened part
   of its own tubing curve (`worstCliffViolation()`) cannot operate there;
   the worst is stopped and the tree rebalanced, one well per pass.
   **B's settle loop** (`settleGroupTree_()`) repeats a pass -- balance, then
   for each network root the network solve and `updateGroupTreeOpenSet()` --
   until the open set stops changing, at most 10 passes. In a pass a root's
   open set changes by: reopening every candidate that flows (first pass
   only, where candidates are offered), stopping all THP wells under a node
   that failed the node test (5), or else stopping the one worst well cliff.
   After the first pass the open set only shrinks. Still changing after 10
   passes: warning `GroupTreeCliffDiagnosisNotStable`, the last
   configuration is kept.
2. **Reopen candidates:** stopped wells with a usable trial IPR enter as THP
   wells starting at their trial bhp (`trial_ipr.md`). A solve with
   candidates that fails is solved again without them.
3. **Rebalancing on the active mode:** a group well's guide rate is looked up
   in its node's active mode (`b8593f0e9`); a well on the balancer's fallback
   is pinned at its balanced rate. With the preferred mode's guide rates,
   Norne's oil-rich wells got 4.7 times their gas share and pushed their
   template past what they could lift.
4. **Rounds:** after the wells' own solves (A3), the IPR mismatch decides
   whether to rebalance and solve again with fresh IPRs, at most
   `--network-max-outer-iterations` rounds per Newton iteration.

## 7. Tried and dropped

| what | result |
|---|---|
| Warm start (each THP well solved alone at its node pressure, node pressures recomputed from the flows, 2 passes) | Helped without branch flattening (model5 STDW iteration hold 54 -> 16 timesteps), nothing with it (unique solution). Kept on branch `gt-warm-start`; `next_steps.md`. |
| Sweep of node pressure guesses, cap-switch dwell, consistent projection at the cap, slope margins 0.8/0.5 | No effect, worse, or stalls (`next_steps.md`). |
| Line search from the first iteration (`--group-tree-network-step=linesearch`) | Stalled at kinks and before cap switches; with the cap and node tests no better than full steps (4.6). Kept as an option. |
| Keeping the flattened solution when pass 2 fails | Removes give-ups, not chops; superseded by the node test (5). |

## 8. Options

| option | default |
|---|---|
| `--group-tree-branch-flattening` | `true` |
| `--group-tree-real-table-pass` | `true` |
| `--group-tree-cap-test` | `node-pressure` |
| `--network-analytic-jacobian` | `true` |
| `--group-tree-network-step` | `full` |
| `--group-tree-flattened-branch-fallback` | `relaxed` |

Details and measurements: `options.md`.
