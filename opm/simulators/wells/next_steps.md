# Next steps for the group-tree network branch

Status as of branch `group-tree-network` @ 369f31a41. Background and design
are in the other notes in this directory:
- `timestep_workflow.md` — the A/B workflow, changes C1–C7, questions Q1–Q7;
- `well_status_oscillation.md` — open/stop oscillation within a timestep;
- `timestep_initialization.md` — initialization, section 6 the initialization
  procedure and the lift margin;
- `master_issues.md` — problems found in master.

---

## 1. Done in this session (for reference)

- **C1:** the balancer owns production control within NUPCOL.
- **C2:** group-tree node pressures are handed to the wells undamped.
- **The A/B workflow:** rounds, IPR-mismatch criterion, stops held, reopen
  candidates in B's first pass.
- **Analytic Jacobian** for `GroupTreeSystem` (`--network-analytic-jacobian`).
- **Network nodes that are not groups** are handled by the simultaneous
  solvers.
- **Step 0:** the metrics script (`network_metrics.py`).
- **Step 1:** `--group-tree-stop-hold` (`iteration` | `timestep` | `nupcol`).
- **Step 2:** `VFPHelpers::liftMargin()` and `maxFlowingThp()`.
- **Step 3:** the refined maximum flowing THP per well (`computeAnchor()`,
  `--log-well-anchors`).
- **Step 4 (first part):** stopped wells' trial IPR from their maximum flowing
  THP; reopen candidates start at the trial IPR's bhp;
  `--group-tree-reopen-thp-margin`.
- **Multi-root networks:** fixed-pressure nodes below the root,
  `independentRoots()`.
- **Datum-depth correction** in the group-tree solve.

---

## 2. Agreed next steps, in order

### Step 3 — Lift margin and refined maximum flowing THP for a well (done)
`WellBhpThpCalculator::maxFlowingThp()` / `calculateMinimumBhpAtThp()` and
`WellInterface::computeAnchor()`: the refined maximum flowing THP from the best
available start (a solve at the well's own bhp, else at a given THP, else the
table's lowest THP), with a bhp bracket and a final solve
(`timestep_initialization.md`, 6.4). `--log-well-anchors` logs it per producer
at step start, with per-round detail; no effect on the run.

### Step 4 — Use it where the risk is lowest (first part done)
Done (e3c75ab58):
- Stopped wells with a tubing table get their trial IPR from
  `computeAnchor()`, started from a solve at the THP limit (the node pressure).
  No trial IPR, so no reopen candidate, if the well cannot flow or its maximum
  flowing THP is below the THP limit plus `--group-tree-reopen-thp-margin`
  (bar, default 0). Wells without a tubing table keep `estimateOperableBhp()`.
- The bhp the trial IPR was taken at (`SingleWellState::stopped_ipr_bhp`) is
  where a reopen candidate starts in `GroupTreeSystem`. Starting at shut-in, a
  candidate on a tubing curve with a minimum never left its cap, so no well was
  ever reopened by the network before this.

Results (default time-stepping, `--group-tree-stop-hold=timestep`, timesteps /
wasted Newton iterations):

| | before | margin 0 | margin 1 bar | margin 3 bar | margin 5 bar |
|---|---|---|---|---|---|
| FLOW-CGC | 26 / 11.7% | 73 / 46.5% | 35 / 23.3% | 26 / 11.7% | 26 / 11.7% |
| model5 STDW | 31 / 51.8% | 12 / 0% | 12 / 0% | 12 / 0% | 31 / 51.8% |
| model5 MSW | 11 / 0% | 12 / 0% | 12 / 0% | 11 / 0% | 11 / 0% |

With margin 0, FLOW-CGC's PROD2 is reopened about 2 bar inside its lift limit
(18 times) and then cycles in its local solves. Runs with the iteration hold
are unaffected: there, A1 reopens wells before B gets the chance.

Remaining:
1. **The reopen margin:** default kept at 0; 3 bar worked best on these decks
   but is model dependent. Experiment further (the large model; a relative
   margin; tying it to the stop decision).
2. **The cheap lift margin outside B** (no solves): as the stop/reopen margin in
   the local operability check and for fixed-THP wells. Not added inside B for
   flowing wells: the flattened curve already answers there
   (`timestep_initialization.md`, 6.5).

Test cases: FLOW-CGC PROD2; model5 STDW/MSW with the timestep hold; B-1H in
`5_NETWORK_MODEL5_MSW` with a 1-day cap and the iteration hold (494 status
events) -- though model5's anchors for flowing wells are capped at the table's
highest THP, so the table can't resolve their lift limit.

### Step 5 — Initialization (medium–large)
`timestep_initialization.md`, 6.2 and 6.2a. Group-tree mode only; other modes
unchanged.
1. **Initial solve in `beginTimeStep()`** for every well without a valid
   previous solution (fresh start, new wells, wells reopened by events, wells
   whose initial solve failed earlier), replacing today's event-well solve for
   them: solve at the strictest individual limit (deck THP ignored for network
   wells); when only the BHP limit binds, one more solve at the bhp for a THP
   guess (node pressure if known, else deck THP or a reservoir-based guess),
   guarded by the lift margin. A network well without a node pressure then gets
   a dynamic THP limit equal to the THP guess.
2. **`prepareTimeStep()` keeps all its well solves** (decided to limit
   changes): the wells from item 1 converge in one iteration there, since they
   start at their solution under the same constraints. Skip wells held stopped.
   **Drop the pre-step rebalance** in group-tree mode.
3. **First round of the first global iteration: skip A1** (all wells were just
   solved at this reservoir state) and go straight to the IPR refresh, B, B4,
   A3. At a fresh start this is the initial network solve; its node pressures
   replace the THP guesses.
4. **Failure handling** (6.6): physical vs numerical; a failed initial solve is
   held stopped for the whole timestep regardless of `--group-tree-stop-hold`
   (a separate flag), retried at the next step, and the well is shut after N
   consecutive numerical failures.
5. **Start-of-step balance on potentials** (`timestep_initialization.md`,
   6.2b): in `beginTimeStep()` after the potentials and guide rates, before
   the initial solves; committed without rates, and the balancer owns
   production from there. `--group-tree-initial-balance` (default true).
6. **The initial solve with the balancer's group target** (6.2b): a well the
   start-of-step balance put on group control is first solved under GRUP at
   its committed target (deck THP still left out); if that fails, the solve at
   the individual limits is the fallback.

**Item 5 done (6e2602f18):**
`WellInterface::estimateStrictestProductionLimitFromPotentials()`,
`BlackoilWellModel::balanceGroupTreeFromPotentials_()`; the balancer takes the
wells' phase fractions from the potentials (`balanceGroupTree(…, rateOverride)`).
Results against the same build without it (default time-stepping): identical
on all opm-tests network decks with both holds (their wells are individually
controlled, so the start-of-step categorization equals round 1's); FLOW-CGC
timestep hold identical, iteration hold 80 → 84 timesteps, wasted Newton
41.9% → 37.9%, FOPT −0.5% (PROD2's cliff sensitivity).
`--group-tree-initial-balance=false` reproduces the baseline exactly.

On FLOW-CGC's fresh start the balance puts all four wells on GRUP (about 2500
each, manifolds at ORAT 5000); the initial solve, which skips group control,
first converges them at their individual ORAT limits, and `prepareTimeStep()`
then switches them to GRUP -- the detour item 6 removes.

**Item 4 implemented (uncommitted):** `initialSolveForGroupTree()` returns
Flows / NoFlow / NotConverged; a well that doesn't flow is stopped at zero
rate and held for the timestep (`holdStoppedForTimestep()`, independent of
`--group-tree-stop-hold`), gets no trial IPR (never a reopen candidate) and is
skipped in `prepareTimeStep()`. NoFlow: closed or not at the end of the step
by the usual well-test update. NotConverged: kept out of the well-test update,
retried at the next timestep, and shut after
`--group-tree-max-initial-solve-failures` (default 3) consecutive timesteps
(if `--shut-unsolvable-wells`). Counts are kept globally, per successful
timestep. Results identical on all test decks (no failures there); exercised on
NETWORK-01 with PROD2's BHP limit at 600 bar (NoFlow: stopped, closed at step
end) and with `--max-inner-iter-wells=2` (NotConverged: all producers retried
and shut after 3 timesteps). Not covered: a well that stopped with no flow is
an ordinary well at the next step (not retried by the initial solve).

**Item 6 done (9067b3a39):** results identical on all opm-tests
network decks and on FLOW-CGC's totals (both holds); on FLOW-CGC's first step
the four ORAT→GRUP switches are gone (same number of inner iterations before
the first network solve).

**Items 1–3 done (369f31a41), behind `--group-tree-initialization`
(default true):** `WellInterface::initialSolveForGroupTree()` in
`beginTimeStep()` for event producers in group-tree mode
(`groupTreeModeActive_()`); `prepareTimeStep()` keeps their solution (no
target reset) and skips the pre-step rebalance; the first workflow round of a
timestep's first global iteration skips A1. Where a node pressure already
exists, the THP limit is honoured in the initial solve (it is the THP guess);
only the deck THP limit is ignored.

Results, on vs off (same build), default time-stepping: identical on most
opm-tests network decks; `5_NETWORK_MODEL5_MSW` 85 → 16 timesteps (iteration
hold), 12 → 11 (timestep hold); `NETWORK-01-REROUTE` 14 → 20 timesteps
(+1.7% FOPT). The REROUTE difference starts at report step 1, where its wells
open: both runs have the same global-Newton blow-up there (CNV 4.6e3 / 2.4e5 in
iteration 1) and diverge only afterwards -- a sensitivity of that deck, not a
systematic effect of the initialization.

Test cases: a fresh start of the large model, with and without restart. (The
report-step-0 give-ups in the NETWORK-01 decks are not a fresh-start problem:
no producer is open at report step 0 in those decks.)

Scope (`timestep_initialization.md` 6.3): network wells need 6.2 in full (the
deck THP limit is ignored when a network THP is available; all other limits
apply). Fixed-THP wells use the existing `estimateOperableBhp()` chain. Wells
without a THP limit need nothing new.

### Step 5b — Wells stopped at the end of a timestep (implemented, uncommitted)
Observed on FLOW-CGC (`--group-tree-stop-hold=timestep`,
`--solver-max-time-step-in-days=5`): PROD2 at its lift limit ends a step
stopped, is open again at the next one, and so on over many chopped steps.
Cause: a dynamic stop lives only on the (recreated) well object, and at the end
of the step `checkWellOperability()` returns early for a stopped well unless
`changed_to_stopped_this_step_` is set -- which the group-tree hold explicitly
clears -- so the well never reaches `updateWellTestStatePhysical()`. In master
the outcome also depends on how the well was stopped (fragile).

Agreed rule: a well stopped cleanly during a timestep and still stopped at the
end of the accepted (converged) step goes to `updateWellTestStatePhysical()`:
shut or stopped per its shut-in instructions until WTEST or a schedule event,
for the rest of the run without WTEST. Accepting the step accepts that the
well can't flow, so **no end-of-step re-check** (that would bring back the
question the reservoir-aware IPR addresses, one level up).
- Record a **stop reason** on the well wherever it is stopped dynamically
  (operability check, no flow in a converged local solve, reopen limit,
  network/balancer hold, numerical).
- Clean reasons close at step end; numerical keeps its handling (unsolvable
  wells as today; initial-solve failures retried, shut after N).
- Not affected: STOP from the schedule, wells open at zero rate under a zero
  target.
- Injectors: same rule; drop the injection-network exception in
  `WellInterfaceGeneric::updateWellTestState()` (kept stopped, not shut).
- Behind a switch first (affects all modes and regression results).

Implemented: `WellInterfaceGeneric::StopReason` set by `stopWell(reason)` at
every dynamic stop (restored with the status where a solve saves and restores
it; `solveWellWithZeroRate()` no longer touches it), `stoppedCleanly()`, and in
`BlackoilWellModel::updateWellTestState()` a cleanly stopped well (not under a
zero group target, not STOP from the schedule) goes straight to
`updateWellTestStatePhysical()`. Switch: `--close-stopped-wells` (default
false). The injection-network exception is bypassed for cleanly stopped
injectors by this path (not removed).

Results, on vs off:
- Off: identical to 3e9d0e05f on all decks.
- On, iteration hold: identical on all decks except FLOW-CGC with a 5-day cap
  (183 vs 184 timesteps) -- master's operability path already closes these
  wells.
- On, timestep hold: model5 STDW/MSW FOPT −8.7% (C-1H, B-3H, C-2H shut; no
  WTEST in model5). With the switch off the timestep hold had kept C-1H from
  being closed (the hold clears `changed_to_stopped_this_step_`) and it was
  retried until it flowed again at day 91; the iteration-hold runs close it
  with or without the switch.
- On, timestep hold, FLOW-CGC: 73 → 39 timesteps (FOPT −0.1%); with
  `--solver-max-time-step-in-days=5 --group-tree-reopen-thp-margin=3`:
  239 → 166 timesteps, wasted Newton 45.7% → 0%, chops 27 → 0, PROD2 status
  events 412 → 3, FOPT unchanged. PROD2 is closed once (reason: network) and
  reopened by WTEST at day 227.

### Step 6 — Reservoir-aware IPR (large, research)
The linearized near-well response from the assembled Jacobian
(`well_status_oscillation.md`, strategy D), including the measured response
from a first oscillation. Only after Steps 3–5, when measurements show how much
oscillation is left. It plugs into the same places as Step 4.

### Along the way, driven by measurements on the large model
- **Stop-hold default:** run the large model with each setting.
  - On FLOW-CGC `timestep` is clearly best.
  - On `5_NETWORK_MODEL5_STDW` (default time-stepping) it is worse than
    `iteration`.
- **Analytic Jacobian:** turn it on in large-model runs; if no worse, make it
  the default for group-tree.
- **Eliminating the Thp-well bhp unknowns** from the dense solve: only if
  profiling shows the dense solve matters.
- **Iterating the injection networks together with the A rounds:** when a case
  needs it.
- **Branch hygiene:** rebase on master regularly; once Steps 3–4 have settled,
  split into reviewable PRs (slope-limited VFP lookup; `GroupTreeSystem` with
  analytic Jacobian; C1; the A/B workflow).

---

## 3. Open question: the group-tree / network solution procedure in B

How the balancer and the network solve interact inside B
(`timestep_workflow.md`, Q7) is deliberately kept open. In particular: should
the balancer ever see a solution in which a well sits on the flattened
(slope-limited) part of its tubing curve?

Variants:
- **B-i:** the stop loop runs with the categorization fixed. The balancer runs
  only once the open set is settled, so it only ever reads a valid state. Not
  implemented; it needs the network solve to handle a well stopped in a fixed
  tree (pinned at zero).
- **B-ii (current, `settleGroupTree_()`):** rebalance after every change of the
  open set. Each stop decision is then made against a categorization that
  accounts for earlier stops; but each intermediate balance reads well data
  from a solution where a well may be on the flattened curve.
- **Flattened wells kept in the balancer with an honest limit.** Instead of
  removing them, give the balancer the rate the well can actually sustain: the
  minimum stable rate at the maximum flowing THP (Step 3's anchor), or zero if
  the node pressure is above the maximum flowing THP. The stop decision then
  follows from the balance rather than preceding it. Possible once Step 3
  exists.
- **Decide before solving.** With the maximum flowing THP known, a well whose
  node pressure is above it (minus a margin) can be treated as stopped before
  the network solve, so that the network solve rarely lands on a flattened
  curve at all. This reduces how often the question arises but doesn't answer
  it.

How to decide: keep the variant switchable, and compare on the cases with
wells at their lift limit (FLOW-CGC PROD2, model5 B-1H, and whatever the large
model shows), on global iterations, status changes and agreement with a
small-timestep reference.

---

## 4. Issues found, not yet addressed

- **Report step 0 on the opm-tests NETWORK-01 decks:** the network solve gives
  up with an empty balanced tree. Not a problem: no producer is open at report
  step 0 in those decks (earlier read as the fresh-start problem; corrected).
- **Group-tree fails earlier than fixed-point with `--enable-tuning=true`** on
  model5 (aborts at day 60/91 against 152). Worth a look once Steps 3–4 are in.
- **`NETWORK-01_STANDARD`:** group-tree needs more iterations than fixed-point
  (18 vs 14 timesteps), unlike the other NETWORK-01 decks.
- **WVFPDP:** the THP-dependent bhp adjustment (`getVfpBhpAdjustment()`) is not
  applied in `GroupTreeSystem`; only the datum-depth correction is.
- **The Newton network solver** still gives up on fixed-pressure nodes below the
  root (only the group-tree solve was extended).
- **GCONPROD actions other than RATE** are not taken while the balancer owns
  production (`timestep_workflow.md`, known limitations). To be handled
  outside the balancer/network work.
- **Fallbacks:** autochoke nodes, gas-lift optimisation and reservoir coupling
  fall back to the old path.
- **Debug log:** the fixed-point "Network pressure computation … Node
  pressures" dump is written directly to `OpmLog` and appears out of order
  with the deferred output; it also shows the fixed-point pressures rather
  than the group-tree result.
- **Chopped timesteps:** whether a chopped timestep restores the node
  pressures of the step start is unchecked (`timestep_initialization.md`,
  section 5).
- **Master:** the fixed-point multi-root bug (`master_issues.md`, item 1), left
  for now.

## 5. Open questions carried over

- **Q3 / Q6:** open/stop oscillation within a timestep; the reopen counter's
  semantics. `--group-tree-stop-hold` is the current handle.
- **Q4:** after NUPCOL categorization and network pressures freeze; the long
  term goal is to remove NUPCOL.
- **Q7:** see section 3.
