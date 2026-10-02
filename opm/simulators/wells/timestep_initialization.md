# How wells, groups and networks are initialized at the start of a timestep

Scope: what the well model starts each timestep from, up to the first global
Newton iteration's well/network update. Two cases:
- a **fresh start**: the first timestep of a run with no restart, where nothing
  from a previous solve exists;
- a **continuation**: a previous timestep's solution is available, either
  within a report step, at a new report step, or from a restart file.

Functions are the stable reference; line numbers are as of branch
`group-tree-network` @ ece86a9ac.

---

## 1. Where the state lives

| State | Holder | Carried between timesteps? |
|---|---|---|
| Well rates, bhp, thp, control mode, group target, potentials, IPRs, persistent status (OPEN/STOP/SHUT) | `WellState` (`SingleWellState` per well) | Yes: committed at the end of each successful step (`commitWGState()`), restored by `resetWGState()` at the next `beginTimeStep()`. |
| Group rates, control modes, targets | `GroupState` | Yes, together with the well state. |
| Network node pressures | `BlackoilWellModelNetworkGeneric::node_pressures_` / `domain_node_pressures_` | Yes: a member of the network object, not part of the WG state. Only the secant/damping history is cleared each step (`network_.beginTimeStep()`). |
| Well objects: dynamic status (`wellIsStopped()`), dynamic THP limit, operability flags, reopen counter | `WellInterface` in `well_container_` | **No**: recreated every timestep by `createWellContainer()`. |
| Guide rates | `GuideRate` | Updated at every `beginTimeStep()` from the current potentials. |

The last row matters: anything the network or the balancer set on a well object
(a dynamic stop, a dynamic THP limit, a hold) is gone at the next timestep
unless something sets it again. The network sets the dynamic THP limit again
(section 3); nothing sets a dynamic stop again.

---

## 2. Sequence

### 2.1 Once per report step: `beginReportStep()` ([BlackoilWellModel_impl.hpp](BlackoilWellModel_impl.hpp))

1. `initializeLocalWellStructure()` → `initializeWellState()` →
   `WellState::init()` ([WellState.cpp:273](WellState.cpp#L273)). It builds a new
   well state from the schedule and the current cell pressures, then copies
   what it can from the previous well state (section 4).
2. `initializeGroupStructure()`, VFP properties, `commitWGState()`.

### 2.2 Every timestep: `beginTimeStep()` ([BlackoilWellModel_impl.hpp:335](BlackoilWellModel_impl.hpp#L335))

In order:

1. `network_.beginTimeStep()`: clears the per-node secant/damping history.
   Node pressures are kept.
2. `resetWGState()`: well and group state back to the last committed state.
3. `wellTesting()` (WTEST). Each tested well is a scratch well object, which
   gets the network node pressure as its THP limit (`network_.initializeWell()`)
   before it is tested.
4. `createWellContainer()` ([:824](BlackoilWellModel_impl.hpp#L824)): new well
   objects. Wells closed by the well-test state become SHUT (not added) or STOP
   (added, `stopWell()`) according to their auto-shut-in setting. It ends with
   **`network_.initialize()`** ([:1069](BlackoilWellModel_impl.hpp#L1069)), which
   gives every well in a network its node's current pressure as the dynamic THP
   limit (section 3).
5. `updateAndCommunicateGroupData(false)`, `initWellContainer()`, efficiency
   factors, `closeCompletions()`.
6. Only with `alternative_well_rate_init_`:
   `initializeProducerWellState()` gives producers multi-phase rates
   proportional to their rates at zero bhp. Its purpose is to avoid
   single-phase or zero rates, and with them 0/0 phase fractions in the VFP
   lookup.
7. `setPrevSurfaceRates()` for wells with a VFP table.
8. **`updateWellPotentials(onlyAfterEvent = true)`**
   ([BlackoilWellModelGeneric.cpp:1669](BlackoilWellModelGeneric.cpp#L1669)):
   potentials are recomputed only for wells with an event this step (new
   well, status, completion, PI, efficiency or control change, ACTIONX). Every
   other well keeps the potentials computed at the end of the previous step
   (section 2.4).
9. `updateGuideRates()`, from those potentials.
10. GPMAINT targets. Then `setBalancerOwnsProduction(false)` and
    **`updateAndCommunicateGroupData(true)`**: the standard group-target
    allocation into every well's `ws.group_target`.
11. For wells with a status/type event, or whose status differs from the
    previous well state: `updateWellStateWithTarget()` and a full
    `solveWellEquation()`.

### 2.3 First global iteration: `assemble()` → `prepareTimeStep()` ([:2615](BlackoilWellModel_impl.hpp#L2615))

1. Wells with control events: `updateWellStateWithTarget()`.
2. If `solve_welleq_initially_`: `solveWellEquation()` for every operable well,
   at the dynamic THP limit set in 2.2 step 4.
3. `resetWellOperability()`.
4. **Pre-step network rebalance** if due (`shouldDoPreStepNetworkRebalance_()`;
   the network-side trigger `needPreStepRebalance()` is a WELL_STATUS_CHANGE
   event on a network well). It runs the whole well/network update of section
   1.3 in `timestep_workflow.md` with `mandatory_network_balance = true`, which
   with the group-tree solver is the A/B workflow, and then one
   `solveEqAndUpdateWellState()` per well.

Then the ordinary per-iteration update (`updateWellControlsAndNetwork()`).

### 2.4 End of every successful timestep: `timeStepSucceeded()` ([:661](BlackoilWellModel_impl.hpp#L661))

`updateWellPotentials(onlyAfterEvent = false)` for all predicting wells, well
test state (economic/physical shut-ins), then `commitWGState()`. **The
potentials the next timestep starts from are computed here**, at the end of
the step, with the well objects' THP limits as they were then (section 3).

---

## 3. Network node pressures and the THP the wells see

`network_.initialize()` → `initializeWell()`
([BlackoilWellModelNetworkGeneric.cpp:1835](BlackoilWellModelNetworkGeneric.cpp#L1835)):
for each well whose group is a node of an active network, **if the network has
node pressures at all**, the well gets `setDynamicThpLimit(node pressure)`.
Injectors have it clamped to their VFP table's THP axis.
`WellInterfaceGeneric::getTHPConstraint()` then returns this dynamic limit in
place of the deck's WCONPROD/WCONINJE THP limit. Everything that uses the
THP-limit reads it through `getTHPConstraint()`: the local well solve, the
operability check, and **the potential calculation**.

What the node pressures are at that point:

| Case | `node_pressures_` at `beginTimeStep()` | THP the wells (and potentials) use |
|---|---|---|
| Fresh start, no restart | **empty** (never solved) | The deck THP limit. `initializeWell()` does nothing, so the network's back-pressure is not seen at all. |
| Restart | From the restart file (`setFromRestart()`, called in the network's constructor) | Restart node pressures. |
| Later timestep | The last pressures the previous timestep's network update produced | Those pressures. Within NUPCOL they are the balanced ones; after NUPCOL (NETBALAN NUPCOL mode) or after the first iteration (TimeStepStart mode) they are whatever was last computed, i.e. frozen, not necessarily balanced against the final reservoir state. |

The first network update of a fresh start is also different: `updatePressures()`
has no previous pressures, so it takes the computed ones as they are, with no
damping or secant step. The group-tree solve's own starting guess falls back
to the terminal pressure for any node without a previous value
(`solveGroupTree()`).

---

## 4. What each piece starts from

### 4.1 Well state

**Fresh well state** (`WellState::base_init()` → `initSingleProducer()` /
`initSingleInjector()`, then `SingleWellState::update_producer_targets()`,
[SingleWellState.cpp:254](SingleWellState.cpp#L254)):
- *status:* OPEN/STOP/SHUT from the schedule;
- *control mode:* the deck's (GRUP becomes BHP if the well isn't available for
  group control);
- *rates:* the target phase set to its target for ORAT/WRAT/GRAT control;
  **zero** for GRUP, THP and BHP control. The zeros are filled in later by
  `initializeProducerWellState()` if enabled, otherwise by the first well
  solve.
- *bhp:* the BHP limit if on BHP control, otherwise 0.99 × the pressure of the
  first connection (producers). A stopped well gets the first-connection
  pressure.
- *thp:* the THP limit if on THP control;
- *group target, IPRs, potentials:* unset or zero.

**Continuation, same wells** (`WellState::init()` with the previous state,
[WellState.cpp:382](WellState.cpp#L382)), for each well that existed before and
wasn't shut:
- bhp, thp, temperature (`init_timestep()`), status;
- control mode, unless the schedule sets new controls this step (WCONPROD
  etc.);
- surface and reservoir rates, perforation rates (if the connections are
  unchanged), well potentials, `group_target`, productivity index.

Within a report step the whole committed well state is restored
(`resetWGState()`), which also brings back the IPRs, the trial IPRs, the ALQ
and so on.

**Restart:** the well state is initialized from the restart file
(`initFromRestartFile()`): rates, bhp, thp and control modes, as the restart
format provides them.

### 4.2 Group state

- *Control modes* come from the committed group state within a run. At a
  fresh start they're the deck's GCONPROD/GCONINJE modes.
- *Group rates* are recomputed from the well rates in every
  `updateAndCommunicateGroupData()`.
- *Group targets per well* (`ws.group_target`) are recomputed by the standard
  allocation at the end of `beginTimeStep()` (2.2 step 10) from those rates
  and the guide rates.

### 4.3 Well potentials and guide rates

| Case | Potentials at `beginTimeStep()` |
|---|---|
| Fresh start | Computed for every new well (NEW_WELL / status events), **with the deck THP limit**, since no node pressures exist yet (section 3). |
| Continuation, well without event | Those of the end of the previous timestep, computed with the node pressures of that step. |
| Continuation, well with event | Recomputed, with the previous step's node pressures. |

Guide rates follow the potentials, so at a fresh start the guide rates are
network-blind too.

---

## 5. Consequences and open points for a start-of-step balancer

On `group-tree-balancing`, the balancer also runs in `beginTimeStep()`, on
limits estimated from freshly evaluated potentials (not IPRs, which don't
exist yet at that point). Against the above:

- **Continuation:** the potentials are consistent with the network as last
  solved, which is a sensible starting point. The one caveat is after NUPCOL,
  when the last node pressures may be frozen rather than balanced.
- **Fresh start:** the potentials, the guide rates, and therefore a
  potentials-based balance all assume the deck THP limit. With a network whose
  back-pressure is well above that limit, or a THP limit defaulted in the deck,
  every network well's potential is overestimated. The balancer then starts
  from a categorization (too many wells on group control, too few on THP) that
  the first network solve will overturn. Section 6 proposes a procedure that
  avoids depending on the deck THP limit at all.
- **Dynamic state is lost at every timestep:** a well the network stopped in
  the previous step starts the next one open, unless the end-of-step checks
  closed it through the well test state. That's consistent with the
  semantics of "stopped at the end of a step" (see `well_status_oscillation.md`),
  but a start-of-step balancer will see it as an open candidate again.
- **Worth checking:** whether a chopped timestep restores the node pressures of
  the step start. The WG state is restored by `resetWGState()`, but the node
  pressures live in the network object; there are `last_valid_*` copies there
  whose use I haven't traced.

---

## 6. Proposed initialization procedure

### 6.1 Principle

Make each well solvable on its own first, then do everything network-related
on IPRs. Only one assumption is needed, and it's the user's responsibility to
meet it through reasonable deck limits:

> With the reservoir as at the start of the step, a well's local solve
> converges at its strictest individual limit (a rate limit, or the BHP limit
> when no rate limit binds).

No THP-controlled solve is needed to initialize a well. THP solves are where
the trouble lives: multiple IPR/tubing-curve crossings, the lift cliff, no
solution at all. And for a well in a network the deck THP limit is ignored
anyway, because the node pressure replaces it.

What matters for the initial solve is that it **converges within the well's
limits** and that its THP is **not too far from the THP the network will
eventually impose**, which is typically not known yet. The initial solve only
has to give a decent IPR in roughly the right region; the network solve and the
A loop do the rest.

The procedure applies to **every well without a valid previous solution**: at
a fresh start, new wells, wells reopened by events, and wells whose initial
solve failed in an earlier step. Every other well starts from its previous
solution and IPR, which are closer to the answer. The network solve only needs
an IPR per well, wherever it came from, so the two kinds mix freely.

### 6.2 The initial solve

1. **Solve at the strictest individual limit.** When a rate limit binds, that
   is close to where the well will run.
2. **When only the BHP limit binds, one more solve at a THP guess.** The deck
   BHP limit can be set very low, far from where the well will operate, and an
   IPR taken there extrapolates badly. Take the bhp the first solve's IPR gives
   for a THP guess, and solve again there. THP guess, best first: the node
   pressure (previous step or restart); otherwise the deck THP limit or a guess
   based on reservoir pressure.
3. **Guard with the lift margin** (section 6.4, no solves): if the THP guess is
   above the maximum flowing THP of the current IPR, lower it to just below,
   or report that the well cannot flow there.
4. **Let the network decide the THP.** The group-tree solve on these IPRs gives
   the node pressures (at a fresh start from a trivial guess: terminal pressure
   everywhere, or a zero-flow pass down the tree); A1/A3 re-solve the wells
   there, and the IPR-mismatch criterion refines. The circularity at a fresh
   start (potentials need a THP, the THP needs node pressures, node pressures
   need rates) is broken, and the balancer can use the same IPRs instead of
   deck-THP potentials.

The maximum flowing THP is deliberately **not** the target of the initial
solve: it is the most marginal state a well has (at the lift limit, at the
lowest stable rate, or at zero rate when the touching point is at shut-in), the
worst place to start a well and to take the IPR the network will extrapolate.

### 6.2a Where it happens in the timestep (group-tree mode)

Today, in the first global iteration of a step, the same wells can be solved
four or five times at the same reservoir state: the event-well solve in
`beginTimeStep()`, `prepareTimeStep()`'s solve of every operable well, the
pre-step rebalance (a full A/B workflow), and the regular update's A1 and A3.
Agreed order, with the group-tree workflow active (other modes unchanged):

0. **`beginTimeStep()`, after the potentials and guide rates:** balance the
   group tree on the wells' potentials and commit its categorization and
   targets (6.2b). The balancer owns production from here, so the standard
   group-target update no longer overwrites them.
1. **`beginTimeStep()`:** the initial solve (6.2, steps 1–3) for wells without
   a valid previous solution, replacing today's event-well solve for them. A
   network well whose node has no pressure yet then gets a **dynamic THP limit
   equal to its THP guess**: the deck THP limit is ignored from then on by the
   existing rule, and step 2 below solves the same problem the initial solve
   did.
2. **`prepareTimeStep()`:** unchanged (control-event target updates, operability
   reset, `solveWellEquation()` for every operable well). Wells already solved
   in step 1 start from their converged state under the same constraints, so
   for them this is a single iteration (formally THP-controlled, but started at
   its own solution). Wells held stopped (6.6) are skipped. **No pre-step
   rebalance:** it would run the same workflow at the same reservoir state as
   step 3.
3. **First round of the regular update in the first global iteration:** skip
   A1 (every well was just solved at this reservoir state); refresh the IPRs,
   then B, B4 (node pressures become THP limits, replacing the guesses) and A3.
   This round is the pre-step rebalance, done once, and at a fresh start the
   initial network solve; a network well's first solve at a real node pressure
   is A3.
4. **Later global iterations:** unchanged (A1 as now, the reservoir has moved).

Wells whose initial solve failed are held stopped for the whole timestep,
independent of `--group-tree-stop-hold` (which may release holds every global
iteration): a separate flag for "initial solve failed".

### 6.2b Start-of-step balance on potentials

The group targets the standard update computes at step start are not trusted
(one reason the balancer exists), so the balancer sets them from the start of
the timestep. Its limits come from the potentials, not from IPRs: no well has
been solved at this reservoir state yet, and the potentials are already
available.

- **The common case, a few new wells among existing ones:** potentials are
  recomputed for every well at the end of each successful step, at the THP
  limits the wells had (the node pressures), and for event wells at
  `beginTimeStep()`, at their node's pressure if it has one. The group totals
  are dominated by the existing wells, so the targets a new well gets are close
  to final, without a single extra well solve.
- **A fresh start** (rarer, allowed to be less accurate): potentials are at the
  deck THP, or at the BHP limit when a network well has no deck THP limit (not
  uncommon), and overestimate. The balance still gives each well a
  categorization and a target; wells that can't reach their target settle at
  their pressure limit, and the first workflow round corrects with IPRs and real
  node pressures.

Each well's limit (`estimateStrictestProductionLimitFromPotentials()`): the
strictest rate limit if the potentials exceed it, else the pressure limit with
the total potential as its rate (THP if the well has a THP limit, else BHP). The
balancer takes the wells' phase fractions from the potentials too, since a new
well has no meaningful rates yet. Committed without rates. Switch:
`--group-tree-initial-balance` (default true).

**The initial solve with the group target** (step 1 above, done): a well the
balance put on group control is solved under GRUP at its committed target
first, so it starts at about its final operating point; the solve at the
individual limits is the fallback.

Possible follow-ups, only if measurements call for them:
- **Fresh start only:** recompute potentials and guide rates once real node
  pressures exist (after the first B4), removing the deck-THP/BHP-limit bias
  from the guide rates for the rest of the first step.
- **A full B (balance and network solve) before `prepareTimeStep()`:** more
  accurate at a fresh start, but a network solve per step that largely
  duplicates the first round.

### 6.3 Which wells need what

| Well | THP | Initialization |
|---|---|---|
| In a network | Unknown: the node pressure | 6.2 in full |
| Fixed deck THP limit, no network | Known: the deck value | The existing chain below |
| No THP limit (no tubing table) | None | Standard solve at the strictest limit; nothing new |

**Network wells and the deck THP limit (settled).** When a network THP is
available, the deck THP limit is ignored; all other deck limits (rates, BHP)
still apply. This is what `WellInterfaceGeneric::getTHPConstraint()` already
does.

**Fixed-THP wells.** The THP is known, so the operating point at it is found by
`estimateOperableBhp()` → `calculateMinimumBhpFromThp()` → `solveWellWithBhp()`
→ `estimateStableBhp()`, which already avoids a direct THP-controlled solve.
They need no new code, only the failure handling of 6.6. The lift margin at the
deck THP is their cheap stop/reopen test.

### 6.4 Lift margin and maximum flowing THP

For a given IPR, a well can flow at a THP if its IPR crosses the tubing curve
there. As the THP rises the crossing moves, and at the **maximum flowing THP**
it disappears:
- **Tubing curve with a minimum** (the lift-cliff case): the IPR line just
  touches the curve.
- **Tubing curve rising with rate** (friction-dominated): the crossing goes to
  zero rate, at the IPR's shut-in bhp.

**Lift margin, no solves** (`VFPHelpers::liftMargin()`): with the IPR, the
phase fractions and the ALQ fixed,

`F(thp) = min over FLO of [VFP(thp, FLO) − bhp_IPR(FLO)]`

over the FLO reachable above the BHP limit. The table is piecewise linear in
FLO and the IPR linear, so the minimum is at a FLO knot, zero rate or the BHP
limit: exact, table lookups only. F ≤ 0: the well can flow at that THP.
`VFPHelpers::maxFlowingThp()` finds the largest such THP in the table's range,
scanning the THP knots from the top (tables need not be monotone in THP) and
bisecting in the bracket; it reports whether the well can flow at all and
whether the result is capped at the table's highest THP (a lower bound).
`WellBhpThpCalculator::maxFlowingThp()` applies it to a well's IPR with the
datum and WVFPDP adjustments. The table is extrapolated below its first FLO
point (accepted).

**Refined, with solves** (`WellInterface::computeAnchor()`): the IPR is only
good near where it was taken, so the cheap estimate is refined: solve at the
touching bhp, take the IPR there, repeat until the touching point is where the
IPR was taken (or, at a zero-rate touching point, until the THP settles). A
bracket of the highest bhp seen to flow and the lowest seen not to (or not to
converge) keeps the steps safe. A final solve gives the well's actual rates,
IPR and THP there. Start, best first: a solve at the well's own bhp if it is
flowing (a solve is needed even there: fresh well objects at step start have no
assembled equations, and a multisegment well's IPR cannot be linearized
without them); a solve on the tubing curve at a given THP (e.g. the node
pressure); the table's lowest THP as last resort.

On FLOW-CGC (`--log-well-anchors`): from the well's own bhp, typically 2
rounds; final solve succeeded in 102 of 104 cases; search and final THP agree
within 0.1 bar in 89.

**Uses:**
- *Stopped wells (reopen):* no current operating point, so the flattened curve
  has nothing to work with. The refined anchor, started from the node pressure,
  replaces the fragile trial IPR (`updateStoppedWellTrialIpr()`).
- *Decisions outside B:* the local operability check, fixed-THP wells, and the
  stop/reopen margin (`well_status_oscillation.md`, strategy F) -- with the
  cheap zero-solve margin on the current IPR.
- *Guard in initialization* (6.2, step 3).
- *Diagnostics* of wells near their lift limit.

### 6.5 Relation to the flattened tubing curves

For a network well inside B, the slope-limited curve and the maximum flowing
THP answer the same question: the well lands on the flattened part when its
node pressure is above its maximum flowing THP (approximately; the slope
limit's 0.95 margin makes the flattening slightly conservative). The slope
limit gives B a unique root and a stop signal without computing anything
extra, so the maximum flowing THP is not added inside B for flowing wells. For
a tubing curve rising with rate there is no cliff for either mechanism: the
network solve drives the rate to zero and the shut-in cap handles it.

### 6.6 Wells that fail

(Implemented in group-tree mode: see `next_steps.md`, Step 5 item 4.)

Two different failures, treated separately.

- **No flow at its initial solve** (the solve converges to zero, or the well is
  inoperable): a physical result. The well is stopped.
- **The initial solve doesn't converge**: numerical, the user's responsibility.
  The well is stopped with a warning naming it and the limit it failed at.

**Within the timestep**, either way, the well is stopped at zero rate and held
for the whole step:
- The balancer leaves it out, so any group-control share goes to the other
  wells. The network gather marks it `has_ipr = false`, and `GroupTreeSystem`
  pins it at zero: it still counts, as zero, in every node and target sum. It
  is never a reopen candidate, since it has no IPR.
- In the reservoir it is assembled like any stopped well (zero-rate solve, no
  source term).
- No retries within the step: a stopped well barely changes the reservoir state
  that made it fail, each retry is a local solve per global iteration, and a
  retry that sometimes succeeds brings back the open/stop oscillation. Held
  from the start, the open set within the step only shrinks.

**At the end of the timestep:**
- *Physical:* stopped or shut per its shut-in instructions, as today; a later
  well test can reopen it.
- *Numerical:* not closed in the well-test state (the well might flow fine, and
  a well test may never be scheduled). At the next timestep it is still a well
  without a valid solution, so the initial solve is attempted again against a
  slightly different reservoir.
- *Numerical, N consecutive steps:* shut with a warning, matching what OPM
  already does for wells that keep failing to converge
  (`forceShutWellByName()`, `--shut-unsolvable-wells`).
