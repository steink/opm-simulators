# Well / group / network workflow within one timestep

Scope: everything the well model does from the start of a timestep up to the
well-side work that precedes each global (reservoir) Newton linear solve, with
`--enable-group-tree-balancer=true --network-solver=group-tree`. Line numbers
are as of branch `group-tree-network` @ 56b59882e plus the uncommitted Part 3/4g
and slope-limit work; function names are the stable reference.

Terminology used below:

| Level | Name | Driven by |
|---|---|---|
| G | global Newton iteration | reservoir nonlinear solver, calls `assemble()` once per iteration |
| O | outer network iteration | `updateWellControlsAndNetwork()`'s `while` loop |
| C | cliff-diagnosis round | new do-while in `updateWellControlsAndNetworkIteration()` |
| S | inner network sub-iteration | `BlackoilWellModelNetwork::update()`'s `for` loop |
| L | local well solve | `WellInterface::prepareWellBeforeAssembling()` / `solveWellEquation()` |

---

## 1. Current workflow

### 1.1 Once per timestep — `BlackoilWellModel::beginTimeStep()` ([BlackoilWellModel_impl.hpp:334](BlackoilWellModel_impl.hpp#L334))

1. `resetWGState()`, default ALQ, gas-lift init.
2. `wellTesting()` — WTEST reopen attempts for closed wells (on a scratch well state).
3. `createWellContainer()` / `initWellContainer()`, efficiency factors.
4. `closeCompletions()` — economic completion closures.
5. `updateWellPotentials(onlyAfterEvent=true)` → guide-rate update.
6. `updateAndCommunicateGroupData(update_wellgrouptarget=true)` — **standard OPM**
   group-target allocation into `ws.group_target`.
7. For wells with a status/type event: `solveWellEquation()` (L).

### 1.2 First global iteration only — `assemble()` → `prepareTimeStep()` ([:2283](BlackoilWellModel_impl.hpp#L2283))

1. Wells with control events: `updateWellStateWithTarget()`.
2. If `solve_welleq_initially_`: `solveWellEquation()` for every operable well (L).
3. `resetWellOperability()`.
4. If a pre-step rebalance is due: `doPreStepRebalance()`
   ([BlackoilWellModelNetwork_impl.hpp:60](BlackoilWellModelNetwork_impl.hpp#L60)),
   which runs the **whole of 1.3 below** with `mandatory_network_balance=true`,
   then one `solveEqAndUpdateWellState()` per well.

### 1.3 Every global iteration — `assemble()` ([:1151](BlackoilWellModel_impl.hpp#L1151))

```
G  assemble()
   └─ updateWellControlsAndNetwork()                                   [:1225]
      O  while (do_network_update && iter < network_max_outer_iterations_)
         └─ updateWellControlsAndNetworkIteration()                    [:1294]
            a. updateAndCommunicateGroupData(true)          ← std. group targets
            b. updateWellControls()                          ← std. OPM switching
                 loop ≤ well_group_constraints_max_iterations_:
                   updateGroupControls(FIELD)
                   per well: updateWellControl(Group), updateWellControl(Individual)
                 → well_group_control_changed
            c. [NEW] if balancer && well_group_control_changed && withinNupcol
                     && production network && group-tree:
               C  do {
                    prepareWellsForBalancing_()   (reads ws.surface_rates, cmode)
                    balanceGroupTree()
                    setBalancedGroupTree(tree)
                    per root: stopWorstGroupTreeCliffViolation()
                              = solveGroupTree() + worstCliffViolation()
                              → well->stopWell() on the worst offender
                  } while (stopped && rounds < 10)
            d. network_.update()                             [Network_impl.hpp:94]
                 computeWellGroupThp()                       (autochoke)
                 per predicting well: updateIPRImplicit(),
                                      updateStoppedWellTrialIpr(), VFP datum dp
               S for i < network_max_sub_iterations_:
                   updatePressures()                         [NetworkGeneric.cpp:1300]
                     fixed-point computePressures()          (always, baseline)
                     group-tree: per root solveGroupTree() → node pressures only
                     per node: secant / damped step toward computed pressure
                     imbalance = size of that step
                     per well: setDynamicThpLimit(node pressure)
                   if imbalance ≤ tol: break
                   resolve(): prepareWellBeforeAssembling() for every predicting
                              producer, then injector                     (L)
                   updateAndCommunicateGroupData(true)
                 → more_network_update = sub-loop unconverged || group thp moved
            e. gas-lift optimisation (if enabled)
            f. prepareWellsBeforeAssembling()  → every well                (L)
            g. [OLD] if balancer && changed && withinNupcol && !(c ran):
                 balanceGroupTree(); setBalancedGroupTree() or applyTreeToState()
            h. guide-rate update if needed
            → do_network_update = shouldBalance && (more_network_update || alq)
   assembleWellEqWithoutIteration(); updateCellRates()
   ⟶ global linear solve
```

### 1.4 Inside a local well solve — `prepareWellBeforeAssembling()` ([WellInterface_impl.hpp:1183](WellInterface_impl.hpp#L1183))

- `checkWellOperability()` (if enabled).
- Only while `shouldRunInnerWellIterations(max_niter_inner_well_iter_)`:
  - if `number_of_well_reopenings_ ≥ max_well_status_switch_`: `stopWell()`,
    zero-rate solve, return (the existing anti-oscillation net).
  - `iterateWellEquations()` → for an OPEN producer with implicit IPR,
    `solveWellWithOperabilityCheck()` ([:728](WellInterface_impl.hpp#L728)):
    - if `wellIsStopped()`: **`openWell()`**, `estimateOperableBhp()`; none found →
      zero-rate solve + `stopWell()`, otherwise re-solve from that bhp;
    - `iterateWellEqWithSwitching()` (local Newton with control switching);
    - THP control: `updateIPRImplicit()`, stability check, `estimateStableBhp()`
      re-solve if unstable;
    - unconverged: explicit fractions, `estimateOperableBhp()` again.
  - converged but stopped although it had rate before: `openWell()`,
    `number_of_well_reopenings_++`.
- not operable → `stopWell()`, zero-rate solve.

So one global iteration performs local well solves in at least three places (d/S
`resolve()`, f, plus 1.2 on the first iteration), each of which may stop or
reopen a well and refresh its implicit IPR.

---

## 2. What makes the current workflow fragile

**P1 — Two independent control categorizations.** Standard OPM switching (1.3 b,
`updateGroupControls` / `updateWellControl`) decides each well's `production_cmode`
and GRUP target (via the guide-rate target calculator, 1.3 a). The balancer decides
its own Group / Individual / THP categorization and its own allocation. With a
production network, `applyTreeToState()` is never called, so the balancer's
categorization only shapes the `GroupTreeSystem`; the wells themselves keep what
standard switching decided. The network solve and the well solves can therefore
disagree on which wells are GRUP, individually limited or on THP, and on what the
group targets are.

**P2 — The coupled solve is reduced to node pressures and then relaxed again.**
`updatePressures()` keeps only the node pressures from `solveGroupTree()`; the
solved well rates, bhps and λ are discarded. The pressures are then secant- /
damped-stepped toward (lines 1487–1521 — "a node a simultaneous solve has placed …
still need[s] the step to it bounded") before becoming THP limits, and the wells
re-derive their rates in their own local solves. A fully converged coupled solve is
thus fed back into the same pressure ↔ well-solve fixed point that the plain
fixed-point solver relies on. There are two sources of truth for well rates.

**P3 — IPRs refresh at uncontrolled points.** Implicit IPRs are refreshed once per
`network_.update()` (d), but also inside every local well solve on THP wells
(1.4), which happens in S `resolve()`, in f, and in 1.2. The balancer (c) and the
network solve (S) therefore see IPRs from different stages of the same outer
iteration.

**P4 — Stop / reopen has no single owner.** Stopping is decided inside local solves
(operability, zero-rate targets), with a reopen counter as the only damping, and
`solveWellWithOperabilityCheck()` reopens a stopped well and retries it on its very
next solve. The cliff diagnosis (c) uses the same `well->stopWell()`, so a
diagnosis-stopped well is retried at the next local solve — in d/S or f of the
*same* outer iteration. A network-level decision cannot persist even for one outer
iteration, and the well-level and network-level views of "stopped" fight.

**P5 — Loop criteria don't match the work each loop does.**
- The outer loop O continues only on network imbalance or ALQ change — not on
  `well_group_control_changed` or on a cliff stop.
- The cliff loop C runs only when standard switching (b) happened to move
  something, and only within NUPCOL; after NUPCOL no cliff violation is ever
  diagnosed.
- With the group-tree solver, the network system is already solved to convergence
  in each S sub-iteration; what S actually iterates is the pressure ↔ well-solve
  fixed point of P2.

**P6 — Balancer input freshness.** `prepareWellsForBalancing_()` reads
`ws.surface_rates` / `production_cmode` from whatever local solve ran last. In c
that is after (b) but before any well is re-solved in this iteration — which is why
early iterations with all-zero rates produced empty trees ("FIELD is not a node in
the balanced group tree").

**P7 — Three definitions of a participating well.**
`prepareWellsForBalancing_()` excludes zero-rate and `wellIsStopped()` wells;
`gatherWellNetworkDataForGroupTree()` keeps STOP wells with `has_ipr=false` and
drops SHUT ones; `updatePressures()` sets THP limits on every predicting well with
a node.

---

## 3. Proposed target workflow

Principle: **one owner per decision, one source of truth per quantity, and loops
nested by cost** — the cheap network/group loop settles fully for a fixed set of
IPRs before the expensive well solves run again.

```
G  per global Newton iteration (within NUPCOL: full; after NUPCOL: see Q4)
   A1  local well solves at the current constraints        (once per G: the
         → rates, bhp, implicit IPR (flowing wells),        reservoir moved, so
           trial IPR (stopped wells, Part 1b)               last IPRs are stale)
   A   well-side loop (≤ N_A rounds)
       A2  network/group settle loop, IPRs frozen
           O := wells open after the last well solve (A1 / A3)
           R := stopped wells with a valid trial IPR    (reopen candidates)
           B   pass loop (≤ N_B passes)
             B1  balance with O — the balancer owns Group / Individual / THP;
                 stopped wells are omitted. (Whether it may also run on a
                 solution with a well on the flattened curve is open — Q7.)
             B2/B3  settle the open set for this balance
                 first pass only: O := O ∪ R — each candidate enters as a Thp
                   well on its trial IPR, starting at rate 0 / bhp at shut-in
                 B2  solve GroupTreeSystem with slope-limited Λ̃
                 B3  any well on the flattened curve (SlopeLimit ≠ Unflattened)?
                       → stop the worst one (remove from O), back to B2
                     opened candidates that settled at their shut-in cap (no
                       flow) stay stopped — not a cliff stop, no re-solve needed
                     else stable: every open well on the real curve
                 [variant B-ii, Q7: rebalance (B1) after each stop instead]
             O unchanged since B1? → B4
             else → B1 with the new O. R is not offered again this A round,
                    so after the first pass O only shrinks and the loop ends.
           B4  commit: open / stop status, node pressures as THP limits, and the
               balancer's categorization + targets (applyTreeToState-equivalent,
               also with a network); the coupled solve's bhp / rates only as a
               seed for A3
       A3  local well solves at the committed constraints → fresh rates, bhp,
           implicit IPR, trial IPR (same operation as A1)
           repeat from A2 with A3's IPRs (no re-solve: A3 already produced
           them, and its fresh trial IPRs give stopped wells a new chance as R)
           if either
             - the IPRs B used are off: per well, the rate mismatch at the new
               operating point |q_solved − (ipr_a_B + ipr_b_B · bhp_solved)|
               exceeds tol, or
             - a well B opened fails to operate in its own local solve
           else done — the wells are already solved at the committed
           constraints, so no separate pre-assembly solve (today's step f)
   injection network update (fixed-point path), after production has settled
   assemble well equations ⟶ global linear solve
```

### Changes this implies

- **C1 — Single categorization authority.** The balancer's tree covers every
  production group and well from FIELD down, networked or not, so with the
  balancer enabled its categorization and targets are committed to all of them
  (B4; the no-network case already does this via `applyTreeToState()`).
  Standard `updateWellControls()` / `updateGroupControls()` switching is bypassed
  for all production wells and groups; injectors and injection group controls
  keep standard switching, since the balancer is production-only. (P1)
  Consequences to handle:
  - `updateGroupControls()` checks a group's production and injection
    constraints together, so the bypass has to separate the two.
  - `well_group_control_changed` currently triggers the balancer (1.3 c); with
    switching bypassed it never fires, so B needs its own trigger.
  - The same flag tells the global Newton "controls moved, not converged" (via
    `last_report_.well_group_control_changed`); the A/B loop must report its own
    categorization changes there instead.
- **C2 — Hand the wells the network's answer, keep their THP solve.** The global
  solve runs Thp wells at fixed THP, so A1 keeps a local THP solve: B4 commits
  node pressures as THP limits, without secant/damping (damping kept only as an
  explicit safeguard, e.g. when B did not converge), and the coupled solve's
  bhp/rates are only a seed for A3. B's bhp need not match the local/global bhp
  within a timestep; A3 measures that mismatch (the IPR linearisation error).
  Agreement at the end of the timestep only follows if balancing continues until
  global convergence — not the case after NUPCOL (Q4). (P2)
- **C3 — Open / stop decided in B, verified in A3.** Within A2 the open set is
  the network's own decision: reopen candidates (stopped wells with a valid trial
  IPR) are offered once per A round, and wells are stopped one at a time until
  every open well sits on the real tubing curve. B4 commits the result to the
  wells' own status; `solveWellWithOperabilityCheck()` must not undo it on its
  own. A3's local solves then verify it — a well B opened that cannot operate
  on its own is a reason for another A round — and supply fresh trial IPRs, so
  a well stopped this round is reconsidered next round rather than retried
  immediately. Open until Q3 is settled: how long a stop persists beyond the
  current global iteration. (P4)
- **Reopen needs the deferred Part 1b mechanism.** A reopen candidate enters
  `GroupTreeSystem` as a Thp well on its trial IPR (`stopped_ipr_a/b`),
  starting at rate 0 / bhp at shut-in, so the solve itself decides whether the
  node pressure lets it flow.
- **C4 — Well solves in two places only.** Wells are solved, and their IPRs
  refreshed, in A1 (once per global iteration) and A3 (once per A round), and
  IPRs are frozen throughout A2. For group-tree, today's re-solves inside S
  (`resolve()`) and in step f are replaced by A3. (P3, P5)
- **C5 — Loop criteria match their work.** (Written for variant B-i; with
  B-ii the same bounds hold, since every rebalance follows a stop.) The B2/B3
  loop ends when no open
  well is on the flattened curve — bounded by the size of the open set, since
  each round removes a well. The B pass loop ends when a pass leaves the open set
  unchanged — bounded because candidates are offered only in the first pass,
  after which the set only shrinks (cap N_B as a safeguard, warn on hitting it).
  A ends when the IPR rate mismatch at the new operating points is below
  tolerance and every well B opened operates in its own solve (cap N_A). The
  global Newton sees non-convergence through the categorization-changed /
  network-balance flags. (P5)
- **C6 — Balancer inputs from A1/A3.** B1 always runs right after a well solve
  (A1 or A3), so its inputs are never stale; an empty candidate set becomes
  "nothing to balance, keep the previous tree" rather than a give-up. (P6)
- **C7 — One participation rule.** OPEN producers in prediction mode participate
  normally; STOP (well-level or network-level) participate pinned at zero, with a
  trial IPR for reopening; SHUT never participate. Applied identically in the
  balancer candidate list, the network gather and the THP-limit write-back. (P7)

---

## 4. Discussion outcomes

- **Q1 — Bypass standard switching?** Yes. It was kept during testing only to
  check agreement. It applies to *all* production wells and groups, not just
  networked ones, since the balancer covers the whole production tree;
  injection keeps standard switching. See C1 for the consequences.
- **Q2 — Thp well bhp from B, or a local THP solve?** Local THP solve, seeded
  from B, because the global solve itself runs at fixed THP. B's bhp and the
  global solve's bhp need not agree within the timestep, only at its end. See
  C2.
- **Q3 — Scope of a network-level stop.** Within a global iteration this is now
  settled by C3 (B decides, A3 verifies, stopped wells are reconsidered each A
  round). Postponed: whether a stop should persist across global iterations.
  Context: a well still stopped at the end of a timestep is then shut or
  stopped for good according to its shut-in instructions, and can reopen at a
  later step only through e.g. a scheduled well test (WTEST) — how PROD2
  comes back on FLOW-CGC. So the open question is only the open/stop
  oscillation *within* a timestep, over global Newton iterations; the
  whole-timestep hold experiment in step 4's status is one data point for it.
- **Q4 — After NUPCOL.** Currently both the categorization and the network
  pressures freeze: `shouldBalance()` balances only while `withinNupcol` in
  NETBALAN's NUPCOL mode (only on the first global iteration in `TimeStepStart`
  mode), and the balancer is gated on `withinNupcol` too. Stop/open oscillations
  after that are handled by the ad-hoc reopen counter. Long term NUPCOL should
  go away altogether: it exists because the network/group system is not solved
  fully coupled with the reservoir, and this work is one step towards that.
- **Q5 — Injection network.** With the reservoir frozen, VREP/REIN is one-way,
  production → injection, with no feedback; it only requires the production
  network to be solved first. So the injection update runs after the production
  A loop has settled (diagram in section 3), instead of being iterated jointly
  with it as in today's outer loop. (Today's sub-iteration `resolve()` already
  re-solves producers before injectors.)
- **Q6 — `number_of_well_reopenings_` / `max_well_status_switch_`.** Deferred on
  purpose: both are ad-hoc safeguards whose role may change once the target
  workflow is running, so they are kept in mind rather than redesigned now.
  Facts to remember when revisiting:
  - The counter is never reset explicitly; it restarts each timestep only
    because `beginTimeStep()` recreates the well objects
    (`createWellContainer()`).
  - It counts reopen *events* in `prepareWellBeforeAssembling()`, which runs
    several times per global iteration (each S `resolve()`, plus step f), and
    line 1228 increments it on every call once the cap is reached. So one
    global iteration can use up several of `max_well_status_switch_`.
  - Intended meaning (a count over global iterations, enforced in A1 by keeping
    a closed well closed once the limit is reached) would need at most one
    increment per global iteration — otherwise several A rounds per iteration
    use up the limit as fast as today.

- **Q7 — Balancer inside or outside the B2/B3 stop loop? Open.** Both
  variants are to be kept possible and compared on real decks later; the
  A/B restructure should not hard-wire either.
  - *B-i (as drawn in section 3):* the stop loop runs with the categorization
    fixed; the balancer runs only once the open set is settled, so it never
    sees a solution with a well on the flattened curve. Argument for: the
    balancer only ever reads a physically valid state.
  - *B-ii:* rebalance after each stop, i.e. stop → B1 → B2. Argument for:
    stopping a well changes what its group can deliver, so the categorization
    (and hence which other wells are on THP, and whether they are past their
    cliff) may change after every stop; with B-i the next stop decision is
    made against a categorization that is already stale. This is what the
    current cliff-diagnosis loop does.
  - Either way, a pass ending with O unchanged since the last balance ends B.

### Suggested order of implementation

1. C1 — commit the balancer's categorization and targets to all production
   wells and groups, and bypass standard production switching (removes the
   biggest inconsistency).
2. C2 (first half) — no damping for group-tree-solved nodes.
3. C4 / C5 / C6 — restructure into the A / B loops, with B's stop-one-at-a-time
   loop (the current cliff diagnosis, moved inside B) and B4 committing the
   open / stop status. Keep the balancer placement in B switchable (Q7).
   *Status:* first version in place (`useGroupTreeWorkflow_()`,
   `updateGroupTreeNetworkIteration_()`): the outer loop's rounds are the A
   rounds (A1 on the first, then IPR refresh, B = `settleGroupTree_()` in
   variant B-ii, B4 commit + production-only `updatePressures()`, A3); a round
   repeats while node pressures move or the IPR mismatch exceeds
   `--group-tree-ipr-tolerance`; wells B stops are held stopped for the rest
   of the global iteration; injection networks run after the rounds. Not yet:
   B-i, per-round injection coupling, autochoke / gas lift / reservoir
   coupling (these fall back to the old path).
4. C3 reopen — reopen candidates entering B2 as Thp wells on their trial IPR
   (the deferred Part 1b consumption).
   *Status:* in place. Trial IPRs are computed for dynamically stopped wells
   (well temporarily opened for the trial). In B's first pass, dynamically
   stopped producers (not persistent STOP) with a trial IPR enter
   `GroupTreeSystem` as Thp candidates at shut-in, counted in their nearest
   Active ancestor's sum (`updateGroupTreeOpenSet()`); one that flows on the
   real tubing curve is reopened and seeded with the network's operating
   point. A Thp well's shut-in cap is now released when its tubing curve at
   zero rate needs less than the shut-in bhp. After B every stopped producer
   is held; a reopened well its own A3 solve stops again counts as IPR
   mismatch 1. On FLOW-CGC candidates rarely arise, because A1 at the next
   global iteration reopens B's stops first (Q3). Holding B's stops for the
   whole timestep instead (tested, not kept) gave 27 timesteps, 0 wasted
   iterations, no switch limit needed.

### Known limitations (left for later)

- **GCONPROD actions other than RATE.** The balancer only handles the RATE
  action on a violated group limit. Since C1, standard switching is bypassed
  for production while the balancer owns it (within NUPCOL), so the other
  actions (WELL, CON, CON+, PLUG, ...) are not taken at all during that
  window. These actions belong outside the balancer/network work and should
  be handled separately. Parked for now.
