# Open/stop oscillation of THP-limited wells within a timestep

Status: discussion note. Context: the group-tree network workflow
(`timestep_workflow.md`), but the problem is general: it affects any THP-limited
well, with or without a network.

---

## 1. The problem

Wells and the network are solved between global Newton iterations with the
reservoir held fixed. For a THP-limited well close to its lift limit this gives a
feedback loop that the fixed-reservoir well solve can't see:

```
global iteration k    : reservoir pressure p_k high enough → well can lift → OPEN, flows
global Newton update  : well withdraws fluid → near-well pressure drops to p_{k+1}
global iteration k+1  : at p_{k+1} the IPR no longer crosses the tubing curve → STOP
global Newton update  : no withdrawal → near-well pressure recovers toward p_k
global iteration k+2  : well can lift again → OPEN … and so on
```

The local well solve and the network solve both decide "can this well operate?"
against a reservoir state that the decision itself will change. Near the lift
limit this is expected, not a bug. With a network the effect spreads to other
wells: stopping one well lowers the flow and back-pressure at its node, which
helps its neighbours, and reopening it hurts them. So several wells can move
together.

**Why it matters.** Real models have hundreds of wells, complex group trees and
networks, and around 10⁶ cells. A global Newton iteration is by far the most
expensive unit of work, so every iteration spent flip-flopping is wasted, and a
non-converging step turns into a timestep chop, which costs even more.

**What "stopped" means across timesteps (settled).** A well still stopped at the
end of a timestep is then shut or stopped according to its shut-in instructions.
It can reopen at a later step only through, for example, a scheduled well test
(WTEST). The open question is only what to do about open/stop switching *within*
a timestep, across global Newton iterations.

### 1.1 Current safeguard

`max_well_status_switch_` (`--max-well-status-switch-for-wells`): after N
re-openings in a timestep the well is kept stopped (`WellInterface::prepareWellBeforeAssembling()`).
It works, but:

- it spends up to N round trips, each costing global Newton iterations, before
  it acts;
- the final state is whatever the well was when the cap was hit, not a decision
  based on whether the well can actually operate;
- it counts reopen *events*, not global iterations. `prepareWellBeforeAssembling()`
  runs several times per global iteration, so the budget can be used up faster
  than intended (see `timestep_workflow.md`, Q6);
- in the group-tree workflow a network-level stop decided in B is undone by the
  next global iteration's A1 local solve, so two mechanisms fight over the same
  decision.

### 1.2 Evidence from the toy model (FLOW-CGC, 4 wells)

| Variant | Timesteps | Wasted Newton | Chops |
|---|---|---|---|
| Stops held one global iteration, no switch limit | 39 | 36 % | 5 |
| Same, `--max-well-status-switch-for-wells=2` | 29 | 0 % | 0 |
| Stops held for the whole timestep (experiment) | 27 | 0 % | 0 |

All the trouble is one well (PROD2) at one report step. Cumulative production
differs by less than 1 % between the variants. The toy model shows the mechanism
but says little about how a rule scales; see section 5.

---

## 2. What a good rule needs

1. **Few extra global iterations.** Decide a well's status at most a small,
   bounded number of times per timestep, ideally once.
2. **Cheap compared with a global iteration.** Per-well or per-network work is
   fine; anything that touches the whole reservoir is not.
3. **Based on physics, not on counting.** The decision should answer "can this
   well operate over this timestep?", not "has it switched too often?".
4. **Monotone or bounded within a step.** Each well's status should change a
   bounded number of times, so convergence is guaranteed.
5. **Consistent with the end-of-step fate.** A well stopped at the end of the
   step gets shut or stopped by its instructions, so wrongly stopping a well
   that could flow costs production until the next well test. Wrongly keeping
   open a well that cannot flow costs convergence.
6. **Scales to many wells and networks.** Works per well, in parallel, and
   combines with the network's own open-set decision (B in the workflow).

---

## 3. Levels of "can this well operate?"

The strategies below differ mainly in which reservoir response the operability
test accounts for:

| Level | Reservoir in the test | Cost | Notes |
|---|---|---|---|
| L0 | Fixed (cell pressures of the current iterate) | Tiny | Today's operability check and trial IPR. Blind to the well's own drawdown. |
| L1 | Linearised near-well response | Small | Uses the already-assembled reservoir Jacobian around the perforated cells (section 4.D). |
| L2 | Nonlinear near-well region | Moderate | Well equations coupled to perforated cells plus a few neighbour layers, far field fixed (section 4.C). |
| L3 | Whole reservoir | A global iteration | What the oscillation effectively does now, several times over. |

The oscillation happens because decisions are taken at L0 and checked at L3.
Moving the decision to L1 or L2 predicts what L3 will say.

---

## 4. Strategies

### A. Status-switch cap (current)

Described in section 1.1. Worth keeping as a last-resort safety net whatever
else is done, but with two cheap fixes:
- count per global iteration, not per call;
- count only stop→open, and only after the first global iterations.

### B. One authority and a monotone rule within the timestep

- One owner for status decisions within a step. In the group-tree workflow
  that's B: the network decides stops, and reopens only via the trial IPR.
  Local solves may still *stop* a well that cannot converge, but never reopen one.
- **Monotone within the timestep:** open→stop is allowed at any time;
  stop→open only during the first k global iterations (or within NUPCOL), and
  only via an explicit reopen test. After that, the open set can only shrink,
  which guarantees the status settles.
- Tested variant: holding B's stops for the whole timestep gave the table's
  last row.
- **Risk:** a well stopped early in the step, from an unconverged reservoir
  state, may stay stopped although it could flow at the converged state. Short
  steps and the end-of-step check (strategy E) limit this.

Cheap, simple, and close to what's already implemented. A good first step,
and the frame the other strategies plug into.

### C. Near-well coupled operability solve (L2)

Solve the well equations coupled to the perforated cells and a few layers of
neighbours, over the current dt, with the far field held fixed at the current
iterate as a boundary condition. The well is operable if this local problem
converges to a flowing state at its THP limit, which here is the network node
pressure.

- **Building blocks already present:** the NLDD nonlinear domain-decomposition
  solver (`NonlinearSystemNldd.hpp`) does local reservoir-plus-well domain
  solves, and its partitioning can keep neighbour layers with a well
  (`--local-domains-partition-well-neighbor-levels`). A per-well domain
  (well cells plus n layers) can probably be built with the same machinery,
  independent of whether NLDD is the nonlinear solver in use. Not yet tested
  for this purpose.
- **When:**
  - *End of step*, before committing a stop or shut: "would it really fail to
    operate?" This is the most valuable place, because it guards the costly
    decision and runs once per step for the few wells concerned.
  - *Start of step*, for wells stopped during the step, as a cheap in-step
    alternative to a well test: a well that the near-well solve says can flow
    is offered for reopening once.
  - *Within the step*, only as the reopen test in strategy B, not every
    iteration.
- **Open points:**
  - the boundary condition on the region's edge (fixed pressure vs no-flow,
    which is optimistic vs pessimistic);
  - how many layers are enough;
  - wells whose regions overlap;
  - the network node pressure is taken as fixed during the test, which is
    consistent with B deciding the network separately.

### D. Reservoir-aware IPR (L1)

*Suggestion.* A cheaper cousin of C that uses what the global Newton has
already assembled: correct the well's IPR for the drawdown its own rate causes
over dt.

**Construction.** Take the well's perforated cells P and a small region
Ω = P plus a few neighbour layers. At the current global iterate:

1. *Reservoir response.* From the assembled reservoir Jacobian, take the block
   for Ω. The CPR pressure matrix is the natural choice: one scalar pressure per
   cell, already built for the preconditioner. Hold the far field fixed
   (δp = 0 on the boundary of Ω) and solve for the pressure change in the
   perforated cells caused by the perforation rates: `δp_P = −S q_P`, with S of
   size n_perf × n_perf. This is one small sparse solve per well, with one
   right-hand side per perforation.
2. *Perforation inflow.* `q_i = M_i (p_i + δp_i − bhp − h_i)`, where
   `M_i = T_i λ_i` is the connection transmissibility times mobility.
3. *Eliminate δp.* `q = (I + M S)⁻¹ M (p − h − bhp·1)`, so the total rate is
   still affine in bhp: `q_w = a′ − b′·bhp`, with `b′ = 1ᵀ(I + M S)⁻¹ M 1`.
   b′ is smaller than the ordinary slope `b = 1ᵀ M 1`, so the IPR is flatter,
   and the flattening is exactly the drawdown that makes the next global
   iteration stop the well.
4. *General form.* For any well type, including multisegment, this is the
   Schur complement taken the other way round. The global solve eliminates the
   well unknowns into the reservoir system (A − C D⁻¹ B). Here the local
   reservoir unknowns are eliminated into the well system,
   `D_eff = D − B A_Ω⁻¹ C`, and the existing implicit-IPR computation
   (`updateIPRImplicit()`) runs with D_eff instead of D.

The accumulation term, pore volume × compressibility / dt, sits on the diagonal
of A_Ω. So S grows with dt: a short step gives almost the ordinary IPR, a long
step gives more drawdown. That's the dependence on step size the oscillation
itself has.

**Stopped wells: no lag needed for the reservoir part.** A reopen test needs
this IPR for a well that isn't flowing.
- A_Ω is the cells' accumulation and flux terms. OPM assembles it for every
  cell every iteration, and the well terms are applied separately (unless
  `--matrix-add-well-contributions` is on, in which case they must be left
  out). A stopped well's cells have a valid A_Ω.
- What does depend on the well flowing is M, the connection mobilities and
  phase fractions. These are degenerate at zero rate. But the trial solve
  (`updateStoppedWellTrialIpr()`) already linearises the well at a flowing
  trial state against the current cells, and its well blocks (D, B, C) in the
  scratch copy are what the formula needs.
- The current iterate is the right linearisation point. The question is "if
  this well reopens, what will the next global Newton step do to its cell
  pressures?" That step linearises at the current iterate, so A_Ω there
  predicts it to first order. A Jacobian lagged from the last iteration where
  the well flowed belongs to a reservoir state the Newton step won't start from.

**Where it falls short: nonlinearity.** When a well stops, the near-well state
relaxes: pressure recovers, and coning or saturation build-up near the well
reverses. A_Ω and the trial mobilities, linearised in that relaxed state, can
be optimistic about the well once it flows. A first-order method can't capture
that whatever point it linearises at; a Jacobian lagged from the flowing state
would err the other way (pessimistic, and stale). Two handles:
- the near-well nonlinear solve (C), the principled fix, used where the
  decision matters most (end of step);
- *a measured response from the oscillation itself.* Once a well has flipped
  once, the full reservoir response has been observed for free: cell pressures
  p_open at iteration k with rate q, then p_stopped at k+1. That gives, per
  perforation (diagonal), `S_emp,i ≈ (p_stopped,i − p_open,i) / q_i`. It
  includes the far field and the nonlinearity along that path. Use the
  Jacobian-based S by default, and switch to the measured one (or the larger of
  the two) once a well has oscillated. The first flip then becomes information
  rather than waste, and the second shouldn't happen.

**Where to use it.** In the operability check (the stop decision) and the trial
IPR (the reopen test), and optionally in B's `GroupTreeSystem` for THP wells
near their limit, so B decides against the pressure the well will actually
see. Stop and reopen must use the same IPR: if only the reopen test used it,
the asymmetry would itself produce oscillation.

**Practical points.**
- Other wells perforating Ω: include their current linearisations in A_Ω, or
  ignore them in a first version.
- MPI: Ω may cross rank boundaries; a first version can be limited to owned
  cells plus the overlap layer.
- Cost: a region of tens to a few hundred cells per THP well near its limit,
  solved only when a status decision is at stake, from an already-assembled
  matrix. No extra assembly, no extra global iteration. It combines well with
  C: D inside the step, C at the end.

### E. Decide at step start, verify at step end

*Suggestion.* Treat status as explicit within the step:

- **At the start of the step:** decide each THP well's status once, from the
  converged state of the previous step (optionally checked with L1 or L2). Keep
  it fixed for all global iterations of the step. Inside the step a well that
  cannot operate stays in the equations at the minimum stable rate or on its
  zero-rate branch; its status doesn't flip.
- **At the end of the step:** verify with the converged reservoir state. A well
  that can't operate is stopped or shut per its instructions (that's already
  the end-of-step semantics). If the mismatch is large, e.g. a well that was
  held open actually produced nothing useful, optionally repeat the step with
  the corrected status, which still costs far less than oscillating.

This removes in-step oscillation entirely: at most one extra step solve, only
when the end-of-step check disagrees. The price is a first-order lag in the
status decision, which is small when steps are short relative to how fast
wells die. It pairs naturally with C at the end-of-step check.

### F. Margin (hysteresis) around the operability boundary

*Suggestion.* Stop only if the well fails with a margin, and reopen only if it
passes with a margin. For example, the IPR/tubing crossing must exist with a
rate or bhp margin δ, or the trial rate must exceed a minimum economic or
stable rate. Wells right at the boundary then stop flip-flopping. Very cheap
and easy to add to the existing checks and to B's reopen outcome, but δ needs
tuning and does not by itself guarantee termination. Best as an add-on to B.

### G. Detect oscillation and damp the well's source term

*Suggestion, more speculative.* When a well's status pattern within a step reads
open/stop/open, stop switching it and instead limit the change of its source
term between global iterations (relax its rate toward the new value), so the
reservoir sees a gradual change. This mimics what a coupled solve would do,
at the cost of Newton consistency: the global residual converges only once the
relaxation is released. Useful as a fallback for the rare wells none of the
above settles.

### H. Move the status into the global system

*Longer term.* Make the open/stop choice part of the global Newton as a
complementarity condition, e.g. rate ≥ 0, lift margin ≥ 0 and their product
= 0, solved with a semi-smooth Newton, with the slope-limited tubing curve
keeping the lifting branch unique. The reservoir response is then implicit and
no oscillation can occur. The catch is the one already noted: the flattened
curve is only an approximation, so a well converging onto the flattened part
needs a re-solve, and robustness of a semi-smooth Newton on the full
reservoir system is a project in itself. Worth keeping in mind as the fully
coupled end state, together with removing NUPCOL (`timestep_workflow.md`, Q4).

---

## 5. Suggested path

1. **Now, cheap (B + A fixes).** One authority, and a monotone rule within the
   timestep: stops allowed, reopening only through B's trial-IPR test and only
   in the first k global iterations or within NUPCOL. Fix the switch counter to
   count per global iteration. Keep the switch cap as a safety net. Mostly
   already implemented; the whole-timestep hold experiment is a variant of it.
2. **Next, reservoir-aware decisions (D, optionally F).** Build the
   linearised near-well response from the assembled Jacobian and use it in the
   operability check and the trial IPR, so decisions within the step anticipate
   the drawdown; switch to the measured response for wells that have already
   oscillated once.
3. **Then, end-of-step verification (C at the end of step, possibly E).**
   Before a stop or shut is committed at the end of a step, run the near-well
   coupled solve, using the NLDD machinery. Decide whether the explicit
   step-start variant (E) is worth adopting based on the results of steps 1–2.
4. **Longer term (H).** Implicit status in the global system, alongside
   removing NUPCOL.

## 6. How to evaluate

- **Metrics per run:**
  - global Newton iterations, both total and wasted;
  - timestep chops;
  - status changes per well per timestep;
  - lost production against a reference run with small fixed timesteps, where
    the reservoir barely moves between iterations and the status decisions
    are close to exact.
- **Test cases:** FLOW-CGC shows the mechanism with one well. Needed as well:
  - a model with many THP wells near their lift limit on a shared network,
    where neighbours interact;
  - a case where a well genuinely dies during a step and must be stopped;
  - a case where a well is only transiently unable to lift, and stopping it
    would be wrong.
- **Scaling:** report the cost of the operability tests (L1/L2) per timestep
  against the cost of one global iteration, as well count grows.
