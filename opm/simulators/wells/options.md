# Command-line options for the group-tree network work

Status as of branch `group-tree-network` @ 3e9d0e05f (plus `--close-stopped-wells`, uncommitted). Defaults are as
registered in `opm/simulators/flow/BlackoilModelParameters.{hpp,cpp}`.

---

## 1. Introduced for the group-tree workflow

The workflow runs only with both of these set:

```
--enable-group-tree-balancer=true --network-solver=group-tree
```

It also needs a production network with no autochoke nodes, no gas-lift
optimisation and no reservoir coupling. Otherwise the old path runs and the
options below that are marked "group-tree only" have no effect.

| Option | Default | What it does | Findings |
|---|---|---|---|
| `--enable-group-tree-balancer` | `false` | Runs the group-tree balancer: it decides which wells are on individual or group control and distributes the group targets. Within NUPCOL the balancer owns production control. | Required. |
| `--network-solver` | `fixedpoint` | Takes `fixedpoint`, `newton` or `group-tree`. `group-tree` solves the production network for the balancer's configuration. The fixed point still runs alongside it for comparison. | Required (`group-tree`). |
| `--group-tree-balancer-tolerance` | `1e-4` | Convergence tolerance of the balancer's guide-rate distribution. | Not varied. |
| `--group-tree-ipr-tolerance` | `1e-2` | Relative rate mismatch between the IPRs the network solve used and the wells' own solves afterwards. Above it, another round (re-balance with fresh IPRs) runs. | Not varied. |
| `--group-tree-stop-hold` | `iteration` | How long a stop decided by the network solve is held in the wells' own solves: `iteration` (one global iteration), `timestep`, or `nupcol` (while within NUPCOL). | Varies by case. On FLOW-CGC, `timestep` is clearly best. On model5 STDW with default time-stepping, `iteration` was better. Still to measure on the large model. |
| `--group-tree-reopen-thp-margin` | `0` (bar) | A stopped producer is offered to the network solve for reopening only if its maximum flowing THP exceeds the node pressure by this margin. | 3 bar was best on the test decks (it removes FLOW-CGC's PROD2 cycling). The default stays 0 for now; to be experimented with. |
| `--group-tree-branch-flattening` | `true` | Group-tree only. Solves the network with every branch table flattened so a branch's pressure drop never falls with flow (slope limit 0), which gives a unique solution. If some branch ends on its flattened part, solves again from there with those branches' real tables. If that fails, warns (`NETWORK_FLATTENED_BRANCH`) and uses the relaxed update. | Removes all failures of the first solve on the test decks. The second solve runs on most solves (MANI-A in FLOW-CGC, C1 or B1+C1 in model5), usually in 2 iterations. It fails on model5 STDW at a few report steps (C1). See `next_steps.md`, "Network convergence". |
| `--group-tree-flattened-branch-fallback` | `relaxed` | Group-tree only, experimental. When the solve with the real branch tables after a flattened solution does not converge: `relaxed` falls back to the relaxed fixed-point update, `flattened` keeps the flattened solution. Both warn (`NETWORK_FLATTENED_BRANCH`). | Norne: `flattened` removes the fallbacks but not the chops (-6 to +17 per run); not a fix. |
| `--group-tree-network-step` | `full` | Group-tree only, experimental. Step control of the network Newton solve: `full` steps, or `linesearch` (halve the step until the residual 2-norm drops, at most 12 times). | Removes the full-step cycles (Norne E-1H) but stalls elsewhere: with the finite-difference Jacobian on a direction that is no descent direction, with the analytic one at table knots and short of a shut-in cap (Norne B-3H). Worse overall (model5 STDW timestep hold: 91 give-ups against 7). |
| `--group-tree-initialization` | `true` | Group-tree only. Producers without a valid previous solution get an initial solve at their own limits in `beginTimeStep()`. The pre-step network rebalance is skipped, and so are the well re-solves (A1) in a timestep's first round. | model5 MSW (iteration hold): 85 → 16 timesteps. Identical on most other decks. |
| `--group-tree-initial-balance` | `true` | Group-tree only. At the start of each timestep, balances the group tree on the wells' current rates, with their potentials as pressure limits, and commits its targets before the first well solves. A well put on group control is first solved at its balancer target. | Identical on the opm-tests decks. On FLOW-CGC, removes the first step's ORAT→GRUP switches. |
| `--group-tree-max-initial-solve-failures` | `3` | Group-tree only. A producer whose initial solve does not converge is stopped for that timestep and retried at the next. After this many consecutive failures it is shut, if `--shut-unsolvable-wells` is true. | Tested on a deck set up to fail. |
| `--close-stopped-wells` | `false` | All modes. A well stopped cleanly during a timestep (no flow, not operable, re-open limit, network) and still stopped at its end is closed in the well-test state (shut or stopped per its shut-in setting until WTEST or a schedule event), without re-checking its operability. | With the timestep hold: FLOW-CGC 73 → 39 timesteps; FLOW-CGC 5-day cap, margin 3: 239 → 166, no chops, PROD2 412 → 3 status events; model5 −8.7% FOPT (no WTEST). Iteration hold: no change. |
| `--log-well-anchors` | `false` | Diagnostic only. Computes and debug-logs every predicting producer's maximum flowing THP (its "anchor") at step start. No effect on the run. | Used for the step 3 analysis. |
| `--network-analytic-jacobian` | `true` (was `false`) | Assembles the network Jacobian from VFP-table derivatives instead of finite differences (`newton` and `group-tree`). | Used on the large model. The finite-difference Jacobian is wrong where its step (0.01 bar of bhp) moves a well with a steep IPR across a tubing-table knot (Norne B-3H: 100 m3/d per 0.01 bar, derivative 9.4 instead of 19.6); the Newton direction is then no descent direction and a line search stalls. |

### Only on `unstable-group-targets` and `integration`

Not on `group-tree-network`: it comes from the port of
`force-bhp-control-from-intersection-in-problematic-vfp-situations`
(`unstable-group-targets`, `fd0d83b4d`), and reaches the `integration`
worktree (`/home/steink/OPM/opm-simulators-int`) by its merge.

| Option | Default | What it does | Findings |
|---|---|---|---|
| `--stabilize-problematic-wells` | `false` | A well whose local solve converges to an unstable solution, or does not converge, is solved from the intersection of its IPR and tubing curve, and its THP constraint is replaced by the bhp found there until its next local solve (`SingleWellState::forced_bhp_from_thp`). Off: `integration` behaves as `group-tree-network`. | No measurements recorded here yet. |

---

## 2. Network options already on the branch (from `pr_network`)

Not introduced in this work, but on the branch and relevant to network runs.
Most apply to `--network-solver=newton` only.

| Option | Default | What it does |
|---|---|---|
| `--network-pressure-update-secant` | `injection` | Networks whose node pressures use the bracketing/secant update instead of the damped update: `injection`, `all` or `none`. |
| `--network-group-control` | `false` | (newton only) The network holds a group's injection total and decides the split itself. |
| `--network-autochoke` | `false` | (newton only) Solves autochoke nodes inside the simultaneous network solve. |
| `--network-autochoke-bracket-samples` | `300` | Number of samples in the legacy autochoke bracket search. |
| `--network-complementarity` | `false` | (newton with analytic Jacobian) Complementarity rows for well limits instead of an active set. |
| `--gas-lift-network-response` | `false` | (newton only) Gas-lift trial evaluations from the simultaneous network solve. |
| `--network-dump-failures` | `""` | Path prefix for dumping network systems that fail to converge. |
| `--convert-to-multisegment-well` | `none` | `per-connection` converts standard wells to multisegment wells. |

---

## 3. Existing OPM options we have used in experiments

| Option | Default | Use |
|---|---|---|
| `--output-dir` | the deck's own directory | Always set for our runs, so the output does not land next to the deck. |
| `--solver-max-time-step-in-days` | `365` | Caps the time step. The model5 decks are very sensitive to time stepping, so compare there with a cap (e.g. 1 or 5 days). |
| `--enable-tuning` | `false` | Honours the TUNING keyword, as in `regressionTests.cmake` for model5. With it, group-tree aborted earlier than fixed-point on model5 (open issue). |
| `--max-well-status-switch-for-wells` | `99` | Maximum open/stop switches per well per timestep. We tried `2` against FLOW-CGC's PROD2 oscillation. `2` also removes the open/stop cycling seen at NUPCOL 99 (FLOW-CGC 47–58 → 26 timesteps; model5 STDW 11–17), but costs 7–8% FOPT on model5 (wells kept stopped). See `next_steps.md`, "Network convergence". |
| `--max-well-status-switch-in-inner-iter-wells` | `99` | The same limit, within inner (local) well iterations. |
| `--max-inner-iter-wells` | `50` | Maximum inner iterations for standard wells. `2` forces initial-solve failures (test of the failure handling). |
| `--shut-unsolvable-wells` | `true` | Shuts wells that keep failing to converge. Also gates the shut after `--group-tree-max-initial-solve-failures`. |
| `--network-pressure-update-damping-factor` | `0.1` | Damping of the fixed-point node pressure update. Group-tree node pressures are handed to the wells undamped (C2). |
| `--network-max-pressure-update-in-bars` | `5` | Cap on the fixed-point pressure update. |
| `--network-max-outer-iterations` | `3` | Bound on the network outer iterations per Newton iteration. In group-tree mode these are the workflow's rounds (B, B4, A3): at most this many per Newton iteration, after which the loop stops with the last round's result. Within a round, B's settle loop (rebalance after a cliff stop or a reopen) is capped at 10 passes in the code (`kMaxCliffDiagnosisIters`), with no option. |
| `--network-max-sub-iterations` | `100` | Bound on the fixed-point network sub-iterations. In group-tree mode it applies only to the injection networks, which still use the fixed point. |
| `--check-group-constraints-inner-well-iterations` | `true` | Group constraints are checked in the local well solves. The initial solve under group control relies on this. |
| `--use-implicit-ipr` | `true` | Implicit IPRs, used by the balancer and the network solve. |
| `--solve-welleq-initially` | `true` | Wells are solved in `prepareTimeStep()`. Also a condition for skipping A1 in the first round. |

---

## 4. Typical command lines

```
# group-tree workflow, default settings
flow_blackoil CASE.DATA --output-dir=OUT \
    --enable-group-tree-balancer=true --network-solver=group-tree

# with the timestep hold and a reopen margin (best on FLOW-CGC)
    ... --group-tree-stop-hold=timestep --group-tree-reopen-thp-margin=3

# the step-5 initialization switched off, for comparison
    ... --group-tree-initialization=false --group-tree-initial-balance=false

# model5 decks: compare with a capped time step
    ... --solver-max-time-step-in-days=1
```

Run metrics are compared with `network_metrics.py RUN [RUN ...] --ref RUN`
(in this directory).
