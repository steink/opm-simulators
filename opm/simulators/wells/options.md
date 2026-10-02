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
| `--group-tree-initialization` | `true` | Group-tree only. Producers without a valid previous solution get an initial solve at their own limits in `beginTimeStep()`. The pre-step network rebalance is skipped, and so are the well re-solves (A1) in a timestep's first round. | model5 MSW (iteration hold): 85 → 16 timesteps. Identical on most other decks. |
| `--group-tree-initial-balance` | `true` | Group-tree only. At the start of each timestep, balances the group tree on the wells' potentials and commits its targets before the first well solves. A well put on group control is first solved at its balancer target. | Identical on the opm-tests decks. On FLOW-CGC, removes the first step's ORAT→GRUP switches. |
| `--group-tree-max-initial-solve-failures` | `3` | Group-tree only. A producer whose initial solve does not converge is stopped for that timestep and retried at the next. After this many consecutive failures it is shut, if `--shut-unsolvable-wells` is true. | Tested on a deck set up to fail. |
| `--close-stopped-wells` | `false` | All modes. A well stopped cleanly during a timestep (no flow, not operable, re-open limit, network) and still stopped at its end is closed in the well-test state (shut or stopped per its shut-in setting until WTEST or a schedule event), without re-checking its operability. | With the timestep hold: FLOW-CGC 73 → 39 timesteps; FLOW-CGC 5-day cap, margin 3: 239 → 166, no chops, PROD2 412 → 3 status events; model5 −8.7% FOPT (no WTEST). Iteration hold: no change. |
| `--log-well-anchors` | `false` | Diagnostic only. Computes and debug-logs every predicting producer's maximum flowing THP (its "anchor") at step start. No effect on the run. | Used for the step 3 analysis. |
| `--network-analytic-jacobian` | `false` | Assembles the network Jacobian from VFP-table derivatives instead of finite differences. Also used by `group-tree`, despite its help text saying newton only. | To try on the large model. If it is no worse there, make it the group-tree default. |

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
| `--max-well-status-switch-for-wells` | `99` | Maximum open/stop switches per well per timestep. We tried `2` against FLOW-CGC's PROD2 oscillation. |
| `--max-well-status-switch-in-inner-iter-wells` | `99` | The same limit, within inner (local) well iterations. |
| `--max-inner-iter-wells` | `50` | Maximum inner iterations for standard wells. `2` forces initial-solve failures (test of the failure handling). |
| `--shut-unsolvable-wells` | `true` | Shuts wells that keep failing to converge. Also gates the shut after `--group-tree-max-initial-solve-failures`. |
| `--network-pressure-update-damping-factor` | `0.1` | Damping of the fixed-point node pressure update. Group-tree node pressures are handed to the wells undamped (C2). |
| `--network-max-pressure-update-in-bars` | `5` | Cap on the fixed-point pressure update. |
| `--network-max-outer-iterations`, `--network-max-sub-iterations` | `3`, `100` | Bounds on the network iterations. In group-tree mode the sub-iterations are the workflow's rounds. |
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
