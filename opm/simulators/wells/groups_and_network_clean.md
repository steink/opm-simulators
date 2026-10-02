# Groups and Network

## Approach

A rapid re-solve strategy, holding fixed:
- the reservoir state,
- each well's IPR,
- the well fractions within a group,
- the initial state of the wells, groups, and network.

**Prerequisites:**
- Every well has a valid IPR computed at the converged solution.
- A stopped well (solved with zero rate) has a valid IPR computed at a flowing solution, preferably close to its maximum THP.

## Design decisions and scope

These are the deliberate simplifications this formulation makes, so the investigation isn't cluttered by everything at once:

- **Constant phase fractions.** The phase fractions used in the VFP-table lookup (water cut, GOR) are held constant. This is exact for a standard well, whose per-phase IPR is affine in a single flow variable with fixed slopes, so total flow is affine in $bhp$ with one well-defined slope. It is only an approximation for multisegment wells, where fractions genuinely vary with rate. Handling that variation is explicitly deferred.
- **All switching lives in the balancer.** Which wells are on Group vs. Individual control, and which groups are Active, is decided entirely by the group-tree balancer, run at the start of each outer iteration. The Newton system never re-derives this categorization itself — it only assembles the continuous equations for whatever configuration it was handed, and does not need to encode the tree's nesting structure internally (see "Group equations" below).
- **Rates are eliminated, not unknowns.** The unknowns are node/well pressures (including $bhp$/$thp$) and the group multipliers $\lambda$. Rates appear throughout the equations, but always as simple (affine) functions of a well pressure, so they are substituted out rather than carried as separate unknowns.
- **Oscillation of the active-set/configuration across outer iterations is a known risk and is explicitly postponed.** For now the only safeguard in scope is a hard cap on outer iterations with a diagnostic on what failed to stabilize; a real cycle-breaker (comparable to what a discrete active-set method needs elsewhere) is future work.

---

## Straightforward approach using the group balancer

**Group-tree balancer.** The balancer provides well rates $q_w$, bottom-hole pressures $p_w$, and the active well/group controls, consistent with the rates/IPR and the well THPs. Assuming a linear IPR, $q_w = a\,p_w + b$, gives $dq_w/dp_w = a$.

**Network equations**
- Branch flow $q$ is obtained by summing over the wells feeding it; $\partial q_i/\partial p$ is non-zero (equal to $a$ from the IPR) for each well included in the sum.
- Node equation ($i$): $p_i - \Lambda(p_j, q_i) = 0$
- Well equations:
  - Individual control (rate/BHP given): $p_w - bhp = 0$
  - THP control: $bhp - \widetilde\Lambda(thp, q_i) = 0$, with $q_i = q_i(p_w)$ substituted from the well's own IPR. This is a plain scalar equation in pressures — nothing needs to be inverted for $q$.
  - Group control: $q_w - g_w \lambda = 0$
- Group equations, written only for a group the balancer has marked Active — see below for what "Active" means once the tree is flattened:
  - Active (group limit binding): $\sum q_w - T = 0$
  - Inactive: $\lambda - T / \sum q_w = 0$

**Notes**
- Linearization needs $dq/dp$, which is not available if the IPR is only used to eliminate the well equation (valid only where that residual is exactly zero). So, to solve fully coupled, the rates $q$ must remain symbolically in the formulation and be substituted algebraically, not dropped — which is exactly what the elimination above does.
- This gives a system of size $n_p + n_g$ (node/well pressures plus group multipliers $\lambda$), $N(p, q(p), \lambda) = 0$, with Jacobian blocks:
  - $J_{pp} = \dfrac{\partial N_p}{\partial p} + \dfrac{\partial N_p}{\partial q}\dfrac{\partial q}{\partial p}$
  - $J_{pg} = \dfrac{\partial N_p}{\partial \lambda}$ (non-zero for group-controlled wells)
  - $J_{gp} = \dfrac{\partial N_g}{\partial q}\dfrac{\partial q}{\partial p}$ (non-zero for group-controlled wells)
  - $J_{gg} = \dfrac{\partial N_g}{\partial \lambda}$ (non-zero for inactive groups)
- The group equation enforces the *active* control mode, while the well equation uses the well's *preferred* mode (via a single preferred-mode guide rate).
- **Group equations are flattened.** For a group the balancer marks Active, the sum $\sum q_w$ runs directly over every well leaf in that group's subtree, each weighted by the cumulative efficiency factor along its path up to the active group — bypassing any intermediate groups rather than composing their individual share equations layer by layer. An intermediate group that is itself Inactive gets no equation of its own in the Newton system at all; its rate is only recovered by summing its children *after* convergence, as a reporting step, not as part of the nonlinear solve. *(Open item: confirm the cumulative-efficiency weighting is part of what the balancer hands off per well, since two wells under the same active group can sit behind different chains of intermediate GEFACs.)*

### $\widetilde\Lambda$ and the cliff diagnosis

The purpose of $\widetilde\Lambda$ is not to approximate the tubing curve for accuracy — it is to remove the jumps a real VFP/IPR crossing can have (a well with both a stable flowing branch and an unstable/dead branch) so that the *whole* network-and-well system is smooth and can be solved to full Newton convergence with no discrete branching inside the solve. $q$ (FLO) is taken positive for production throughout this section. $\widetilde\Lambda$ is built from $\Lambda$ by discarding flow intervals until $\partial\Lambda/\partial FLO > \partial(ipr)/\partial FLO + \epsilon$ holds everywhere — flatter than the IPR (i.e. $\Lambda$ rises at least $\epsilon$ faster per unit flow than the IPR's own, already-negative slope, matching the stable branch's own behavior: $\Lambda - IPR$ strictly increasing, so at most one root), so the two curves cross at most once. Because this comparison uses a single fixed IPR slope, $\widetilde\Lambda$ depends only on the well's tables and its constant fractions, not on the current iterate, and can be built once per well up front rather than inside the Newton loop. *(Corrected 2026-09-18: the inequality's direction was flipped in the original draft.)*

Once converged, the diagnosis step compares the *unflattened* $\Lambda$ against $\widetilde\Lambda$ at the well's converged operating point:
- If they agree, the well's rate is trustworthy as computed.
- If they disagree, the converged point fell in a discarded interval — the well cannot really operate there — and it is a candidate to be stopped (see "Outer solution loop" below).

*(Open item: flattening guarantees at most one crossing, not that one exists. Whether a crossing is guaranteed to exist from the boundary behavior of $\Lambda$ vs. the IPR at $FLO = 0$ and at the table's largest flow, or whether "no crossing at all" is itself a case the diagnosis step must recognize separately from "crossing existed but was discarded," is still to be checked.)*

## Iteration

Leaning toward solving each frozen configuration to full Newton convergence before re-running the balancer, rather than taking a single Newton step per balancer call — the balancer's categorization should only be trusted on a physically consistent state, which a partially converged iterate is not. (Not fully settled; single-stepping is cheaper per outer pass and most steps change nothing, so it remains an option, especially since the difficult case is a step where *a lot* changes rather than the common case where nothing does.)

Given an initial guess for the THPs:
1. Run the group-tree balancer at the current THPs to get well rates and the control configuration (which wells are Group/Individual, which groups are Active, and — per Active group — the flattened, efficiency-weighted well list and target $T$).
2. Assemble the network equations for that configuration.
3. Solve to convergence (or take one Newton step — see above) and update the well THPs.

## Outer solution loop

1. Iterate the above until $\lVert N(p, q(p), \lambda)\rVert < \epsilon$ and the control configuration has stopped changing.
2. If any well has $\Lambda \neq \widetilde\Lambda$ at its converged point, set the worst offender to stopped and return to step 1. (One well at a time, worst first, so a single restart doesn't have to absorb several simultaneous control changes.)
3. Re-solve the wells and recompute their IPRs. If any IPR has changed by more than some tolerance, return to step 1. A well stopped in step 2 is reconsidered here — it is retried whenever its IPR changes, rather than staying stopped for the rest of the run.

## Extensions

- **Autochoke** (common THP, fixed total rate):
  - No special treatment needed in the balancer.
  - In the network solve:
    - One extra unknown: the common THP.
    - The choke is active when $thp > p$ (the node pressure).
    - One extra equation: $\sum q_w - T = 0$ when active, $thp - p = 0$ otherwise.
- **Injection network** — open question: does this need any different treatment?

## Open items to verify

- Cumulative efficiency-factor weighting in the balancer's flattened per-group well lists.
- Existence (not just uniqueness) of the $\Lambda$/IPR crossing under $\widetilde\Lambda$.
- Newton-step granularity: full convergence per frozen configuration vs. one step per balancer call.
- Oscillation/cycling of the configuration across outer iterations: deferred; only a hard iteration cap plus a diagnostic is in scope for now.
- Multisegment wells: constant-fraction assumption is an approximation there; not yet addressed.
