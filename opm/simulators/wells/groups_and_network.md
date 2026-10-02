Groups and Network
------

Rapid solve using
- constant reservoir
- well IPR
- constant well fractions
- initial state wells, groups and network

Prerequisits:
- All wells have valid IPR computed at converged solution
- Stopped wells (solved with zero rate), has a valid IPR computed at a flowing solution, preferably close to maximal thp

---
Straight forward approach using group balancer:
---

Group tree balancer:
The balancer provides well-rates qw, bhps pw, active well and group controls (consistent with rates/ipr and well thps)
Assume iprs qw = a pw + b, then dqw/dpw = a

Network equations:
- Full q obtained by summing branches, $\frac{\partial q_i}{\partial p}$ is non-zero (a from ipr) for wells in sum
- Node eq i: $p_i - \Lambda(p_j, q_i) = 0$
- Well eq:
    - Individual: $pw - bhp = 0$ (given)
    - Thp: $bhp - \widetilde{\Lambda}(thp, q_i) = 0$
    - Group: $q_w - g_w\lambda = 0$
- Group eq:
    - Active: $\sum q_w - T = 0$
    - Inactive: $\lambda - T/\sum q_w$ = 0

Notes:
- Need dq/dp for linearization which is not available if ipr is exchanged with well equations (only if well eq residuals are zero). Hence, to solve fully coupled we need to include rates in formulation
- Get a system of size np+ng (node/well-pressures + group lambdas) unknwons $N(p, q(p), \lambda) = 0$ with Jacobian

    - $J_{pp} = \frac{\partial N_p}{\partial p} + \frac{\partial N_p}{\partial q}\frac{\partial q}{\partial p}$

    - $J_{pg} = \frac{\partial N_p}{\partial \lambda}$ (non-zero for group controlled wells)

    - $J_{gp} = \frac{\partial N_g}{\partial q}\frac{\partial q}{\partial p}$ (non-zero for group controlled wells)

    - $J_{gg} = \frac{\partial N_g}{\partial \lambda}$ (non-zero for inactive groups)

- Group equation is for active mode, but well equations is for *preferred* mode (use single preffered mode guide-rate)
- Special treatment of thp-control: $\widetilde{\Lambda}$ is the vfp-function with all entries removed satisfying $\frac{\partial\Lambda}{\partial FLO}  + \epsilon < \frac{\partial ipr}{\partial FLO}$ (less is flatter). This implies $ipr$ and $\widetilde\Lambda$ will always have an intersection. In practice $\widetilde\Lambda(p, FLO.,..)$ can be computed starting with $\Lambda(p, FLO.,..)$ and $\frac{\partial\Lambda}{\partial FLO}  + \epsilon \geq \frac{\partial ipr}{\partial FLO}$, jumping one flo-interval at a time until inequality no longer holds.

Iteration given initial guess for thps:
1. Run group-tree balancer for given thps -> well rates and control config
2. Setup network equations for given well rates and config
3. Solve Newton-step and update well thps

Solution:
1. Iterate until $N(p, q(p), \lambda) < \epsilon$ and control config is fixed
2. If any well has $\Lambda\neq\widetilde\Lambda$, set worst offending to stopped and go to 1
3. Resolve wells and recompute iprs. If iprs have changed more than some tollerance, go to 1 

Extensions:
- autochoke (common thp, fixed total rate)
    - no special treatment in balancer
    - network solve:
        - one extra unknown (common thp) choke is active if thp > node-pressure (p)
        - one extra equation, sum q_w - target = 0 if active, thp - p = 0 otherwise
- injecton network (any difference?) 
