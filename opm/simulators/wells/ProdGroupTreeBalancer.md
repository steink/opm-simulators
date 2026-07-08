# ProdGroupTreeBalancer — Algorithm Description

## Overview

The balancer distributes production targets across a group tree so that each
sub-tree respects the most restrictive of its active rate limits.  The tree
is traversed **bottom-up**: leaf subtrees are balanced first, then their
parent groups, and so on up to the FIELD node.

The key simplifying assumption is that **well phase-rate fractions are fixed**
at the current operating point.  This means that scaling a well's total rate
is equivalent to scaling each phase independently, and the strictest per-phase
constraint on a well can be estimated either from

* the well's own rate limits combined with its pressure-derived potential
  (BHP / THP), or
* the well's rate limits combined with a converged IPR at the current BHP.

The result of this pre-processing step is a single `(mode, limit)` pair per
well that is stored in the tree before balancing begins.

---

## Recursive balancing of a single subtree

`balanceGroupTree` is called with a target mode and target rate for the root
of the subtree.  It recurses through three distinct cases.

### 1. Sorted distribution to direct children

The direct guide-rate children of the current node are collected and sorted
by their **guide-rate-to-limit ratio**

```
ratio_k = min over active limits { (limit - current_sum) / guide_rate_slope }
```

A small ratio means the child will hit its limit early (before the group's
budget is exhausted); a large ratio means the child can absorb a large share
before becoming constrained.  Distributing in increasing-ratio order
guarantees that each child gets its full guide-rate share of the *remaining*
budget after earlier children have been capped at their limits.

For the common case where all direct children are individually constrained
wells this single sorted pass always produces the correct answer.  When some
children are themselves groups (item 2), the pass may need to be repeated.

### 2. Resorting when group-controlled children are present

After a group child is recursively balanced, its effective mode category can
change: it may end up at a limit (*Individual*) or still consuming a
guide-rate share (*Group*).  Because the category of earlier children
determines the remaining budget available to later children, the sort order
can become invalid once a group child has been processed.

Whenever a newly balanced child turns out to be *Individual* but a later
child is still *Group*, the algorithm re-sorts and restarts the distribution
pass (up to a configurable maximum number of re-sorts).  The goal of the
sorting is precisely to ensure that all *Individual*-destined children are
served before *Group*-controlled ones, so that the final budget split is
consistent.

### 3. Dynamic limit checking for transparent groups (no guide rate)

Some groups sit between the current node and its guide-rate children but
carry their own rate limits; because they have no guide rate of their own
they appear in the `c_trans` list rather than the sorted `c` list.

The guide-rate-to-limit ratio of a transparent group changes dynamically as
each guide-rate child is processed, because the numerator (remaining budget
at the transparent group) grows with each child allocation.  The algorithm
therefore **checks the transparent-group ratios between each child
distribution step**: if any transparent group's ratio has dropped below the
ratio of the next child to be served, that group is balanced first (its limit
becomes binding) before continuing with the remaining children.

---

## Post-balancing cleanup

After all subtrees are balanced, groups that have no limits of their own and
therefore never entered the main loop (pure reporting groups, typically above
the topmost balanced node) are updated with rates equal to the
efficiency-weighted sum of their children.

