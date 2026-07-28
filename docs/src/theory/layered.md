# The layered solution and filtering

This page defines the semantics `resolve` implements — the *layered
solution* — and proves the two guarantees that make the implementation
trustworthy:

1. **filtering preserves satisfiability** — a solvable problem never
   becomes unsolvable by filtering; and
2. **filtering preserves the answer** — the solution computed on the
   filtered problem equals the layered solution of the *raw* problem, so
   the semantics is independent of filtering details.

The second result is the important one: it means the layered solution has a
canonical, filter-free definition, and the filter is licensed by a theorem
rather than trusted as part of the spec. It also yields a general contract
under which *future, more aggressive* reductions are automatically safe.

## The layered solution

Fix a problem, requirements, and a strict total priority order ``\prec`` on
packages (from the `by` argument). The **layered trace** is the following
deterministic procedure over the problem's models (= valid solutions
covering the requirements, junk allowed):

1. If no model exists, the answer is *unsatisfiable* (`nothing`).
2. Otherwise set the current layer ``L := \mathrm{reqs}`` and start with no
   pins.
3. For each ``q \in L`` in ``\prec``-order: pin ``q`` at the **best rank
   feasible given the pins so far** — the minimum
   ``\mathrm{rank}_M(q)`` over models ``M`` that extend all current pins.
4. Let the next layer be the not-yet-visited dependencies of the versions
   just pinned; if nonempty, go to 3.
5. The answer ``\mathrm{Lay}(D, \mathrm{reqs}, \prec)`` is the set of pins.

!!! note "Theorem 1 (well-definedness)"
    The layered trace is deterministic — it depends only on the problem,
    the requirements, and the priority order, not on which models a solver
    happens to produce — and its answer is a tight solution.

    *Proof.* Each step-3 choice is a minimum over a set determined by the
    problem and the current pins; layers are determined by pinned versions'
    dependency sets; so the trace is a deterministic recursion.
    Feasibility is maintained at every pin (the minimizing model extends
    the pins), so the final pin set extends to a model; the pins are
    exactly the dependency closure of the requirements under the pinned
    versions (layer construction), hence a tight solution. ∎

The layered solution is *pointwise optimal*: each package, taken in layer
order and priority order, is at the best version compatible with the
choices already made. It is also Pareto-optimal in the sense made precise
on the [next page](front.md) — no tight solution beats it on every package
— though Pareto-optimality is a strictly weaker statement than the trace
itself, which pins the unique answer.

## Filtering

Real problems are shrunk before solving by two passes (`pkg_info` with
`filter=true`, the default):

**Reachability** (`find_reachable`) computes, per package, a version
*prefix* to keep, by a fixpoint of *conflict-driven degradation*:

- requirements' first versions are reachable;
- the first version of each dependency of a reachable version is reachable;
- if a reachable version conflicts with a reachable version, both packages'
  prefixes extend one further;
- if every version of a package has been pushed past (the package
  *saturates* — meaning conflicts could make it uninstallable), packages
  depending on it have their prefixes pushed past the depending versions —
  and a saturated package keeps **all** its versions.

**Redundancy elimination** (`mark_necessary!`) then deletes any version
``q@j`` for which some better version ``q@i`` (``i < j``) has a *subset* of
its constraints — every dependency and every (surviving) conflict of
``q@i`` is also one of ``q@j`` — iterated to a fixpoint.

The intuition behind reachability is: *an optimal solution only settles for
a worse version for a reason, and the reason is a conflict.* That intuition
is exactly right for the layered (pointwise) semantics, and the theorems
below make it precise. (It is *not* right for Pareto-front reasoning — see
the [next page](front.md#Why-the-filter-cannot-serve-the-front) — which is
one of the sharpest ways to see the difference between the two notions of
optimality.)

## Theorem 2: the reach fixpoint invariant

Everything rests on one observation about what the reachability fixpoint
guarantees. Call a kept version *last* if it is the highest-rank (worst)
kept version of an unsaturated package.

!!! note "Theorem 2 (fixpoint invariant)"
    At the reachability fixpoint, the last kept version ``p@r(p)`` of every
    unsaturated package ``p``:

    1. conflicts with **no kept version of any package** (saturated
       packages' versions included); and
    2. depends on **no saturated package**.

    *Proof.* If ``p@r(p)`` conflicted with a kept ``q@k``, then whichever
    of the two was processed later scans the other (already reachable) and
    fires the degradation rule, extending ``p``'s prefix past ``r(p)`` —
    contradicting that ``r(p)`` is last. If ``p@r(p)`` depended on a
    saturated package, the saturation rule pushes dependents past the
    depending version, likewise extending the prefix. ∎

In other words: **the cascade pushes prefixes exactly until the
"all-worst-kept" assignment becomes self-consistent.** That immediately
gives the first correctness result.

## Theorem 3: filtering preserves satisfiability

!!! note "Theorem 3"
    If the raw problem has a solution covering the requirements, so does
    the filtered problem. (The converse is trivial, since filtered models
    are a subset of raw models.)

    *Proof.* Let ``S`` be a raw solution. First suppose no requirement's
    dependency chain touches a saturated package. Build ``T`` as the
    dependency closure from the requirements choosing every package's
    *last kept* version. By invariant (2) this closure never enters a
    saturated package; by invariant (1) every chosen version conflicts
    with nothing kept, so ``T`` is conflict-free; it is dependency-closed
    by construction and covers the requirements. ``T`` is a filtered
    solution.

    In general, splice: for saturated packages use ``S``'s own versions —
    all kept, since saturation drops nothing, and ``S`` is
    dependency-closed so saturated dependencies of ``S``-members are
    assigned in ``S`` — and for unsaturated packages use last kept
    versions. Conflicts between an unsaturated choice and anything vanish
    by invariant (1); conflicts among the saturated choices are absent
    because they come from the valid ``S``. Dependency edges are satisfied
    branch by branch: unsaturated members' dependencies are unsaturated by
    invariant (2) and get last-kept versions; saturated members'
    dependencies are resolved inside ``S`` or by last-kept substitutes
    (dependency edges are package-level, so any assigned version satisfies
    them). ∎

For the redundancy pass, satisfiability is preserved by direct
substitution: any solution using a deleted ``q@j`` maps to a solution using
its dominating ``q@i`` (fewer dependencies to satisfy, fewer conflicts to
avoid), and domination chains terminate at a surviving version because
reachability keeps prefixes — a better version is never deleted while a
worse one survives.

## Theorem 4: reduction invariance — the optimization license

The next theorem is the general principle; the filter's correctness is an
instance of it.

!!! note "Theorem 4 (reduction invariance)"
    Let ``D'`` be obtained from ``D`` by deleting versions (so
    ``\mathrm{models}(D') \subseteq \mathrm{models}(D)``), and suppose the
    layered solution ``\mathrm{Lay}(D)`` survives in ``D'`` — every version
    it uses is kept. Then
    ``\mathrm{Lay}(D') = \mathrm{Lay}(D)``.

    *Proof.* By induction along the trace, both runs make identical pins.
    At each step the choice for ``q`` is the best rank feasible given the
    (identical) pins. Deleting models can only make ranks *harder* to
    achieve, so the ``D'``-best is no better than the ``D``-best. And the
    ``D``-best is witnessed *in* ``D'`` by ``\mathrm{Lay}(D)`` itself: it
    extends all current pins (its values on pinned packages are exactly
    the pins, by induction), assigns ``q`` its ``D``-best rank (that is
    what the raw trace pinned), and survives in ``D'`` by hypothesis. So
    the choices agree; identical pins induce identical layers. The
    satisfiable/unsatisfiable branch also agrees, since
    ``\mathrm{Lay}(D)`` witnesses satisfiability in both. ∎

This is the "decoupling" result: the layered semantics is defined once, on
the raw problem, and *any* reduction whatsoever is licensed by exhibiting a
single surviving solution — the answer. Two consequences:

- **Certification.** Answer preservation can be checked *post hoc*,
  cheaply: resolve the filtered problem, then verify each pin against the
  raw problem — feasibility is free (the answer is itself a raw model), and
  optimality of each pin is one unsatisfiable solve ("no raw model with
  these pins and a better version of ``q``") per pinned package. A resolver
  can offer machine-checked raw-canonical answers for the cost of building
  the raw instance plus one probe per package.
- **Headroom.** Future reductions can be far more aggressive than the
  current filter — anything that keeps the answer is safe, and the
  certification above detects violations per-resolve rather than trusting
  the reduction globally.

The deep reason a theorem this convenient is available is that every SAT
query the layered trace makes is *existential* ("is a better version
feasible?") and is witnessed by the final answer, which survives the
reduction. Universal queries — "is this package in *every* solution?" — do
not enjoy this: deleting models changes their truth value. That asymmetry
is the fundamental obstruction explored on the [next page](front.md).

## Theorem 5: the filter preserves the answer

It remains to show the actual filter satisfies Theorem 4's hypothesis.

!!! note "Theorem 5 (reach keeps the layered answer)"
    Every version chosen by the raw layered trace is kept by reachability
    filtering.

    *Proof.* By induction along the trace; assume all pins so far are kept
    versions. Consider the step for package ``q``, and suppose first ``q``
    is unsaturated. We claim ``q``'s last kept version ``q@r(q)`` is
    feasible given the pins, so the best feasible rank — the trace's choice
    — is at most ``r(q)``, i.e. a kept version. To see feasibility,
    construct a model: take the pins (jointly feasible by the trace's
    invariant, and kept by induction), add ``q@r(q)``, and complete with
    the last-kept dependency closure of everything present. By invariant
    (1) of Theorem 2, ``q@r(q)`` and every last-kept version conflict with
    nothing kept — in particular not with the pins or each other; by
    invariant (2) the completion never needs a saturated package; pins'
    dependencies are satisfied by pinned or last-kept packages. The result
    is a model containing the pins and ``q@r(q)``. If instead ``q`` is
    saturated, all its versions are kept and there is nothing to prove. ∎

!!! note "Theorem 6 (redundancy keeps the layered answer)"
    The layered trace never chooses a version deleted by redundancy
    elimination.

    *Proof.* Suppose the trace chose ``q@j`` and ``q@j`` is dominated by a
    kept ``q@i`` with ``i < j``. Any model witnessing ``q@j``'s
    feasibility given the pins converts to one witnessing ``q@i``'s: swap
    the version; ``q@i``'s dependencies are a subset of ``q@j``'s (already
    present) and its conflicts a subset (already avoided). So ``q@i`` was
    feasible and better, contradicting that the trace chose ``j``.
    Domination chains terminate at a kept version, and the argument
    applies to whichever kept dominator terminates the chain. ∎

Composing: the raw answer survives reachability (Theorem 5), so by
Theorem 4 the reach-filtered problem has the same layered answer; that
answer survives redundancy elimination (Theorem 6 applied on the
reach-filtered problem), so by Theorem 4 again the fully filtered problem
has the same layered answer. **The implementation's answer is the raw
layered solution.**

!!! warning "Formalization boundary"
    Theorems 2, 3, 5, and 6 are proved against the rule set as stated
    above, which is a faithful but hand-made reading of `find_reachable`
    and `mark_necessary!`; the correspondence of code to rules is checked
    by review and by the test suite, not mechanically. The package's tests
    compare the implementation against a brute-force reference that
    implements the layered trace *on raw data* across hundreds of
    thousands of random instances — a direct, ongoing empirical check of
    Theorems 4–6 — and, whenever the implementation reports
    unsatisfiable, verify by enumeration that no raw solution covers the
    requirements (Theorem 3).
