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

Real problems are shrunk before solving by two passes
(`filter_pkg_info!`, run per resolve from `prepare_pkg_info`):

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

## Iterating the filter

It is natural to ask whether a single reachability pass followed by a
single redundancy pass leaves anything on the table — whether iterating
the two would prune more. The answer illuminates what each pass can and
cannot see.

Reachability alone is a closure operation:

!!! note "Proposition (reachability is idempotent)"
    Re-running reachability on its own output keeps everything.

    *Proof.* The reach set is a least fixpoint, built by a derivation
    whose steps mention only versions that end up kept: requirements'
    first versions, first versions of kept dependencies, and
    conflict/saturation extensions among reachable versions. Reachability
    deletes only version *suffixes* and whole packages, none of which
    appear in any derivation step, so the entire derivation replays on
    the filtered problem and re-derives every kept version. ∎

The composite filter, however, is **not** idempotent: redundancy
elimination can invalidate reachability's reasons. A minimal example:

!!! note "Example (redundancy changes reachability)"
    Take requirements ``\{p, q, r\}`` with two versions of each of ``p``,
    ``q``, ``r``, one version of ``s``, the single dependency
    ``p@1 \to \{s\}`` (``p@2`` depends on nothing), and conflicts
    ``q@1 \otimes r@1``, ``q@2 \otimes r@1``, ``q@2 \otimes p@1``.

    Reachability: ``q@1 \otimes r@1`` extends both packages, making
    ``q@2`` reachable; then ``q@2 \otimes p@1`` extends ``p``, so ``p@2``
    is kept — that conflict is ``p@2``'s *only* justification.
    Redundancy: ``q@2``'s constraints strictly contain ``q@1``'s, so
    ``q@2`` is deleted. But the now-orphaned ``p@2`` survives: domination
    by ``p@1`` fails because ``p@1``'s dependency on ``s`` has no
    counterpart in ``p@2``. Re-running reachability on the filtered
    problem would drop ``p@2``, whose reason to exist is gone.

The general shape is: a deleted version was another package's only
reason to degrade, and the orphaned version escapes the domination test.
Two things can shield an orphan. One is a dependency-set difference
between it and its better sibling, as above — version-varying dependency
sets again, the same dividing line as the
[fragment lemma](front.md#Why-the-filter-cannot-serve-the-front).
The other is saturation: saturation pushes are dependency-driven and
invisible to the domination test, and deleting versions can unsaturate a
package on a re-run, vaporizing the pushes that justified its
dependents' worse versions.

(A note on the domination test: it ignores conflicts whose partner
version was itself dropped by reachability — such conflicts constrain no
model over the kept versions, so ignoring them is licensed by exactly
the swap argument of Theorem 6. Counting them would spuriously block
domination and orphan far more versions than the two genuine mechanisms
above.)

None of this threatens correctness. By Theorem 4, *any* number of
additional passes preserves the layered answer, so iterating the filter
to a mutual fixpoint is perfectly legal — and since the passes are
cheap, that is what the implementation does: redundancy elimination
first (it needs no reachability marks, and on registry-scale problems
it does most of the shrinking), then reachability and redundancy
alternating until neither deletes anything. The returns diminish fast:
a further round prunes anything at all in only about 0.007% of
deliberately conflict-dense random instances with version-independent
dependencies and about 1% with version-varying ones, and what it prunes
is marginal. Essentially all of the reduction is captured by one
reachability pass and one redundancy pass; the remaining rounds just
collect the stragglers because they cost next to nothing.

## Preprocessing for caching

Everything above concerns one resolve. When many resolves share a
registry state, the question becomes which artifacts can be computed
once and reused — and the answer is organized by what each artifact
depends on: the registry content, the version ordering, and the
requirements. Three results pin down the tiers.

!!! note "Proposition (filtering for a larger requirement set)"
    Let ``\mathrm{reqs} \subseteq \mathrm{reqs}'``. The filter run with
    requirements ``\mathrm{reqs}'`` keeps every version of
    ``\mathrm{Lay}(D, \mathrm{reqs})``, so by Theorem 4
    ``\mathrm{Lay}(F_{\mathrm{reqs}'}(D), \mathrm{reqs}) =
    \mathrm{Lay}(D, \mathrm{reqs})``.

    *Proof.* Reachability is monotone in its seed set on a fixed
    problem — every derivation from ``\mathrm{reqs}`` is a derivation
    from ``\mathrm{reqs}'`` — so the ``\mathrm{reqs}'``-reach keeps a
    superset of the ``\mathrm{reqs}``-reach, which keeps the answer
    (Theorem 5). Redundancy elimination is requirement-agnostic
    outright: Theorem 6's swap argument converts any witness model
    containing a dominated version into one containing its dominator
    without touching coverage, so on whatever problem it runs, the
    layered trace of *every* requirement set avoids the versions it
    deletes. Each round of the iterated filter therefore preserves the
    answer, and Theorem 4 applies after each. ∎

!!! note "Corollary (whole-registry prefiltering)"
    Filtering once with *every* package required yields a reduced
    problem on which any requirement set resolves to its raw answer,
    with satisfiability verdicts preserved (unsatisfiability trivially,
    since models only shrink; satisfiability because the answer
    survives). Filtering is, in this sense, a monotone function of the
    requirements, and the all-requirements filter is its safe upper
    envelope.

The prefiltered problem still depends on the version *ordering* — both
reachability's prefixes and redundancy's better-dominates-worse are
stated in terms of it. How much preprocessing can be shared across
orderings has an exact answer:

!!! note "Proposition (order-free deletion is exactly model-freeness)"
    A version can be deleted without changing the layered answer for
    every version ordering and every requirement set iff it appears in
    no model.

    *Proof.* If a version appears in no model, deleting it leaves the
    model set — and hence every trace — untouched. Conversely, if
    ``p@v`` appears in a model ``M``, order ``p``'s versions with ``v``
    first and require ``p``: the trace pins ``p`` at the best feasible
    rank, and ``M`` witnesses rank 1, so the answer contains ``p@v``
    and deleting it changes the answer. ∎

So the strongest ordering-independent filter is exactly the removal of
model-free versions (of which the arc-consistency test — a version one
of whose dependencies has no compatible version left — is a cheap
sound approximation, and per-version SAT probes are the complete one).
One more piece of ordering-independent structure exists:

!!! note "Lemma (interchangeable versions)"
    Call two versions of a package *interchangeable* when they have the
    same dependency set and every compatibility entry in the problem —
    their own and every other package's — constrains them equally
    (equivalently: their constraint rows are identical). Interchangeable
    versions can be collapsed to whichever of them the active ordering
    ranks best, preserving that ordering's layered answers: the
    collapsed versions are dominated by the representative (identical
    constraints are in particular a subset), so Theorem 6's argument
    applies on the raw problem and Theorem 4 finishes. The classes are
    ordering-independent; only the representative choice is not.

    Beware that a package's *own* dependencies and compat matching is
    not sufficient: another package's compatibility bound can admit one
    of the two versions and not the other, and then they are genuinely
    distinguishable — solutions needing the excluded pairing must pick
    the right one. The equal-incoming-constraints half of the
    definition is what rules this out.

Together these give a three-tier cache: per registry state (parsed
data, the model-free set, interchangeability classes), per ordering
(the whole-registry prefilter), and per resolve (requirement-specific
filtering, which is fast and not worth caching). A non-standard
ordering — say, preferring versions closest to an installed manifest —
invalidates only the middle tier, which is orders of magnitude cheaper
to rebuild than the first.

The implementation draws the line between the first tier and the rest
at `pkg_info`, whose output — conflict matrices indexed by the classes
(`version_classes`), with the versions that cannot be installed
whatever else is chosen already dropped — is exactly the registry-only
artifact, and it merges the second and third tiers into one per-resolve
pass, `prepare_pkg_info`, since the middle tier is cheap enough not to
be worth storing.

The artifact is indexed by those classes throughout: a `PkgInfo`'s
conflict matrix has one row per class and one column per partner class.
Building it takes one row per version first, since that is what the
partition is computed from, and then merges them (`collapse_classes!`).
Nothing is lost in the merge: members of a class have equal rows by
definition, and — since conflicts are recorded symmetrically — members
of a partner class index equal columns, so the merged matrix states
exactly the same conflicts.

What the per-resolve pass then does is choose, per class, which member
it stands for. The user's constraints act only by forbidding members
(`exclusion_masks`): a class keeps whichever of its admissible members
the query's ordering ranks best, and competes at that member's rank,
which is why the classes have to be laid out in the order those
representatives put them in (`class_ranking`). A class the query
admits no member of is *deactivated* — nothing can choose it — and
that, not deletion, is the whole of a constraint's effect. Then the
result is filtered (`filter_pkg_info!`).

Deactivation is deliberately not deletion. A deactivated class stays in
the matrices: reachability continues its prefix past one, since a class
that cannot be selected is not somewhere a package can be installed,
but still follows its dependencies, and redundancy elimination neither
deletes one nor lets one license a deletion. What that buys is that the
universe a later question about relaxing a constraint would need is
still present — and, one level down, that such a question is even well
posed, since a constraint cannot make two versions' rows identical: it
does not touch the rows, and versions whose rows *are* identical are
one class sharing one column, with no pair of separately deletable
objects to fall between.

The ordering is a `resolve` parameter rather than part of the problem,
which is what makes one artifact serve every query: it selects among
the valid solutions instead of changing which solutions are valid. What
*does* change the model set — the user's compat bounds and pins, and
the admission of prerelease versions — is a `Problem`, and enters as
constraints on the artifact's versions rather than as deletions from it.
So there is one T1 artifact per registry state, full stop: every
version the registry state offers is in it, in canonical order and in
its class, and each query says which of them it will accept and in what
order it prefers them.

"Every version the registry state offers" is meant strictly, and it is
where yanked versions come in. Yankedness is not a query fact but part
of what the registry state says — a yanked version is one the registry
has struck, so it is not on offer — and the default artifact is built
with those versions filtered out at the provider (see
`bin/Registries.jl`). That is a deletion, not a constraint, and it is
the right shape: a resolve then cannot produce a yanked version or
suggest one as an alternative, with nothing left to enforce per query.
The cost is that `--allow-yanked`, the query that wants them back, is
the one query the baked artifact does not serve: it asks for a universe
the registry state does not describe, and so rebuilds T1. That is rare
and cheap next to the alternative, which is carrying struck versions
through every layer and every diagnostic in order to forbid them again
at the end.

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
