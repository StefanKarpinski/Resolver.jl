# Diagnosing unsatisfiability

When a resolve fails, the resolver explains why: which requirements clash,
which constraints are implicated, and what relaxations would make the
problem solvable. Every one of those questions is answered by *more SAT
queries on the very instance that just failed* — the prepared production
instance, with the user's constraints turned back into assumptions.

That is a strong claim, because the instance was built for one purpose (to
resolve one problem) and is now being asked a family of different
questions. This page works out exactly how far it is allowed. The answer is
not "always": preparation is stable under relaxations of a certain
*granularity* and demonstrably unstable below it, and that boundary is what
determines the units a diagnosis is allowed to blame.

Three results:

1. **Requirement-subset queries** are already licensed by the existing
   [requirement-monotonicity proposition](layered.md#Preprocessing-for-caching).
2. **Relaxation stability** (new here) covers switching off *constraint
   groups* — but only column-closed ones. A companion proposition shows the
   finer granularity really does fail, with a four-version counterexample.
3. **Closure exactness** licenses answering a cluster's questions on a
   sub-instance restricted to that cluster's dependency closure.

## Relaxations

Fix a problem ``D`` as in the [setting](setting.md), with requirements
``\mathrm{reqs}``. Designate a finite family ``G`` of **relaxable conflict
groups**: each ``g \in G`` is a set of conflicting version pairs. Conflicts
in no group are **permanent**.

A **relaxation** is a pair ``\sigma = (\mathrm{reqs}_0, G_0)`` with
``\mathrm{reqs}_0 \subseteq \mathrm{reqs}`` and ``G_0 \subseteq G``. The
**relaxed problem** ``D|\sigma`` has the same packages, versions and
dependency sets as ``D``; its requirements are ``\mathrm{reqs}_0``; and its
conflict relation is the permanent conflicts together with
``\bigcup_{g \in G_0} g``. Groups may overlap: a pair asserted by two groups
survives until both are switched off, which is exactly what the
selector-guarded encoding does.

Two families of group matter in practice:

- **User-constraint groups.** A compat bound or a pin is read as a *virtual
  package* with one always-present version whose conflict rows are the
  versions the constraint forbids (see `Problem`). Switching one off is
  "what if the user relaxed this?".
- **Third-party bound groups.** A registry compatibility declaration
  contributes a rectangle of conflicting pairs. Switching one off is "what
  if upstream released a version that loosened this?".

Write ``\mathrm{Lay}(D, \sigma)`` for the layered solution of ``D|\sigma``
(undefined when ``D|\sigma`` is unsatisfiable), and ``F`` for **preparation**
— everything between the raw universe and the SAT instance, run **once**, with
the full requirement set and every group active. Following
[`pkg_info`](@ref Resolver.pkg_info) and
[`prepare_pkg_info`](@ref Resolver.prepare_pkg_info), that is four steps:

1. **arc consistency**, deleting versions one of whose dependencies has no
   versions at all;
2. the **interchangeability collapse**, keeping one representative of each
   class of versions nothing can tell apart (refined by ``\sigma``'s own
   exclusions — see below);
3. **reachability**, keeping a prefix of each package's versions; and
4. **redundancy elimination**, deleting dominated versions,

with 3 and 4 alternating to a fixpoint. Each deletes versions and nothing
else, which is what lets them be treated uniformly below.

## Requirement subsets

The existing
[proposition on filtering for a larger requirement set](layered.md#Preprocessing-for-caching)
and its corollary already say that filtering seeded with ``\mathrm{reqs}``
preserves the layered answer — and hence the satisfiability verdict — of
every ``\mathrm{reqs}_0 \subseteq \mathrm{reqs}``. That licenses the
requirement-level half of a diagnosis: clustering the requirements into
independent conflicts, and testing whether dropping some of them helps.
Nothing new is needed.

What it does not cover is the other half — lifting a compat bound, a pin or
a registry bound. Those change the *conflict relation*, which is what the
preparation's reachability prefixes and domination tests are computed
against — and, in the case of the collapse, what its class representatives are
chosen from.

## Reachability is stable; redundancy is the delicate half

!!! note "Lemma R (reachability is monotone in seeds and conflicts)"
    Let ``D`` and ``D'`` have the same packages, versions and dependency
    sets, with conflict relation ``\otimes' \subseteq \otimes`` and
    requirement seeds ``\mathrm{reqs}' \subseteq \mathrm{reqs}``. Then
    ``r'(p) \le r(p)`` for every package ``p``, where ``r`` is the reach
    prefix map and saturation counts as ``+\infty``.

    *Proof.* The reach set is the least fixpoint of the derivation system of
    [Filtering](layered.md#Filtering): seed requirements' first versions;
    first versions of kept dependencies; conflict-driven degradation of both
    sides of a conflict between reachable versions; and the saturation rule
    pushing dependents of a saturated package. Every rule instance of the
    primed system is a rule instance of the unprimed one — the seeds are a
    subset, the conflict premises are a subset, and the dependency premises
    are identical. So every primed derivation replays verbatim, and the
    primed least fixpoint is contained in the unprimed one. ∎

The exclusion rule for user constraints ("a reachable excluded version can
never be a resting point") is the ordinary conflict rule applied to the
virtual package's always-present version, so the same argument covers it.
Lemma R asks nothing of the *shape* of the removed conflicts: it holds for
**arbitrary** group removal.

Redundancy elimination does not. It deletes ``q@j`` when some better kept
``q@i`` has a *subset* of its constraints, and a relaxation that clears a
constraint of ``q@i`` without clearing the corresponding constraint of
``q@j`` destroys the inclusion. Whether that can happen depends entirely on
the shape of the group:

!!! note "Definition (column-closed group)"
    Redundancy compares two versions of a package column by column of that
    package's constraint matrix: one column per dependency, one per
    ``(\text{partner package}, \text{partner version})`` pair, and one for
    the user-constraint virtual package. A group ``g`` is **column-closed**
    if switching it off clears, in every package's matrix, a set of *whole
    columns* — never part of one.

Two natural granularities are column-closed, and they are the coarsest
sensible ones:

- **Per package, for user constraints.** All of a package's user
  constraints together — its compat bound *and* its pin — are exactly its
  exclusion column.
- **Per interacting pair, for third-party bounds.** All conflicts between
  packages ``p`` and ``q`` are exactly ``q``'s whole column block in ``p``'s
  matrix and ``p``'s whole block in ``q``'s.

!!! note "Theorem D1 (relaxation stability)"
    Suppose every group in ``G`` is column-closed. Then for every relaxation
    ``\sigma``,

    ```math
    \mathrm{Lay}(F(D), \sigma) \;=\; \mathrm{Lay}(D, \sigma),
    ```

    both sides being undefined together. In particular ``F(D)|\sigma`` is
    satisfiable iff ``D|\sigma`` is.

    *Proof.* If ``D|\sigma`` is unsatisfiable, so is ``F(D)|\sigma``: its
    models are a subset. So assume ``D|\sigma`` is satisfiable. We show each
    pass of preparation — each of which deletes versions from ``D``, and
    hence induces the same deletion on ``D|\sigma`` — keeps every version of
    the current ``\sigma``-answer, and conclude with
    [Theorem 4](layered.md#Theorem-4:-reduction-invariance-—-the-optimization-license)
    applied to ``D|\sigma`` after each pass.

    *Arc consistency.* The pass deactivates a version exactly when one
    of its dependencies has no version at all, iterated to a fixpoint. The
    criterion mentions no conflicts and no requirements, so it is
    ``\sigma``-independent, and such a version appears in no model of
    ``D|\sigma`` either — dependency edges must be satisfied whatever the
    conflict relation.

    *The collapse.* The pass keeps, of each interchangeability class, the
    member the active order ranks best — where the classes are the registry's
    row-equality classes
    ([`version_classes`](@ref Resolver.version_classes)) *refined by the
    exclusion mask*, so that a class splits into its forbidden and its
    admissible members and each half keeps a representative
    ([`class_representatives`](@ref Resolver.class_representatives)).

    Two things have to hold. First, the deletion is answer-preserving on the
    unrelaxed problem: class members have identical constraint rows, hence
    dominate each other, so Theorem 6 applies — this is the
    [interchangeability lemma](layered.md#Preprocessing-for-caching), and it is
    why the pass is legitimate at all.

    Second, it survives ``\sigma``. Let ``L`` be a class of ``D|\sigma``'s
    partition — the partition of the *relaxed* problem, refined by ``\sigma``'s
    exclusion mask. Two claims:

    * ``L`` is a union of classes of the unrelaxed refined partition. The
      registry-level classes coarsen under ``\sigma`` (removing conflict
      columns can only merge row-equality classes, never split them), and the
      exclusion refinement coarsens too, because ``\sigma`` clears the mask
      for whole packages at a time — the mask is one column, and column-closure
      is exactly the hypothesis that ``\sigma`` takes all of it or none of it.
    * Hence the best member of ``L`` in the active order is the best member of
      one of those finer classes, and is therefore kept.

    So every ``\sigma``-class retains its best member, and the members that
    are gone are dominated by it *under ``\sigma``* (they are in the same
    ``\sigma``-class, so their rows agree there too). Theorem 6 applied to
    ``D|\sigma`` says the ``\sigma``-layered trace never chooses them.

    Note where column-closure did the work: were ``\sigma`` allowed to clear
    *part* of a package's exclusion column, the refinement could split a
    ``\sigma``-class between a forbidden and an admissible half whose
    representative is the wrong one, and the best ``\sigma``-admissible member
    could be gone. That is the collapse's version of Proposition D1′ below, and
    it has the same cause.

    *Reachability.* Let ``r`` be the prefix map the pass computes with the
    full seed set and every group active, and ``r_\sigma`` the map it would
    compute on ``D|\sigma``. By Lemma R, ``r_\sigma \le r``. By
    [Theorem 5](layered.md#Theorem-5:-the-filter-preserves-the-answer)
    applied to ``D|\sigma``, every version of the ``\sigma``-answer lies
    within ``r_\sigma``, hence within ``r``.

    *Redundancy elimination.* The pass deletes ``q@j`` when some kept
    ``q@i`` with ``i < j`` satisfies ``X[i,c] \Rightarrow X[j,c]`` for every
    active column ``c``. Because every group is column-closed, passing to
    ``\sigma`` clears a set ``C_\sigma`` of whole columns — *the same*
    columns for ``i`` and for ``j``. Restricting an implication that holds
    columnwise to a subset of the columns leaves it holding, so ``q@j`` is
    dominated in ``D|\sigma`` too, and
    [Theorem 6](layered.md#Theorem-5:-the-filter-preserves-the-answer)
    applied to ``D|\sigma`` says the ``\sigma``-layered trace never chooses
    it. The test's dropped-partner exemption — conflicts whose partner
    version reachability already deleted are ignored — is
    ``\sigma``-independent, since which versions were deleted is decided by
    the full-constraint filter and not by ``\sigma``; it drops the same
    columns on both sides.

    *Iteration.* Each pass preserves the ``\sigma``-answer of the problem it
    runs on, so Theorem 4 applies after each and the answers chain — through
    the collapse as much as through the filter proper, the collapse being
    handed to the filter as marks and materialized by its first deletion round
    rather than as a separate rebuild — exactly as they do in the unrelaxed
    composition at the end of
    [Theorem 6](layered.md#Theorem-5:-the-filter-preserves-the-answer). ∎

### Below column-closure it is false

It is tempting to relax at the granularity of the *sources* a report would
like to name: this compat bound as against this pin, this one compat
declaration as against another. Neither is column-closed — a constraint
source clears part of a package's exclusion column, and a single compat
declaration clears part of a partner's block — and neither is safe.

!!! warning "Proposition D1′ (finer relaxations are not stable)"
    Relaxing user constraints per *source*, or third-party bounds per
    *declaration*, does not preserve satisfiability verdicts. One package
    with two versions suffices for the first.

    *Counterexample (source granularity).* A single package ``p`` with two
    versions and no dependencies. The user's compat bound on ``p`` admits
    only ``p@2``; the user's pin holds ``p`` at ``p@1``. Requirement
    ``\{p\}``.

    Together the two sources forbid everything, so ``p``'s exclusion column
    is *full*: ``p@1`` and ``p@2`` have identical constraint rows, ``p@1``
    dominates, and redundancy deletes ``p@2``. The filtered instance carries
    only ``p@1``.

    Now take ``\sigma``: drop the pin, keep the compat. The raw problem is
    satisfiable, at ``p@2``. The filtered instance no longer has ``p@2`` and
    reports unsatisfiable. ∎

    *Counterexample (declaration granularity).* Packages ``p`` and ``q``,
    two versions each, ``p@2`` depending on ``q``; ``q@1`` declares compat
    on ``p`` admitting only ``p@1``, and ``q@2`` declares compat on ``p``
    admitting nothing. Those are two declarations, hence two groups.
    Dropping only ``q@2``'s leaves ``q@1``'s in force, clearing part of
    ``p``'s column block in ``q``'s matrix and not the rest — and the same
    collapse follows.

    The general shape is what the proof needs: ``\sigma`` cleared a cell of
    the *dominating* version without clearing the corresponding cell of the
    dominated one, and the domination evaporated.

    Neither is a rare corner. Sweeping small random instances over **every**
    relaxation:

    | preparation | group granularity | σ swept | wrong |
    |:--|:--|--:|--:|
    | full | source + declaration | 212,664 | 1.358% |
    | full | package + declaration | 916,020 | 0.401% |
    | full, collapsed | package + pair | 629,740 | **0** |
    | full, uncollapsed | package + pair | 629,740 | **0** |
    | full, collapsed, reordered | package + pair | 629,740 | **0** |
    | reach only | source + declaration | 3,177,516 | **0** |

    The operative rows are the middle three and the last: column-closure is
    what makes the difference; the collapse neither helps nor hurts, in the
    canonical version order or a reversed one; and redundancy elimination is
    the entire source of the instability.

The two results localize the damage precisely: by Lemma R, reachability is
stable at *every* granularity, so the entire instability is redundancy
elimination — the last row of the table above is the empirical face of
Lemma R. A diagnosis that needed source-level relaxations could therefore
buy them back, at the cost of a second, redundancy-free filtered universe.
The implementation does not: it stays at column-closed granularity and
keeps one instance.

### What transfers, and at what granularity

!!! note "Corollary D2 (cores, corrections and enumeration transfer)"
    Let ``\Sigma`` be the lattice of relaxations over a column-closed family
    ``G``, and ``\mathrm{sat} : \Sigma \to \{0,1\}`` the satisfiability
    verdict. Theorem D1 says the verdict functions of ``D`` and ``F(D)`` are
    the same function on ``\Sigma``. Therefore, computed on ``F(D)``: a
    minimal unsatisfiable relaxation is one of ``D``; a minimal correction
    set is one of ``D``; and an enumeration of either is complete for ``D``
    exactly when it is complete for ``F(D)``.

    *Proof.* All three notions are defined purely in terms of
    ``\mathrm{sat}`` on ``\Sigma``, and the two verdict functions
    coincide. ∎

Because Theorem D1 preserves the layered *answer* and not merely the
verdict, the solution reported alongside a fix — obtained by re-resolving
the filtered instance with the fix's groups switched off — is the answer the
raw problem would give under that fix, not an artifact of filtering.

The price of column-closure is paid in the report's vocabulary, and it is
small:

- A fix relaxes **all** of a package's constraints at once — compat bound,
  pin and admission knobs together. For the overwhelming majority of packages,
  which carry one source, that is the same thing as naming it. When a package
  carries several, the fix names each of them that still forbids a version the
  instance has; the ones that forbid nothing there are exactly the ones the
  relaxation does not depend on, so leaving them unnamed keeps the advice both
  minimal and honest.
- A bound-level conflict is blamed on an **incompatible pair of packages**
  rather than on one compat declaration, and an upstream suggestion names
  both sides. That was going to be true anyway: package info records
  compatibility symmetrically, so which side declared a bound is not
  recoverable from the filtered instance (see the boundary note below).

!!! warning "Relaxations only"
    Every diagnostic query must be a relaxation over a column-closed family.
    Assuming a subset of the requirements and a subset of the selector
    variables is one; adding a *clause* is not.

    This is why fix enumeration is done MARCO-style with assumption-subset
    queries — force-keeping a previously dropped group and re-running the
    greedy pass — rather than by adding a blocking clause
    ``\bigvee_{g \in \mathrm{MCS}} s_g`` over the selectors after each fix.
    A blocking clause is an added constraint; the resulting instance is
    outside the ``\sigma``-family and Theorem D1 licenses nothing about it.

## Closure exactness

Bound-level questions are asked per cluster, on a sub-instance covering only
the packages that cluster can reach. That restriction is exact.

!!! note "Lemma D3 (closure exactness)"
    Fix requirements ``R_0`` and let ``C`` be the least set containing
    ``R_0`` and closed under

    ```math
    p \in C \;\Longrightarrow\; \bigcup_{v} \mathrm{deps}(p, v) \subseteq C
    ```

    — the union over **all** of ``p``'s versions, not just a chosen one. Let
    ``D[C]`` be ``D`` restricted to ``C``: packages outside ``C`` deleted,
    dependencies and conflicts among ``C`` unchanged. Then ``D`` has a
    solution covering ``R_0`` iff ``D[C]`` does, and the tight solutions
    covering ``R_0`` are literally the same set.

    *Proof.* (⇒) Let ``s`` be a solution of ``D`` covering ``R_0`` and
    ``\pi(s)`` its pruning, a tight solution with support ``C(s)``
    ([pruning lemma](setting.md#Tightness-and-pruning)). Then
    ``C(s) \subseteq C``: ``R_0 \subseteq C``, and if ``p \in C(s)`` then
    ``\mathrm{deps}(p, s(p)) \subseteq \bigcup_v \mathrm{deps}(p, v)
    \subseteq C``. So ``\pi(s)`` lives inside ``C`` and is a solution of
    ``D[C]`` — coverage, closure and conflict-freedom all inherited.

    (⇐) Let ``t`` be a solution of ``D[C]`` covering ``R_0``. Read as a
    partial assignment on ``D`` that assigns nothing outside ``C``, it is
    still a solution: coverage is unchanged; dependency closure holds
    because every dependency of an assigned version lies in ``C`` and is
    assigned by ``t``; conflict freedom is inherited, since ``t`` assigns no
    version outside ``C`` to conflict with. ∎

    At the level of the SAT encoding the extension is even more transparent:
    every clause an outside package occurs in is guarded by one of its own
    literals — dependency clauses ``\lnot p@i \lor q``, conflict clauses
    ``\lnot p@i \lor \lnot q@j`` — so setting the whole package false
    satisfies them vacuously. Outside packages are never *forced*, because
    only dependency edges force packages and every edge out of ``C`` stays
    inside ``C``.

Composing with Theorem D1: a per-cluster sub-instance built from the
filtered universe ``F(D)``, restricted to the closure of the cluster's
requirements, and equipped with selectors for the column-closed groups it
contains, answers every relaxation query for that cluster exactly as the raw
problem would.

One honest caveat: ``C`` unions over *all* versions' dependency sets, so it
is generally larger than the support of any single solution — and for a
cluster rooted at a hub package it can be a substantial fraction of the
universe. Exactness is not the same as smallness.

## What the implementation does with this

- The production instance asserts the user-constraint selectors as unit
  clauses inside a single `sat_push` frame. One `sat_pop` turns every user
  constraint into an assumable selector, converting the failed instance into
  the diagnostic instance in place — no rebuild, no second copy of the
  universe.
- Diagnosis relaxes user constraints **per package**, never per source, and
  third-party bounds **per interacting pair**, never per declaration. That
  is Theorem D1's hypothesis, and Proposition D1′ is what happens if it is
  ignored. "Per package" covers every source at once — the compat bound, the
  pin, and every admission knob (`:prerelease`, `:yanked`, …) that touches the
  package — because they all write into the one exclusion column.
- Clustering and fix enumeration run on the converted instance with
  assumption-only queries.
- Each fix's reported solution comes from re-resolving the same instance
  with the fix's groups switched off inside a temporary frame; the answer
  half of Theorem D1 makes those versions the real ones.

!!! warning "Formalization boundary"
    Theorem D1 and Lemma R are proved against the same hand-made rule set as
    Theorems 2–6 on the [previous page](layered.md) — a faithful reading of
    `find_reachable` and `mark_necessary!`, checked against the code by
    review and by the test suite rather than mechanically. Lemma D3 is
    proved against the abstract [setting](setting.md), and its SAT-level
    restatement against the encoding described there.

    Everything on this page has a direct empirical check, in the same spirit
    as the brute-force reference-resolve tests that back Theorems 4–6. The
    suite sweeps small random instances and, for **every** relaxation
    ``\sigma`` — every subset of the requirements crossed with every subset
    of the constraint groups, not a sample — asserts that resolving the
    filtered instance under ``\sigma`` gives exactly what a brute-force
    reference resolve of the raw data with those groups deleted gives,
    verdict and answer alike. It asserts stability at package-and-pair
    granularity, reproduces Proposition D1′'s counterexample at source
    granularity — with the collapse on and off, and in a non-canonical version
    order — and confirms Lemma R by sweeping a reach-only universe at the finer
    granularity where full preparation fails. A separate testset pins the
    property the collapse case of the proof turns on: that refining classes by
    the exclusion mask leaves a forbidden version in the instance for every
    kind of source, so switching that source off can bring it back.

    One representational caveat, orthogonal to the theorems: package info
    records compatibility *symmetrically*, so by the time a problem reaches
    the filtered instance the resolver no longer knows which side of an
    incompatibility declared it. Bound facts in a diagnosis therefore name
    an incompatible pair of version groups rather than a directed "``X``
    requires ``Y`` at …" declaration, and upstream suggestions name both
    sides.
