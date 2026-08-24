# Relaxation-stable filtering

The previous page proves that filtering preserves the answer to *the query
it was run for*. This page proves something stronger, which is what lets a
failed resolve be explained rather than merely reported: the filtered
universe answers every **relaxation** of that query too.

A relaxation is the same query with some of its demands relaxed — a
requirement dropped, a compat entry ignored, a pin released. Those are
exactly the repairs a user can make to an unsatisfiable problem, so
"filtered for ``Q``, valid for every ``R \le Q``" is precisely the property
under which a diagnosis may test candidate fixes on the universe the failed
resolve already built, instead of re-filtering per fix.

The cost of the stronger property is that the filter deletes slightly less.
Nothing else changes: the passes, their order, and the fixpoint they iterate
to are the ones the previous page describes.

## What a relaxation is

Fix an artifact ``D``. A **query** ``Q`` is a requirement set together with a
set of constraints — one per package a compat entry or a pin of the user's
problem names, plus one per admission kind (prereleases, yanked versions),
which names no package in particular.

!!! note "Definition (the relaxation order)"
    ``R \le Q`` — "``R`` is a relaxation of ``Q``" — iff
    ``\mathrm{reqs}_R \subseteq \mathrm{reqs}_Q`` and ``R``'s constraints are
    a subset of ``Q``'s.

This is a lattice with ``2^{|\mathrm{reqs}_Q| + |\mathrm{cons}_Q|}``
points, ``Q`` at the top and the empty query at the bottom. Note ``Q \le Q``:
every statement below about "every relaxation" includes the query itself,
which is what makes the results here strengthen rather than replace the
previous page's.

Write ``\mathrm{Lay}(D, R)`` for the layered answer to ``R`` on ``D``.

## Two keys per class, and a third

A class competes at the rank of the best member the query admits. This is
the whole source of the difficulty: a compat bound that forbids a class's
best member does not remove the class, it *moves* it, possibly behind a class
it would otherwise outrank. Relaxing the bound moves it back.

Three quantities per class ``c``, all of them version ranks in the query's
version order:

- ``\mathrm{adm}_Q(c)`` — the members of ``c`` that ``Q`` admits. ``c`` is
  **deactivated** under ``Q`` when this is empty.
- ``\mathrm{key}_Q(c) = \min \mathrm{adm}_Q(c)``, and ``\infty`` when ``c`` is
  deactivated. This is what the walk minimises over.
- ``\mathrm{key}_\emptyset(c)`` — the rank of ``c``'s first member, ignoring
  every constraint. Its key under the empty query.
- ``\kappa_Q(c) = \mathrm{key}_Q(c)`` for an active class, and
  ``\mathrm{key}_\emptyset(c)`` for a deactivated one. This is the **layout
  key**: classes are laid out in ascending ``\kappa_Q``.

The distinction between ``\mathrm{key}_Q`` and ``\kappa_Q`` matters and is
easy to lose. A deactivated class has no admissible member, so
``\mathrm{key}_Q`` is ``\infty`` — but it is not laid out last. It is laid
out where its best member would put it, because that is where a relaxation
reviving it would put it. `class_ranking_keys` computes exactly this: the
first loop fills each class's key from its best admissible member, and the
second fills the classes that got none from their best member, period.

!!! note "Fact 0 (the keys are ordered)"
    ``\mathrm{key}_\emptyset(c) \le \kappa_Q(c) \le \mathrm{key}_Q(c)`` for
    every class ``c``.

    *Proof.* If ``c`` is active, ``\mathrm{adm}_Q(c)`` is a subset of ``c``'s
    members, so its minimum is no smaller:
    ``\mathrm{key}_\emptyset(c) \le \mathrm{key}_Q(c)``, and
    ``\kappa_Q(c) = \mathrm{key}_Q(c)``. If ``c`` is deactivated,
    ``\mathrm{key}_\emptyset(c) = \kappa_Q(c)`` by definition and
    ``\mathrm{key}_Q(c) = \infty``. ∎

## The one order fact a relaxation cannot reverse

!!! note "Lemma 1 (monotonicity)"
    For every class ``c`` and every ``R \le Q``:
    ``\mathrm{key}_\emptyset(c) \le \mathrm{key}_R(c) \le \mathrm{key}_Q(c)``.

    *Proof.* Adding constraints only removes admissible members, so
    ``\mathrm{adm}_Q(c) \subseteq \mathrm{adm}_R(c) \subseteq
    \mathrm{adm}_\emptyset(c)``. A minimum over a superset is no larger
    (with ``\min \emptyset = \infty``). ∎

!!! note "Lemma 2 (definite dominance)"
    If ``\mathrm{key}_Q(c') < \mathrm{key}_\emptyset(c)`` then
    ``\mathrm{key}_R(c') < \mathrm{key}_R(c)`` for every ``R \le Q``.

    *Proof.* ``\mathrm{key}_R(c') \le \mathrm{key}_Q(c') <
    \mathrm{key}_\emptyset(c) \le \mathrm{key}_R(c)`` by Lemma 1. ∎

Two consequences the rules below lean on. The hypothesis forces
``\mathrm{key}_Q(c') < \infty``, so ``c'`` is admissible under ``Q`` and
hence under every ``R \le Q`` — **a deactivated class can never license
anything**. And this is the *only* order fact stable under relaxation:
"``c'`` precedes ``c`` in ``Q``'s layout" is not enough, because ``R`` can
reverse it.

!!! note "Definition (possibly-best)"
    For a set ``S`` of classes of one package, ``c \in S`` is
    **possibly-best in ``S``** iff no ``c' \in S`` has
    ``\mathrm{key}_Q(c') < \mathrm{key}_\emptyset(c)``. Write ``PB(S)``.

Since ``\mathrm{key}_\emptyset(c) \le \mathrm{key}_Q(c)``, no class
disqualifies itself, and the test collapses to one scalar: with
``t = \min \{\mathrm{key}_Q(c') : c' \in S\}``,

```math
PB(S) = \{\, c \in S : \mathrm{key}_\emptyset(c) \le t \,\}.
```

!!! note "Lemma 3 (``PB`` covers every relaxation's choice)"
    For every ``R \le Q`` and every non-empty ``S``:
    ``\arg\min_{c \in S} \mathrm{key}_R(c) \in PB(S)``.

    *Proof.* Let ``c^*`` minimise ``\mathrm{key}_R`` over ``S`` and suppose
    ``c^* \notin PB(S)``: some ``c' \in S`` has
    ``\mathrm{key}_Q(c') < \mathrm{key}_\emptyset(c^*)``, so
    ``\mathrm{key}_R(c') < \mathrm{key}_R(c^*)`` by Lemma 2, contradicting
    minimality. ∎

Both filtering rules are the same substitution: wherever the prefix filter
said "the first class", the relaxation-stable filter says "``PB`` of what
remains", and Lemma 3 says that covers whatever any relaxation would pick.

## Rule A — redundancy

!!! note "Rule A"
    Strike class ``c_2`` of package ``p`` only if some class ``c_1`` of ``p``
    has registry rows ``\subseteq`` ``c_2``'s **and**
    ``\mathrm{key}_Q(c_1) < \mathrm{key}_\emptyset(c_2)``.

The prefix rule asked only that ``c_1`` precede ``c_2`` in ``Q``'s layout —
the fact Lemma 2's discussion shows a relaxation can reverse.

!!! note "Theorem A (a struck class is needed by no relaxation)"
    If ``c_2`` is struck by Rule A then ``\mathrm{Lay}(D, R)`` does not use
    ``c_2``, for every ``R \le Q``.

    *Proof.* Suppose the layered trace for ``R`` pins ``p`` at ``c_2``,
    witnessed by a model ``M`` extending the pins so far. Build ``M'`` from
    ``M`` by putting ``p`` at ``c_1`` instead, at ``c_1``'s best
    ``R``-admissible member — which exists, since
    ``\mathrm{key}_R(c_1) \le \mathrm{key}_Q(c_1) < \infty``.

    ``M'`` is a model: ``c_1``'s dependencies are a subset of ``c_2``'s and
    are already satisfied in ``M``; ``c_1``'s conflicts are a subset of
    ``c_2``'s, so a conflict of ``c_1`` with a member of ``M`` would be one
    of ``c_2`` too, contradicting ``M``'s validity. ``M'`` agrees with ``M``
    everywhere else, hence with the pins.

    By Lemma 2, ``\mathrm{key}_R(c_1) < \mathrm{key}_R(c_2)``, so ``M'``
    gives ``p`` a strictly better rank under ``R``. The trace pins the best
    *feasible* class, so it would have pinned ``c_1`` or better, not
    ``c_2``. Contradiction. ∎

A struck class is not thereby forgotten: `mark_necessary!` records its members
in the shadow list of ``c_1``, which is what lets a diagnosis account for them
— see [Shadows](layered.md#Shadows:-what-a-deletion-leaves-behind) on the
previous page. The record is a side table; Theorem A quantifies over the same
classes either way.

!!! warning "Side condition"
    The domination test ignores conflicts whose partner class is absent from
    the universe. That exemption is sound only when evaluated against the
    *conservatively* kept classes — Rule B's output, or all classes — never
    against a single query's reach prefix. A partner absent under Rule B is,
    by Theorem B, needed by no relaxation, so conflicts with it constrain no
    relevant model; a partner that ``Q``'s prefix alone would drop may well
    matter under ``R``. This is why the two rules are adopted together.

## Rule B — reachability

!!! note "Rule B"
    Least fixpoint of sets ``S(p)`` with ``\mathrm{push}(p) \subseteq S(p)``
    recording the classes the walk has been forced past. Writing
    ``t(p) = \min\{\mathrm{key}_Q(c) : c \notin \mathrm{push}(p)\}``, an
    unseeded package has ``S(p) = \emptyset`` and a seeded one has
    ``S(p) = \{c : \mathrm{key}_\emptyset(c) \le t(p)\}``, subject to:

    1. ``p \in \mathrm{reqs}_Q`` ⇒ ``p`` is seeded.
    2. ``c \in S(p)`` and ``c`` depends on ``q`` ⇒ ``q`` is seeded.
    3. ``c \in S(p)`` deactivated ⇒ ``c \in \mathrm{push}(p)``.
    4. ``c \in S(p)``, ``c' \in S(q)``, ``c`` conflicts with ``c'`` ⇒
       ``c \in \mathrm{push}(p)`` and ``c' \in \mathrm{push}(q)``.
    5. ``c \in S(p)`` depends on ``q`` with
       ``\mathrm{push}(q) = \mathrm{classes}(q)`` ⇒
       ``c \in \mathrm{push}(p)``.

!!! note "Theorem B (every relaxation's reachable set is kept)"
    For every ``R \le Q`` and every package ``p``, every class the prefix
    filter would keep for ``R`` lies in ``S(p)``.

    *Proof.* By induction over the derivation of ``R``'s reach fixpoint,
    maintaining the invariant

    > (∗) ``\mathrm{passed}_R(p) \subseteq \mathrm{push}(p)`` — the walk has
    > been forced past everything ``R``'s walk was.

    *Invariant.* ``R`` passes ``c`` for one of three reasons, each firing the
    matching rule: (a) ``c`` conflicts with a class reachable under ``R``,
    both in ``S`` by the induction hypothesis, and conflicts are
    query-independent, so rule 4 fires; (b) ``c`` is not selectable under
    ``R`` — and since ``\mathrm{adm}_Q(c) \subseteq \mathrm{adm}_R(c)``, not
    under ``Q`` either, so rule 3 fires; (c) ``c`` depends on a package
    saturated under ``R``, which by (∗) is saturated here too, so rule 5
    fires.

    *Seeds.* ``p \in \mathrm{reqs}_R \subseteq \mathrm{reqs}_Q``, so rule 1
    seeds ``p``, and ``R``'s first class of ``p`` minimises
    ``\mathrm{key}_R`` over all of ``p``'s classes, which Lemma 3 places in
    ``S(p)``. Dependencies are the same argument via rule 2.

    *Frontier.* Let ``c^*`` be ``R``'s next class after its passed set:
    ``c^* = \arg\min \mathrm{key}_R`` over
    ``\mathrm{classes}(p) \setminus \mathrm{passed}_R(p)``. If
    ``c^* \in \mathrm{push}(p)`` then ``c^* \in S(p)``, since pushed classes
    are kept. Otherwise ``c^*`` lies in
    ``\mathrm{classes}(p) \setminus \mathrm{push}(p)``, which by (∗) is a
    subset of ``\mathrm{classes}(p) \setminus \mathrm{passed}_R(p)``;
    ``c^*`` minimises ``\mathrm{key}_R`` over the larger set and belongs to
    the smaller, so it minimises over the smaller too, and Lemma 3 puts it in
    ``S(p)``. ∎

## The query's own answer, without a second walk

Theorem B with ``R = Q`` already says ``S(p)`` contains what ``Q``'s own
prefix walk would keep. That is worth stating concretely, because it is what
licenses the implementation to run *one* walk: the possibly-best walk is not
unioned with a prefix walk, and no prefix walk is run at all.

Number each package's classes ``1, \dots, m_p`` in layout order, so
``\kappa_Q(p, \cdot)`` is non-decreasing. The prefix walk ``W`` is the least
fixpoint of ``r(p) \in \{0, \dots, m_p+1\}`` under

1. ``p \in \mathrm{reqs}`` ⇒ ``r(p) \ge 1``;
2. ``j \le r(p)`` and ``p@j`` depends on ``q`` ⇒ ``r(q) \ge 1``;
3. ``j \le r(p)`` and ``p@j`` deactivated ⇒ ``r(p) \ge j+1``;
4. ``j \le r(p)``, ``k \le r(q)``, ``p@j`` conflicts with ``q@k`` ⇒
   ``r(p) \ge j+1`` and ``r(q) \ge k+1``;
5. ``j \le r(p)``, ``p@j`` depends on ``q``, ``r(q) \ge m_q + 1`` ⇒
   ``r(p) \ge j+1``.

with saturation written ``r(p) = m_p + 1``. This is `find_reachable` as the
previous page states it, under the precondition that every class of the
universe is active on entry — which `drop_unmarked!` establishes, since it
ends by setting every surviving class's flag.

!!! note "Lemma D1 (a pushed prefix keeps a prefix)"
    If ``p`` is seeded and ``\{1, \dots, j-1\} \subseteq \mathrm{push}(p)``
    for some ``1 \le j \le m_p + 1``, then
    ``\{1, \dots, \min(j, m_p)\} \subseteq S(p)``.

    *Proof.* If ``j = m_p + 1`` then every class is pushed, so
    ``t(p) = \min \emptyset = \infty`` and ``S(p)`` is all of ``p``'s
    classes.

    Otherwise ``j \le m_p``. Every ``c < j`` is pushed, so ``t(p)`` is a
    minimum of ``\mathrm{key}_Q`` over a subset of ``\{c : c \ge j\}``, hence

    ```math
    t(p) \;\ge\; \min_{c \ge j} \mathrm{key}_Q(c)
          \;\ge\; \min_{c \ge j} \kappa_Q(c)
          \;=\; \kappa_Q(j),
    ```

    using Fact 0 for the middle step and monotonicity of ``\kappa_Q`` in the
    layout position for the last. Then for every ``i \le j``,
    ``\mathrm{key}_\emptyset(i) \le \kappa_Q(i) \le \kappa_Q(j) \le t(p)``,
    again by Fact 0 and monotonicity, so ``i \in S(p)``. ∎

!!! note "Lemma D2 (the prefix walk's pushes are matched)"
    Every fact ``r(p) \ge j`` derivable by ``W`` implies, at Rule B's
    fixpoint, that ``p`` is seeded and
    ``\{1, \dots, j-1\} \subseteq \mathrm{push}(p)``.

    *Proof.* Induction on the derivation of ``r(p) \ge j``.

    *(1)* ``j = 1``: rule 1 seeds ``p``, and the claimed inclusion is empty.

    *(2)* The premise ``j \le r(p)`` is a derived fact, so by the induction
    hypothesis ``p`` is seeded with ``\{1, \dots, j-1\} \subseteq
    \mathrm{push}(p)``, whence ``j \in S(p)`` by Lemma D1. ``p@j`` depends on
    ``q``, so rule 2 seeds ``q``; the inclusion for ``q`` at ``j = 1`` is
    empty.

    *(3)* As in (2), ``j \in S(p)``. ``p@j`` is deactivated, so rule 3 puts
    ``j \in \mathrm{push}(p)``; with the induction hypothesis this gives
    ``\{1, \dots, j\} \subseteq \mathrm{push}(p)``.

    *(4)* As in (2), ``j \in S(p)`` and ``k \in S(q)``. Conflicts are
    query-independent, so rule 4 pushes both, giving
    ``\{1, \dots, j\} \subseteq \mathrm{push}(p)`` and
    ``\{1, \dots, k\} \subseteq \mathrm{push}(q)``.

    *(5)* By the induction hypothesis applied to ``r(q) \ge m_q + 1``,
    ``\{1, \dots, m_q\} \subseteq \mathrm{push}(q)``, i.e. ``q`` is saturated
    in Rule B's sense. As in (2), ``j \in S(p)``, so rule 5 pushes it. ∎

!!! note "Theorem D (containment: the union is never necessary)"
    At the two fixpoints, for every package ``p``:

    ```math
    \{1, \dots, \min(r(p),\, m_p)\} \;\subseteq\; S(p).
    ```

    *Proof.* If ``r(p) = 0`` the claim is empty. Otherwise ``r(p) \ge r(p)``
    is derivable, so Lemma D2 gives ``p`` seeded with
    ``\{1, \dots, r(p)-1\} \subseteq \mathrm{push}(p)``, and Lemma D1 with
    ``j = r(p)`` gives the conclusion. ∎

So Theorems 5 and 6 of the previous page transfer unchanged: whatever the
prefix walk kept for ``Q`` is still kept, and the layered answer to ``Q``
survives filtering exactly as before. Running a prefix walk alongside and
keeping the union would add nothing.

!!! note "Remark (why one integer suffices)"
    ``\mathrm{push}(p)`` only grows, so ``t(p)`` — a minimum over its
    complement — only rises, so ``S(p) = \{c : \mathrm{key}_\emptyset(c) \le
    t(p)\}`` only grows and is **downward closed in ``\mathrm{key}_\emptyset``
    order**. Every rule that pushes a class has that class already in
    ``S(p)`` as a premise, so ``\mathrm{push}(p) \subseteq S(p)`` needs no
    separate enforcement and the kept set is exactly the final ``PB``.

    A set-valued frontier is therefore not required: the kept set is a
    prefix, in ``\mathrm{key}_\emptyset`` order rather than layout order, and
    one monotone pointer along each of the two orders represents it. That is
    what `find_reachable` implements.

## Composition

!!! note "Theorem C (relaxation-stable filtering)"
    Let ``F_Q`` be the filter under ``Q`` using Rules A and B, iterated to a
    fixpoint. Then for every ``R \le Q``:
    ``\mathrm{Lay}(F_Q(D), R) = \mathrm{Lay}(D, R)``.

    *Proof.* By Theorems A and B no pass ever deletes a class that
    ``\mathrm{Lay}(D, R)`` uses, for any ``R \le Q``; deletions compose, so
    ``\mathrm{Lay}(D, R)`` survives in ``F_Q(D)``. Theorem 4 of the previous
    page (reduction invariance) then gives equality. ∎

!!! note "Corollary C1 (the point)"
    A diagnosis may answer any relaxation of a failed query on the universe
    filtered for that query. Re-filtering per candidate fix is unnecessary.

[`relax`](@ref Resolver.relax) builds that question and `resolve` answers it.

Because ``\le`` is a lattice with the empty query at the bottom, the same
argument run at ``Q = \emptyset`` says a universe filtered for no constraints
at all serves every query over the same requirements — which is what makes a
filtered artifact worth caching.

## What this does not give

The guarantee is over relaxations of one query, not over arbitrary queries.
A *different* query — one that adds a requirement, or tightens a bound —
is not below ``Q`` in the order and is not covered; it needs its own filter
run. In particular this says nothing about upstream fixes, which change the
registry rather than the query, and therefore change ``D``.

Nor does it make the class *order* stable. Relaxing genuinely reorders
classes, and the layered answer depends on that order. What Theorem C says
is that the surviving universe still contains everything the reordered
answer needs — so a relaxation is answered by re-ranking and re-solving the
filtered universe, not by re-filtering it.

!!! warning "Formalization boundary"
    As on the previous page, the rule sets above are a faithful but
    hand-made reading of `class_ranking_keys`, `mark_necessary!` and
    `find_reachable`; the correspondence is checked by review and by the
    test suite rather than mechanically. `test/relaxation_stable.jl` checks the
    conclusions directly: for random universes and random queries it
    enumerates relaxations, re-ranks and re-solves the filtered universe
    under each, and compares against resolving that relaxation from the
    unfiltered artifact — an empirical check of Theorem C, and through it
    of the rules it rests on.
