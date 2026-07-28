# Pareto front optimality

This page documents a fully worked-out *alternative* semantics for the
resolver — necessity on the Pareto front — including its theorems, its
pure-SAT implementation, and the reasons it was ultimately set aside. It is
a digression from what the package implements, but a productive one: its
impossibility results explain *why* the layered semantics and the filter
fit together, its proof techniques power the theorems on the
[previous page](layered.md), and it closes the question of a purely
declarative optimality criterion rather than leaving it open. A reference
implementation lives on the `front-necessity` branch.

## Motivation

The layered solution is defined by a procedure. One can ask for something
more declarative: a criterion picking out *the* optimal solution using only
the solution set — no layers, no reference to the dependency walk, ideally
nothing but satisfiability queries. The natural candidate is a *necessity*
descent: repeatedly pin the highest-priority package that is **necessary**
— present in every remaining solution — at the best version available.
Necessity, however, is a *universal* property of a solution set, so its
meaning depends on which set is quantified over:

- **All tight solutions** (raw necessity): too weak — a package needed by
  every *sensible* solution counts as avoidable if one strictly-dominated
  solution technically avoids it, e.g. by taking a requirement's worst
  version merely because that version happens to drop a dependency edge.
  This distorts the descent, and it also forbids problem reduction
  entirely, which is ruinous at registry scale.
- **Models of the filtered problem**: fast, but then the *answer* depends
  on filtering details — sharpening or dulling the filter changes the
  result, which is unacceptable as a semantics.

The canonical middle ground is necessity over the **Pareto front**: quantify
over the solutions that are not strictly worse than others.

## The dominance order

For tight solutions ``s, t``:

```math
t \trianglerighteq s \quad\iff\quad
\mathrm{supp}(t) \supseteq \mathrm{supp}(s), \quad
\mathrm{rank}_t(p) \le \mathrm{rank}_s(p) \;\; \forall p \in \mathrm{supp}(s),
```

with strict improvement at some ``p \in \mathrm{supp}(s)`` (*coverage-
monotone version dominance*). This is a strict partial order on the finite
set ``R^*`` of tight solutions, so every solution is either maximal or
strictly dominated by a maximal one. The **front** ``F`` is the set of
maximal elements.

Both design choices in ``\trianglerighteq`` are forced:

- **Orientation.** An order in which *leaner* solutions dominate cannot
  work: the solutions that falsely break necessity are support-minimal
  (a requirement alone at a bad version), and a support-minimal solution
  has no leaner competitor, so no leaner-dominates order can remove one.
  Also, since a tight solution's support is a function of its versions,
  "leaner at equal quality" never occurs — leanness always trades against
  version quality, and the order must price the trade. ``\trianglerighteq``
  prices it as: a bigger solution wins only by covering everything at
  least as well and strictly better somewhere.
- **Junk is excluded by the space, not the order.** With the ``\supseteq``
  orientation, unjustified packages would ride along for free, so the
  candidate space must be the *tight* solutions. Tightness (well-founded
  justification, the same notion as stable models in answer-set
  programming) is not expressible by any dominance relation over versions
  — it is not monotone in the version lattice. Empirically, defining the
  front over junk-tolerant solutions changes answers on 92 of 7,464 small
  instances.

**The front-necessity semantics** is the forced descent over ``F``:
repeatedly pin the highest-priority package present in *every* front member
consistent with the pins so far, at the best version among them, until no
unpinned package is common.

## Theorem A: pin exactness

The linchpin of everything that follows.

!!! note "Theorem A"
    Let `pins` be the pin state after any number of descent rounds. If
    ``f \in F`` contains every pinned package at rank ≤ its pinned rank,
    then ``f`` agrees *exactly* with the pins.

    *Proof.* Suppose not; among pinned packages where ``f`` is strictly
    better, let ``\hat q`` be the one pinned earliest, at rank ``\hat r``.
    On packages pinned before ``\hat q``, ``f`` is ≤ with no strict
    improvement, hence equal — so ``f`` was in the remaining front set when
    ``\hat q`` was pinned. But then
    ``\hat r = \min\{\mathrm{rank}_g(\hat q)\} \le \mathrm{rank}_f(\hat q)
    < \hat r``. ∎

!!! note "Corollary A1"
    If ``t \in F`` and ``t \trianglerighteq s`` (or ``t = s``) for any
    ``s`` agreeing with the pins, then ``t`` agrees exactly with the pins:
    *dominating a pin-consistent solution can never buy pin-violating
    quality.*

## Theorem B: the answer is a front member

!!! note "Theorem B"
    The descent terminates with the remaining front set equal to exactly
    ``\{\mathrm{pins}\}``; in particular the answer is itself a valid,
    tight, undominated solution.

    *Proof.* Requirements are in every front member, so they get pinned.
    If ``q`` is pinned, every remaining ``f`` has ``q`` at the pinned
    version, so ``\mathrm{deps}(q, \mathrm{pins}(q))`` is common and gets
    pinned — the pins are dependency-closed. Now take any final survivor
    ``f`` and walk its closure: it starts at requirements (pinned), and
    whenever a closure node ``u`` is pinned, ``f(u) = \mathrm{pins}(u)``,
    so the next closure nodes are pinned too. Hence
    ``C(f) \subseteq \mathrm{dom}(\mathrm{pins})``, and by tightness
    ``\mathrm{supp}(f) = C(f) \subseteq \mathrm{dom}(\mathrm{pins})
    \subseteq \mathrm{supp}(f)`` — so ``f = \mathrm{pins}``. ∎

The same argument shape, one level down, shows the *layered* answer is
``\trianglerighteq``-undominated: a dominator would have to agree with
every layered pin before its earliest strict improvement, contradicting
that pin's optimality. So the layered solution always lies *on* the front;
front-necessity differs only in *which* front point it selects (230 of
7,464 small instances) — front-necessity makes priority the only ordering
principle among necessary packages, where the layered order lets
requirements outrank higher-priority dependencies.

## Theorem C: a pure-SAT implementation exists

Front-necessity is implementable against a SAT instance that encodes
validity only — no tightness encoding, which matters because eager
tightness encodings blow up (rooted-chain/level constructions), the naive
support encoding admits unrooted cycles, and lazily this is exactly
answer-set programming's unfounded-set problem. Instead, two cheap
operations run *outside* the solver: pruning ``\pi`` (a closure
computation) and support comparison. The key sub-results:

- **Witness transfer.** A model lacking ``p`` prunes to a tight solution
  lacking ``p``; conversely every tight solution is a model. So junk never
  distorts the existential side of a necessity query.
- **Dominator transfer.** A tight solution ``s`` has a tight strict
  dominator iff some model ``d`` covers ``\mathrm{supp}(s)`` at ranks ≤
  with a strict improvement *and* passes the **prune-check**
  ``\pi(d) \supseteq \mathrm{supp}(s)`` — junk can fake domination whose
  tight core doesn't actually cover ``s``, and the prune-check filters
  exactly those fakes.
- **Chain-top admissibility.** Climbing dominators terminates at a front
  member, which by Corollary A1 agrees exactly with the pins even though
  intermediate steps may not.
- **Cone blocking.** When a chain-top ``t`` contains ``p``, everything
  ``\trianglelefteq t`` is either dominated or is ``t`` itself, so one
  clause excludes the entire cone without losing any witness.
- **Greedy construction.** In the implementation, candidates are
  *constructed* rather than searched for: the layered greedy runs under
  per-query assumptions, and a junk-free greedy answer is provably
  Pareto-maximal within its query space (the Theorem A argument again).
  Restricting all queries to *supported* models — one clause per package:
  present implies depended-upon — loses no tight solution and reduces junk
  to dependency cycles.

With these, the descent's two oracle types (is ``p`` in every pin-exact
front member? what is the best rank among them?) are decided exactly, by
sequences of ordinary SAT solves. Correctness of the composed pipeline was
additionally verified against a brute-force implementation of the
definition on roughly 19,500 instance/order pairs with zero deviations.

## Theorem D: no version-level reduction is faithful

Here is the decisive negative result. One would like to run the necessity
descent on a *reduced* problem — that is how the layered semantics gets its
speed. For front-necessity this fails not just for the actual filter but
**in principle**:

!!! note "Theorem D"
    There is a problem on which the forced descent over models of the
    *exact front-support reduction* — keep a version iff some front member
    uses it, the strongest version-level reduction that could possibly be
    faithful — differs from the forced descent over the front itself.

    *Proof.* Requirements ``\{a\}``, priority ``p > a > d > b``, version 1
    best:

    | package | dependencies | conflicts |
    |---|---|---|
    | ``a`` | ``a@1 \to \{d\}``, ``a@2 \to \{b\}`` | |
    | ``b`` | ``b@1 \to \{p\}`` | |
    | ``d`` | ``d@1 \to \{p\}``, ``d@2 \to \{\}`` | ``d@1 \otimes p@1`` |
    | ``p`` | ``p@1 \to \{d\}``, ``p@2 \to \{\}`` | |

    ``R^*`` has four members and the front is
    ``F = \{f_1, f_2\}`` with ``f_1 = \{a@2, b@1, p@1, d@2\}`` and
    ``f_2 = \{a@1, d@1, p@2\}`` — both contain ``p``, so the front descent
    pins ``p`` first at rank 1 and lands on ``f_1``. But every version
    above appears in some front member, so the exact front-support
    reduction keeps everything — including the **recombinant**
    ``c = \{a@1, d@2\}`` (``f_2``'s ``a`` with ``f_1``'s ``d``), a tight
    solution missing ``p``. In the reduced problem ``p`` is no longer
    common; the descent pins ``a`` first, and rank-minimization at ``a``
    eliminates ``f_1``: the reduced answer is ``f_2 \ne f_1``. (Verified
    computationally.) ∎

The moral: necessity is a property of the solution *set*, and a version
universe only bounds that set from above — recombinations of individually
front-supported versions can erase necessity. Universal queries survive
only when asked against the front itself. Note the contrast with
Theorem 4 of the [layered page](layered.md): the layered semantics needs a
reduction to preserve *one solution*; front-necessity would need it to
preserve the *entire front and nothing that changes commonality* — and
even the latter is unattainable at the version level. Such counterexamples
require engineered structure (two dependency routes, a recombination
orphaning the necessity package, a conflict pinching the diagonal): none
arose in 12,000 random 4-package instances, which is also why smaller
empirical sweeps had suggested the reduction was safe.

## Why the filter cannot serve the front

The actual reachability filter fails front-preservation in an even more
basic way: it **drops front members outright**. On the Theorem D instance
it keeps only ``\{a@1, d@1, d@2, p@1, p@2\}`` — never reaching ``a@2``,
deleting ``b`` — yet ``f_1`` uses both. The filter's designed invariant is
*"a worse version is only chosen under conflict pressure against the
better version,"* and ``a@1`` is never conflicted. But ``f_1`` takes
``a@2`` for a different kind of reason: under ``a@1``, the only dependency
edge that can justify ``p`` is ``d@1 \to p``, and ``p@1 \otimes d@1`` — so
on the ``a@1`` route the better ``p@1`` is unreachable-as-justified.
Switching to ``a@2`` reroutes justification through ``b`` and frees
``p@1``. The downgrade at ``a`` is paid for by an upgrade at ``p`` — an
*opportunity-cost trade*, mediated not by any conflict against ``a`` but by
**version-varying dependency sets** rewiring the justification topology.

That is precisely the boundary of the filter's intuition:

!!! note "Fragment lemma"
    If every version of each package has the same dependency set, then on
    the front every non-best version is blocked by a conflict with a
    co-installed package — the filter's invariant, as a theorem.

    *Proof (one step).* Let ``s \in F`` with ``q`` at a non-best version
    and ``v'`` better. The in-place upgrade ``s' = s[q \mapsto v']`` has
    the same support, and — dependency sets being version-independent —
    the same justification edges, hence remains tight. If ``s'`` were
    conflict-free it would strictly dominate ``s``. So some conflict
    between ``q@v'`` and a member of ``s`` blocks it. ∎

Pointwise (layered) optimality never makes opportunity-cost trades, which
is why conflict-driven degradation covers everything it can ever choose
(Theorem 5 of the layered page); the front additionally contains the
trade-mediated points, and front-*necessity* quantifies over exactly those.

## Why it was set aside

| | layered | front-necessity |
|---|---|---|
| canonical, filter-free definition | yes (Theorems 4–6) | yes (by construction) |
| declarative (procedure-free) statement | no — defined by the trace | yes |
| SAT query character | existential, answer-witnessed | universal (necessity, dominance) |
| admissible reductions | any preserving *one* solution | none at the version level (Theorem D) |
| registry performance | milliseconds–seconds | small: ms; `DataFrames`: ~85 ms; `DifferentialEquations` (26k versions, unreduced): intractable |

The costs are intrinsic, not implementation slack: universal queries must
run against the unreduced problem (Theorem D), every semantically-avoidable
candidate needs a witness construction plus a dominance query, and witness
certificates are invalidated by later pins. Even under optimistic
assumptions the descent is an order of magnitude or more behind the layered
implementation — and end users do not observe the difference between the
two Pareto points the semantics select. What the exploration contributes is
what this page records: the impossibility theorems that close the
declarative route, the proof machinery reused for the layered guarantees,
and a precise understanding of what the filter does and does not preserve —
satisfiability (provably), the layered answer (provably), the front
(provably not).
