# Necessity on the Pareto front: canonical single-solution semantics

Status: design notes for the `front-necessity` branch. Builds on the
single-solution API (`resolve` returns a package ⇒ version `Dict`, or
`nothing` when unsatisfiable). Complements issue #33.

## 1. Motivation

The single-solution resolver pins packages one at a time, best version
first, in a caller-supplied priority order. Two questions determine its
semantics: *which package gets pinned next*, and *relative to what set of
solutions "best" and "next" are judged*.

The dependency-layered descent answers the first question syntactically
(requirements first, then dependencies of chosen versions, layer by
layer), which subordinates the priority order to dependency depth and ties
the resolver to the dependency structure. The natural alternative —
repeatedly pin the highest-priority package that is *necessary*, i.e.
present in every remaining solution — makes priority the only ordering
principle and needs nothing but a SAT oracle. But necessity is a
*universal* property of a solution set, so its meaning depends on which
set you ask about:

- **All valid solutions** (raw necessity): a package needed by every
  *sensible* solution still counts as avoidable if some strictly-dominated
  solution technically avoids it (e.g. by taking a requirement's worst
  version just to drop a dependency edge). This both distorts pin order
  and is prohibitively slow at registry scale, since it forbids problem
  reduction: on a `DifferentialEquations` resolve the unreduced problem
  takes ~41 s vs ~0.8 s reduced.
- **Models of the reduced problem** (what `pkg_info` filtering yields):
  fast, but the answer then depends on the reduction — resolving the same
  requirements with a sharper or duller filter can change the result,
  which is not acceptable as a semantics.

The dependency-layered descent is immune to this because all of its SAT
queries are *existential* and witnessed by its own final answer, which
survives any front-preserving reduction; universal queries are inherently
sensitive to removing models. So a necessity-based resolver needs a
canonical solution set. These notes develop that set — the Pareto front
under a specific dominance order — prove that the resulting semantics is
implementable with plain SAT queries plus two cheap outside-the-solver
operations, and prove that no version-level problem reduction (not even a
perfect one) could serve instead.

## 2. Definitions

Fix package data (each package a finite, best-first ordered version list;
each version a dependency set and conflict constraints) and requirements
`reqs`. Write `rank_s(p)` for the rank (1 = best) of `p`'s version in a
solution `s`, and `supp(s)` for the set of packages `s` assigns.

- A **solution** is a partial version assignment that covers `reqs`, is
  dependency-closed (deps of every assigned version are assigned), and
  conflict-free.
- The **closure** `C(s)` is the least set containing `reqs` and closed
  under following dependencies of `s`'s chosen versions. A solution is
  **tight** if `supp(s) = C(s)`; `R*` denotes the tight solutions.
  **Pruning** `π(s)` restricts `s` to `C(s)`.

  *Lemma 1.* `π(s)` is a tight solution with `s`'s versions, and
  `C(π(s)) = C(s)`. (The closure computation only ever inspects versions
  of packages already in the closure.)

- **Dominance** (coverage-monotone version dominance):

  > `t ⊵ s` iff `supp(t) ⊇ supp(s)`, `rank_t(p) ≤ rank_s(p)` for every
  > `p ∈ supp(s)`, and strictly for at least one such `p`.

  `⊵` is a strict partial order, so every member of the finite set `R*`
  is either maximal or strictly dominated by a maximal member
  (*Lemma 2*). The **front** `F` is the set of maximal members.

- **Forced descent** `FD(C)` over a candidate set `C`: starting from
  empty pins, repeatedly take the highest-priority unpinned package
  present in *every* remaining candidate, pin it at the minimum rank any
  remaining candidate gives it, and drop candidates that disagree; stop
  when no unpinned package is common. The pins are the answer.

**The semantics adopted here is `FD(F)`**: repeatedly pin the
highest-priority package that every Pareto-optimal tight solution
consistent with the pins must include, at the best version among them.

### Why this dominance order

Two design decisions in `⊵` are load-bearing, and both were forced by
failed alternatives:

- **Orientation.** Making *leaner* solutions dominate (`supp(t) ⊆
  supp(s)`) cannot work: the solutions that falsely break necessity are
  support-minimal (a requirement alone at a version that drops its
  dependency edges), and a support-minimal solution has no leaner
  competitor, so no leaner-dominates order can ever remove one. Removing
  them requires the `⊇` orientation — a bigger solution may dominate by
  covering everything at least as well and strictly better somewhere.
  Note also that for tight solutions support is a *function* of the
  version assignment (the closure), so "leaner at equal quality" never
  occurs between distinct tight solutions; leanness differences always
  come bundled with version differences that the order must price.
- **Junk is excluded by the space, not the order.** With the `⊇`
  orientation, unjustified extra packages would ride along for free, so
  the candidate space must be the *tight* solutions. Tightness
  (well-founded justification from the requirements — the same notion as
  ASP's stable-model semantics) is not expressible by any dominance
  relation over versions: it is not monotone in the version lattice.
  Empirically, defining the front over junk-tolerant solutions changes
  the answer on 92 of 7,464 tiny instances.

## 3. Theorem A (pin exactness)

*Let `pins` be the pin state after any number of `FD(F)` rounds. If
`f ∈ F` contains every pinned package at rank ≤ its pinned rank, then `f`
agrees exactly with the pins.*

**Proof.** Suppose not. Among pinned packages where `f` is strictly
better, let `q̂` be the one pinned earliest, at round `j` with rank `r̂`.
On packages pinned before round `j`, `f` is ≤ with no strict improvement
(`q̂` was earliest), hence equal — so `f` agrees exactly with the
round-`(j−1)` pins and is in the remaining front set at round `j`. But
then `r̂ = min{rank_g(q̂) : g remaining} ≤ rank_f(q̂) < r̂`. ∎

**Corollary A1.** If `t ∈ F` and `t ⊵ s` (or `t = s`) for any `s` that
agrees with the pins, then `t` agrees exactly with the pins — domination
gives containment and rank-≤ on `supp(s) ⊇ dom(pins)`, and Theorem A does
the rest. *Dominating a pin-consistent solution can never buy
pin-violating quality.*

## 4. Theorem B (the answer is a front member)

*`FD(F)` terminates with the remaining front set equal to `{pins}`; in
particular the answer is itself a valid, tight, undominated solution.*

**Proof.** (i) Requirements lie in every front member, so they stay
common until pinned. (ii) If `q` is pinned, every remaining `f` has `q`
at `pins(q)`, so `deps(q, pins(q)) ⊆ supp(f)` for every remaining `f`;
those packages are common and get pinned before termination — the pins
are dependency-closed. (iii) Take any final survivor `f` and walk its
closure `C(f)`: it starts at `reqs` (pinned, by (i)); whenever a closure
node `u` is pinned, `f(u) = pins(u)`, so the next closure nodes
`deps(u, pins(u))` are pinned by (ii). Inductively every node of `C(f)`
is pinned. (iv) `f` is tight, so
`supp(f) = C(f) ⊆ dom(pins) ⊆ supp(f)`; hence `supp(f) = dom(pins)` and,
by agreement, `f = pins`. ∎

## 5. Theorem C (a pure-SAT implementation)

The SAT instance encodes validity only (packages, versions, deps,
conflicts, requirement units) — **no tightness encoding**. Models may
therefore carry junk. Pins are carried as *assumptions* (they must not
constrain dominance queries, which range over all tight solutions). Two
outside-the-solver operations are used: pruning `π` (a closure
computation over the chosen versions) and support-comparison.

The oracle `EXISTS-FRONT(Π, φ)` — is there a front member agreeing
exactly with pins `Π` and satisfying `φ` (either "`p` absent" or "`p` at
version `v`") — is computed as:

1. Find a model `m` of the instance under assumptions `Π ∧ φ` (minus
   blocked assignments). UNSAT ⟹ **no**.
2. `s := π(m)`. If `s` dropped a pinned package (or `p@v`), block `m`
   and go to 1 — junk carried the constraint.
3. **Climb**: while `s` has a strict tight dominator, replace `s` by one.
   A single climb step is: solve (without pin assumptions) for a model
   `d` containing every `q ∈ supp(s)` at rank ≤ `rank_s(q)`, strictly
   better somewhere; check `π(d) ⊇ supp(s)` (the **prune-check**); on
   failure block `d` and re-solve. The chain-top `s*` has no tight
   dominator, so `s* ∈ F`.
4. By Corollary A1, `s*` agrees exactly with `Π`. If `s*` satisfies `φ`,
   **yes** (with `s*` as a certificate). Otherwise add the *cone block*
   ¬(⊴ `s*`) — one clause: some package outside `supp(s*)` is installed,
   or some package beats `s*`'s rank — plus an exact block for `m`, and
   go to 1.

Correctness rests on:

- **C1 (witness transfer).** Pruning preserves validity, tightness
  (Lemma 1) and "`p` absent"; models whose pruned form drops a
  constraint are safely skipped because a true witness `f ∈ F` is its own
  model with `π(f) = f`, so completeness never depends on skipped models.
- **C2 (dominator transfer).** `s` has a strict tight dominator iff some
  loose model `d` passes the prune-check: `π(d)` keeps `d`'s versions,
  the strict package lies in `supp(s) ⊆ supp(π(d))`, and conversely a
  tight dominator passes as itself.
- **C3 (climb termination).** `⊵` is a strict order on a finite set.
- **C4 (chain-top admissibility).** `s* ⊵ ⋯ ⊵ π(m) ⊇ Π` and transitivity
  put `s*` under Corollary A1: chain-tops never violate pins. This is
  where Theorem A carries the construction: intermediate climb elements
  may improve pinned packages; the top provably cannot.
- **C5 (cone-block soundness).** Everything ⊴ `s*` is dominated or is
  `s*` itself, which failed `φ`; no true witness is excluded. Exact
  blocks guarantee progress against junky models outside the cone;
  termination is finite. (The exact block is the simple, complete
  fallback; ASP-style unfounded-set clauses for the orphaned component
  are the sharper alternative if this loop ever matters at scale.)
- **C6 (version probing).** After `p` is known forced, probe its ranks
  best-first with `φ = p@v`. A chain-top at probe rank `r` contains `p`
  at rank ≤ `r`; strictly better would have been a certificate for an
  earlier probe that answered no — contradiction. So chain-tops of
  version probes are immediate certificates, and probes advance only on
  UNSAT, which is conclusive.

With both oracles exact, the descent reproduces `FD(F)` round by round.
∎ (Empirically: zero deviations over ~19,500 instances at ≤ 4 packages.)

## 6. Theorem D: no version-level reduction is faithful

It is tempting to implement `FD(F)` as "reduce the problem to the
versions that appear in some front member, then run the ordinary
model-necessity descent" — the reduction the current `pkg_info` filter
approximates. Over 7,464 instances with ≤ 3 packages this agreed with
`FD(F)` every time. It is nevertheless **wrong**, and not because of
filter slack: even the *exact* front-support reduction fails.

Counterexample (4 packages; reqs = `{a}`; priority `p > a > d > b`;
version 1 best):

| package | deps | conflicts |
|---|---|---|
| `a` | `a@1 → {d}`, `a@2 → {b}` | |
| `b` | `b@1 → {p}` | |
| `d` | `d@1 → {p}`, `d@2 → {}` | `d@1 ⊗ p@1` |
| `p` | `p@1 → {d}`, `p@2 → {}` | |

`R*` has four members; the front is `F = {f₁, f₂}` with
`f₁ = {a@2, b@1, p@1, d@2}` and `f₂ = {a@1, d@1, p@2}` — **both contain
`p`**, so `FD(F)` pins `p` first at rank 1 and lands on `f₁`. But every
version above appears in some front member, so the exact front-support
reduction keeps everything — including the **recombinant**
`c = {a@1, d@2}` (`f₂`'s `a` with `f₁`'s `d`), a tight solution missing
`p`. In the reduced problem `p` is no longer common; the descent pins `a`
first, and rank-minimization at `a` eliminates `f₁`: the reduced answer
is `f₂ ≠ f₁`.

The moral: necessity is a property of the solution *set*, and a version
universe only bounds that set from above — recombinations of
individually-front-supported versions can erase necessity. Universal
queries survive only when asked against the front itself. Reductions
(including the existing filter) remain sound *performance* devices for
finding models and witness candidates; they can never adjudicate
necessity. (Such counterexamples require engineered structure — two
disjoint dependency routes, a recombination orphaning the necessity
package, and a conflict pinching the diagonal; none arose in 12,000
random 4-package instances.)

## 7. Empirical summary

Brute-force studies over exhaustive 2-package and random ≤ 4-package
instances, both priority orders (~19,500 instance/order pairs):

| claim | outcome |
|---|---|
| filtered-model necessity = `FD(F)` | fails (54 / 7,464 at ≤ 3 pkgs) |
| priority-lex-max over `F` = `FD(F)` | fails (61 / 7,464) |
| `FD(F)` = dependency-layered answer | differs (230 / 7,464, by design) |
| front over junk-tolerant solutions = `FD(F)` | fails (92 / 7,464) |
| exact front-support reduction = `FD(F)` | fails (constructed, §6) |
| lazy pipeline (§5) = `FD(F)` | exact everywhere tested |

## 8. Implementation notes (this branch)

- `resolve_core` pins by front-necessity: scan, in priority order, the
  intersection of the supports of *certified front members* (any package
  outside it is certifiably avoidable); the necessity test runs the §5
  oracle. Certified front members are cached and reused as avoidability
  certificates; the cache is filtered on each pin (keep members agreeing
  with the pin) and reseeded by the certifying chain-top.
- Cheap positive fast paths, all sound for front-necessity since
  raw-forced ⟹ front-forced: requirements; dependencies of pinned
  versions (maintained as a set outside the solver, replacing the
  top-level-propagation trick, which pins-as-assumptions would defeat);
  top-level forcings from the requirement units.
- Version selection: optimize in model space first (the existing
  improvement loop, pins assumed per solve) to a lower bound `r₀`, then
  certify by climbing the pruned optimum; by C6, the first certified rank
  is the answer, and probes advance only on UNSAT.
- Dominance queries open a temporary clause scope (containment + bound +
  strictness clauses) and never assume pins.
- The dependency structure is consulted only by the closure oracle
  (pruning, prune-checks, dep-forced fast path) — a checking device, not
  a search structure. Tightness is inherently a dependency notion; it
  cannot be recovered from the SAT abstraction (and eager SAT encodings
  of it blow up; lazily it is ASP's unfounded-set problem).
- The brute-force reference in the test harness implements §2 verbatim on
  raw data: enumerate tight solutions, form the front, run the descent.
  The reference no longer involves `pkg_info` at all.

## 9. Open questions

- Registry-scale cost: climb lengths, dominator-query cost, cone-block
  frequency. (Benchmarks tracked on this branch.)
- Proof-hardening: Theorems A–C are believed complete as stated; a
  mechanized or adversarially-reviewed pass would be prudent before this
  replaces the dependency-layered descent.
- Whether requirements should outrank higher-priority dependencies is a
  policy choice orthogonal to this document: it can be recovered by
  putting requirements first in the priority order.
