# Explaining an unsatisfiable resolve — the complete theory

This document is the whole theory of diagnosing and explaining an
unsatisfiable resolver query: what a report must contain, what every part of
it means, and the proof of every claim the design rests on. It is written to
be implemented from. It names no functions and assumes no code; what it
assumes of the surrounding resolver is collected in one section
(Section 10) as a short list of
interfaces.

The contract, up front. A report for an unsatisfiable query owes the reader
three things:

1. **The fixes.** How many independent things are wrong, and for each, a menu
   of alternatives — such that *any* combination of choices, one per menu,
   repairs the query, giving up as little as any repair can.
2. **The explanations.** For each conflict, an accurate account of where it
   comes from: which of the user's own facts collide, at which package, and
   through which of the registry's statements.
3. **The witnesses.** For each fix, what taking it yields, with enough on the
   page that the reader can see why the fix works.

What a report does **not** owe is a hand-checkable proof. Every sentence
printed must be machine-verified true, and the contradiction must be visible
on the page — but the reader's job is to *understand*, not to re-derive. It
also does not owe completeness of every enumeration; it owes **honesty about
every gap**, at three named places
(Section 11).

Sections 1–2 build the logic. Sections 3–5 are the fix half: repairs, menus,
and which reasons each conflict answers for. Sections 6–7 are the explanation
half: the pivot form and what it guarantees about the fixes. Sections 8–11
say what to verify, how to word things, what to assume, and where the edges
are.

## 1. Setting

There is a finite set of **packages**; each package `p` has a finite nonempty
set of **versions** `V(p)`, totally ordered ("newer"). Write
`V⊥(p) = V(p) ∪ {⊥}`, where `⊥` means *not installed*. An **installation**
`ι` assigns every package an element of `V⊥(p)`.

All constraints — the registry's and the query's — are **clauses** in the
logic of Section 2. Two actors supply them:

* **The registry** is a finite, fixed, satisfiable set `ℛ` of clauses. Each
  clause of `ℛ` must be one independently intelligible statement — a
  dependency edge over a range of the depender's versions, or one
  compatibility bound between a range of one package's versions and a range
  of another's. This granularity is a requirement, not a convenience: a
  clause of `ℛ` is the unit of blame *and* the unit of print, so each must be
  checkable against a single registry entry.
* **The query** is a finite set `Q` of **facts**, of two kinds:
    * `Requirement(p)` — the clause `⟨p : V(p)⟩`: *p is installed*.
    * `Constraint(p)` — the clause `⟨p : L(p) ∪ {⊥}⟩`, where `L(p)` is the set
      of versions the query's own constraints (compat, pins, admission
      policy) leave of `p`: *if p is installed, it is one of these*.

Each fact must be a single clause **and** a single user-performable action:
dropping the requirement on `p`, or lifting one's own constraints on `p`.
Package granularity is the finest grain that stays a user action — no user
edit gives back one excluded version of a package while keeping its siblings
excluded — and that is why facts are per-package, not per-version.

A **withdrawal** of `W ⊆ Q` keeps `ℛ ∪ (Q ∖ W)`. Withdrawing more facts only
removes clauses, so satisfiability is monotone in withdrawal: every model
survives. We assume throughout that `ℛ ∪ Q` is unsatisfiable (otherwise there
is nothing to report) and that `ℛ` alone is satisfiable (otherwise the
registry is broken, which is a different report). The second assumption
guarantees that repairs exist: withdrawing all of `Q` always works.

## 2. The clause logic

A **literal** on package `p` is a subset of `V⊥(p)`. A **clause** `C` assigns
a literal `C(p)` to each package, with `C(p) = ∅` for all but finitely many;
its **support** is the set of packages with nonempty literal. An installation
satisfies `C`, written `ι ⊨ C`, when `ι(p) ∈ C(p)` for some `p`. A clause
with empty support is unsatisfiable; a clause with `C(p) = V⊥(p)` is a
tautology. Clause `C` **subsumes** `D` when `C(p) ⊆ D(p)` for every `p`: then
every model of `C` models `D`.

Six readings of one shape cover everything a report says:

| clause | reads as |
|---|---|
| `⟨A : V(A)⟩` | A is installed |
| `⟨A : L ∪ {⊥}⟩` | the query (or the registry) leaves only `L` of A |
| `⟨A : (V(A)∖R) ∪ {⊥}⟩` | A cannot be any of `R` |
| `⟨A : (V(A)∖R) ∪ {⊥}, B : S⟩` | `A@R` **requires** `B@S` |
| `⟨A : (V(A)∖R) ∪ {⊥}, B : S ∪ {⊥}⟩` | `A@R` **constrains** `B@S` |
| `⟨A : (V(A)∖R) ∪ {⊥}, B : {⊥}⟩` | `A@R` **leaves no version of** `B` |

The difference between the verbs is only what the consequent admits:
*requires* brings the package in, *constrains* permits its absence, and a
consequent of `{⊥}` alone rules the package out. Whether a statement forces
is data in the clause, never a flag beside it — which is why either package
may be taken as the consequent, and the verb follows (Section 9).

**One rule.** Clauses compose by resolution on a package. It is stated for
any finite set of clauses at once, because literals are sets and more than
two clauses can be needed to empty one:

> **Rule (resolution on `q`).** For clauses `C₁, …, Cₙ`, the resolvent
> `Res_q(C₁,…,Cₙ)` has literal `⋂ᵢ Cᵢ(q)` at `q` and `⋃ᵢ Cᵢ(p)` at every
> other `p`.

**Lemma 1 (soundness).** If `ι ⊨ Cᵢ` for every `i`, then
`ι ⊨ Res_q(C₁,…,Cₙ)`.

*Proof.* If some `Cᵢ` is satisfied at a package `p ≠ q`, then
`ι(p) ∈ Cᵢ(p) ⊆ ⋃ⱼ Cⱼ(p)`. Otherwise every `Cᵢ` is satisfied at `q`, so
`ι(q) ∈ ⋂ᵢ Cᵢ(q)`. ∎

In particular every resolvent is entailed by the clauses it came from —
soundness of everything derived is inherited, never re-established.

**Elimination.** To eliminate package `q` from a clause set `Γ`: let `Γ_q` be
the clauses whose support contains `q` and `Γ'` the rest, and put

```
Elim_q(Γ) = Γ' ∪ { Res_q(S) : S ⊆ Γ_q minimal with ⋂_{C∈S} C(q) = ∅ }.
```

Each such resolvent has `∅` at `q`, i.e. does not mention `q`. Restricting to
minimal emptying subsets loses nothing: the resolvent of a larger emptying
subset is subsumed by that of a minimal one inside it. This is the coverage
condition of the design, stated as the rule itself: eliminating a package
resolves **every** minimal set of clauses that leaves it nothing to be —
pairs do not suffice, because three set-valued literals can have empty
intersection while every two overlap.

**Theorem 2 (elimination preserves satisfiability).** `Γ` is satisfiable if
and only if `Elim_q(Γ)` is.

*Proof.* (⇒) A model of `Γ` satisfies every resolvent by Lemma 1, and since
a resolvent's `q`-literal is empty it is satisfied at some other package; so
the model, unchanged, models `Elim_q(Γ)`.

(⇐) Let `ι` model `Elim_q(Γ)`; we choose a value at `q`. Let
`U = { C ∈ Γ_q : ι(p) ∉ C(p) for every p ≠ q }` — the clauses that must be
satisfied at `q` itself. If `⋂_{C∈U} C(q) = ∅`, pick a minimal `S ⊆ U` with
empty intersection; its resolvent is in `Elim_q(Γ)`, so `ι` satisfies it at
some `p ≠ q`, i.e. `ι(p) ∈ C(p)` for some `C ∈ S ⊆ U` — contradicting
`C ∈ U`. So the intersection is nonempty; set `ι(q)` to any member (if `U` is
empty, any value at all). Then `ι ⊨ Γ`. ∎

**Projection.** For a keep-set `K` of packages, `proj_K(Γ)` eliminates the
packages outside `K` one at a time, in any order. By induction on Theorem 2 it
preserves satisfiability, and every clause of it has support inside `K`.
Different elimination orders may produce different (mutually subsuming) clause
sets; the theory licenses any of them, and the preferences of Section 6 choose
among what they print.

## 3. Reasons and repairs

Throughout, `ℛ` is fixed and satisfiability of a set `X ⊆ Q` means
satisfiability of `ℛ ∪ X`.

* A **reason** is a minimal unsatisfiable subset of `Q`: it cannot hold, and
  every proper subset can. Reasons say *why* the query fails.
* A **correction set** is a `W ⊆ Q` whose withdrawal is satisfiable. A
  **repair** is a minimal correction set; `F` is the set of repairs, and
  `F_min ⊆ F` the repairs of smallest cardinality `k`. Repairs say *what to
  change*; `F_min` are the ones that give up as few of the user's facts as
  any repair can.

Since `ℛ` is satisfiable, `Q` is a correction set, so `F ≠ ∅`. Since `ℛ ∪ Q`
is unsatisfiable and everything is finite, reasons exist: every unsatisfiable
subset contains one.

**Theorem 3 (hitting duality).** `W` is a correction set if and only if `W`
intersects every reason. Hence the repairs are exactly the minimal hitting
sets of the reasons.

*Proof.* If some reason `R` is disjoint from `W`, then `R ⊆ Q ∖ W`, so the
withdrawal of `W` keeps an unsatisfiable subset and is unsatisfiable.
Conversely, if `W` meets every reason, then `Q ∖ W` contains no reason; were
`Q ∖ W` unsatisfiable it would contain a minimal unsatisfiable subset — a
reason — by finiteness. So it is satisfiable and `W` is a correction set.
Minimality on both sides is the same condition on the same sets. ∎

Two distinct repairs are incomparable (both minimal), a fact used silently
below.

**Finding `F_min`.** One oracle suffices: *is there an installation
satisfying `ℛ` that violates at most `k` facts of `Q`, while for each blocked
set `M` satisfying at least one fact of `M`?* (In a solver this is the failed
instance plus a cardinality bound over the fact indicators plus one clause
per blocked set; but any oracle with this signature will do.)

**Theorem 4 (bounded enumeration).** (a) The oracle with bound `k` and no
blocked sets succeeds iff some correction set of size ≤ `k` exists; the least
such `k` is the size of the cheapest repairs. (b) At that `k`, iterating with
each found repair blocked enumerates exactly `F_min` and terminates.

*Proof.* (a) A model violating `U`, `|U| ≤ k`, witnesses that `U` is a
correction set; a correction set `C` with `|C| ≤ k` has a model of its
withdrawal, which violates a subset of `C`. (b) Let `k` be minimal. A found
model's violated set `U` is a correction set with `|U| ≤ k`; it contains a
repair, which has size ≥ `k`; so `U` *is* that repair and `|U| = k`, i.e.
`U ∈ F_min`. The blocking clauses exclude exactly the already-found members
(distinct repairs are incomparable, so no unfound member contains a found
one... more precisely: `U` satisfies "some fact of each blocked `M` holds"
exactly when `M ⊄ U`, and for distinct members of `F_min` that is automatic),
so each iteration finds a new member, and while an unfound `M ∈ F_min`
remains, the model witnessing `M` violates exactly `M`, contains no other
member, and satisfies every block — the oracle succeeds. Finiteness
terminates it. ∎

**Theorem 5 (one question decides whether anything larger exists).** With no
cardinality bound and all of `F_min` blocked, the oracle succeeds iff
`F ≠ F_min`.

*Proof.* A model's violated set `U` contains no member of `F_min`; the repair
inside `U` is therefore not in `F_min`, so `F ∖ F_min ≠ ∅`. Conversely for
`N ∈ F ∖ F_min`, the model witnessing `N` violates exactly `N`, which
contains no member of `F_min` by incomparability, so every block is
satisfied. ∎

## 4. Menus

A **presentation** is a family of pairwise disjoint nonempty fact-sets
`D₁, …, D_m` — the menus — one per **conflict**. A **selection** picks one
entry from each menu; its pick-set is the set of picked facts. The
presentation contract has three clauses:

* **(C1)** every selection's pick-set is a repair — the offer is sound;
* **(C2)** every selection's pick-set is in `F_min` — the offer is as cheap
  as any repair;
* **(C3)** if some member of `F_min` is no selection, the report says so.

Everything a menu page claims follows from these three, and Section 5's
theorems need only (C1). Two constructions meet the contract.

**The product.** Sometimes `F_min` *is* a product, and the check is cheap:

**Theorem 6 (a family that counts right is a product).** Let `C₁, …, C_k`
partition the facts `F_min` uses, suppose every `M ∈ F_min` contains exactly
one fact of each `Cᵢ`, and `|F_min| = ∏ᵢ |Cᵢ|`. Then
`F_min = C₁ × ⋯ × C_k` (as pick-sets).

*Proof.* Sending each `M` to the tuple of its facts is well defined by the
first hypothesis and injective because a set is determined by its elements.
An injection between finite sets of equal size is a bijection, so every tuple
is some `M`. ∎

**Lemma 7 (the partition is forced).** If `F_min = C₁ × ⋯ × C_k`, then two
facts share a repair in `F_min` iff they lie in different factors. Hence the
factors are the connected components of the "never share a repair" relation
on the used facts, and the candidate partition for Theorem 6 is unique and
computable without search.

*Proof.* Each member has exactly one fact per factor, so same-factor facts
never share. For `c ∈ Cᵢ`, `c′ ∈ Cⱼ`, `i ≠ j`, the product contains a member
picking both. So "never share" holds exactly within factors; within a factor
every pair is directly related, across factors no pair is, and the components
are the factors. ∎

Offering `Cᵢ` as the `i`-th menu satisfies (C1) and (C2) with every member of
`F_min` reached, so (C3) has nothing to say. The conflict count `k` is not an
artifact of any search: it is the size of the smallest fix. And the shape is
not a choice:

**Theorem 8 (canonicity).** Suppose `F_min = C₁ × ⋯ × C_k` and let
`D₁, …, D_m` be any presentation in which every selection's pick-set is in
`F_min` and every member of `F_min` is some selection's pick-set. Then
`m = k` and the `D`s are the `C`s, up to order.

*Proof.* First, the factors `Cᵢ` are pairwise disjoint: a fact in two factors
would let a tuple pick it twice, giving a pick-set of size `k−1` in `F_min`,
which has members of size `k` only. Menus are disjoint by definition, so every
pick-set has exactly `m` elements; pick-sets lie in `F_min`, so `m = k`.

No two distinct facts of one factor `Cᵢ` are pickable at different menus: a
selection taking both would have a pick-set in `F_min` with two facts of
`Cᵢ`. Every fact of `Cᵢ` occurs in some member of `F_min`, hence in some
menu; so all of `Cᵢ` lies in a single menu `D_{σ(i)}`, and no fact lies in
two menus (disjointness). `σ` is injective: if `Cᵢ, Cⱼ ⊆ D_p` with `i ≠ j`,
a member of `F_min` needs one fact of each, but both are available only as
the single pick from `D_p`. So `σ` is a bijection of `{1..k}`. Finally
`D_{σ(i)}` holds nothing beyond `Cᵢ`: an extra fact `d` would enter some
pick-set, hence lie in `⋃ⱼ Cⱼ`, hence in `Cⱼ` for some `j ≠ i`, hence in
`D_{σ(j)}` — against disjointness. ∎

So in the product case the report's shape is forced by the repair family
itself. This is the menu-side twin of the pivot theorem of Section 6: there
the choice is proved entirely free, here it is proved entirely absent, and
each is the strongest possible statement of its kind.

**Entries may be sets.** Nothing above needs an entry to be a single fact.
With entries as fact-sets and a selection taking the union of one entry per
menu, (C1) and (C2) read verbatim, and Theorem 9 generalizes with
"contains some whole menu" read as "meets every entry of some menu" — the
proof is the same, picking per menu an entry disjoint from the reason. A
menu whose entries are whole repairs is the degenerate case: one choice,
compound entries, trivially a product. The primary presentation below keeps
entries small — single facts wherever the family lets it — so everything in
Sections 5–7 applies to it unchanged.

**Theorem 23 (the finest factorization).** Call a partition `P` of the used
facts *factoring* when `F_min` is exactly the unions of one trace per block
— `F_min|B = {M ∩ B : M ∈ F_min}` — with every combination realized. Then:
the trivial partition factors; if `P` and `Q` factor, so does their common
refinement; hence a unique finest factoring partition exists. Moreover each
block's traces have one size.

*Proof.* The trivial partition factors by definition. For the refinement:
it suffices that any one block `A ∩ B` of the refinement swaps freely, and
the general combination follows block by block. Given members `M, M′`, form
`M₁` as `M′` with its `A`-trace replaced by `M`'s — a member, since `P`
factors — and then `M₂` as `M′` with its `B`-trace replaced by `M₁`'s — a
member, since `Q` factors. `M₂` is `M′` everywhere except `A ∩ B`, where it
is `M`. Factoring partitions are thus closed under common refinement, and
the refinement of all of them is the finest. Sizes: were two traces of one
block of different sizes, combining every block's smallest trace would give
a member smaller than `k`. ∎

The finest factorization is the canonical presentation *when it is fine*:
menus are its blocks, entries its traces, compound exactly where a trace is
not a singleton, and the plain product is the case where every trace is.
Measured on entangled registry queries, however, the families are bimodal:
either the plain product holds, or the family is **prime** — no nontrivial
factoring partition at all — and the finest factorization degenerates to
the one-menu list of whole repairs. So structure, where it exists, is found
by Theorem 23; where it does not, it is found inside `F_min` a piece at a
time, and coverage is completed by layers.

### Covering the cheapest fixes by products

A **generalized rectangle** in `F_min` is a presentation `M₁, …, M_j` over
pairwise disjoint groups of facts — the group of `Mᵢ` being the union of its
entries, an entry a set of facts — such that **every** selection's pick-set
is a member of `F_min`. It is *exact*: it offers nothing the family does not
hold, so (C1) and (C2) come for free and only (C3) is left to answer. Its
**coverage** is the number of members its selections reach, which is `∏ᵢ|Mᵢ|`
exactly, since the groups are disjoint and a set determines its elements.
Three things follow immediately from the definition. Every member `M` alone
is a rectangle, of coverage one — one menu, one entry. Any subfamily
`S ⊆ F_min` that factors in the sense of Theorem 23 is a rectangle, its menus
the blocks' traces. And `F_min` itself is a rectangle of one menu whose
entries are whole members, which is the flat list.

A **product cover** of `F_min` is a sequence of generalized rectangles — the
**layers** — whose selections' pick-sets are together exactly `F_min`. Build
one greedily: take a rectangle, emit it, delete the members it reaches, and
recurse on what is left, which is again a family of size-`k` fact sets over
the same facts.

**Theorem 29 (a product cover exists, is finite, and presents the family
exactly).** The construction terminates, and its layers present `F_min`: no
selection of any layer is outside `F_min`, and no member of `F_min` is the
selection of no layer.

*Proof.* Exactness is the definition of a rectangle, and a rectangle in the
family that is left is a rectangle in `F_min`, since what is left is a subset
of it. For termination, each layer reaches at least one member — a single
member is always available as a rectangle — so the remainder strictly
shrinks and the construction stops with it empty. Completeness is then the
stopping condition: a member never reached would still be there. ∎

Which rectangle to take is an explainability preference among sound forms,
not a soundness question — Section 5 needs (C1) alone — and the preference
is what a layer costs to read against how much of the family it saves the
reader from reading: a line for each menu, a line for each entry under it,
and a phrase for each action an entry names. That price is what keeps the
flat list from winning by default: it reaches everything, and it costs an
entry, and every action of that entry, per member.

The first layer is constrained beyond that, and Corollary 25 says why: it is
the one the proofs answer for. It is required to be a product of two or more
menus wherever one exists, so the reader is shown the family's structure and
not a list of it.

**So a conflict is a menu of the leading layer.** The leading layer satisfies
(C1), so by Theorem 9 each of its menus owns a reason, and by Theorem 10 one of
its very own; the layers after it own none (Corollary 25). A conflict is
whatever a page can argue for, and that is exactly a leading menu: it gets a
heading, the reason it owns as its chain, its own menu of fixes, and the
verdicts on the actions those make tempting. A menu of compound entries earns a
further owned reason wherever the first leaves a package the menu offers
unspoken of — a page that asks the reader to change something its argument
never mentioned is talking past them, which is what V3 of Section 8 checks.

The block's later layers are then not conflicts but **alternatives** to all of
its conflicts at once, printed after them. Each is a product of menus like any
layer, and every one of its selections is a cheapest repair the block's
conflicts do not offer between them. What an alternative replaces is its own
block and nothing else: by Theorem 30 the blocks are independent, so a repair
of the query takes one selection from each block, and where the reader takes an
alternative for one block the other blocks stand exactly as printed. That is
the whole of what "without any of the fixes for Conflicts 1–3" means, and the
label carries it.

**Lemma 31 (an unreached member avoids a whole menu).** Let `M ∈ F_min` be no
selection of the leading layer. Then `M` is disjoint from some menu of that
layer.

*Proof.* Suppose `M` meets every menu, say in the entry `eᵢ ⊆ M` of `Mᵢ` — a
menu's entries need not be disjoint from one another, so fix one such `eᵢ` per
menu. The selection `S = ⋃ᵢ eᵢ` is a member of `F_min` by exactness, and
`S ⊆ M`. Members of `F_min` are minimal correction sets, hence pairwise
incomparable, so `S = M` — and `M` is a selection after all. ∎

So every selection of a later layer declines a whole menu, and what the label
may say the layer as a whole is instead of is the menus the layer shares no
fact with at all — those its every selection declines. Where each of them
offers a single fix, that fix is a single thing the reader is doing without and
the label says it; where one of them offers a choice there is no single thing
to name, and the label names those conflicts instead. Lemma 31 is per member,
not per layer, so the set can be empty — the layer's selections declining
different menus and none of them all — and then the label says only that the
layer is another way.

**Lemma 24 (what a layer leaves is a blocked instance's cheapest repairs).**
Add to the instance, for each member `M` a layer reaches, the clause that
some fact of `M` holds. The blocked instance's minimal correction sets of
size `k` are exactly the members left.

*Proof.* As in Theorem 5: a model of the blocked instance violates a set
containing no reached member, and if that set has size `k` it is a minimal
correction set of the original instance — a member of `F_min` — that is not
reached; conversely a left member's witnessing model satisfies every
blocking clause, since distinct members of `F_min` are incomparable. ∎

So what is left is defined by the family and not by any choice of fixes, and
it is itself a repair family over the same facts — which is what licenses the
recursion. Covers are not unique and least covers need not be; the order of
construction is a stated preference, like the pivot tie-break of Section 6.

**Corollary 25 (reasons do not layer).** The first layer satisfies (C1), so
by Theorem 9 every reason meets every entry of one of *its* menus. The
reasons — and with them the explanations of Section 6 — attach entirely to
the first layer; later layers repartition repairs, not reasons, and owe the
reader menus and witnesses, never proofs. (The blocked instance of Lemma 24
does have minimal unsatisfiable sets of its own, but they mix the query's
facts with the refusal clauses — "given that none of the fixes above is
taken" — and a proof from refusals is contingent on them. The page never
prints one.)

**Theorem 30 (factoring composes, so nothing is left across factors).** Let
`P = {B₁, …, B_r}` be a factoring partition of `F_min` (Theorem 23) and give
each block `Bₜ` a product cover of its own trace family `F_min|Bₜ`. Reading
each block as one conflict, and a repair of the whole query as one selection
from each, the presentation is exact and complete: the selections are
exactly `F_min`.

*Proof.* `P` factoring says `F_min` is the set of unions of one trace per
block, every combination realized. Each block's cover presents `F_min|Bₜ`
exactly (Theorem 29), so the selections available at block `t` are exactly
its traces; the combinations of them are therefore exactly the unions of one
trace per block, which is `F_min`. ∎

Two consequences. Completeness of the whole presentation is the *product* of
the factors' completeness, so where the finest factorization is nontrivial
there is no cross-factor remainder to print: a page that covers each factor
covers the query, and a residue enumerating combinations across factors is
enumerating repairs the factors have each already offered. And a block that
is prime is covered on its own facts alone: the search runs over that block's
trace family, not over the whole of `F_min`.

**What a one-entry menu may say about itself.** With `F_min` covered
(Theorems 29 and 30) and Theorem 5 answered, the gap a menu's wording must
confess has two sources: costlier fixes, and an enumeration cut short (the
one way coverage can silently fail — and one further blocked solve at the cap
decides even that, so it is known, not guessed). "Only" is a claim about the
world:

| what exists beyond this entry | a menu of one says |
|---|---|
| nothing | *The only fix:* |
| only costlier fixes | *The only minimal fix:* |
| an alternative to this conflict's block, or a cut enumeration | *One fix:* |

Never derive "only" from the length of a vector; derive it from the decided
questions. The alternative row is block-local: an alternative to some *other*
block is another repair of the query and no repair of this conflict, so it
takes nothing away from what this entry may say about itself, while an
alternative to this conflict's own block is a cheapest repair that settles the
block without this entry and denies "only" outright.

## 5. Ownership: which reasons a conflict answers for

Conflict `i` **owns** a reason `R` when `Dᵢ ⊆ R`. Its page explains the
reasons it owns. Everything below needs only contract clause (C1).

**Theorem 9 (ownership).** In any presentation where every selection's
pick-set is a repair, every reason contains some whole menu.

*Proof.* Suppose reason `R` misses an element of every menu: pick
`dᵢ ∈ Dᵢ ∖ R` for each `i`. The pick-set `S = {d₁, …, d_m}` is disjoint from
`R`, so after withdrawing `S` every fact of `R` is still held, and `R` is
unsatisfiable — so the withdrawal is unsatisfiable, contradicting (C1). ∎

Note what the hypothesis is *not*: no product, no equal sizes, no minimality
of the pick-sets. Any presentation meeting the contract — product, rectangle,
or anything future — inherits the theorem.

**Theorem 10 (privacy).** For every conflict `i` and every transversal `T`
of the other menus (one entry from each `Dⱼ`, `j ≠ i`), there is a reason
disjoint from `T` that contains the whole of `Dᵢ`; and any reason disjoint
from such a `T` is owned by conflict `i` alone.

*Proof.* Menus are disjoint, so `|T| = m − 1 = k − 1 < k`, the minimum
repair size; hence `T` is not a correction set and the withdrawal of `T` is
unsatisfiable. Take any minimal unsatisfiable subset `R` of `Q ∖ T` — a
genuine reason, since minimality and unsatisfiability of a set do not depend
on the pool it was found in. By Theorem 9, `R` contains some whole menu; for
each `j ≠ i` it misses `T`'s entry of `Dⱼ`, so the menu it contains is `Dᵢ`,
and it contains no other whole menu. ∎

**Corollary 11.** Every conflict owns at least one reason, and indeed one of
its very own — containing no other conflict's whole menu — for *every* way
the other conflicts might be settled.

**Corollary 12 (what is left is settled).** Fix a transversal `T` of the other
menus. Every reason of the query that avoids `T` — in particular, every
reason of the `T`-withdrawn query — is owned by conflict `i`.

*Proof.* As in Theorem 10: it contains a whole menu and misses every other
menu's `T`-entry. ∎

Corollary 12 is why no further machinery is needed for requirements that
fail "for the same reason" as one another: with the other conflicts settled,
*everything* still wrong is this conflict's to explain, and each still-broken
requirement sits in its own reason with this conflict's facts. Explaining one
owned reason per page section, every such requirement gets its account.

**Ownership is not a partition.** It is tempting to read Theorem 9 as "the
reasons divide among the conflicts". They do not:

> **Example (a reason owned twice).** Facts `a₁ a₂ b₁ b₂ x y z v`; reasons
> exactly
>
> ```
> {a₁ a₂ b₁ b₂}   {a₁ a₂ x}   {a₁ a₂ y}   {b₁ b₂ z}   {b₁ b₂ v}
> ```
>
> (realizable: give `ℛ` one clause per line forbidding that conjunction —
> then the unsatisfiable subsets are exactly the supersets of these five).
> A hitting set must cover `{a₁a₂x}` and `{a₁a₂y}` — with one element only
> via `a₁` or `a₂` — and `{b₁b₂z}`, `{b₁b₂v}` likewise via `b₁` or `b₂`;
> both together also cover `{a₁a₂b₁b₂}`. Any pair using `x, y, z, v` misses
> one of the five. So `F_min = {a₁,a₂} × {b₁,b₂}` exactly: a perfect
> product, two conflicts. Yet `{a₁ a₂ b₁ b₂}` is a genuine minimal reason
> containing **both** menus entire.

So a reason may be owned by several conflicts, and either conflict's action
settles a shared one. The sound policy is to print a shared reason under
every conflict that owns it — redundant, never wrong — and Corollary 11
guarantees the redundancy is never load-bearing: each conflict also has
private grounds. (A deduplicating presentation is possible but buys little;
measure the frequency of sharing before spending design on it.)

**Costlier fixes and second reasons are the same phenomenon.** Two results
pin down when repairs beyond the cheapest exist, and they meet exactly at
the blocked-fix entries of Section 9.

**Theorem 26 (a second owned reason forces a costlier fix).** In the
product case, if a conflict owns two reasons, then `F ≠ F_min`.

*Proof.* Let conflict `i` own `R₁ ≠ R₂` and let `T` be a transversal of the
other menus. Every reason not owned by `i` contains some other menu whole
and so contains its `T`-pick. Every reason owned by `i` contains `Dᵢ`
properly — were `Dᵢ` itself a reason, the antichain property would leave it
the only owned one — so the owned reasons' parts outside `Dᵢ` are nonempty,
and some finite set `S` of facts outside `Dᵢ` meets them all. Then `T ∪ S`
is a correction set avoiding `Dᵢ` entirely, and the minimal correction set
inside it is no selection; since `F_min` is exactly the selections
(Theorem 8), it is not of size `k`, so it is strictly larger. ∎

**Theorem 27 (a costlier fix forces a second owned reason).** If some
minimal correction set has size greater than `k`, then some conflict owns
two reasons.

*Proof.* Let `N` be minimal with `|N| > k`. Minimality gives every `e ∈ N`
a witness reason meeting `N` in `{e}` alone, and distinct elements get
distinct witnesses — so more than `k` reasons exist. Each contains a whole
menu (Theorem 9) and there are `k` menus, so some menu lies inside two of
them. ∎

The consequence for the page: a blocked entry of the *unless* form is a
costlier fix, named (Section 7) — so where one printed, the footer that
would announce costlier fixes in the abstract says nothing the page has not
already said better, and is omitted. Where only idle verdicts printed, or
none, the footer speaks for what remains unshown.

**Finding the owned reasons.** Reasons are enumerated by removal, over the
full fact set — never over a restricted pool, since a reason can contain this
conflict's menu *and* facts other menus offer, and a pool with those held out
cannot contain it:

**Theorem 13 (the removal walk is complete).** Define
`walk(Π)`: if `Π` is satisfiable, stop; otherwise take any one minimal
unsatisfiable `m ⊆ Π`, record it, and recurse on `Π ∖ {x}` for each
`x ∈ m`. Started at `Π = Q`, the walk records every reason.

*Proof.* Induction on `|Π|`: every reason `R′ ⊆ Π` is recorded by
`walk(Π)`. If the recorded `m` is `R′`, done. Otherwise `m` and `R′` are
distinct minimal sets, so neither contains the other, and some `x ∈ m` has
`x ∉ R′`; then `R′ ⊆ Π ∖ {x}`, and the recursion records it by the induction
hypothesis. ∎

The walk is exponential in the worst case and must be run with a budget.
Truncation is sound — nothing below depends on the enumeration having
finished; a conflict simply explains fewer of the reasons it owns — but it is
**recorded**: when the budget stops the walk the diagnosis says so, though
the page does not (Section 11).

## 6. The explanation of a reason

Each owned reason `R` gets one explanation. Its shape is fixed, and this
section proves the shape always exists and always carries what the fixes
need.

### Blame before proof

The facts of `R` were determined against the registry *as it is* — `ℛ`
inviolable, only `Q` withdrawable. This order matters: ask which facts are
at fault while the registry's statements are also withdrawable, and a
minimal set is free to drop a registry statement and blame a requirement in
its place — "you require C" is no part of why a query fails when something
else requires C anyway. Minimal, and a lie about what the user can do. So:
blame first, on the real instance; proof second, with the blame held.

The proof half asks the same kind of question one level down. Fix `R`; a
**core** for `R` is a minimal `C_R ⊆ ℛ` with `C_R ∪ R` unsatisfiable. (To
extract one, the registry's statements must be individually assumable — see
Substrate obligations. Note a core may be *empty*: a requirement whose
package the query's own constraint leaves nothing of contradicts it with no
help from the registry. Empty is a value here, not a failure.)

**Lemma 14 (the combined set is minimal).** `Γ = C_R ∪ R` is a minimal
unsatisfiable clause set.

*Proof.* Unsatisfiable by construction. Drop a registry clause: satisfiable
by minimality of `C_R` given `R`. Drop a fact `f ∈ R`: were
`C_R ∪ (R ∖ {f})` unsatisfiable then, since `C_R ⊆ ℛ`, the set
`ℛ ∪ (R ∖ {f})` would be too, contradicting the minimality of `R` as a
reason. ∎

### The pivot

**Theorem 15 (the pivot theorem).** Let `Γ` be a minimal unsatisfiable
clause set and `P` a package in the support of some clause of `Γ`. Then
`proj_{{P}}(Γ)` is a set of **bounds** `P ∈ S₁, …, P ∈ Sₙ` — unit clauses,
each entailed by the subset of `Γ` it was derived from — with every `Sᵢ`
nonempty and `⋂ᵢ Sᵢ = ∅`.

*Proof.* By Theorem 2 (applied per elimination) the projection is
unsatisfiable, and every clause of it has support inside `{P}` — a unit
bound, or the empty clause. We rule the empty clause out. Since `P` is never
eliminated, every resolution step unions the participants' `P`-literals, so
by induction any derived clause `d` has `d(P) ⊇ C(P)` for every original
clause `C` in its derivation. A derived empty clause `e` has `e(P) = ∅`, so
its derivation uses only clauses with empty `P`-literal — a subset of `Γ`
omitting the clause that mentions `P`, hence proper. By Lemma 1 that proper
subset entails the empty clause, i.e. is unsatisfiable — contradicting
minimality of `Γ`. So every projected clause is a bound with `Sᵢ ≠ ∅`; a set
of nonempty unit bounds is unsatisfiable exactly when no value lies in all
of them, i.e. `⋂ Sᵢ = ∅`. ∎

**Corollary 16.** (a) The choice of pivot is free: *any* package the
argument touches supports the whole contradiction. (b) At least two bounds
are always needed (`⋂` of one nonempty set is itself). (c) At least one
bound excludes `⊥` — otherwise `⊥` lies in the intersection. Read on the
page: some line must *force* — a meet of sides that all admit "not
installed" is no contradiction, and a report all of whose registry sides
admit `⊥` has elided the registry's half of the argument.

The **explanation of `R` at pivot `P`** is: the facts of `R`, printed as the
query's own lines, and a subfamily of the projection's bounds with empty
intersection — the **sides** of the **meet** — each side carrying its
support and route (below). Which subfamily: any with empty intersection is
sound (each side is entailed; a set of true bounds with empty meet is a
contradiction); prefer the smallest, then the preferences below. Theorem 18
will show no such minimisation can lose anything that matters.

### What the meet can look like

Call a bound's **version-part** `Sᵢ ∖ {⊥}`, and call it **convex** when it
is an interval of the version order. Registry compatibility bounds are
usually convex; disconnected version-parts arise from multi-range compat
entries, and from composition — a side derived along a chain can come out
with a hole even from convex inputs.

**Lemma 17 (Helly, on a line).** If finitely many intervals of a total order
pairwise intersect, they all share a point.

*Proof.* Let `m` be the greatest left endpoint, attained by interval `I`,
and `M` the least right endpoint, attained by `J`. `I ∩ J ≠ ∅` forces
`m ≤ M`, and every interval contains `[m, M]`. ∎

Say a family of bounds is **irredundant** when its intersection is empty but
every proper subfamily's is not — the meets a report prints, since a
redundant side would be dropped.

**Proposition 18 (the shape of convex meets).** In an irredundant family of
bounds with every version-part convex: (a) there are at most three bounds;
(b) if there are three, exactly one excludes `⊥`, and the other two overlap
only at `⊥` — their version-parts are disjoint.

*Proof.* At least one bound excludes `⊥`, else `⊥` is in the full meet.

*At most one excludes `⊥`, given three or more bounds.* Suppose `S₁, S₂`
both exclude `⊥` in an irredundant family of `n ≥ 3`. For any two indices
`a, b` choose `j ∉ {a, b}`; the subfamily dropping `j` has a common element
`x_j` (irredundancy), and since the subfamily retains `S₁` or `S₂` (it drops
only one set), `x_j ≠ ⊥`; so `x_j` is a common *version* of all
version-parts but possibly `T_j`, in particular `x_j ∈ T_a ∩ T_b`. So the
version-parts pairwise intersect; by Lemma 17 they share a version, which
then lies in every `Sᵢ` — contradicting the empty meet.

(b) With `n = 3` and (by the above) exactly one `⊥`-excluding bound `S₁`:
if the other two version-parts met at `x`, then — irredundancy giving
`S₁∩S₂ ≠ ∅` and `S₁∩S₃ ≠ ∅`, both `⊥`-free since `⊥ ∉ S₁` — all three
version-parts pairwise intersect, and Lemma 17 again contradicts the empty
meet. So `T₂ ∩ T₃ = ∅` and `S₂ ∩ S₃ = {⊥}`.

(a) Let `n ≥ 4`. By the second paragraph at most one bound excludes `⊥`;
exactly one, say `S₁`. Drop any `⊥`-admitting `Sⱼ`: the remaining family
contains `S₁` and at least two other `⊥`-admitting bounds, and has a common
element `x_j ≠ ⊥`. As before, for any pair `a, b ≥ 2` choose `j ∉ {a,b}`:
`x_j ∈ T_a ∩ T_b`, and `x_j ∈ T₁` too. All version-parts pairwise intersect;
Lemma 17 contradicts the empty meet. ∎

So in a convex world a meet has two sides, or the one three-sided shape:
*something must be installed, and the two constraints on it admit no common
version* — and any meet of four or more, or a three-sided meet of any other
shape, is possible only through a disconnected version-part. That is why
wide meets are rare, and where to look for them: only where a compat entry
or a derived side has a hole. (Empirically, in the General registry, a few
hundred of tens of thousands of bounded compat entries are disconnected, and
genuine same-shape triples are countable on one hand. This is a remark, not
a premise; nothing above depends on it.)

This is also what bounds how often a meet has to be *displayed* as one:
two sides always read out as a chain and three need not (Proposition 22).

### Sources: what a side rests on

Say `X ⊨ σ` ("`X` entails `σ`") when every installation satisfying every
clause of `X` satisfies `σ`. For a side `σ` of the meet:

* a **support** of `σ` is an `A ⊆ R` with `ℛ ∪ A ⊨ σ`; minimal supports
  exist by finiteness (and `A = ∅` is possible: a side the registry entails
  alone, e.g. "the registry offers only these versions of `P`");
* `σ` is **rooted in** fact `f` when `ℛ ∪ (R ∖ {f}) ⊭ σ` — no support
  avoids `f`.

Rootedness is a semantic notion, deliberately: a derivation that happened to
use `f` is a *witness* that may over-attribute (the side might also be
entailed without `f`), and implementations may track derivations as an
approximation, but the theorems are about entailment, and entailment is one
oracle query per check.

**Theorem 19 (source coverage — and no minimisation can break it).** Let `Π`
be **any** family of bounds on `P`, each entailed by `ℛ ∪ R`, with
`⋂Π = ∅`. Then for every fact `f ∈ R`, some member of `Π` is rooted in `f`.

*Proof.* `R` is minimal, so `ℛ ∪ (R ∖ {f})` is satisfiable; let `ι` model
it. Every side entailed by `ℛ ∪ (R ∖ {f})` is satisfied by `ι`, and a unit
bound is satisfied only at `P`: `ι(P)` lies in each such side's set. If
every member of `Π` were entailed without `f`, then `ι(P) ∈ ⋂Π`,
contradicting `⋂Π = ∅`. ∎

This holds for *whatever* subfamily the presentation chose to print — the
theorem quantifies over all of them — so trimming the meet to its smallest
irredundant form cannot lose a source. Every fact of the reason is visible in
every printed meet, as the root of at least one side. In particular the
verification of source coverage is a sanity assertion, not a gatekeeper.

A side with a minimal support of size one is **separated**: fully attributed
to a single fact ("requiring `A` forces `P ∈ S`, through …"). A side whose
minimal supports all have several facts is **conditional**, and is printed as
the implication from all of them ("with `A` at its forced version, `B`'s
chain forces `P ∈ S`") — honest, attributed, and always available, since
every side has some minimal support. So the form is total: separated where
possible, conditional where not. A **separating pivot** is one where every
printed side is separated; one need not exist, which is why separation is a
preference and not an obligation.

### Choosing the pivot

Every package the core touches is a lawful pivot (Theorem 15), so the choice
is made for the reader, by ranking what each candidate would print:

1. **(mandatory)** the meet is printed whole — its sides listed together and
   their emptiness visible on its face — and includes a forcing side
   (Corollary 16(c) guarantees one exists);
2. prefer pivots where every side rooted in a **menu fact** of this conflict
   is separated — the fixes are the sides the reader will trace, and a fix
   whose only side is conditional reads as a mechanism blurred;
3. then pivots where every side is separated;
4. then fewest sides;
5. then least version fracture and shortest routes (below).

Ranks 2–5 are preferences among sound pages; nothing below depends on how
ties break, only that they break deterministically.

Two designs this replaces, and why, for the implementor tempted to
rediscover them. *Searching for a chain* — walking from the requirements,
narrowing each package by everything so far, printing one registry edge per
link — cannot be made valid: each link's stated range is narrowed by
packages that are no part of any minimal reason, so either the page states a
narrowing no printed line accounts for, or it attributes one package's bound
to another. The projection has no such gap, because elimination resolves
every minimal emptying set. *Deriving first and choosing afterwards* —
keeping a derivation graph and picking a best node out of it — confines the
choice to statements the derivation happened to make and prices them by the
work done rather than the page printed. Choosing the pivot first inverts
both: the vocabulary is known before anything is derived, and the cost *is*
the page.

### Routes

A side's **route** is a set of packages whose elimination derived it — a
courtesy pointer ("through X and Y") telling the reader where in the registry
the entailment lives. The route makes no claim: the side's truth is its
supported entailment (Verification, V1), and routes carry no verification
obligation of their own. Prefer short routes; order a route from the side's
source toward the pivot; and where a route is long, that is the registry's
length, not the presentation's — say it as a route rather than spending a
line per hop.

### Licensed coarsening

Real cores carry **parallel families**: many statements about one set of
packages, differing only in version thresholds — a staircase of couplings the
argument may never need. Left as they are, the projection pays for their
combinations; the report, if it survives, pays in lines.

The remedy is the resolution rule read as a merge. For a parallel family
`{C₁, …, Cₙ}` over one package set, resolving on a package `q` of it —
intersect at `q`, union everywhere else — yields one clause entailed by the
family (Lemma 1), so printing it in their place keeps every line registry-true
and inherits the union of their supports and routes. On the antecedent package
this is exactly the join of implications: `A@R₁ → B@S₁` with `A@R₂ → B@S₂`
gives `A@(R₁∪R₂) → B@(S₁∪S₂)`, the antecedent literals being stored
complemented.

The join is weaker than the family, so it needs a licence, and the licence is
the whole criterion: **the claim's clause set must still contradict with the
join in the family's place** — one satisfiability check, against that
reason's facts and the rest of its core, per claim and never against a union
of claims. No local rule can stand in for it: whether two thresholds may be
treated as one depends on what happens links away, and a boundary the
contradiction stands on must refuse to join while its neighbours collapse.
Bisecting a family that refuses to join whole gives up exactly the boundaries
the proof is standing on.

Two disciplines make the pass sound as a whole. Families are coarsened
sequentially, each licence asked against the set *as it stands* — earlier
joins included — so the invariant after every accepted join is that the whole
current set still contradicts; licensed against the original set instead, two
joins could each pass and jointly satisfy. And coarsening is not confluent, so
the order of attempts is a stated preference, not an accident.

Where it runs: when no pivot's elimination finishes within budget, the core is
coarsened and the projection tried once more — which is what makes a lockstep
family of nine parallel edges projectable at all — and if the projection still
fails, the coarsened core is the fallback, every join of it true and the set
still contradicting.

## 7. Fixes on the page

The two halves now connect, with no machinery at the joint.

**Theorem 20 (fix visibility).** In any presentation meeting contract (C1),
every entry of every menu is a fact of every reason its conflict owns; hence
in every explanation printed under that conflict — at any pivot, after any
minimisation — the entry roots at least one side; and the sides *not* rooted
in it have nonempty intersection.

*Proof.* An owned reason contains the whole menu (Theorem 9), so the entry's
fact `f` is in `R`; Theorem 19 gives a side of the printed meet rooted in
`f`. For the last part: the sides not rooted in `f` are each entailed by
`ℛ ∪ (R ∖ {f})`, and the model `ι` from Theorem 19's proof puts `ι(P)` in
all of them. ∎

That is the reading a report should support typographically: each fix names a
fact; that fact is the root of at least one visible side; cover those sides
with a thumb and the remaining bounds visibly agree on a version (or on
`⊥`). Why each fix works is a corollary of the form, not a feature to build.

**Theorem 21 (witness coherence).** Let a fix be taken: entry `e` of menu
`Dᵢ`, together with one chosen entry of every other menu — a selection, whose
pick-set `W` is a repair by (C1). Let `ι` be any solution of the withdrawn
query (one exists), e.g. the resolved witness the report shows. Then for
every explanation printed under any conflict, `ι(P)` lies in every side
entailed by `ℛ ∪ (R ∖ W)`; in particular the witness lands inside the meet
opened by the withdrawal.

*Proof.* `ι` models `ℛ ∪ (Q ∖ W)`, and `R ∖ W ⊆ Q ∖ W`, so `ι` satisfies
every side entailed by `ℛ ∪ (R ∖ W)`; a unit bound on `P` is satisfied only
at `P`. ∎

Note the theorem is stated against the *full* withdrawal `W`, not the single
entry `e` — necessarily so: an owned reason can contain other conflicts'
facts (shared reasons contain their whole menus; ordinary ones may contain
other menus' unpicked entries), and the witness respects only the sides whose
supports survive everything withdrawn. The display rule that goes with it:
the versions shown for a fix (**"allows:"**) always include the pivots of the
conflict's printed meets, so the reader sees the witness land where the
opened meet says it can.

**Versions come from a resolve.** No version in a report — witness or
otherwise — is ever read off the failed instance or any model found during
diagnosis. A model witnesses satisfiability; which versions the user would
*get* is the resolver's optimising answer for the withdrawn query, computed
by the substrate (S4). Diagnosis decides what is true; the resolver decides
what is chosen.

### Blocked fixes: which, and how

The page's own lines make some actions tempting: the reader sees "your
compat leaves QuantumLattices 0.15.4" and asks why relaxing it is not
offered. A **tempting action** is an action on a package a conflict's lines
or heading name that appears nowhere in the cover. The blocked-fixes
section answers for exactly these — one sentence per action, so nothing is
said twice and nothing tempting goes unanswered.

What the sentence says is decided by solves, not judged, and the question
must be *well posed*: `x` is **load-bearing** when it lies in **some**
minimal repair, **idle** when it lies in none. Asking of one repair a solver
happened to return is not enough — a tie between equal repairs could call
the same action idle in one sentence and a rescue in the next, which is
exactly the inconsistency this phrasing exists to avoid. So the search
favours `x`: take a repair `W ∋ x` within the bound, shrink its *other*
members while it still repairs, and ask whether `x` too can go. If it
cannot, the shrunk `W` is a minimal repair carrying `x` and `x` is
load-bearing; the sentence names its price, *"«x» would not help unless you
also «W ∖ {x}»."* — in full while the price is a few actions, and as a count,
*"without «n» other changes"*, beyond that: a long completion is a verdict on
the road and not a list anyone will follow, and the exhibit behind it is
checked either way (Section 8). If `x` can go, that witnesses a smaller repair without it —
block `W` and try again, a few rounds, before concluding **idle** and
printing the flat *"«x» does not help."* Because the cover holds every
cheapest repair, a tempting `x` is in none of them, so a load-bearing `x`'s
minimal repair is always strictly costlier than `k` (Lemma 28), and idle
means dead weight in every repair within reach.

**One repair, said once.** Two tempting actions can be load-bearing in the
*same* minimal repair — each insufficient alone, the repair carrying both.
Mirroring it from each end ("relax A unless you also drop B" beside "drop B
unless you also relax A") says one thing twice and reads as a contradiction.
So entries are grouped by the repair they exhibit: the actions that share
one are one entry, said as the choice they are not: *"«x₁» or «x₂» would
only help if you do both."* An
entry of this form still exhibits a costlier fix, so it suppresses the
footer exactly as an *unless*-entry does.

**Lemma 28 (completions never restate the cover).** A load-bearing
completion `W ∖ {x}` contains no printed fix.

*Proof.* A correction set's supersets are correction sets, so if
`W ∖ {x}` contained a fix it would itself repair — and then `x` is idle in
`W`, not load-bearing. ∎

So the unless-sentences and the cover cannot duplicate one another: the
cover lists the cheapest fixes, and each unless-sentence exhibits a
strictly costlier one in which its action genuinely participates. That
sharpens the footer rule of Section 5: an *unless*-entry names a costlier
fix, so where one printed, the costlier-fixes footer says nothing new and
is omitted; idle verdicts name none, and do not suppress it.

None of this consults the reason walk: blocked fixes are questions about
repairs, and the completed cover plus one bounded search per tempting
action decides them. The walk's remaining duty is the conflicts' own
stories.

## 8. Verification

What a checker checks, and what it need not.

* **(V1) Side truth.** For every printed side `σ` with printed support `A`:
  `ℛ ∪ A ⊨ σ`. One entailment query each (equivalently: `ℛ ∪ A` plus the
  negation of `σ` — the assertion that `P` takes a value outside `σ`'s set —
  is unsatisfiable).
* **(V2) Visible closure.** The printed sides of each meet intersect
  emptily. Pure set arithmetic on the page; no solver.
* **(V3) Source coverage.** Every fact of the reason roots a printed side.
  Guaranteed by Theorem 19 once V1 and V2 hold; assert it as a sanity check
  on the rootedness bookkeeping, not as a gate.
* **(V4) Cover contract.** Every selection's pick-set is in `F_min`, and
  every member of `F_min` is some selection. What the page offers is defined
  over blocks: a block's own repairs are one entry from each of its conflicts'
  menus in every combination, together with every selection of each
  alternative to it, and a repair of the query is one of each block's, in
  every combination. By construction (Theorem 6's check; a generalized
  rectangle's selections are members by definition; Theorems 29 and 30 for the
  layers and the blocks) — and asserted per report, since it is a set
  comparison over at most `|F_min|` members and costs nothing.
* **(V5) Witness coherence.** For each fix's witness and each printed meet:
  the witness's pivot value lies in every side whose printed support avoids
  the withdrawn facts. Membership tests; no solver. Silent breakage here is
  invisible to every other check, which is exactly why this one exists.
* **(V6) Disclosures.** The menu wording matches the decided questions
  (Section 4's table), its alternative row read block by block; the
  alternatives appear exactly when some block's cover has more than one layer,
  one per layer after the leading one, and their entries' witnesses are checked
  like any entry's (V5); the headline's *", pick a fix for each"* appears
  exactly when there is more than one conflict; the enumeration-cut
  sentence appears exactly when the repair enumeration was truncated and the
  deciding solve found more; the costlier-fixes footer appears exactly when
  Theorem 5 answered yes and no unless-entry printed. A reason walk cut
  short is recorded on the diagnosis and not announced on the page (Section
  11).

**Every check is per-explanation, never against a union.** Where two
explanations' line-sets are `S₁ ∪ S₂` and `S₂` alone is contradictory, the
union stays contradictory whatever is deleted from `S₁` — a union-level
check would pass a page that silently destroyed the whole account of `S₁`.
This is the one failure mode in this design that is silent by nature; the
rule exists to make it impossible rather than detectable.

## 9. Rendering

Rules, not preferences — each guards a truth-condition or an attribution:

* **A side is stated from its support.** The requirements the conflict
  answers for are said by its heading and nowhere else; the reason's other
  facts print once, as the query's own lines ("your compat leaves A 1.2").
  Each side prints as the implication from its support's packages to its
  bound, with its route in parentheses. A conditional side names all of its
  support ("A 1.2 and B ≥ 2 together …"). This is what makes every fix
  traceable: the fact a fix withdraws is the stated antecedent of a visible
  side.
* **The user's constraints are named as the user's.** Which constraint kinds
  removed which versions is read directly off the query (S5) — no solver, no
  license needed — and only query lines may say "your"; everything else is
  the registry's, and a registry statement never attributes a bound to a
  `Project.toml` it cannot see.
* **The verb is read off the printed consequent** — see the verb table
  below. Never carried beside the clause: a line said the other way round
  re-derives its own verb, which is what makes contraposition and
  composition safe to print.
* **Direction is chosen at print time.** A clause has no direction; state a
  two-package statement from the side that declares the dependency edge
  where exactly one does, and symmetrically ("X and Y are incompatible")
  where neither. Among otherwise equal readings prefer the one whose stated
  consequent has versions to name: a consequent with an empty range prints
  shortest and says least, and choosing it throws away the very bound the
  reader was about to hold against another. A line that *continues* a chain
  overrides all of this: it is said to the package the chain has not
  reached, whichever that is.
* **A range is printed when something narrowed it.** An antecedent that is a
  package's entire offering prints as the bare package name; naming the full
  range reads as a narrowing that never happened.
* **Reading order: a proof is one chain.** A proof prints as a
  polysyllogism — a root fact, then implications each arguing from what the
  line above it left, ending at the fact that contradicts the accumulated
  bound. The root is the heading's own subject, said as the user's compat
  where the reason narrowed it and taken as the heading's premise where it
  did not; a subject the statements only *conclude* about is where the chain
  ends instead, so a source is preferred to it. A line may argue only from
  packages the query named or a line above introduced, and the query's fact
  about a package prints where the package arrives: after the line that
  reaches it, before the line that argues from it. Prefer to close at a
  printed compat line where one contradicts the chain; otherwise the last
  line closes against the heading's implicit *given the rest of the
  requirements*, and the package it leaves nothing of is named in no further
  line. A line continuing through a package reads its antecedent against the
  bound the chain left there — widening its literal by everything already
  ruled out weakens the clause, so the line stays true and says what the
  reader is holding. Only a meet of three or more sides breaks the chain
  (Proposition 22), and that alone prints whole, its sides adjacent.
* **The heading names the primary reason's requirements, and claims
  nothing.** "Conflict 1: QuantumLattices" — the colon form says these
  packages are the conflict's core and stops. A sentence would either claim
  too much ("cannot be satisfied" is false absolutely, true only under the
  unstated *given the rest*) or spell the context out at absurd length; the
  bare form is as clear, shorter, and never wrong. The one sentence heading
  that survives is the absolute truth: "no version of X is available." Where
  two conflicts share a root requirement their headings would read as one
  conflict said twice, so each is extended with the package it contradicts
  at — the package its own lines name and the others' do not, taken where
  the chain closes — which leaves the heading a true bare list and stops it
  reading as a duplicate. The checkers still premise every requirement the
  conflict answers for: what shrinks is the sentence, not the claim set.
* **Blocked fixes answer for actions, one sentence each.** For every
  tempting action — named by the conflict's lines or heading, in no fix of
  the cover — the section prints its solve-decided verdict (Section 7): an
  idle action gets *"«x» does not help."*, a load-bearing one gets
  *"«x» would not help unless you also «W ∖ {x}»."*, or *"«x» would not help
  without «n» other changes."* once the completion is more than a few
  actions long. Action-indexed, so nothing is said twice; solver-licensed,
  so nothing is judged; and no proof prints — why a fix is not offered is
  second-order, and the verdict has already said it.
* **The headline tells the reader what to do.** With more than one conflict
  it reads *"Unsatisfiable — N conflicts, pick a fix for each:"*. That is an
  instruction, and a true one: by hitting-set duality every solution resolves
  every conflict (Theorem 9), and one entry from each menu, in every
  combination, is a cheapest repair (Theorem 6, or the leading layer's
  exactness). It is not a claim that those are the only cheapest repairs, and
  it must not be read as one where an alternative exists — an alternative is
  a cheapest repair taking no entry of some conflict's menu (Lemma 31). The
  page does not leave that to the reader: every alternative opens with
  *"Or, to fix …"*, so the instruction and its exceptions stand side by side,
  the instruction first because it is the route the conflicts argue for.
* **An alternative prints after the conflicts, as fixes and not as a
  conflict.** Where a block's cover has more than one layer, each layer after
  the leading one prints after the last conflict, introduced by what it
  declines: *"Or, to fix without dropping requirement Knet:"* where every
  menu it avoids offers a single fix (their gerunds joined by *or*: without
  either), *"Or, to fix without any of the fixes for Conflicts 1 and 3:"*
  where one of them offers a choice, and *"Or, to fix another way:"* in the
  case Lemma 31 does not rule out, where the layer's selections avoid
  different menus and none of them all. Then what it offers: its menus of one
  entry said as one thing to do, joined by *and*, since they are all to be
  done; a menu with a choice numbered like a conflict's menu, *"…, and one
  of:"*, one witness under each entry, since an entry there completes the
  whole alternative; a single line and one witness where the layer leaves no
  choice at all; and, where several menus offer a choice, *"settle each of
  these:"* with one bullet per menu and its entries joined by *or*. Each
  witness is taken with the layer's other menus at their first entries and
  every other block settled the first way its menus offer. No proofs print
  there: reasons do not layer (Corollary 25), and the conflicts above have
  already explained every one. An alternative replaces the menus of its own
  block and nothing else — the rest of the report stands as printed, which is
  what its label says.
* **A coupled factor prints as a choice of whole options.** Where a menu's
  entries are compound — a factor the family couples, so that no one action
  of it is a choice on its own — each entry names the whole of what it asks
  for: *"relax both, or relax P and drop T, or drop both"*, never a list of
  actions the reader has to recombine. An entry is a thing to do, and a menu
  is a choice between things to do, at every arity.
* **Footers.** *Costlier fixes also exist.* prints only when Theorem 5
  answered yes and no *unless*-entry printed — an unless-entry exhibits a
  costlier fix (Lemma 28 says it never merely restates a cheap one), while
  an idle verdict exhibits none and does not suppress the news. *There are more minimal
  fixes than are shown.* prints only when the repair enumeration hit its
  cap and the deciding solve confirmed more. The reader is told about every
  gap in the page's own vocabulary — fixes, not solutions.
* **Menu wording.** Exactly the three-state table of Section 4.

### The verbs

Four forms, and which one is printed is data in the clause rather than a
flag beside it. Write `S` for the literal on the package a line is said to —
its printed consequent — and read `⊥ ∈ S` as "this statement is satisfied by
the package not being installed".

| form | when | what it says |
| --- | --- | --- |
| `your ‹kinds› leaves X r` | the line is one of the query's own facts | of `X`, the query's own constraints of those kinds left `r` — and only such a line may say *your* |
| `… requires X r` | `⊥ ∉ S` | `X` is installed, at one of `r` |
| `… constrains X r` | `⊥ ∈ S`, some version in `S` | if `X` is installed it is at one of `r`; the statement is silent if it is not |
| `… leaves no version of X` | `S = {⊥}` | `X` is not installed at all |

*Leaves* is the user's; *requires* forces; *constrains* binds only what is
there; *leaves no version of* is forced-empty. The bare package name is a
consequent whose range is the package's whole offering, never a consequent
with nothing in it — the fourth row exists so that a flipped line cannot
print as the package arriving when what it says is the package's going.

A trailing `(through P, Q and R)` is the whole of the difference between a
statement the registry makes directly and one an elimination composed: the
packages the route passed through, carrying no claim of their own (Routes,
Section 6). Its absence says the statement stands as it is.

Because the verb follows the printed consequent, flipping a line re-derives
it. That is what makes the next result safe to print.

### Linearization

**Proposition 22 (two sides linearize; three need not).** A meet of two
sides prints as one chain. A meet of three pairwise-overlapping sides does
not.

*Proof.* A clause has no direction (Section 2), so which of its packages is
the printed consequent is chosen when it is printed, and the contrapositive
is the same clause said the other way. Let the sides at pivot `P` be `S₁`
and `S₂`. Say `S₁` forward, from the packages its support names to the bound
it puts on `P`; the reader now holds `P ∈ S₁`. Say `S₂` to the package *it*
rests on, with `P` as its antecedent: what the chain has left at `P` — `S₁`,
narrowed by the query's own line about `P` where there is one — meets `S₂`
in nothing, since the family closes, so everything still open at `P` denies
`S₂` and the line fires on exactly what the lines before it left. That is
the resolution rule of Section 2, read out one step at a time, and the chain
ends at `S₂`'s own root fact. With three sides the same walk cannot be
arranged: one side opens the chain and one closes it, and the
third has both of its ends already spoken for — its support is introduced by
nothing the chain has said. The page is then two chains meeting at `P`,
which is what the meet display says in one place rather than leaving the
reader to find. And three sides can genuinely be needed: `{1,2}`, `{2,3}`,
`{1,3}` overlap pairwise and have empty intersection, so no two of them
contradict and no pair stands for the family. ∎

The pivot's own facts are not sides in this counting: a query line about `P`
is not an implication needing a root, and prints where the chain reaches `P`
— which is why "the compat leaves `P` `r`, and the one side says `P` `s`"
is a chain of two lines and not a meet at all. So the meet display survives
exactly at three sides and more, and Proposition 18 says how rare that is:
in a convex world an irredundant meet has at most three sides, and three of
them require a package one of whose bounds has a hole in it.

## 10. Substrate obligations

What the surrounding resolver must provide. Each is one interface; none
constrains how.

* **(S1) Assumable query facts.** Satisfiability of `ℛ ∪ X` for arbitrary
  `X ⊆ Q`, with each fact independently assumable, plus the bounded/blocked
  oracle of Theorem 4 and MUS extraction (any minimal unsatisfiable subset
  of an unsatisfiable assumption set).
* **(S2) Assumable registry statements.** A second instance over the same
  universe in which each clause of `ℛ` is individually assumable (e.g. via a
  selector per clause), for extracting cores `C_R`; the query facts are held,
  never selectable, while it is used. Built for diagnosis only — no resolve
  pays for it.
* **(S3) Entailment queries.** `ℛ ∪ A ⊨ σ` for `A ⊆ Q` and unit bounds `σ`
  (reducible to S1/S2-style satisfiability of the negation).
* **(S4) Resolution of withdrawals.** For any repair `W`, the resolver's
  actual solution of the withdrawn query — the versions a user would get —
  with the guarantee that the universe the diagnosis ran on answers for
  every withdrawal of `Q` it is asked about (in the present architecture
  this is the relaxation-stability of the filtered universe; whatever the
  architecture, the guarantee is the obligation).
* **(S5) Attribution data.** For each package, which of the query's
  constraint kinds exclude which versions — readable without any solver.
* **(S6) The version order**, per package, for rendering ranges and for
  Proposition 18's convexity talk.

## 11. Boundaries and honesty

Three enumerations in this design can be cut short, and each is accounted
for — two by a sentence on the page, one by a record on the diagnosis;
nothing else in the design is allowed to be incomplete.

1. **The reason walk** (Theorem 13) is exponential. Truncation costs
   completeness of each conflict's account, never soundness: every conflict
   still has a reason of its very own (Theorem 10), and every fix and every
   verdict on the page is decided by a solve and not by the walk. So the
   diagnosis records it, and the page does not announce it — a sentence
   saying "there may be more to say" gives the reader nothing to do, and the
   page owes only what can be acted on.
2. **Rectangle coverage** can be below `|F_min|`; the report says other
   equally cheap solutions exist (C3), and a one-entry menu then says "One
   fix:", never "The only".
3. **Beyond the cheapest**: whether strictly larger repairs exist is one
   oracle query (Theorem 5) and its answer sets the wording; *exploring*
   those larger repairs is out of scope by design, and the report says they
   exist without pretending to enumerate them.

Costs worth knowing, none load-bearing: finding `k` is a handful of
bounded solves (`k` is small in practice); enumerating `F_min` is one solve
per member; one core per owned reason; one projection per candidate pivot,
each a local elimination over a core's worth of clauses, with a work cap
whose only correct response on overflow is to fall back to printing the core
itself — a true, unshaped account being better than a shaped, partial one.

Finally, the two free-choice results deserve restating together, because
they bound the design space from both ends. The pivot theorem says the
explanation's anchor is *entirely* free — any package the argument touches
can carry it, so choose for the reader. The canonicity theorem says the
menus are *entirely* forced — the repair family admits exactly one
presentation, so there is nothing to choose. Everything else in a report
lives between those two poles, and every preference in this document lives
on the pivot side of the line.
