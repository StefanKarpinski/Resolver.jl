# Setting

This page fixes the objects and notation used throughout the theory
section. Everything downstream — the semantics the resolver implements, the
correctness proofs for its filtering, and the Pareto-front analysis — is
stated in these terms.

## The resolution problem

A **resolution problem** consists of:

- a finite set of **packages** ``\mathcal{P}``;
- for each package ``p`` a finite, totally ordered list of **versions**,
  written by **rank**: rank ``1`` is the best version, rank ``2`` the next,
  and so on (the meaning of "best" is supplied by the caller — newest,
  oldest, closest-to-installed — the theory only uses the order);
- for each version ``p@v`` a **dependency set**
  ``\mathrm{deps}(p, v) \subseteq \mathcal{P}`` — note that dependencies are
  attached to *versions*, and different versions of the same package may
  have different dependency sets, a fact that turns out to carry a lot of
  weight;
- a symmetric **conflict** relation on version pairs (in the package data
  this arises from compatibility constraints: ``p@v`` conflicts with ``q@w``
  when ``w`` is outside ``p@v``'s compatibility bound for ``q``, or vice
  versa);
- a set of **requirements** ``\mathrm{reqs} \subseteq \mathcal{P}``.

A **solution** is a partial assignment ``s`` of versions to a support set
``\mathrm{supp}(s) \subseteq \mathcal{P}`` such that:

1. *coverage*: ``\mathrm{reqs} \subseteq \mathrm{supp}(s)``;
2. *dependency closure*: for every ``p \in \mathrm{supp}(s)``,
   ``\mathrm{deps}(p, s(p)) \subseteq \mathrm{supp}(s)``;
3. *conflict freedom*: no two assigned versions conflict.

Write ``\mathrm{rank}_s(p)`` for the rank of ``s(p)``.

## Tightness and pruning

The **closure** ``C(s)`` of a solution is the least set containing
``\mathrm{reqs}`` and closed under following dependencies of ``s``'s chosen
versions:

```math
C(s) \;=\; \bigcup_k L_k, \qquad
L_0 = \mathrm{reqs}, \quad
L_{k+1} = L_k \cup \bigcup_{q \in L_k} \mathrm{deps}(q, s(q)).
```

A solution is **tight** if ``\mathrm{supp}(s) = C(s)`` — it contains
nothing that isn't reachable from the requirements along its own chosen
dependency edges. ``R^*`` denotes the set of tight solutions. **Pruning**
``\pi(s)`` restricts ``s`` to ``C(s)``.

!!! note "Lemma (pruning)"
    ``\pi(s)`` is a tight solution with ``s``'s versions, and
    ``C(\pi(s)) = C(s)``.

    *Proof.* Coverage and dependency closure hold because the closure
    construction adds exactly the dependencies of its members; conflict
    freedom is inherited from ``s``; and the closure computation only ever
    inspects versions of packages already in the closure, which agree with
    ``s``. ∎

Two structural facts about tight solutions recur constantly:

- **Support is a function of the versions.** For a tight solution, the
  support is determined by the version assignment (it *is* the closure), so
  two distinct tight solutions with nested supports must already differ on
  a version within the smaller support.
- **Junk.** A valid solution need not be tight — it may carry packages
  ("junk") not justified from the requirements. Junk matters because the
  resolver works with a SAT solver, and SAT models are exactly the valid
  solutions, junk included.

## The SAT encoding

The resolver encodes a problem instance with one boolean variable per
package and one per version:

- package implies some version: ``p \Rightarrow \bigvee_i p@i``;
- version implies its package: ``p@i \Rightarrow p``;
- versions are mutually exclusive: ``\lnot p@i \lor \lnot p@j``;
- dependencies: ``p@i \Rightarrow q`` for each ``q \in \mathrm{deps}(p,i)``;
- conflicts: ``\lnot p@i \lor \lnot q@j`` for conflicting pairs.

With requirement unit clauses added, the models of this formula are
precisely the valid solutions covering the requirements (with junk
allowed). All constraints are dependency edges (positive, from versions to
packages) and *pairwise* version conflicts (negative) — a structural fact
several proofs exploit: any infeasibility ultimately traces to pairwise
conflicts, because dependencies alone can always be satisfied by including
more packages.

Before solving, the problem is shrunk by *filtering* —
[reachability analysis and redundancy elimination](layered.md#Filtering) —
whose correctness is the subject of the next page.
