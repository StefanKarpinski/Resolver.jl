# Resolver.jl

A package version resolver built on a real SAT solver (PicoSAT). Given
package data — versions, dependencies, and compatibility constraints — and a
set of required packages, [`resolve`](@ref) returns a single optimal solution
mapping each needed package to a version, or a [`Diagnosis`](@ref) — which
requirements cannot hold together, why, and a menu of verified fixes — when
they cannot be satisfied.

What "optimal" means, precisely, and why the resolver's aggressive problem
filtering provably does not change the answer, is worked out in the Theory
section:

- [Setting](theory/setting.md) — the resolution problem, solutions,
  tightness, and the SAT encoding.
- [The layered solution & filtering](theory/layered.md) — the semantics the
  resolver implements, and proofs that filtering preserves both
  satisfiability and the answer.
- [Pareto front optimality](theory/front.md) — an alternative, purely
  declarative notion of optimality that was fully worked out and ultimately
  set aside; its theorems and counterexamples illuminate why the layered
  semantics and its filter fit together as well as they do.
- [Relaxation-stable filtering](theory/relaxation.md) — the stronger property
  that a universe filtered for one query answers every relaxation of it, which
  is what lets a failed resolve be explained rather than merely reported.
- [Diagnosing an unsatisfiable resolve](theory/diagnostics.md) — what
  licenses the questions a `Diagnosis` is assembled from, and why the
  versions it shows come from a resolve rather than from the instance.

For the package's motivation and project status, see the
[README](https://github.com/StefanKarpinski/Resolver.jl).
