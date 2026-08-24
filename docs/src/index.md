# Resolver.jl

A package version resolver built on a real SAT solver (PicoSAT). Given
package data — versions, dependencies, and compatibility constraints — and a
set of required packages, [`resolve`](@ref) returns a single optimal solution
mapping each needed package to a version, or a [`Diagnosis`](@ref) — which
requirements cannot hold together, why, and a menu of verified fixes for
each of them — when they cannot be satisfied.

A diagnosis reads like this. Two requirements fail for unrelated reasons, so
each gets its own menu, and any one alternative from each menu resolves the
query:

```
Unsatisfiable — 2 conflicts, each of which must be fixed:

Conflict 1: DataFrames cannot be satisfied.
  • you require DataFrames
  • DataFrames: ≤ 0.22.7 excluded by your compat
  • DataFrames 1.3.5–1.8.2 requires PrettyTables
  • no version of PrettyTables is available: 0.1.0–3.4.8 excluded by your
    compat
  Fix it by any one of:
    1. relax your compat on DataFrames
       → allows: DataFrames 0.21.8
    2. relax your compat on PrettyTables
       → allows: DataFrames 1.8.2, PrettyTables 3.4.8
    3. drop requirement DataFrames

Conflict 2: Plots cannot be satisfied.
  • you require Plots
  • Plots 1.0.6–1.41.7 requires RecipesBase
  • no version of RecipesBase is available: 0.4.0–1.3.4 excluded by your
    compat
  Fix it by any one of:
    1. relax your compat on RecipesBase
       → allows: Plots 1.41.0, RecipesBase 1.3.4
    2. drop requirement Plots
```

Every fix is *verified*: the version list after each one is what a resolve of
the fixed query actually returns, not a guess. When the menus do not reach
every minimal fix the report says so — "Larger solutions also exist." when the
ones left out cost more, "Other solutions also exist." when it cannot promise
even that — and says nothing at all when the menus are the whole of it.

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
