# Resolver.jl

A new, experimental package version resolver for Julia, built on a real SAT solver.

Resolver.jl was presented at JuliaCon 2025. For the full story, watch the talk or read the slides:

- 📺 **Video:** [Pkg's new SAT-based version resolver](https://www.youtube.com/watch?v=Vab6en-6o6c) (JuliaCon Global 2025)
- 🖥️ **Slides:** [Markdown](https://github.com/StefanKarpinski/JuliaCon2025Talks/blob/main/New_SAT-based_Resolver.md) · [PDF](https://github.com/StefanKarpinski/JuliaCon2025Talks/blob/main/New_SAT-based_Resolver.pdf)

## Why a new resolver?

Julia's current `Pkg` resolver (a belief-propagation solver by Carlo Baldassi) works well, but it bakes in the assumption that **newer versions are always better**. That assumption runs deep enough that it's hard to support things people increasingly want:

- **Pre-release and build versions**, which don't fit a simple "higher is better" ordering.
- **Version fixing** — changing as little of an existing manifest as possible.
- **Download avoidance** — preferring versions already installed locally.
- **Downgrade resolution** — deliberately choosing the *lowest* compatible versions, e.g. to check that lower compat bounds actually work.

All of these need a resolver whose notion of "better" is configurable rather than hard-wired to version order.

## How it works

Resolver.jl separates *what* to optimize from *how* to solve, handing the solving to PicoSAT, an actual SAT solver:

- Packages and their versions become boolean variables, and compatibility and dependency rules become SAT clauses (a package implies one of its versions, a version implies its package, versions are mutually exclusive, dependencies must be present, and conflicts are forbidden).
- Before solving, the problem is shrunk by **reachability analysis** (keeping only versions that could appear in a solution) and **redundancy elimination** (dropping versions that are strictly dominated by others).
- Solutions are optimized **lexicographically, one package at a time**, against a caller-supplied preference order — so the resolver can enumerate multiple **Pareto-optimal** resolutions, where none is strictly better than another.
- Because the preference order is supplied per package, different strategies — newest, oldest/downgrade, prefer-installed, keep-current — can be mixed **within a single resolve**.

The core resolver in [`src/`](src/) is generic: it knows nothing about Julia versions or registries. The [`bin/resolve.jl`](bin/resolve.jl) command-line tool adapts it to Julia's real registries and can emit a `Manifest.toml`.

## Status

Experimental and research-stage. Resolver.jl is a standalone package rather than a part of `Pkg`, but the core resolver is complete, and `bin/resolve.jl` resolves real projects against the General registry today — it already backs [`julia-downgrade-compat`](https://github.com/julia-actions/julia-downgrade-compat), which uses it to compute minimum-version manifests in CI.

The main open problem is giving good, actionable feedback when a resolve turns out to be **unsatisfiable**. Ideas welcome.
