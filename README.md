# Resolver.jl

A new, experimental package version resolver for Julia, built on a real SAT solver.

Resolver.jl was presented at JuliaCon 2025. For the full story, watch the talk or read the slides:

- 📺 **Video:** [Pkg's new SAT-based version resolver](https://www.youtube.com/watch?v=Vab6en-6o6c) (JuliaCon Global 2025)
- 🖥️ **Slides:** [Markdown](https://github.com/StefanKarpinski/JuliaCon2025Talks/blob/main/New_SAT-based_Resolver.md) · [PDF](https://github.com/StefanKarpinski/JuliaCon2025Talks/blob/main/New_SAT-based_Resolver.pdf)

## Why a new resolver?

Julia's current `Pkg` resolver, a belief-propagation solver by Carlo Baldassi, is a **heuristic**: it's fast and usually works well, but it isn't guaranteed to find a valid resolution even when one exists, and it makes no guarantees about *which* resolution it returns when several are possible. It also bakes in the assumption that **newer versions are always better** — an assumption that runs deep enough to make several things people increasingly want hard to support:

- **Pre-release and build versions**, which don't fit a simple "higher is better" ordering.
- **Version fixing** — changing as little of an existing manifest as possible.
- **Download avoidance** — preferring versions already installed locally.
- **Downgrade resolution** — deliberately choosing the *lowest* compatible versions, e.g. to check that lower compat bounds actually work.

All of these need a resolver whose notion of "better" is configurable rather than hard-wired to version order.

## How it works

Resolver.jl separates *what* to optimize from *how* to solve, handing the solving to PicoSAT, an actual SAT solver — which, unlike a heuristic, is guaranteed to find a solution whenever one exists:

- Packages and their versions become boolean variables, and compatibility and dependency rules become SAT clauses:
    - a package implies one of its versions,
    - a version implies its package,
    - a package's versions are mutually exclusive,
    - a version's dependencies must be present, and
    - conflicting versions cannot both be chosen.
- Before solving, the problem is shrunk by **reachability analysis** (keeping only versions that could appear in a solution) and **redundancy elimination** (dropping versions that are strictly dominated by others).
- The solution is optimized **lexicographically, one package at a time**, against a caller-supplied priority order — the resolver returns the unique optimal resolution that order determines (`nothing` if the requirements are unsatisfiable), and it is **Pareto-optimal**: no valid resolution is strictly better.
- Because the preference order is supplied per package, different strategies — newest, oldest/downgrade, prefer-installed, keep-current — can be mixed **within a single resolve**.

The core resolver in [`src/`](src/) is generic: it knows nothing about Julia versions or registries. The [`bin/resolve.jl`](bin/resolve.jl) command-line tool adapts it to Julia's real registries and can emit a `Manifest.toml`.

## Status

Experimental and research-stage. Resolver.jl is a standalone package rather than a part of `Pkg`, but the core resolver is complete, and `bin/resolve.jl` resolves real projects against the General registry today — it already backs [`julia-downgrade-compat`](https://github.com/julia-actions/julia-downgrade-compat), which uses it to compute minimum-version manifests in CI.

A major in-progress goal is providing much better diagnostics when a manifest is **unresolvable** (unsatisfiable). The rich theory behind SAT problems helps here: rather than a bare failure, the resolver can compute

- **MUSes (minimal unsatisfiable subsets)** — smallest sets of requirements that already conflict on their own, explaining *why* resolution failed; and
- **MCSes (minimal correction sets)** — smallest sets of requirements whose removal would make the manifest resolvable, pointing at *what to change* to fix it.
