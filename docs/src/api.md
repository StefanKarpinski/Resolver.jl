# API

```@docs
resolve
Problem
```

## Diagnosing unsatisfiability

An unsatisfiable resolve returns a [`Diagnosis`](@ref) rather than a bare
"no". Its parts are plain data, so a caller can render its own report — see the
[guide](diagnostics.md) — and the theory that licenses computing them is on the
[diagnostics theory](theory/diagnostics.md) page.

```@docs
Diagnosis
Resolver.Conflict
Resolver.Fix
Resolver.UpstreamFix
Resolver.Requirement
Resolver.Uninstallable
Resolver.UserCompat
Resolver.Pin
Resolver.Admission
Resolver.Bound
```

## Internals

```@docs
Resolver.pkg_info
Resolver.version_classes
Resolver.prepare_pkg_info
Resolver.version_permutations
Resolver.class_representatives
Resolver.find_reachable
Resolver.exclusion_masks
Resolver.diagnose
Resolver.Relaxation
Resolver.bound_story
Resolver.rank_upstream!
Resolver.superseded
```
