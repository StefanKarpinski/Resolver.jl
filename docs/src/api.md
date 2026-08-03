# API

```@docs
resolve
issatisfiable
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
Resolver.Dependency
```

## Explaining a surprising success

A resolve can succeed and still leave a package below the best version its
universe admits. [`holdbacks`](@ref Resolver.holdbacks) says what held it there.

```@docs
Resolver.holdbacks
Holdback
Resolver.Holdbacks
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
Resolver.prepare_goal_info
Resolver.Goal
Resolver.diagnose
Resolver.Relaxation
Resolver.bound_story
Resolver.rank_upstream!
Resolver.superseded
```
