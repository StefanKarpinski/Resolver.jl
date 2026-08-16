# API

```@docs
resolve
issatisfiable
Problem
Diagnosis
```

## Diagnostics

What a `Diagnosis` is made of. All of it is plain data, so a caller can render
its own report without asking the resolver anything further.

```@docs
Resolver.Diagnostics
Resolver.Diagnostics.Conflict
Resolver.Diagnostics.Fix
Resolver.Diagnostics.Action
Resolver.Diagnostics.action_phrase
Resolver.Diagnostics.Fact
Resolver.Diagnostics.Requirement
Resolver.Diagnostics.Availability
Resolver.Diagnostics.unavailable
Resolver.Diagnostics.diagnose
```

## Internals

```@docs
Resolver.PkgInfo
Resolver.pkg_info
Resolver.version_classes
Resolver.class_members
Resolver.collapse_classes!
Resolver.Universe
Resolver.prepare_pkg_info
Resolver.version_permutations
Resolver.class_ranking
Resolver.relaxations
Resolver.relax
Resolver.RelaxedProblem
Resolver.SparsePerm
Resolver.find_reachable
Resolver.exclusion_masks
```
