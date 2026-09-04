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
Resolver.Diagnostics.Line
Resolver.Diagnostics.diagnose
```

### The statements a report is made of

```@docs
Resolver.Clauses
Resolver.Clauses.Clause
Resolver.Clauses.literal
Resolver.Clauses.resolve_on
Resolver.Clauses.resolve_raw
Resolver.Clauses.resolve_all
Resolver.Clauses.subsumes
Resolver.Clauses.Lit
Resolver.Clauses.clause
Resolver.Clauses.isbottom
Resolver.Clauses.version_order
Resolver.Clauses.range_phrase
Resolver.Clauses.clause_phrase
Resolver.Clauses.letters
```

### How a report is arrived at

```@docs
Resolver.Diagnostics.VarMap
Resolver.Diagnostics.with_emptied_packages
Resolver.Diagnostics.project
Resolver.Diagnostics.clause_of
Resolver.Diagnostics.clause_versions
Resolver.Diagnostics.clauses_satisfiable
Resolver.Relation
```

### What a report may be checked against

```@docs
Resolver.Diagnostics.print_conflict
Resolver.Diagnostics.report_problems
```

## Internals

```@docs
Resolver.PkgInfo
Resolver.pkg_info
Resolver.version_classes
Resolver.class_members
Resolver.no_shadows
Resolver.collapse_shadows
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
