# Diagnostics

A resolve that fails does not just say so. `resolve` returns the
[`Diagnosis`](@ref) itself, and printing it gives a report:

```julia-repl
julia> resolve(data, ["A", "B"])
Unsatisfiable — 1 conflict, 2 fixes:
Conflict 1: A, B cannot be satisfied together.
  • you require A
  • you require B
  • A (all versions) works with C only at 1
  • B (all versions) works with C only at 2

Verified fixes:
  1. drop requirement B
     → allows: A 1, C 1
  2. drop requirement A
     → allows: B 1, C 2

Upstream fixes:
  • a release of A or C relaxing their compat on each other
    → allows: A 1, B 1, C 2
  • a release of B or C relaxing their compat on each other
    → allows: A 1, B 1, C 1
```

Three sections, in decreasing order of how much you can do about them.

**Conflicts** are independent: each is a minimal set of requirements that
cannot be satisfied together, and fixing one does nothing for the others. The
bullets tell that conflict's story, ordered from the requirements outward
along the dependency graph.

**Verified fixes** are global — each one repairs *every* conflict at once —
and minimal: none of them relaxes a strict superset of what another does.
Each was checked by re-resolving, and the versions after `→ allows:` are what
you would actually get, restricted to the packages the conflicts implicate.
(The complete assignment is on the `Fix`, for programmatic use.)

**Upstream fixes** are the ones you cannot make yourself: a single new
release loosening one incompatibility would resolve that conflict. A tangled
conflict can produce a dozen of these, so they are ranked — by how current the
blamed versions are, then by how good the result is — and at most three are
printed, with a count of the rest. The structured `Diagnosis` keeps them all.

A story line may also be followed by an indented `(likewise for …)`. Those are
the facts a *minimal* conflict needs in order to close off versions your own
constraints already exclude: real, but not what you are reading for.

## The API

```julia
sol = resolve(data, reqs; compat, pins)
if sol isa Diagnosis
    show(stdout, MIME("text/plain"), sol)  # the report above
    for fix in sol.fixes
        # fix.actions  :: Vector{Fact}
        # fix.solution :: Dict{P,V}
    end
end
```

Everything is plain data, so a caller — Pkg, say — can render its own report
instead. The pieces:

| type | what it is |
|:--|:--|
| [`Diagnosis`](@ref) | requirements, conflicts, fixes, and the version lists to render against |
| [`Resolver.Conflict`](@ref) | one independent conflict: its requirement cluster, its fact chain, its upstream fixes |
| [`Resolver.Fix`](@ref) | a set of `Fact`s to relax, and the resolution that follows |
| [`Resolver.UpstreamFix`](@ref) | an incompatibility a new release could dissolve, and what would follow |

and the `Fact` kinds a chain or a fix is built from:
[`Resolver.Requirement`](@ref), [`Resolver.Uninstallable`](@ref),
[`Resolver.Pin`](@ref), [`Resolver.UserCompat`](@ref),
[`Resolver.Admission`](@ref), [`Resolver.Bound`](@ref).

An [`Resolver.Admission`](@ref) is an *admission knob* — the `excludes` of a
[`Problem`](@ref), which is where "no prereleases" lives. Those are constraints
like any other, so they get facts and fixes like any other:

```
Conflict 1: Wine_jll cannot be satisfied.
  • you require Wine_jll
  • every version of Wine_jll is a prerelease, which is not allowed

Verified fixes:
  1. allow prereleases of Wine_jll
     → allows: julia 1.11.9, Wine_jll 7.0.0-rc3+1
  2. drop requirement Wine_jll
     → allows: julia 1.11.9
```

Yankedness is deliberately *not* one of them: see
[below](#Modelling-withdrawn-versions).

Computing a diagnosis costs several extra solves. When you only need the
verdict — a satisfiability probe in a loop, say — pass `diagnose = false` and
get back `nothing`.

## Modelling platform packages

Some packages are not the user's choice: julia is present in every resolution
whatever anyone wants, and "drop requirement julia" is not advice. That is a
*modelling* question, not a policy knob — model such a package as a
**dependency edge**, never as a requirement. `bin/resolve.jl` does: the
registry provider gives every package version a dependency on `julia`, so
julia is in the solution because things need it, and cannot be dropped because
nothing "requires" it in the sense a fix could undo. Witness lists stay quiet
about it for the same reason: a report names the requirements and the packages
some story blames, and a dependency edge is neither.

The one knob left is wording. `Problem(…; labels)` gives a compat entry a
source flavour, so a bound from `Pkg.add(name, version)` reads "you requested
DataFrames at 1.7" and "relax the version you requested for DataFrames"
instead of calling it project compat. Unknown labels render neutrally.

One suggestion is filtered without being asked for: admitting prereleases of a
package is only proposed when the prerelease you would end up with has not been
*superseded* by a release of the same base version — nobody should be told to
install `1.2.3-alpha1` when `1.2.3` is out.

## Modelling withdrawn versions

A yanked version is the other thing no report should ever propose — and the way
to guarantee that is not a policy flag either. Model a withdrawn version by
**leaving it out of the universe**: the provider filters yanked versions before
the resolver sees them, so no solution contains one and no fix can name one, by
construction rather than by rule. `bin/resolve.jl` does exactly this, with
`--allow-yanked[=<pkgs>]` rebuilding the universe with them back in for anyone
who wants them.

What pre-filtering costs is a *story*, and that is recovered at the rendering
layer rather than by putting the constraint back. A bound admitting only yanked
versions admits nothing at all, and "no versions" is not the answer anyone
wants, so the provider keeps the set it dropped as explanation data:

```
Conflict 1: libsodium_jll cannot be satisfied.
  • you require libsodium_jll
  • your compat restricts libsodium_jll to no versions

Verified fixes:
  1. relax your compat on libsodium_jll
     → allows: julia 1.11.9, libsodium_jll 1.0.21+0
Note: the versions your compat on libsodium_jll admits are yanked (1.0.22+0);
--allow-yanked accepts them anyway.
```

## Two things the report will not do

**It relaxes a package's constraints together, not one at a time.** A compat
bound, a pin and an admission knob on the same package are one relaxable unit.
That is not conservatism: the three share one row in the resolver's
bookkeeping, and relaxing part of it is exactly the counterexample of
[Proposition D1′](theory/diagnostics.md#Below-column-closure-it-is-false) —
the answer would be wrong. Packages with a single constraint source, which is
almost all of them, are unaffected.

A fix does still name only the sources that rule out a version the instance
still has; one that forbids nothing there had no part in the relaxation, so
asking you to change it would be noise. That is why the Wine_jll example above
names nothing besides the prerelease knob.

**It blames incompatible *pairs*, not compat declarations.** A report says
"A (all versions) works with C only at 1", not "A declares compat C = 1", and
an upstream suggestion names both sides. Two reasons, and either alone would
be decisive: package info records compatibility symmetrically, so by the time
a problem reaches the resolver's instance which side declared a bound is no
longer recoverable; and per-declaration relaxation is not
[column-closed](theory/diagnostics.md#Relaxation-stability) either.

## Why this is trustworthy

The diagnosis runs on the very SAT instance the resolve failed on — the
prepared one, collapsed to interchangeability-class representatives and then
filtered, with the user's constraints popped back into assumptions. That
instance was built for a different question, so the
[theory page](theory/diagnostics.md) proves it answers this family of
questions correctly: satisfiability verdicts *and* layered answers are
preserved under every relaxation of the kind a diagnosis performs. The
verified solutions attached to fixes are a consequence of the stronger,
answer-preserving half.

The collapse is part of what has to be justified, not a complication to be
waved away: it deletes versions, and diagnosis then asks questions in which
some of the constraints that motivated those deletions are switched off.
Refining the classes by the exclusion mask before choosing representatives is
what keeps a forbidden version in the instance, so that relaxing its constraint
can bring it back.

Bound-level detail needs to relax conflicts the production instance keeps as
hard clauses, so it runs on a small sub-instance restricted to the conflict's
dependency closure. That restriction is exact, not an approximation
([Lemma D3](theory/diagnostics.md#Closure-exactness)). It costs a solve per
incompatibility, so very large closures are skipped rather than paid for: the
report then stops at requirement level, which for a hub package's conflict is
usually the whole story anyway.

## From the command line

`bin/resolve.jl` diagnoses automatically. An unresolvable project prints the
report to stderr, with uuids replaced by package names, and exits nonzero.
