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

## Asking for something: goals

`resolve` takes two more keywords, and they change the *question* rather than
the problem:

```julia
resolve(info, prob; with = "AnthropicClient")          # can I have this?
resolve(info, prob; with = "DataFrames" => v"1.8.2")   # ... at this version?
resolve(info, prob; without = "FillArrays")            # can I avoid this?
resolve(info, prob; with = ["CUDA", "MKL"], without = "FillArrays")
```

The return contract does not change: the optimal solution satisfying the
problem *and* the goal, or a `Diagnosis` of what it would take.

```
resolve(info, prob; with = "AnthropicClient")     # project already has BulkSMS
```
```
Unsatisfiable — 1 conflict, 1 fix:
Conflict 1: the goal and BulkSMS cannot be satisfied together.
  • you require BulkSMS
  • AnthropicClient (all versions) works with HTTP only at 1.0.1–1.11.0
  • BulkSMS (all versions) works with HTTP only at 0.8.19
Verified fixes:
  1. drop requirement BulkSMS
     → allows: AnthropicClient 0.1.0, HTTP 1.11.0
```

The whole difference between `with = "AnthropicClient"` and adding it to
`prob.reqs` is one line of that report: the problem's constraints are
*negotiable*, and a fix may propose relaxing any of them; the goal is the
*fixed ask*, so "drop requirement AnthropicClient" can never be offered. The
header says which conflicts the goal is part of.

A `without` goal's story is a chain of forced dependencies rather than of
incompatible pairs — nothing conflicts with anything, the package is simply
unavoidable:

```
Conflict 1: the goal and Makie cannot be satisfied together.
  • you require Makie
  • Makie (all versions) requires GeometryBasics
  • GeometryBasics (all admissible versions) requires FillArrays
```

Only *forced* edges appear: a path some other version choice could route
around explains nothing, which is what separates this from `pkg> why
FillArrays` printing every path there is. And when the package is avoidable,
there is no conflict at all — `resolve` returns the optimal FillArrays-free
solution, and [`changes`](@ref Resolver.changes) against the plain resolve is
how you price it.

### Three tiers of effort

Computing a diagnosis costs several extra solves, and so does the optimizing
descent. Three calls answer the same question with more or less work:

| call | cost | answer |
|:--|:--|:--|
| [`issatisfiable(data, prob; with)`](@ref issatisfiable) | one solve | `Bool` |
| `resolve(data, prob; with, diagnose = false)` | one solve + descent | the solution, or `nothing` |
| `resolve(data, prob; with)` | + explanation | the solution, or a `Diagnosis` |

Nothing but the effort separates them; the verdicts agree. A UI greying out an
upgrade button wants the first, a caller that needs the versions the second,
and only the one that will actually print a report pays for the third.

### Which universe answers a goal

A goal query runs on a sub-instance built from the T1 artifact, not on the
per-resolve one, because the per-resolve filter prunes exactly the escape
routes a goal needs. `A→B`, `B@3.1→X`, `B@2.3→∅`, no conflicts: reachability
keeps only `B@3.1`, so the prepared instance calls X unavoidable while the raw
problem avoids it via `B@2.3`. The [theory
page](theory/diagnostics.md#Goal-safety:-which-universe-answers-which-question)
proves what survives — arc consistency and the interchangeability collapse are
goal-safe, reachability and redundancy are not — and the counterexample above
is a test.

## Rendering your own reports

The `Diagnosis` retains *everything*. Every opinionated choice in the output
above is a documented function, and the default `show` is merely its first
client — no policy parameter was passed into the resolver to produce it, which
is the point of the design. Three tiers of control:

**Tier 1 — the default.** `show(io, MIME("text/plain"), d)`. Zero
configuration; what `bin/resolve.jl` prints.

**Tier 2 — the knobs.** [`report`](@ref) is what `show` calls, with its choices
as arguments:

```julia
report(io, d;
       max_upstream = 3,          # cap + "(N more …)"; typemax(Int) for all
       demote_incidental = true,  # fold the "likewise …" facts into an aside
       trim_witnesses = true,     # "→ allows:" lists the contested packages
       labels = Dict("DataFrames" => :requested))   # source wording
```

**Tier 3 — your own layout, our sentences** (or neither). The helpers `show`
is built from are public:

```julia
result = resolve(info, prob)
if result isa Diagnosis
    for c in result.conflicts
        print_header(c.reqs, c.goal)
        for f in c.chain
            f isa Resolver.Bound && f.incidental && continue  # drop, not demote
            println(io, "  • ", sprint(Resolver.render_fact, f, result))
        end
    end
    fixes = result.fixes                       # already verified, already optimal
    ups = filter(u -> maintained(u.bound.pkg), result.upstream)  # your knowledge
    Resolver.rank_upstream!(ups)
end
```

  * [`render_fact`](@ref Resolver.render_fact) writes one bullet's sentence;
    [`render_action`](@ref Resolver.render_action) writes a fix's imperative
    ("relax your compat on B") and [`blame_phrase`](@ref Resolver.blame_phrase)
    the same fact as a noun ("your compat on B"). A client that wants its own
    sentences reads the `Fact` fields instead — they are everything the
    sentences are generated from.
  * [`rank_upstream!`](@ref Resolver.rank_upstream!) sorts suggestions by how
    *current* the blamed versions are: a new release relaxing its latest bound
    is plausible, resurrecting an ancient range is not. Every `UpstreamFix`
    carries its `bound` and a verified `solution`, so a client with better
    knowledge — Pkg knows which packages are maintained — re-sorts or filters
    by its own criteria and ignores ours.
  * `Bound.incidental` is set during MUS construction when every version on the
    bound's own side is already excluded by the user's constraints: the fact's
    only job is closing off versions they cannot have anyway. Minimality is
    preserved by the remaining facts *for the versions that matter*, so a
    renderer may drop them outright.
  * [`superseded`](@ref Resolver.superseded) is the prerelease-supersession
    predicate fix enumeration already consults; it is public so a client
    filtering its own suggestion lists applies the same policy.
  * [`changes`](@ref) diffs two solutions into `Change(pkg, from, to)`, which
    is how a feasible goal query gets its price tag.

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

## When the resolve *succeeds* and still surprises

A resolution can be correct and still not be what you expected. A user compat
bound on a non-upgradable stdlib steers the julia choice — a julia release is
compatible with exactly the stdlib versions it bundles — so a leftover
`LinearAlgebra = "~1.10"` quietly pins the whole toolchain to julia 1.10, and
nothing in a successful resolve says so.

[`holdbacks`](@ref Resolver.holdbacks) explains that:

```julia
info = prepare_pkg_info(pkg_info(data, prob), prob)
sol  = resolve_prepared(info, prob)
for h in holdbacks(info, prob, sol, ["julia"])
    show(stdout, MIME("text/plain"), h)
end
```

```
julia resolved to 1.10.11; 1.12.6 is available.
  • you require LinearAlgebra
  • your compat restricts LinearAlgebra to 1.10.0–1.10.11
  • julia 1.12.6 works with LinearAlgebra only at 1.12.0
  ⇒ held back by your compat on LinearAlgebra.
  relax your compat on LinearAlgebra → allows: julia 1.12.6, LinearAlgebra 1.12.0
```

It is the satisfiable-side sibling of everything above, and reuses it: assume
the better version on the instance, get an unsatisfiable answer, and the same
biased group-MUS names the responsible constraints. The fix is verified by the
same optimizing descent, so its versions are the ones you would really get —
`Resolver.summarize` gives the same thing as one line.

Three things it will not do. It does not report a package whose version is
merely what the priority order settled on — if the better version was feasible
all along there is no constraint to blame, and inventing one would be worse
than silence. It measures "best" against the versions the problem *admits*, not
against the raw universe: a prerelease this query would never accept is not an
option being missed, and advertising it would bury the bound that is the real
story. And when nothing the user set is to blame — the version is held back by a
bound only an upstream release could move — it says `nothing you can change
would move it` rather than manufacturing a fix.

From the command line, `bin/resolve.jl` explains julia automatically whenever
it lands below its best admissible version — that is the surprise nobody went
looking for — as a single note:

```
Note: julia is held back to 1.10.11 by your compat on LinearAlgebra; relaxing it allows julia 1.12.6.
```

`--explain` gives the reasoning for everything that landed below its best
version. Explaining one costs several solves, so `max_probes` bounds how many
are probed; the result carries how many candidates that left over and the
report ends with `(N more packages resolved below their best version, not
examined)` rather than quietly showing a partial answer.

`--explain=<pkg>` answers whichever question is not already answered. If the
package is below its best version, that is the holdback above. If it is not —
because it is at its best, or because it is not in the solution at all — the
only question left is whether you can have it, which is the goal
`with = pkg`, and the answer is a price:

```
$ resolve.jl --explain=Statistics
You can have Statistics, at this price:
  • Statistics 1.11.1 added
```

`--explain=<pkg>@<version>` asks the version goal instead, and an impossible
one comes back as a full diagnosis rather than a bare refusal.

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
