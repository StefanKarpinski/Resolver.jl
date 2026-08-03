# An empty dictionary that costs nothing to make: a singleton, since it has no
# fields. Constructing an unconstrained `Problem` — which every convenience
# `resolve` entry point does — must not allocate a pair of dictionaries just to
# say "no constraints", and the filter must not allocate a mask dictionary per
# round just to say "no masks".
struct EmptyDict{K,V} <: AbstractDict{K,V} end

Base.length(::EmptyDict) = 0
Base.iterate(::EmptyDict, state...) = nothing
Base.haskey(::EmptyDict, key) = false
Base.get(::EmptyDict, key, default) = default
Base.getindex(::EmptyDict, key) = throw(KeyError(key))

# The kinds with no versions to forbid: the shared empty value for a problem's
# `excludes`, for the same reason `EmptyDict` exists.
const NoKinds = Pair{Symbol,Any}[]

"""
    Problem(reqs; compat = ..., pins = ..., excludes = ...)

A resolution problem: the required packages `reqs`, plus the constraints that say
which versions are admissible. A `Problem` carries everything that determines
*satisfiability*; orderings (the `by` priority order, the `order` version rank
order) are `resolve` parameters instead.

Three kinds of constraint, in two shapes:

  * `compat`: per package, the set of allowed versions (queried with `in`).
  * `pins`: per package, the one version it is held at.
  * `excludes`: the *admission* knobs — "no prereleases", "no yanked versions" —
    which say something about versions rather than about particular packages, and
    so come as `kind => predicate` pairs, where `predicate(p, v)` is true for the
    versions that source forbids and `kind` is a symbol naming it.

Constraints for packages that don't exist, and constraints that exclude nothing,
are allowed and have no effect.

Every constraint here is a constraint on a *shared* package universe rather than
a deletion from a private one: that is what lets a single T1 artifact (see
[`pkg_info`](@ref Resolver.pkg_info)) serve queries that admit different things,
and what lets diagnostics eventually name the kind that ruled a version out.
"""
struct Problem{P,V,S}
    reqs   :: Vector{P}
    # allowed-version sets, queried via `in(v, s)`, and exact-version pins.
    # abstractly typed so that the unconstrained case can share one immutable
    # empty dictionary instead of allocating two per problem
    compat :: AbstractDict{P,S}
    pins   :: AbstractDict{P,V}
    # admission knobs: per kind, a predicate `(p, v) -> Bool` that is true for
    # the versions that kind forbids
    excludes :: Vector{Pair{Symbol,Any}}
    # how to describe a package's compat bound in a report, when "your compat"
    # would be wrong
    labels :: AbstractDict{P,Symbol}
end

function Problem(
    reqs   :: SetOrVec{P};
    compat :: AbstractDict{P} = EmptyDict{P,Any}(),
    pins   :: AbstractDict{P} = EmptyDict{P,Any}(),
    excludes = NoKinds,
    labels :: AbstractDict{P} = EmptyDict{P,Symbol}(),
) where {P}
    Problem{P, valtype(pins), valtype(compat)}(
        P[p for p in reqs], adopt(compat), adopt(pins), adopt_kinds(excludes),
        adopt(labels))
end

# a caller's dictionary is copied, so later mutation can't change the problem;
# the shared empties are immutable and shared on purpose
adopt(d::AbstractDict) = Dict(d)
adopt(d::EmptyDict) = d


adopt_kinds(kinds) = isempty(kinds) ? NoKinds :
    Pair{Symbol,Any}[Symbol(kind) => forbids for (kind, forbids) in kinds]

# does the problem constrain anything at all? the fast paths below lean on this:
# an unconstrained resolve must not do any per-version work
is_constrained(prob::Problem) =
    !isempty(prob.compat) || !isempty(prob.pins) || !isempty(prob.excludes)


# The user constraints are read as a *virtual package*: one always-required
# version, no dependencies, whose conflict rows are exactly the exclusions
# below. Under that reading a constrained problem is an ordinary resolution
# problem, so the filtering theorems (see the manual's Theory section) apply
# verbatim. The implementation special-cases the virtual package as the
# per-package exclusion masks computed here — one extra degradation trigger in
# `find_reachable`, one extra conflict column in `mark_necessary!`, and one
# selector-guarded clause group per constraint source in `SAT`.

# does the problem forbid version v of package p?
function is_excluded(prob::Problem{P}, p::P, v) where {P}
    s = get(prob.compat, p, nothing)
    s === nothing || v ∈ s || return true
    w = get(prob.pins, p, nothing)
    w === nothing || v == w || return true
    for (_, forbids) in prob.excludes
        forbids(p, v)::Bool && return true
    end
    return false
end

# the packages a constraint could forbid a version of. Compat bounds and pins
# name their packages, and only a handful of packages are named, which is what
# makes constraining cheap; an admission knob is stated about versions instead, so
# it reaches every package in the universe.
exclusion_candidates(info::AbstractDict{P}, prob::Problem{P}) where {P} =
    isempty(prob.excludes) ?
        union(keys(prob.compat), keys(prob.pins)) : keys(info)

# ... the same set in a deterministic order, for the places that number things
constrained_packages(info::AbstractDict{P}, prob::Problem{P}) where {P} =
    sort!(P[p for p in exclusion_candidates(info, prob)])

# The constraints are read as a *virtual package*: one always-required version,
# no dependencies, whose conflict rows are exactly the exclusions. This lists,
# per source that says anything about `p`, the kind and a predicate on `p`'s
# versions — one selector per source, so diagnostics can tell a compat bound, a
# pin and an admission knob on the same package apart.
function exclusion_sources(prob::Problem{P}, p::P) where {P}
    srcs = Pair{Symbol,Any}[]
    s = get(prob.compat, p, nothing)
    s === nothing || push!(srcs, :compat => (v -> v ∉ s))
    w = get(prob.pins, p, nothing)
    w === nothing || push!(srcs, :pin => (v -> v != w))
    for (kind, forbids) in prob.excludes
        push!(srcs, kind => (v -> forbids(p, v)::Bool))
    end
    return srcs
end

"""
    exclusion_masks(info, prob) :: AbstractDict{P, BitVector}

The versions of each package that `prob`'s constraints forbid, as a mask over
that package's version list in `info`. Only packages with a constraint that
actually excludes something get an entry, so packages nothing constrains — the
overwhelming majority, absent an admission knob — cost nothing downstream.

Masks are index-based, so they must be recomputed whenever versions are
dropped and renumbered (see `filter_pkg_info!`).
"""
function exclusion_masks(
    info :: AbstractDict{P}, # a PkgInfo dict (declared before PkgInfo)
    prob :: Problem{P},
) where {P}
    is_constrained(prob) || return EmptyDict{P,BitVector}()
    excl = Dict{P,BitVector}()
    for p in exclusion_candidates(info, prob)
        haskey(info, p) || continue # constraint on an absent package
        vers = info[p].versions
        e = falses(length(vers))
        any = false
        for (i, v) in enumerate(vers)
            is_excluded(prob, p, v) || continue
            e[i] = true
            any = true
        end
        any && (excl[p] = e) # allow-everything entries get no mask
    end
    return excl
end
