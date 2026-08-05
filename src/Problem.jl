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
  * `excludes`: the *admission* knobs — "no prereleases" is the one the resolver's
    own tooling uses — which say something about versions rather than about
    particular packages, and so come as `kind => predicate` pairs, where
    `predicate(p, v)` is true for the versions that source forbids and `kind` is a
    symbol naming it.

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
    # the versions that kind forbids ("no prereleases", say)
    excludes :: Vector{Pair{Symbol,Any}}
end

function Problem(
    reqs   :: SetOrVec{P};
    compat :: AbstractDict{P} = EmptyDict{P,Any}(),
    pins   :: AbstractDict{P} = EmptyDict{P,Any}(),
    excludes = NoKinds,
) where {P}
    Problem{P, valtype(pins), valtype(compat)}(
        P[p for p in reqs], adopt(compat), adopt(pins), adopt_kinds(excludes))
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

# What a user constraint does, and the only thing it does, is forbid *members*
# of interchangeability classes. Members are indistinguishable to everything in
# the registry, so forbidding some of them changes nothing a class can be asked
# to do — unless it forbids all of them, and then the class is empty and cannot
# be selected. Nothing of a constraint therefore appears in a conflict matrix.
# Constraints are evaluated once, here, into the per-version masks below;
# `class_ranking` turns those into each class's representative; and a class left
# without one is deactivated in the universe (see `Universe`), which is a single
# relaxable unit clause in the SAT instance and one degradation trigger in the
# filter.
#
# The consequence worth naming: a constraint cannot make two versions'
# constraint rows identical, because it does not touch the rows — and versions
# whose rows *are* identical are one class sharing one column, so there is no
# pair of separately-deletable objects for a partial relaxation to fall
# between.

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

# Which of the problem's constraints forbid version `v` of package `p`: one
# symbol per source that does — `:compat` for an allowed-version set the version
# is not in, `:pin` for a pin at another version, and its own symbol for each
# admission kind whose predicate holds — in that order, and empty when the
# problem admits the version. Every source is named, where `is_excluded` stops
# at the first: lifting one source leaves the others in force, so what a
# question about lifting one needs is all of them.
function exclusion_sources(prob::Problem{P}, p::P, v) where {P}
    srcs = Symbol[]
    s = get(prob.compat, p, nothing)
    s === nothing || v ∈ s || push!(srcs, :compat)
    w = get(prob.pins, p, nothing)
    w === nothing || v == w || push!(srcs, :pin)
    for (kind, forbids) in prob.excludes
        forbids(p, v)::Bool && push!(srcs, kind)
    end
    return srcs
end

# Why class `c` of package `p` holds nothing this problem admits: the class's
# members in version order, each paired with the sources that exclude it
# (`exclusion_sources`). A constraint forbids versions, and a class is empty
# when every member of it is forbidden — so the class is empty exactly when no
# pair here has an empty source list, and these pairs are then the whole of what
# emptied it: which versions the class holds, and what rules each one out.
#
# Constraints and the partition are all this reads; no solver is involved.
function class_exclusions(
    prob   :: Problem{P},
    p      :: P,
    info_p,        # `p`'s PkgInfo (declared after this file)
    c      :: Integer,
) where {P}
    vers = info_p.versions
    return [vers[i] => exclusion_sources(prob, p, vers[i])
            for i in info_p.members[c]]
end

# the packages a constraint could forbid a version of. Compat bounds and pins
# name their packages, and only a handful of packages are named, which is what
# makes constraining cheap; an admission knob is stated about versions instead, so
# it reaches every package in the universe.
exclusion_candidates(info::AbstractDict{P}, prob::Problem{P}) where {P} =
    isempty(prob.excludes) ?
        union(keys(prob.compat), keys(prob.pins)) : keys(info)

"""
    exclusion_masks(info, prob) :: AbstractDict{P, BitVector}

The versions of each package that `prob`'s constraints forbid, as a mask over
that package's version list in `info`. Only packages with a constraint that
actually excludes something get an entry, so packages nothing constrains — the
overwhelming majority, absent an admission knob — cost nothing downstream.

This is the only place a constraint is ever evaluated, and what reads it is
[`class_ranking`](@ref Resolver.class_ranking), which turns it into
representatives. Masks are index-based, so they are read against the version
list they were built from.
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
