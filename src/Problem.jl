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

"""
    Problem(reqs; compat = ..., pins = ...)

A resolution problem: the required packages `reqs`, plus the user constraints
that say which versions are admissible — `compat` bounds (allowed-version sets,
queried with `in`) and `pins` (hold a package at exactly one version). A
`Problem` carries everything that determines *satisfiability*; orderings (the
`by` priority order, the version rank order) are `resolve` parameters instead.

Constraints for packages that don't exist, and constraints that exclude
nothing, are allowed and have no effect.
"""
struct Problem{P,V,S}
    reqs   :: Vector{P}
    # allowed-version sets, queried via `in(v, s)`, and exact-version pins.
    # abstractly typed so that the unconstrained case can share one immutable
    # empty dictionary instead of allocating two per problem
    compat :: AbstractDict{P,S}
    pins   :: AbstractDict{P,V}
end

function Problem(
    reqs   :: SetOrVec{P};
    compat :: AbstractDict{P} = EmptyDict{P,Any}(),
    pins   :: AbstractDict{P} = EmptyDict{P,Any}(),
) where {P}
    Problem{P, valtype(pins), valtype(compat)}(
        P[p for p in reqs], adopt(compat), adopt(pins))
end

# a caller's dictionary is copied, so later mutation can't change the problem;
# the shared empties are immutable and shared on purpose
adopt(d::AbstractDict) = Dict(d)
adopt(d::EmptyDict) = d

# does the problem constrain nothing at all? the fast paths below lean on this:
# an unconstrained resolve must not do any per-version work
is_constrained(prob::Problem) = !isempty(prob.compat) || !isempty(prob.pins)

# The user constraints are read as a *virtual package*: one always-required
# version, no dependencies, whose conflict rows are exactly the exclusions
# below. Under that reading a constrained problem is an ordinary resolution
# problem, so the filtering theorems (see the manual's Theory section) apply
# verbatim. The implementation special-cases the virtual package as the
# per-package exclusion masks computed here — one extra degradation trigger in
# `find_reachable`, one extra conflict column in `mark_necessary!`, and one
# selector-guarded clause group per constraint source in `SAT`.

# does the user forbid version v of package p?
function is_excluded(prob::Problem{P}, p::P, v) where {P}
    s = get(prob.compat, p, nothing)
    s === nothing || v ∈ s || return true
    w = get(prob.pins, p, nothing)
    w === nothing || v == w || return true
    return false
end

# the packages the user constrained, in a deterministic order
constrained_packages(prob::Problem{P}) where {P} =
    sort!(collect(union(keys(prob.compat), keys(prob.pins))))

"""
    exclusion_masks(info, prob) :: AbstractDict{P, BitVector}

The versions of each package that `prob`'s user constraints forbid, as a mask
over that package's version list in `info`. Only packages with a constraint
that actually excludes something get an entry, so packages the user did not
constrain — the overwhelming majority — cost nothing downstream.

Masks are index-based, so they must be recomputed whenever versions are
dropped and renumbered (see `filter_pkg_info!`).
"""
function exclusion_masks(
    info :: AbstractDict{P}, # a PkgInfo dict (declared before PkgInfo)
    prob :: Problem{P},
) where {P}
    is_constrained(prob) || return EmptyDict{P,BitVector}()
    excl = Dict{P,BitVector}()
    for p in constrained_packages(prob)
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
