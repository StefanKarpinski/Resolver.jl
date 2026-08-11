# An empty dictionary that costs nothing to make: a singleton, since it has no
# fields. Constructing an unconstrained `Problem` — which every convenience
# `resolve` entry point does — must not allocate a dictionary just to say "no
# constraints", and the filter must not allocate a mask dictionary per round
# just to say "no masks".
struct EmptyDict{K,V} <: AbstractDict{K,V} end

Base.length(::EmptyDict) = 0
Base.iterate(::EmptyDict, state...) = nothing
Base.haskey(::EmptyDict, key) = false
Base.get(::EmptyDict, key, default) = default
Base.getindex(::EmptyDict, key) = throw(KeyError(key))

## constraints
#
# A problem's constraints are one namespace: a map from *kind* — the symbol
# naming what a constraint is — to the constraint itself. `:compat` and `:pin`
# are the two well-known kinds; any other kind is a predicate its caller
# supplied. They are stored alike because every consumer treats them alike: a
# constraint answers three questions and nothing else.
#
#   excludes(c, p, v)   does it forbid version `v` of package `p`?
#   named(c)            which packages does it speak about, or `nothing` for
#                       every package — a predicate cannot enumerate them
#   relax(c, pkgs)      the same constraint, no longer applying to `pkgs`
#                       (`nothing` again being every package), or `nothing`
#                       when nothing of it is left
#
# The map's value type is a closed union, so `is_excluded` — which
# `exclusion_masks` calls once per package per version — splits over three
# branches instead of dispatching dynamically.

# `:compat`: per package, the set of versions it admits, queried with `in`
struct Compat{P,S}
    allowed :: Dict{P,S}
end

# `:pin`: per package, the one version it is held at
struct Pins{P,V}
    at :: Dict{P,V}
end

# any other kind: a rule about versions — "no prereleases", say — so it reaches
# every package, less the ones it has been relaxed for
struct Predicate{P}
    forbids           # (p, v) -> Bool, true for the versions this forbids
    except :: Set{P}  # the packages it has been relaxed for
end

const Constraint{P,V,S} = Union{Compat{P,S}, Pins{P,V}, Predicate{P}}

# does this constraint forbid version `v` of package `p`?
excludes(c::Compat{P}, p::P, v) where {P} =
    (s = get(c.allowed, p, nothing); s !== nothing && v ∉ s)
excludes(c::Pins{P}, p::P, v) where {P} =
    (w = get(c.at, p, nothing); w !== nothing && v != w)
excludes(c::Predicate{P}, p::P, v) where {P} =
    p ∉ c.except && c.forbids(p, v)::Bool

# the packages this constraint speaks about, or `nothing` for every package
named(c::Compat{P}) where {P} = Set{P}(keys(c.allowed))
named(c::Pins{P}) where {P} = Set{P}(keys(c.at))
named(::Predicate) = nothing

# this constraint no longer applying to `pkgs` (`nothing` being every package),
# or `nothing` when nothing of it is left — so `relax(c, named(c))` is always
# `nothing`, whatever the shape
relax(::Compat, ::Nothing) = nothing
relax(c::Compat{P,S}, pkgs) where {P,S} =
    (d = filter(kv -> first(kv) ∉ pkgs, c.allowed);
     isempty(d) ? nothing : Compat{P,S}(d))
relax(::Pins, ::Nothing) = nothing
relax(c::Pins{P,V}, pkgs) where {P,V} =
    (d = filter(kv -> first(kv) ∉ pkgs, c.at);
     isempty(d) ? nothing : Pins{P,V}(d))
relax(::Predicate, ::Nothing) = nothing
# a predicate speaks about every package, so no finite set of them is all of it
relax(c::Predicate{P}, pkgs) where {P} =
    Predicate{P}(c.forbids, union(c.except, pkgs))

"""
    Problem(reqs; compat = ..., pin = ..., kind = predicate, ...)

A resolution problem: the required packages `reqs`, plus the constraints that say
which versions are admissible. A `Problem` carries everything that determines
*satisfiability*; orderings (the `by` priority order, the `order` version rank
order) are `resolve` parameters instead.

A constraint has a *kind* — a symbol naming it — and the kinds are one
namespace: every keyword here is a kind, and the keyword's name is the kind's
name, so no two constraints can share one. Two kinds are well known, and want a
dictionary keyed by package:

  * `compat`: per package, the set of allowed versions (queried with `in`).
  * `pin`: per package, the one version it is held at.

Any other keyword names a kind of its own, and wants a predicate: the
*admission* knobs — "no prereleases" is the one the resolver's own tooling uses
— say something about versions rather than about particular packages, so they
are stated as `predicate(p, v)`, true for the versions that kind forbids.

    Problem(reqs;
        compat = Dict(:B => [v"1.0", v"1.1"]),
        pre = (p, v) -> !isempty(v.prerelease))

Constraints for packages that don't exist, and constraints that exclude nothing,
are allowed and have no effect.

Every constraint here is a constraint on a *shared* package universe rather than
a deletion from a private one: that is what lets a single T1 artifact (see
[`pkg_info`](@ref Resolver.pkg_info)) serve queries that admit different things,
and what lets diagnostics eventually name the kind that ruled a version out.
"""
struct Problem{P, C<:AbstractDict{Symbol}}
    reqs :: Vector{P}
    # the constraints, by kind. The type is a parameter so the unconstrained
    # case can share one immutable empty dictionary, and so the value type
    # stays the concrete union `is_excluded` splits on
    constraints :: C
end

# A keyword's name is its kind, so nothing tells a typo apart from a kind of
# its own: `compt = Dict(...)` names a predicate perfectly well, and only the
# shape catches it. Hence checking before anything is built, and saying both
# the shape wanted and the shape given.
function check_shape(::Type{P}, kind::Symbol, value) where {P}
    if kind === :compat || kind === :pin
        value isa AbstractDict{P} || wrong_shape(P, kind, value)
    else
        # a predicate is anything callable; what is worth catching is a
        # collection handed to a kind that never wanted one
        isempty(methods(value)) && wrong_shape(P, kind, value)
    end
    return
end

# out of line: naming a type allocates enough to show up in `Problem` if the
# message is built on the success path
@noinline function wrong_shape(::Type{P}, kind::Symbol, value) where {P}
    want = kind in (:compat, :pin) ?
        "a dictionary keyed by package ($P)" :
        "a predicate, `(p, v) -> Bool`, true for the versions it forbids"
    got = isempty(methods(value)) ? "a $(typeof(value))" : "a callable"
    throw(ArgumentError("constraint kind $(repr(kind)) wants $want; got $got"))
end

# The constraint a keyword names, or `nothing` when it names none: an empty
# dictionary excludes nothing, so it makes nothing. A caller's dictionary is
# copied, so later mutation can't change the problem.
function constraint(::Type{P}, kind::Symbol, value) where {P}
    kind === :compat && return isempty(value) ? nothing : Compat(Dict(value))
    kind === :pin    && return isempty(value) ? nothing : Pins(Dict(value))
    return Predicate{P}(value, Set{P}())
end

# Every keyword is a kind, and its name is the kind's name — which is what
# makes a repeated kind a syntax error rather than something to check for.
Problem(reqs::SetOrVec{P}; kinds...) where {P} = problem(reqs, values(kinds))

function problem(reqs::SetOrVec{P}, kws::NamedTuple) where {P}
    for (kind, value) in pairs(kws)
        check_shape(P, kind, value)
    end
    # `:pin` and `:compat` supply the version and version-set types; a problem
    # without them is typed by neither
    pin, compat = get(kws, :pin, nothing), get(kws, :compat, nothing)
    C = Constraint{P,
        pin    === nothing ? Any : valtype(pin),
        compat === nothing ? Any : valtype(compat)}
    r = P[p for p in reqs]
    # made only once there is something to put in it, so a problem that names
    # no constraint costs no dictionary — whether it named no kind at all or
    # only empty ones
    cons = nothing
    for (kind, value) in pairs(kws)
        c = constraint(P, kind, value)
        c === nothing && continue
        cons === nothing && (cons = Dict{Symbol,C}())
        cons[kind] = c
    end
    cons === nothing && return Problem(r, EmptyDict{Symbol,C}())
    return Problem(r, cons)
end

# does the problem constrain anything at all? the fast paths below lean on this:
# an unconstrained resolve must not do any per-version work
is_constrained(prob::Problem) = !isempty(prob.constraints)

# `prob` with the requirements `drop_reqs` relaxed and each constraint named
# in `drops` relaxed for the packages `drops` maps it to — `nothing` there
# being every package, which relaxes that kind outright. Kinds `prob` does not
# carry are ignored, a constraint that isn't there being already relaxed; and
# what the map is consulted with is `haskey`, so *absent* and *relaxed
# outright* cannot collide.
function relax(
    prob      :: Problem{P},
    drop_reqs :: SetOrVec{P},
    drops     :: AbstractDict{Symbol},
) where {P}
    dropped = Set{P}(drop_reqs)
    reqs = P[p for p in prob.reqs if p ∉ dropped]
    cons = prob.constraints
    isempty(drops) && return Problem(reqs, cons)
    kept = Dict{Symbol, valtype(cons)}()
    for (kind, c) in cons
        if haskey(drops, kind)
            c = relax(c, drops[kind])
            c === nothing && continue
        end
        kept[kind] = c
    end
    isempty(kept) && return Problem(reqs, EmptyDict{Symbol, valtype(cons)}())
    return Problem(reqs, kept)
end

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
    for (_, c) in prob.constraints
        excludes(c, p, v) && return true
    end
    return false
end

# Which of the problem's constraints forbid version `v` of package `p`: the kind
# of each one that does, in kind order, and empty when the problem admits the
# version. Every one is named, where `is_excluded` stops at the first: lifting
# one constraint leaves the others in force, so what a question about lifting one
# needs is all of them.
function exclusion_kinds(prob::Problem{P}, p::P, v) where {P}
    kinds = Symbol[]
    for (kind, c) in prob.constraints
        excludes(c, p, v) && push!(kinds, kind)
    end
    return sort!(kinds)
end

# Why class `c` of package `p` holds nothing this problem admits: the class's
# members in version order, each paired with the kinds that exclude it
# (`exclusion_kinds`). A constraint forbids versions, and a class is empty
# when every member of it is forbidden — so the class is empty exactly when no
# pair here has an empty kind list, and these pairs are then the whole of what
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
    return [vers[i] => exclusion_kinds(prob, p, vers[i])
            for i in info_p.members[c]]
end

# the packages a constraint could forbid a version of: what the problem's
# constraints name between them, or every package in the universe as soon as one
# of them names none in particular
function exclusion_candidates(
    info :: AbstractDict{P},
    prob :: Problem{P},
) where {P}
    pkgs = Set{P}()
    for (_, c) in prob.constraints
        names = named(c)
        names === nothing && return keys(info)
        union!(pkgs, names)
    end
    return pkgs
end

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
