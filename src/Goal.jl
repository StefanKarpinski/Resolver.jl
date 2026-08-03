# Goals: what a resolve is *for*, beyond satisfying the problem.
#
# `prob`'s constraints are negotiable -- a diagnosis may propose relaxing any of
# them -- and the goal is the fixed ask. That one sentence is the whole semantic
# difference between `resolve(info, prob; with = "CUDA")` and adding CUDA to
# `prob.reqs`: "drop requirement CUDA" can never appear as a fix.
#
# In the instance a goal is a hard assumption: one selector-guarded clause per
# term, held on by every probe a diagnosis runs, and absent from the keep-order
# fix enumeration works over. Held on rather than asserted at level 0 only so
# that a conflict can be asked whether it *needed* the goal, which is what tells
# the report's header apart from an ordinary one.

"""
    Goal(with, without)

A resolve's fixed ask: `with` names packages that must be present — each as
`pkg => versions` where `versions` is `nothing` for "any version" and otherwise
anything `in` answers for — and `without` names packages that must be absent.

Built from `resolve`'s `with` / `without` keywords rather than constructed
directly; nothing about it is exported.
"""
struct Goal{P}
    with    :: Vector{Pair{P,Any}}
    without :: Vector{P}
end

Base.isempty(g::Goal) = isempty(g.with) && isempty(g.without)

# the packages a goal is about: what the universe has to cover for it to be
# answerable at all
goal_packages(g::Goal{P}) where {P} =
    sort!(unique!(P[P[p for (p, _) in g.with]; g.without]))

# does this version satisfy a version goal? the target is any `in`-supporting
# set -- and a bare version of the universe's own version type reads as the
# singleton set containing it, so `with = "DataFrames" => v"1.8.2"` means what
# it looks like
matches_goal(v::V, s::V) where {V} = v == s
matches_goal(v, s) = v ∈ s

# Normalize the keyword arguments. Each accepts a package, a `pkg => versions`
# pair, or any collection of those; `nothing` (the default) means no goal, and a
# goal with nothing in it is `nothing` again, so the unconstrained path never
# pays for the machinery below.
function make_goal(::Type{P}, with, without) where {P}
    with === nothing && without === nothing && return nothing
    g = Goal{P}(Pair{P,Any}[], P[])
    add_with!(g.with, P, with)
    add_without!(g.without, P, without)
    isempty(g) ? nothing : g
end

add_with!(v::Vector{Pair{P,Any}}, ::Type{P}, ::Nothing) where {P} = v
add_with!(v::Vector{Pair{P,Any}}, ::Type{P}, p::P) where {P} =
    push!(v, p => nothing)
add_with!(v::Vector{Pair{P,Any}}, ::Type{P}, x::Pair) where {P} =
    push!(v, convert(P, first(x)) => last(x))
add_with!(v::Vector{Pair{P,Any}}, ::Type{P}, xs) where {P} =
    (foreach(x -> add_with!(v, P, x), xs); v)

add_without!(v::Vector{P}, ::Type{P}, ::Nothing) where {P} = v
add_without!(v::Vector{P}, ::Type{P}, p::P) where {P} = push!(v, p)
add_without!(v::Vector{P}, ::Type{P}, xs) where {P} =
    (foreach(x -> add_without!(v, P, x), xs); v)

# May a sub-instance still be filtered for reachability and redundancy under
# this goal? Yes with no goal, and yes for a pure presence ask -- whose packages
# the caller has already made requirements of the sub-instance, so reachability
# keeps their prefix. Anything that names versions, or bans a package, is
# exactly what those passes prune away.
goal_reach_safe(::Nothing) = true
goal_reach_safe(g::Goal) =
    isempty(g.without) && all(s === nothing for (_, s) in g.with)

# The versions a goal names by value, per package, as a mask over the universe's
# version list. The class collapse must keep these as representatives of their
# own: two versions in one interchangeability class are indistinguishable to
# every constraint, but a goal is stated on the version *value*, so collapsing
# the one the caller asked for into a sibling would answer a different question.
function goal_masks(
    info :: AbstractDict{P,PkgInfo{P,V}},
    goal :: Goal{P},
) where {P,V}
    masks = Dict{P,BitVector}()
    for (p, s) in goal.with
        s === nothing && continue
        haskey(info, p) || continue
        vers = info[p].versions
        m = falses(length(vers))
        for (i, v) in enumerate(vers)
            m[i] = matches_goal(v, s)
        end
        any(m) && (masks[p] = m)
    end
    return masks
end

"""
    prepare_goal_info(info, prob, goal; order = nothing) :: Dict{P,PkgInfo{P,V}}

The **goal sub-instance**: the T1 artifact `info`, collapsed to
interchangeability-class representatives and arc-consistent, but *not* filtered
for reachability or redundancy.

Those two passes are what the per-resolve preparation
([`prepare_pkg_info`](@ref Resolver.prepare_pkg_info)) adds, and they prune
exactly the escape routes a goal needs: reachability keeps only the versions an
*optimal* solution could use, and redundancy drops a version a better one
dominates — either can delete the only version of a package that avoids (or
reaches) the package the goal is about. See the manual's diagnostics theory page
on goal safety.

What survives is licensed: arc consistency deletes only versions no model
contains, which can witness no goal; and interchangeability classes never merge
across different dependency sets, so a package that avoids `X` can never hide
inside a class whose representative needs it. Versions a goal names by value are
kept out of the collapse besides, since a class sibling is not the version that
was asked for.
"""
function prepare_goal_info(
    info  :: AbstractDict{P,PkgInfo{P,V}},
    prob  :: Problem{P},
    goal  :: Goal{P};
    group :: Bool = true,
    order = nothing,
) where {P,V}
    perms = version_permutations(info, order)
    keep = group ?
        class_representatives(info, prob, perms, goal_masks(info, goal)) :
        nothing
    work = Dict{P,PkgInfo{P,V}}()
    copy_marked!(work, info, keep, perms)
    drop_unmarked!(work)      # materialize the collapse
    mark_installable!(work)   # ... which can leave a dependency unsatisfiable
    drop_unmarked!(work)
    return work
end
