# Diagnostics for resolutions that *succeeded* and still surprised someone.
#
# A user compat bound on a non-upgradable stdlib steers the julia choice: a
# julia version is compatible with exactly the stdlib versions it bundles, so a
# stale `LinearAlgebra = "~1.10"` quietly pins the whole toolchain to julia
# 1.10 -- correctly, and without a word. The resolution is right; what is
# missing is anyone saying *why*.
#
# This is the satisfiable-side sibling of the UNSAT machinery, and it reuses
# all of it. For a package resolved below the best version its universe would
# admit, assume the better version on the instance: that is unsatisfiable, and
# the same biased group-MUS that explains an unsatisfiable resolve explains
# this one. One probe with the blamed constraints relaxed verifies the fix and
# produces its versions, by the same optimizing descent, so a holdback's
# witness is as honest as an UNSAT fix's.

"""
    Holdback(pkg, resolved, best, chain, fix, versions)

Why `pkg` resolved to `resolved` rather than to `best`, the newest version of it
the problem's admission settings allow.

`chain` tells the story in the same [`Fact`](@ref Resolver.Requirement)
vocabulary an unsatisfiable diagnosis uses, and `fix` — when there is one — is a
verified correction whose `solution` says what `pkg` would resolve to instead.
`fix` is `nothing` when nothing the caller is allowed to relax would help, which
is a real answer and not a failure: a version held back by bounds only an
upstream release could move is held back for a reason the report should state
rather than paper over.

`versions` records the version lists the report needs, as in
[`Diagnosis`](@ref).
"""
struct Holdback{P,V}
    pkg      :: P
    resolved :: V
    best     :: V
    chain    :: Vector{Fact}
    fix      :: Union{Nothing,Fix{P,V}}
    versions :: Dict{P,Vector{V}}
end

# what `pkg` would resolve to under the holdback's fix, or `nothing`
improved(h::Holdback{P,V}) where {P,V} =
    h.fix === nothing ? nothing : get(h.fix.solution, h.pkg, nothing)

"""
    Holdbacks(held, unexamined)

What [`holdbacks`](@ref Resolver.holdbacks) returns: the explanations, as an
ordinary vector of [`Holdback`](@ref) — iterate it, index it, `isempty` it —
plus `unexamined`, the number of candidate packages the probe budget stopped it
from looking at.

Explaining a package costs several solves over the whole instance, so the budget
is real and a caller that silently dropped the rest would be misreporting. The
count is carried rather than implied, and `show` says it out loud.

`unexamined` counts packages that resolved below their best version and were
never probed — not packages known to be held back. Whether a constraint is
actually to blame is exactly what a probe decides, so claiming more than that
would be guessing.
"""
struct Holdbacks{P,V} <: AbstractVector{Holdback{P,V}}
    held       :: Vector{Holdback{P,V}}
    unexamined :: Int
end

Base.size(hs::Holdbacks) = size(hs.held)
Base.getindex(hs::Holdbacks, i::Int) = hs.held[i]
Base.IndexStyle(::Type{<:Holdbacks}) = IndexLinear()

function Base.show(io::IO, ::MIME"text/plain", hs::Holdbacks)
    for h in hs
        show(io, MIME("text/plain"), h)
    end
    render_unexamined(io, hs)
end

# the line that keeps a truncated report honest
function render_unexamined(io::IO, hs::Holdbacks)
    hs.unexamined > 0 || return
    println(io, "(", hs.unexamined, " more package",
            hs.unexamined == 1 ? "" : "s",
            " resolved below their best version, not examined)")
end

"""
    holdbacks(info, prob, sol, [pkgs]; by = identity, max_probes = 8)

Explain the packages in `pkgs` that `sol` resolved below their best *admissible*
version, as a `Holdback` each; packages that are at their best, or whose version
is merely the one the priority order settled on, produce nothing.

"Best" is the newest version `prob`'s admission knobs allow: compat bounds and
pins are what a fix relaxes, so they do not bound it, but a prerelease no query
would accept is not an option being missed. A package with no admissible version
at all is an unsatisfiable problem's business and produces nothing here.

`info` must be the *prepared* universe the solution came from (see
`prepare_pkg_info`) and `prob` its problem. The instance this builds is its own,
so nothing the caller holds is disturbed.

A holdback whose only remedy is something no fix can propose — an
incompatibility that only an upstream release could dissolve — is reported with
no fix rather than with an invented one.

`max_probes` bounds how many packages are *probed*, since each costs a handful
of solves over the whole instance; the result's `unexamined` says how many
candidates that left over, and rendering it says so.
"""
function holdbacks(
    info :: Dict{P,PkgInfo{P,V}},
    prob :: Problem{P},
    sol  :: AbstractDict{P,V},
    pkgs = keys(sol);
    by   :: Function = identity,
    max_probes :: Integer = 8,
) where {P,V}
    out = Holdback{P,V}[]
    none = Holdbacks{P,V}(out, 0)
    # the candidates, cheapest test first: a package already at the best version
    # the problem admits needs no solving to rule out
    cands = Pair{P,V}[]
    for p in pkgs
        haskey(info, p) && haskey(sol, p) || continue
        best = reachable_best(prob, p, info[p].versions)
        best === nothing && continue
        best == sol[p] && continue
        push!(cands, p => best)
    end
    isempty(cands) && return none
    sort!(cands)
    unexamined = max(0, length(cands) - max_probes)

    sat = SAT(info, prob)
    try
        with_relaxed_selectors(sat) do
            R = relaxation(sat, prob.reqs)
            # the requirement and constraint groups every probe holds on, and
            # the order the story prefers to blame: requirements are deleted
            # first so that a holdback tellable through the user's own compat
            # is blamed on the compat, not on "you require Foo"
            on = Int[R.reqs; R.user]
            shrink = Int[R.reqs; R.user]
            softof = Dict{P,Int}(R.rel[i].pkg => i for i in R.user)
            for (p, best) in Iterators.take(cands, max_probes)
                h = holdback(sat, R, p, best, sol, on, shrink, softof, by)
                h === nothing || push!(out, h)
            end
        end
    finally
        finalize(sat)
    end
    return Holdbacks{P,V}(out, unexamined)
end

# The best version of `p` a fix could actually reach: the newest one the
# problem's own *admission* settings allow. Compat bounds and pins do not bound
# it -- those are exactly what a fix proposes relaxing -- but an admission kind
# does: advertising a prerelease as "available" to a query that does not admit
# prereleases offers something that is not on the table, and worse, it makes the
# probe ask why *that* version was missed and answer with the knob rather than
# with the bound that is the real story. Yanked versions are pre-filtered out of
# the universe (see the manual's diagnostics page), so they never come up here.
#
# `nothing` when no version is admissible at all: that is an unsatisfiable
# problem's business, not a holdback's.
function reachable_best(prob::Problem{P}, p::P, vers::Vector{V}) where {P,V}
    isempty(prob.excludes) && return isempty(vers) ? nothing : first(vers)
    i = findfirst(vers) do v
        !any(forbids(p, v)::Bool for (_, forbids) in prob.excludes)
    end
    i === nothing ? nothing : vers[i]
end

function holdback(
    sat    :: SAT{P,V},
    R      :: Relaxation{P},
    p      :: P,
    best   :: V,
    sol    :: AbstractDict{P,V},
    on     :: Vector{Int},
    shrink :: Vector{Int},
    softof :: Dict{P,Int},
    by     :: Function,
) where {P,V}
    vers = sat.info[p].versions
    b = findfirst(==(best), vers)::Int

    # Could `p` be at `best` at all? If so its version is what the priority
    # order settled on rather than what a constraint forced, and there is no
    # constraint to blame -- reporting one would be inventing a culprit.
    with_best(g) = (sat_assume(sat, p, b); relax_sat(sat, R, g))
    with_best(on) && return nothing

    # what a minimal explanation needs, with `p@best` assumed throughout
    mus = group_mus_with(sat, R, on, shrink) do
        sat_assume(sat, p, b)
    end
    isempty(mus) && return nothing

    facts = Fact[]
    blamed = P[]
    for i in mus
        r = R.rel[i]
        append!(facts, group_facts(sat, r))
        r.kind === :user && push!(blamed, r.pkg)
    end
    # the incompatibility that actually stops `p@best`, read straight off the
    # conflict matrix: true by construction, and what makes the story concrete
    for q in blamed
        append!(facts, blocking_bound(sat.info, p, best, q))
    end

    # the fix: relax the blamed packages' constraints and see what `p` gets
    drop = Set{Int}(softof[q] for q in blamed if haskey(softof, q))
    fix = nothing
    if !isempty(drop)
        kept = Int[i for i in Int[R.reqs; R.user] if i ∉ drop]
        got = fix_solution(sat, R, kept, by)
        # only a fix if it actually moves `p` forward
        if haskey(got, p) && better(vers, got[p], sol[p])
            actions = Fact[]
            for i in drop
                append!(actions, group_facts(sat, R.rel[i]))
            end
            isempty(actions) ||
                (fix = Fix{P,V}(order_facts(actions), got))
        end
    end

    versions = Dict{P,Vector{V}}()
    record(q) = haskey(sat.info, q) &&
        (versions[q] = copy(sat.info[q].versions))
    record(p)
    for f in facts
        record(fact_pkg(f))
        f isa Bound && record(f.dep)
    end
    return Holdback{P,V}(p, sol[p], best, order_chain(facts, P, V),
                         fix, versions)
end

# is `u` a better version of this package than `w`? versions are listed best
# first, so the earlier index wins
function better(vers::Vector{V}, u::V, w::V) where {V}
    i = findfirst(==(u), vers)
    j = findfirst(==(w), vers)
    i !== nothing && j !== nothing && i < j
end

# the conflicts that stop `p@v` from coexisting with `q`, as one `Bound` naming
# the versions of `q` that would still work. Read off the matrix rather than
# probed, so it costs nothing and cannot be wrong -- it is a restatement of the
# instance, not a claim about minimality.
function blocking_bound(
    info :: Dict{P,PkgInfo{P,V}},
    p    :: P,
    v    :: V,
    q    :: P,
) where {P,V}
    facts = Fact[]
    p == q && return facts
    info_p = info[p]
    haskey(info_p.interacts, q) || return facts
    i = findfirst(==(v), info_p.versions)
    i === nothing && return facts
    b = info_p.interacts[q]
    vers_q = info[q].versions
    bad = [j for j = 1:length(vers_q) if info_p.conflicts[i, b+j]]
    isempty(bad) && return facts
    badset = Set(bad)
    push!(facts, Bound{P,V}(p, V[v], q,
        V[vers_q[j] for j = 1:length(vers_q) if j ∉ badset]))
    return facts
end

# `group_mus`, with extra assumptions re-applied before every solve. picosat
# drops assumptions after each call, so a query that fixes a version for the
# whole shrink has to say so each time.
function group_mus_with(
    setup  :: Function,
    sat    :: SAT,
    R      :: Relaxation,
    on     :: AbstractVector{Int},
    shrink :: AbstractVector{Int},
)
    probe(g) = (setup(); relax_sat(sat, R, g))
    probe(on) && return Int[]
    mus = collect(on)
    for g in shrink
        k = findfirst(==(g), mus)
        k === nothing && continue
        deleteat!(mus, k)
        probe(mus) && insert!(mus, k, g) # needed: keep it
    end
    return sort!(mus)
end

## rendering

function Base.show(io::IO, h::Holdback)
    print(io, "Holdback(", h.pkg, " ", h.resolved, " < ", h.best,
          h.fix === nothing ? ", no fix)" : ", fixable)")
end

function Base.show(io::IO, ::MIME"text/plain", h::Holdback)
    println(io, h.pkg, " resolved to ", h.resolved,
            "; ", h.best, " is available.")
    for f in h.chain
        f isa Bound && f.incidental && continue
        print(io, "  • ")
        render_fact(io, f, h)
        println(io)
    end
    if h.fix === nothing
        println(io, "  ⇒ nothing you can change would move it.")
        return
    end
    println(io, "  ⇒ held back by ",
            join(map(blame_phrase, h.fix.actions), " and "), ".")
    println(io, "  ", join(map(render_action, h.fix.actions), " and "),
            " → allows: ", holdback_solution(h))
end

# "your compat on Statistics", the noun form of the action a fix would take
blame_phrase(f::UserCompat) =
    f.label === :requested ? string("the version you requested for ", f.pkg) :
    f.label === :compat    ? string("your compat on ", f.pkg) :
                             string("the ", f.label, " on ", f.pkg)
blame_phrase(f::Pin)         = string("your pin on ", f.pkg)
blame_phrase(f::Admission)   = string("the ban on ", f.kind, " versions of ", f.pkg)
blame_phrase(f::Requirement) = string("your requirement on ", f.pkg)
blame_phrase(f::Fact)        = string(fact_pkg(f))

# the packages of a holdback's witness worth naming: the one held back, and
# whatever else the story is about
function holdback_solution(h::Holdback{P,V}) where {P,V}
    sol = h.fix.solution
    rel = Set{P}((h.pkg,))
    for f in h.chain
        push!(rel, fact_pkg(f))
        f isa Bound && push!(rel, f.dep)
    end
    pkgs = sort!(P[q for q in keys(sol) if q ∈ rel])
    sort!(pkgs; by = q -> q ≠ h.pkg) # the held-back package first
    isempty(pkgs) ? "nothing" :
        join((string(q, " ", sol[q]) for q in pkgs), ", ")
end

# one line, for a caller that wants the answer and not the reasoning
function summarize(h::Holdback)
    v = improved(h)
    v === nothing && return string(h.pkg, " is held back to ", h.resolved,
        " (", h.best, " available); nothing you can change would move it")
    string(h.pkg, " is held back to ", h.resolved, " by ",
           join(map(blame_phrase, h.fix.actions), " and "),
           "; relaxing ", length(h.fix.actions) == 1 ? "it" : "them",
           " allows ", h.pkg, " ", v)
end
