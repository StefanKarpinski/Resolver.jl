"""
    Resolver.Diagnostics

Why a resolve failed and what would make it succeed. `diagnose` puts further
satisfiability questions to the instance a failed resolve leaves behind and
answers with a [`Diagnosis`](@ref): the independent
[`Conflict`](@ref Resolver.Diagnostics.Conflict)s — sets of requirements that
cannot hold together, each with a chain of
[`Fact`](@ref Resolver.Diagnostics.Fact)s about the packages that says why —
and a menu of [`Fix`](@ref Resolver.Diagnostics.Fix)es, each a set of
[`Action`](@ref Resolver.Diagnostics.Action)s the user could carry out together
with the solution that carrying them out yields.

A `Diagnosis` is plain data, so a caller can render its own report from one
without asking the resolver anything further; `show(io, MIME("text/plain"), d)`
and [`action_phrase`](@ref Resolver.Diagnostics.action_phrase) are the report
this module writes.

These names are this module's API. `Resolver` exports `Diagnosis`, since it is
what `resolve` answers with when the requirements cannot be satisfied, and
re-exports none of the others.
"""
module Diagnostics

export Diagnosis, Conflict, Fix, Action, Fact, Requirement, Availability,
    unavailable, action_phrase, diagnose

using ..Resolver: SAT, Problem, PkgInfo, resolve, installed_lit,
    with_temp_clauses, with_classes_relaxed, sat_new_variable, sat_add_var,
    sat_add, exclusion_sources, class_exclusions
using ..UnsatCores: sat_mus, sat_disjoint_muses, sat_mcses

# What an unsatisfiable resolve says, and how it is found.
#
# A diagnosis answers two questions about a query the resolver could not
# satisfy: which requirements cannot hold together, and what the user could
# change so that they can. Both are questions about *relaxations* of the
# problem — the same problem with some requirements dropped and some of the
# classes this query emptied made choosable again — and both are put to the
# instance the resolve failed on, before it is discarded.
#
# Two things can be relaxed, and nothing else:
#
#   * a requirement, which the instance states as `installed_lit(sat, p)`;
#   * a class this query emptied, which it states as `forbidden_lit(sat, p, c)`.
#
# Everything else is a clause the diagnosis never touches: dependency edges,
# registry conflicts, the structural at-most-one, and the classes that are not
# empty. Every query below is therefore a solve under a subset of one fixed set
# of assumptions — which is what `UnsatCores` answers, and what the licensing
# argument on the manual's Diagnostics page covers.
#
# Candidate order is the one knob those answers take, and both uses here put
# the requirements first, for reasons that happen to agree: a repair is grown
# by keeping candidates in order, so requirements are kept and the menu prefers
# relaxing a constraint to dropping a dependency; a conflict is shrunk by
# deleting candidates in order, so a requirement the story does not need drops
# out of it.
#
# Both are asked a *package* at a time, over one assumable literal per package
# the query left short (`with_emptied_packages`). That is the granularity of
# both answers: a report says what has become of a package's versions, and a
# fix says which of a package's constraints to relax, neither of which is a
# sentence about one class. It is also what keeps the menu short — a package's
# classes usually go the same way, and one repair per class of each of them
# would be the same handful of instructions over and over.
#
# A version named in a report always comes from a real resolve of a real
# problem — never from a model of the failed instance. The instance was
# filtered and its classes were given representatives under constraints a fix
# relaxes, so a version read off it is a version that fix would not produce.
# Fixes are verified by resolving the modified problem against the artifact,
# and the solution that comes back is what the report shows.

## the report

"""
    Resolver.Diagnostics.Fact

One statement in a conflict's chain: something true of the packages that helps
explain why the requirements cannot hold together. The two kinds are
[`Requirement`](@ref Resolver.Diagnostics.Requirement) — "you require this
package" — and [`Availability`](@ref Resolver.Diagnostics.Availability) — "this
is what your constraints leave of the package's versions".

Facts are plain data, so a caller can render its own report from a
[`Diagnosis`](@ref) without asking the resolver anything further.
"""
abstract type Fact end

"""
    Resolver.Diagnostics.Requirement{P}(pkg)

The fact that the query requires `pkg`.
"""
struct Requirement{P} <: Fact
    pkg :: P
end

"""
    Resolver.Diagnostics.Availability{P,V}(pkg, members, excluded)

A package's whole offering, and what the query's own constraints make of it:
the versions of `pkg` the resolve could have chosen, each paired with the
constraints that rule that version out. A conflict's chain carries one of these
for every package the query took versions away from, so that a reader can see
why the package could not supply what the conflict needed.

  * `members` — every version of `pkg` the resolve had to choose between, in
    the order the package lists them (best first). Empty when the package has
    no installable version at all, which is a fact about the registry rather
    than about this query.
  * `excluded` — per member, the constraint sources that forbid it: `:compat`
    for an allowed-version set the version is not in, `:pin` for a pin at
    another version, and its own symbol for each exclusion kind. A member with
    no source is one this query admits.

The whole of the offering is here, so what the fact leaves the resolver is what
it does not exclude: the package is
[unavailable](@ref Resolver.Diagnostics.unavailable) exactly when every member
has a source, and relaxing the sources of any one member is enough to give that
version back.
"""
struct Availability{P,V} <: Fact
    pkg      :: P
    members  :: Vector{V}
    excluded :: Vector{Vector{Symbol}}
end

"""
    Resolver.Diagnostics.unavailable(fact::Availability) :: Bool

Whether the fact leaves its package with no version at all — every version it
offers excluded, or none to offer in the first place.
"""
unavailable(f::Availability) = all(!isempty, f.excluded)

"""
    Resolver.Diagnostics.Action{P}(kind, pkg)

One thing the user could change about `pkg`:

  * `:drop` — stop requiring it;
  * `:compat` — relax the allowed-version set constraining it;
  * `:unpin` — unpin it;
  * anything else — admit the versions that exclusion kind forbids, for this
    package (`:prerelease`, say).

Every kind but `:drop` relaxes one of the [`Problem`](@ref)'s constraint
sources, and carrying out the action means taking that source off that package.
Kinds are named for what the user does, which for a pin is to take it off.
"""
struct Action{P}
    kind :: Symbol
    pkg  :: P
end

Base.:(==)(a::Action, b::Action) = a.kind == b.kind && a.pkg == b.pkg
Base.hash(a::Action, h::UInt) = hash(a.pkg, hash(a.kind, hash(:Action, h)))

Base.:(==)(a::Requirement, b::Requirement) = a.pkg == b.pkg
Base.hash(a::Requirement, h::UInt) = hash(a.pkg, hash(:Requirement, h))

Base.:(==)(a::Availability, b::Availability) =
    a.pkg == b.pkg && a.members == b.members && a.excluded == b.excluded
Base.hash(a::Availability, h::UInt) =
    hash(a.pkg, hash(a.members, hash(a.excluded, hash(:Availability, h))))

"""
    Resolver.Diagnostics.Fix{P,V}(actions, solution)

One repair of the whole problem: a set of
[`Action`](@ref Resolver.Diagnostics.Action)s to carry out, and `solution`,
what resolving with them carried out actually yields. The solution is a real
resolve of the modified problem, so the versions in it are the versions the
user would get.

No fix's actions are a superset of another's, and no two fixes ask for the same
things.
"""
struct Fix{P,V}
    actions  :: Vector{Action{P}}
    solution :: Dict{P,V}
end

"""
    Resolver.Diagnostics.Conflict{P,V}(reqs, chain)

One independent reason the query fails: the requirements `reqs`, which cannot
hold together, and `chain`, a short sequence of
[`Fact`](@ref Resolver.Diagnostics.Fact)s about the packages that says why.
Every fact in the chain is load-bearing — take any one of them away and the
rest of the conflict dissolves.

Conflicts within a [`Diagnosis`](@ref) share no requirement, so each is a
separate thing to fix.
"""
struct Conflict{P,V}
    reqs  :: Vector{P}
    chain :: Vector{Fact}
end

"""
    Diagnosis{P,V}

What [`resolve`](@ref) returns instead of a solution when the requirements
cannot be satisfied: the independent
[`Conflict`](@ref Resolver.Diagnostics.Conflict)s, each with the chain of facts
that explains it, and a menu of verified
[`Fix`](@ref Resolver.Diagnostics.Fix)es, each of which repairs the whole
problem.

`show(io, MIME("text/plain"), d)` prints the report; the fields are plain data,
so a caller that wants a different report can build one from them without
asking the resolver anything further.

Pass `diagnose = false` to [`resolve`](@ref) for the bare `nothing` when only
the verdict is wanted.
"""
struct Diagnosis{P,V}
    conflicts :: Vector{Conflict{P,V}}
    fixes     :: Vector{Fix{P,V}}
end

## rendering
#
# The report states what is true about the packages. Classes, literals,
# assumptions and the solver are how the answer was found, not what it says, so
# none of them is ever named.

# how a constraint source is named to the user: `:compat` and `:pin` are the
# problem's own, anything else is an exclusion kind and goes by its own name
source_phrase(kind::Symbol) =
    kind === :compat ? "your compat" :
    kind === :pin    ? "your pin" :
                       "your $kind setting"

sources_phrase(kinds) = join(map(source_phrase, kinds), " and ")

"""
    Resolver.Diagnostics.action_phrase(action) :: String

An [`Action`](@ref Resolver.Diagnostics.Action) as an imperative the user could
carry out — "drop requirement Foo", "relax your compat on Bar".
"""
action_phrase(a::Action) =
    a.kind === :drop   ? "drop requirement $(a.pkg)" :
    a.kind === :compat ? "relax your compat on $(a.pkg)" :
    a.kind === :unpin  ? "unpin $(a.pkg)" :
                         "allow $(a.kind) versions of $(a.pkg)"

# a run of a package's versions, compressed: packages list their versions best
# first, so a run reads low to high
range_phrase(members, i::Int, j::Int) =
    i == j ? string(members[i]) : "$(members[j])–$(members[i])"

# join a list of phrases into an English list
function list_phrase(parts::Vector{String})
    length(parts) ≤ 1 && return join(parts)
    length(parts) == 2 && return "$(parts[1]) and $(parts[2])"
    return join(parts[1:end-1], ", ") * " and " * parts[end]
end

# An availability fact as one line. The fact runs over the package's whole
# version list, so every maximal run of versions excluded by the same sources is
# a range and gets one clause, and the runs it admits are simply not mentioned;
# the clauses after the first drop the verb, as an English list of them would.
function availability_phrase(f::Availability)
    members, excluded = f.members, f.excluded
    clauses = String[]
    i = 1
    while i ≤ length(members)
        j = i
        while j < length(members) && excluded[j+1] == excluded[i]
            j += 1
        end
        if !isempty(excluded[i])
            verb = i == j ? "is" : "are"
            by = sources_phrase(excluded[i])
            push!(clauses, isempty(clauses) ?
                "$(range_phrase(members, i, j)) $verb excluded by $by" :
                "$(range_phrase(members, i, j)) by $by")
        end
        i = j + 1
    end
    head = unavailable(f) ? "no version of $(f.pkg) is available" : string(f.pkg)
    return isempty(clauses) ? head : "$head: $(join(clauses, ", "))"
end

# print `text` filled to `width` columns, `lead` before the first line and
# `rest` before every later one
function print_wrapped(io::IO, text::AbstractString, lead::String, rest::String;
                       width::Int = 78)
    col = textwidth(lead)
    print(io, lead)
    for (k, word) in enumerate(split(text))
        w = textwidth(word)
        if k > 1
            if col + 1 + w > width
                print(io, "\n", rest)
                col = textwidth(rest)
            else
                print(io, " ")
                col += 1
            end
        end
        print(io, word)
        col += w
    end
    println(io)
end

count_phrase(n::Int, one::String, many::String) =
    n == 0 ? "no $many" : n == 1 ? "1 $one" : "$n $many"

function Base.show(io::IO, d::Diagnosis)
    print(io, "Diagnosis: ",
        count_phrase(length(d.conflicts), "conflict", "conflicts"), ", ",
        count_phrase(length(d.fixes), "fix", "fixes"))
end

function Base.show(io::IO, ::MIME"text/plain", d::Diagnosis{P,V}) where {P,V}
    println(io, "Unsatisfiable — ",
        count_phrase(length(d.conflicts), "conflict", "conflicts"), ", ",
        count_phrase(length(d.fixes), "fix", "fixes"), ":")
    for (k, c) in enumerate(d.conflicts)
        println(io)
        subject = list_phrase(String[string(p) for p in c.reqs])
        println(io, "Conflict $k: $subject cannot be satisfied",
            length(c.reqs) > 1 ? " together." : ".")
        for f in c.chain
            f isa Requirement ?
                println(io, "  • you require $(f.pkg)") :
                print_wrapped(io,
                    availability_phrase(f::Availability), "  • ", "    ")
        end
    end
    isempty(d.fixes) && return
    println(io)
    println(io, "Verified fixes:")
    # the versions worth showing are the ones the report has been talking about
    named = P[]
    for c in d.conflicts, f in c.chain
        p = f isa Requirement ? f.pkg : (f::Availability).pkg
        p in named || push!(named, p)
    end
    sort!(named)
    for (k, fix) in enumerate(d.fixes)
        print_wrapped(io, list_phrase(map(action_phrase, fix.actions)),
            lpad(k, 3) * ". ", "     ")
        # what the fix gets the user, of the packages the report has named; a
        # fix that installs none of them has nothing to show here
        allows = String["$p $(fix.solution[p])" for p in named
                        if haskey(fix.solution, p)]
        isempty(allows) ||
            print_wrapped(io, join(allows, ", "), "     → allows: ", "       ")
    end
end

## the diagnosis

# How many repairs to enumerate before settling for the ones already found.
# Independent conflicts multiply, so there can be exponentially many, and a menu
# is only useful while it is short; deduplication and dominance cut this down
# further. What is enumerated over is the requirements and the packages this
# query left short, of which there are few, so this is a ceiling rather than the
# usual stopping point.
const MAX_FIXES = 32

# The literals a diagnosis assumes, read back as what they say about the
# universe. Package variables are laid out in blocks — the package itself, then
# one variable per class — so the block a literal falls in names the package
# and the offset within it names the class, with offset 0 the package itself.
# The instance hands out variables in package-name order, so `pkgs` is sorted,
# and asking whether the instance has a package at all is a search of it.
struct VarMap{P}
    bases :: Vector{Int}
    pkgs  :: Vector{P}
end

function VarMap(sat::SAT{P}) where {P}
    bases = sort!(collect(values(sat.vars)))
    pkgs = Vector{P}(undef, length(bases))
    for (p, v) in sat.vars
        pkgs[searchsortedfirst(bases, v)] = p
    end
    return VarMap{P}(bases, pkgs)
end

# whether the instance has variables for `p` at all — its universe holds an
# installable version of the package
Base.haskey(vm::VarMap{P}, p::P) where {P} = insorted(p, vm.pkgs)

# `(p, 0)` for "package p is installed", `(p, c)` for "class c of p is forbidden"
function decode(vm::VarMap{P}, l::Integer) where {P}
    v = abs(l)
    i = searchsortedlast(vm.bases, v)
    return vm.pkgs[i], v - vm.bases[i]
end

# What this query's constraints do to `p`'s versions, as a fact. The versions a
# class holds and the sources excluding each of them are `class_exclusions`;
# taking every class at once is what lets the fact say which versions the query
# leaves standing as well as which it takes away, and a package is one thing to
# say a sentence about however many of its classes the story needed.
availability_fact(sat::SAT{P,V}, prob::Problem{P}, p::P) where {P,V} =
    Availability{P,V}(p, copy(sat.info[p].versions),
        Vector{Symbol}[exclusion_sources(prob, p, v)
                       for v in sat.info[p].versions])

# One assumable literal per package this query emptied something of, and the
# packages they stand for, in package-name order.
#
# A story is told a package at a time — what has become of a package's versions
# is one thing to say, however many of its classes the query emptied — so it is
# minimized a package at a time, and that takes a literal per package to
# minimize over. Each is a fresh variable that forbids every class of its
# package the query emptied, defined by one clause per class inside a frame of
# its own, so the instance is exactly as it was found afterwards. The variable
# occurs positively only in what is assumed, so defining it costs the instance
# no model of its own variables.
function with_emptied_packages(body::Function, sat::SAT{P}, vm::VarMap{P}) where {P}
    with_temp_clauses(sat) do
        lits = Int[]
        pkgs = P[]
        for l in sat.deact # sorted by variable, hence by package name
            p, _ = decode(vm, l)
            if isempty(pkgs) || last(pkgs) != p
                push!(pkgs, p)
                push!(lits, sat_new_variable(sat))
            end
            sat_add_var(sat, -last(lits))
            sat_add_var(sat, l)
            sat_add(sat)
        end
        body(lits, pkgs)
    end
end

# The story of one conflict: the smallest set of the query's own facts that is
# still unsatisfiable. Deletion is attempted in the order given, so putting the
# requirements first drops the ones the story does not need and keeps the facts
# the user can act on. The answer is a subsequence of what was asked, so the
# requirements come back in the problem's requirement order and the emptied
# packages in package order.
function conflict_story(
    sat      :: SAT{P,V},
    prob     :: Problem{P},
    vm       :: VarMap{P},
    cluster  :: Vector{Int},
    emptied  :: Dict{Int,P}, # per package literal, the package
    pkg_lits :: Vector{Int},
) where {P,V}
    mus = sat_mus(sat, [cluster; pkg_lits])
    reqs = P[]
    chain = Fact[]
    avails = Fact[]
    for l in mus
        p = get(emptied, l, nothing)
        if p === nothing
            q, _ = decode(vm, l)
            push!(reqs, q)
            push!(chain, Requirement(q))
        else
            push!(avails, availability_fact(sat, prob, p))
        end
    end
    append!(chain, avails)
    return Conflict{P,V}(reqs, chain)
end

# the kind of action that relaxes a constraint source. Kinds are named for what
# the user does, and what a user does with a pin is take it off; every other
# source goes by its own name
action_kind(src::Symbol) = src === :pin ? :unpin : src

# The user actions one correction set asks for. Dropping a requirement is
# already an action; giving a package's versions back is not, since what the
# user relaxes is a *constraint* — so for each class the query emptied, take the
# member that costs the fewest sources (ties going to the better version) and
# name those sources. Admitting one member is enough to make a class choosable
# again, so the union over the package's emptied classes gives all of them back.
function correction_actions(
    sat     :: SAT{P},
    prob    :: Problem{P},
    vm      :: VarMap{P},
    emptied :: Dict{Int,P}, # per package literal, the package
    mcs     :: Vector{Int},
    dropped :: Vector{P},
) where {P}
    actions = Action{P}[]
    for l in mcs
        p = get(emptied, l, nothing)
        if p === nothing
            q, _ = decode(vm, l)
            push!(actions, Action{P}(:drop, q))
            continue
        end
        reps_p = sat.reps[p]
        for c in eachindex(reps_p)
            iszero(reps_p[c]) || continue
            ex = class_exclusions(prob, p, sat.info[p], c)
            best = argmin(i -> (length(ex[i][2]), i), eachindex(ex))
            for src in ex[best][2]
                push!(actions, Action{P}(action_kind(src), p))
            end
        end
    end
    for p in dropped
        push!(actions, Action{P}(:drop, p))
    end
    unique!(actions)
    # relaxations before drops, so the menu reads as an instruction
    sort!(actions; by = a -> (a.kind === :drop, a.pkg, a.kind))
    return actions
end

# `prob` with the actions carried out: the named requirements gone, and every
# named constraint source taken off the package it names. An exclusion kind is
# stated about versions rather than packages, so relaxing one means excepting
# that package from it rather than dropping it.
function relaxed_problem(
    prob    :: Problem{P},
    actions :: Vector{Action{P}},
) where {P}
    drop_pkgs   = Set{P}(a.pkg for a in actions if a.kind === :drop)
    compat_pkgs = Set{P}(a.pkg for a in actions if a.kind === :compat)
    unpin_pkgs  = Set{P}(a.pkg for a in actions if a.kind === :unpin)
    kind_pkgs   = Dict{Symbol,Set{P}}()
    for a in actions
        a.kind in (:drop, :compat, :unpin) && continue
        push!(get!(Set{P}, kind_pkgs, a.kind), a.pkg)
    end
    reqs = P[p for p in prob.reqs if p ∉ drop_pkgs]
    compat = isempty(compat_pkgs) ? prob.compat :
        filter(kv -> first(kv) ∉ compat_pkgs, prob.compat)
    pins = isempty(unpin_pkgs) ? prob.pins :
        filter(kv -> first(kv) ∉ unpin_pkgs, prob.pins)
    excludes = isempty(kind_pkgs) ? prob.excludes :
        Pair{Symbol,Any}[
            kind => excepting(forbids, get(kind_pkgs, kind, nothing))
            for (kind, forbids) in prob.excludes]
    return Problem(reqs; compat, pins, excludes)
end

excepting(forbids, ::Nothing) = forbids
excepting(forbids, pkgs::Set) = (p, v) -> p ∉ pkgs && forbids(p, v)::Bool

"""
    Resolver.Diagnostics.diagnose(sat, prob, info; by, order) :: Diagnosis

Why the resolve of `prob` against the instance `sat` failed, and what would
make it succeed. `info` is the T1 artifact (see
[`pkg_info`](@ref Resolver.pkg_info)) the instance's universe was prepared
from; each fix is verified by resolving the modified problem against it, with
the same `by` and `order` the failed resolve used, and the solution that comes
back is the fix's witness.

`sat` must be the instance exactly as the failed query posed it — the
deactivation frame in place and no clause added since — and is left that way.
"""
function diagnose(
    sat  :: SAT{P,V},
    prob :: Problem{P},
    info :: AbstractDict{P, PkgInfo{P,V}};
    by   :: Function = identity,
    order = nothing,
) where {P,V}
    vm = VarMap(sat)
    # the same package required twice is one thing to assume and one to drop
    reqs = unique(prob.reqs)
    # a requirement the universe holds no installable version of is not
    # something the instance can be asked about: it has no variable, and no
    # constraint of this query is what took it away, since a constraint empties
    # classes rather than deleting them. It is a conflict of its own, and every
    # fix has to drop it
    gone = P[p for p in reqs if !haskey(vm, p)]
    req_lits = Int[installed_lit(sat, p) for p in reqs if haskey(vm, p)]

    conflicts = Conflict{P,V}[
        Conflict{P,V}([p], Fact[Requirement(p),
            Availability{P,V}(p, V[], Vector{Symbol}[])]) for p in gone]

    # partition the requirements into independent conflicts, with the query
    # exactly as production posed it
    clusters = sat_disjoint_muses(sat, req_lits)

    # distinct correction sets can ask for the same actions, and one action set
    # can be a strict superset of another; keep the minimal ones, in
    # enumeration order. Doing this before the witnesses are resolved saves the
    # resolves the surviving fixes do not need
    action_sets = Vector{Action{P}}[]
    seen = Set{Set{Action{P}}}()
    with_classes_relaxed(sat) do
        with_emptied_packages(sat, vm) do pkg_lits, pkgs
            emptied = Dict{Int,P}(zip(pkg_lits, pkgs))
            for cluster in clusters
                push!(conflicts,
                    conflict_story(sat, prob, vm, cluster, emptied, pkg_lits))
            end
            # every minimal repair, requirements first so they are kept and the
            # menu prefers relaxing a constraint to dropping a dependency
            for mcs in sat_mcses(sat, [req_lits; pkg_lits]; limit = MAX_FIXES)
                actions = correction_actions(sat, prob, vm, emptied, mcs, gone)
                key = Set(actions)
                key in seen && continue
                push!(seen, key)
                push!(action_sets, actions)
            end
        end
    end
    filter!(a -> !any(b -> Set(b) ⊊ Set(a), action_sets), action_sets)

    fixes = Fix{P,V}[]
    for actions in action_sets
        sol = resolve(info, relaxed_problem(prob, actions);
            by, order, diagnose = false)
        # relaxing a source admits at least what the correction set asked for
        # and never less, so a repair the instance verified resolves here too
        @assert sol !== nothing
        push!(fixes, Fix{P,V}(actions, sol))
    end
    return Diagnosis{P,V}(conflicts, fixes)
end

end # module
