"""
    Resolver.Diagnostics

Why a resolve failed and what would make it succeed. `diagnose` puts further
questions to the universe a failed resolve was run against and the instance it
failed on, and answers with a [`Diagnosis`](@ref): the independent
[`Conflict`](@ref Resolver.Diagnostics.Conflict)s — sets of requirements that
cannot hold together, each with a chain of
[`Fact`](@ref Resolver.Diagnostics.Fact)s about the packages that says why, and
a menu of [`Fix`](@ref Resolver.Diagnostics.Fix)es for it, each a set of
[`Action`](@ref Resolver.Diagnostics.Action)s the user could carry out together
with the solution that carrying them out yields. One fix from every conflict's
menu, chosen any way at all, repairs the query.

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

using ..Resolver: SAT, Problem, Universe, relax, resolve, installed_lit,
    with_temp_clauses, with_classes_relaxed, sat_new_variable, sat_add_var,
    sat_add, exclusion_kinds, class_exclusions
using ..UnsatCores: sat_mus, sat_mcses

# What an unsatisfiable resolve says, and how it is found.
#
# A diagnosis answers two questions about a query the resolver could not
# satisfy: which requirements cannot hold together, and what the user could
# change so that they can. Both are questions about *relaxations* of the
# problem — the same problem with some requirements dropped and some of the
# classes this query emptied made choosable again — and both are put to the
# universe the resolve was run against and the instance it failed on, before
# they are discarded.
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
# The two questions are answered together, and that is what lets the report
# come apart. The cheapest repairs are enumerated first, and they factor into
# independent choices (`factor_repairs`); each choice is a conflict, its menu
# is what the choice is between, and its story is what is left unsatisfiable
# once every other choice has been made. So a conflict explains the menu under
# it rather than being arrived at separately, and one fix from each menu —
# chosen any way at all — is a repair of the whole query, without anything
# further being asked. See the manual's Diagnostics page.
#
# Both are asked a *package* at a time, over one assumable literal per package
# the query left short (`with_emptied_packages`). That is the granularity of
# both answers: a report says what has become of a package's versions, and a
# fix says which of a package's constraints to relax, neither of which is a
# sentence about one class. It is also what keeps the menu short — a package's
# classes usually go the same way, and one repair per class of each of them
# would be the same handful of instructions over and over.
#
# And it is what keeps every question a *constraint* relaxation, which is what
# the filtered universe is entitled to answer. Giving a package's empty classes
# back is exactly what lifting every kind that excludes anything of it does;
# giving one class back and leaving its siblings empty need be no constraint's
# doing at all. See the manual's Diagnostics page.
#
# A version named in a report always comes from a real resolve of a real
# problem — never from a model of the failed instance. The instance's classes
# were given representatives under constraints a fix relaxes, and were laid out
# in the order those representatives induce, so a version read off it is a
# version that fix would not produce. A fix's witness comes from resolving the
# relaxation `relax` derives from carrying it out — together with the first fix
# of every other conflict, since it is a combination that repairs the query —
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
constraint kinds that rule that version out. A conflict's chain carries one of
these for every package the query took versions away from, so that a reader can
see why the package could not supply what the conflict needed.

  * `members` — every version of `pkg` the resolve had to choose between, in
    the order the package lists them (best first). Empty when the package has
    no installable version at all, which is a fact about the registry rather
    than about this query.
  * `excluded` — per member, the kinds of the [`Problem`](@ref)'s constraints
    that forbid it: `:compat` for an allowed-version set the version is not in,
    `:pin` for a pin at another version, and its own name for each other kind
    whose predicate holds. A member with no kind against it is one this query
    admits.

The whole of the offering is here, so what the fact leaves the resolver is what
it does not exclude: the package is
[unavailable](@ref Resolver.Diagnostics.unavailable) exactly when every member
has a kind against it, and relaxing the kinds against any one member is enough
to give that version back.
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

One thing the user could change about `pkg`. All but one `kind` is the kind of
a [`Problem`](@ref) constraint, and the action is to take that constraint off
that package: `:compat` to widen the allowed-version set, `:pin` to take the
pin off, `:prerelease` to admit what that kind forbids, and so on for whatever
kinds the problem was built with.

The exception is `:drop`, which names a *requirement* rather than a constraint:
the action is to stop requiring the package, and no constraint of the problem
is touched. It is not a kind — it is the one action that is about the
requirements — so a problem that names a constraint kind `:drop` has two things
called by one name here.
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

One way of fixing one [`Conflict`](@ref Resolver.Diagnostics.Conflict): a set
of [`Action`](@ref Resolver.Diagnostics.Action)s to carry out, and `solution`,
what carrying them out actually yields. Fixing every conflict is what repairs
the query, so `solution` is a real resolve of the problem this fix and the
first fix of every *other* conflict relax — the versions in it are the versions
the user would get if that is how the rest are fixed.

No two fixes in a conflict's menu ask for the same things, and none of them
asks for more than another: the menu is a choice, not an order.
"""
struct Fix{P,V}
    actions  :: Vector{Action{P}}
    solution :: Dict{P,V}
end

"""
    Resolver.Diagnostics.Conflict{P,V}(reqs, chain, fixes)

One independent choice the query leaves the user: the requirements `reqs`,
which cannot hold together, `chain`, a short sequence of
[`Fact`](@ref Resolver.Diagnostics.Fact)s about the packages that says why, and
`fixes`, the menu of [`Fix`](@ref Resolver.Diagnostics.Fix)es that settle it.
Every fact in the chain is load-bearing — take any one of them away and the
rest of the conflict dissolves.

The conflicts of a [`Diagnosis`](@ref) are independent in that each has to be
fixed and any fix of one goes with any fix of another. Their *stories* may well
overlap: two conflicts can be about the same requirement without being the same
thing to do about it.
"""
struct Conflict{P,V}
    reqs  :: Vector{P}
    chain :: Vector{Fact}
    fixes :: Vector{Fix{P,V}}
end

"""
    Diagnosis{P,V}

What [`resolve`](@ref) returns instead of a solution when the requirements
cannot be satisfied: the independent
[`Conflict`](@ref Resolver.Diagnostics.Conflict)s, each with the chain of facts
that explains it and its own menu of
[`Fix`](@ref Resolver.Diagnostics.Fix)es. Every conflict has to be fixed, and
one fix from each — chosen any way at all — repairs the query, so the fixes of
the whole problem are the combinations.

`others` says what those combinations leave out:

  * `:none` — nothing at all. They are exactly the minimal fixes of the query,
    and every one of them changes as few things as any fix can.
  * `:larger` — only fixes that change more things than these do. Every
    smallest fix is a combination of these menus.
  * `:some` — fixes that change no more than these do, as well.

`show(io, MIME("text/plain"), d)` prints the report; the fields are plain data,
so a caller that wants a different report can build one from them without
asking the resolver anything further.

Pass `diagnose = false` to [`resolve`](@ref) for the bare `nothing` when only
the verdict is wanted.
"""
struct Diagnosis{P,V}
    conflicts :: Vector{Conflict{P,V}}
    others    :: Symbol
end

## rendering
#
# The report states what is true about the packages. Classes, literals,
# assumptions and the solver are how the answer was found, not what it says, so
# none of them is ever named.

# how a constraint kind is named to the user: `:compat` and `:pin` are what the
# user calls them too, anything else goes by its own name
kind_phrase(kind::Symbol) =
    kind === :compat ? "your compat" :
    kind === :pin    ? "your pin" :
                       "your $kind setting"

kinds_phrase(kinds) = join(map(kind_phrase, kinds), " and ")

"""
    Resolver.Diagnostics.action_phrase(action) :: String

An [`Action`](@ref Resolver.Diagnostics.Action) as an imperative the user could
carry out — "drop requirement Foo", "relax your compat on Bar".
"""
action_phrase(a::Action) =
    a.kind === :drop   ? "drop requirement $(a.pkg)" :
    a.kind === :compat ? "relax your compat on $(a.pkg)" :
    a.kind === :pin    ? "unpin $(a.pkg)" :
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
# version list, so every maximal run of versions excluded by the same kinds is a
# range and gets one clause, and the runs it admits are simply not mentioned;
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
            by = kinds_phrase(excluded[i])
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

# how many ways there are to repair the query: one fix out of every conflict's
# menu, chosen independently
nfixes(d::Diagnosis) =
    isempty(d.conflicts) ? 0 : prod(c -> length(c.fixes), d.conflicts)

function Base.show(io::IO, d::Diagnosis)
    print(io, "Diagnosis: ",
        count_phrase(length(d.conflicts), "conflict", "conflicts"), ", ",
        count_phrase(nfixes(d), "fix", "fixes"))
end

# What the menus leave out, once the reader has made a choice in each of them.
# Nothing to say when they leave out nothing.
others_phrase(others::Symbol) =
    others === :larger ? "Other, larger solutions exist." :
                         "Other solutions exist."

# what the fix gets the user, of the packages this conflict has named; a fix
# that installs none of them has nothing to show here
function print_allows(io::IO, fix::Fix{P}, named::Vector{P},
                      lead::String, rest::String) where {P}
    allows = String["$p $(fix.solution[p])" for p in named
                    if haskey(fix.solution, p)]
    isempty(allows) || print_wrapped(io, join(allows, ", "), lead, rest)
end

function Base.show(io::IO, ::MIME"text/plain", d::Diagnosis{P,V}) where {P,V}
    println(io, "Unsatisfiable — ",
        count_phrase(length(d.conflicts), "conflict", "conflicts"),
        length(d.conflicts) > 1 ? ", every one of which has to be fixed:" : ":")
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
        # the versions worth showing are the ones this conflict is about
        named = P[]
        for f in c.chain
            p = f isa Requirement ? f.pkg : (f::Availability).pkg
            p in named || push!(named, p)
        end
        sort!(named)
        if length(c.fixes) == 1
            # nothing to choose between: this is the one thing that can be done
            fix = only(c.fixes)
            print_wrapped(io, list_phrase(map(action_phrase, fix.actions)),
                "  The only fix: ", "    ")
            print_allows(io, fix, named, "    → allows: ", "      ")
        else
            println(io, "  Fix it by any one of:")
            for (j, fix) in enumerate(c.fixes)
                print_wrapped(io, list_phrase(map(action_phrase, fix.actions)),
                    lpad(j, 5) * ". ", "       ")
                print_allows(io, fix, named, "       → allows: ", "         ")
            end
        end
    end
    d.others === :none && return
    println(io)
    print_wrapped(io, others_phrase(d.others), "", "")
end

## the diagnosis

# How many repairs to enumerate before settling for the ones already found.
# What is enumerated over is the requirements and the packages this query left
# short, of which there are few — 6 to 36 repairs on the registry queries this
# was measured on — so it is a ceiling rather than the usual stopping point.
# Independent conflicts multiply, though, so there can be exponentially many,
# and what a truncated enumeration costs is not a fix but a claim: the smallest
# of the repairs found still factor into conflicts and their menus still repair
# the query, but the report can no longer say that nothing else would.
const MAX_REPAIRS = 256

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
# class holds and the kinds excluding each of them are `class_exclusions`;
# taking every class at once is what lets the fact say which versions the query
# leaves standing as well as which it takes away, and a package is one thing to
# say a sentence about however many of its classes the story needed.
availability_fact(sat::SAT{P,V}, prob::Problem{P}, p::P) where {P,V} =
    Availability{P,V}(p, copy(sat.info[p].versions),
        Vector{Symbol}[exclusion_kinds(prob, p, v)
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
# still unsatisfiable, drawn from the ones it is asked about. Deletion is
# attempted in the order given, so putting the requirements first drops the ones
# the story does not need and keeps the facts the user can act on. The answer is
# a subsequence of what was asked, so the requirements come back in the
# problem's requirement order and the emptied packages in package order.
function conflict_story(
    sat     :: SAT{P,V},
    prob    :: Problem{P},
    vm      :: VarMap{P},
    emptied :: Dict{Int,P}, # per package literal, the package
    cands   :: Vector{Int},
) where {P,V}
    mus = sat_mus(sat, cands)
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
    return reqs, chain
end

# relaxations before drops, so a set of actions reads as an instruction
sort_actions!(actions::Vector{<:Action}) =
    sort!(actions; by = a -> (a.kind === :drop, a.pkg, a.kind))

# The user actions one literal of a repair asks for. Dropping a requirement is
# already an action; giving a package's versions back is not, since what the
# user relaxes is a *constraint* — so for each class the query emptied, take the
# member that costs the fewest kinds (ties going to the better version) and name
# those kinds. Admitting one member is enough to make a class choosable again,
# so the union over the package's emptied classes gives all of them back.
function correction_actions(
    sat     :: SAT{P},
    prob    :: Problem{P},
    vm      :: VarMap{P},
    emptied :: Dict{Int,P}, # per package literal, the package
    lit     :: Int,
) where {P}
    p = get(emptied, lit, nothing)
    if p === nothing
        q, _ = decode(vm, lit)
        return Action{P}[Action{P}(:drop, q)]
    end
    actions = Action{P}[]
    reps_p = sat.reps[p]
    for c in eachindex(reps_p)
        iszero(reps_p[c]) || continue
        ex = class_exclusions(prob, p, sat.info[p], c)
        best = argmin(i -> (length(ex[i][2]), i), eachindex(ex))
        for kind in ex[best][2]
            push!(actions, Action{P}(kind, p))
        end
    end
    return sort_actions!(unique!(actions))
end

# One fix out of every conflict's menu, as the actions to carry out together
combined_actions(parts::Vector{Vector{Action{P}}}) where {P} =
    sort_actions!(unique!(reduce(vcat, parts; init = Action{P}[])))

# The smallest repairs, factored into independent choices.
#
# Two literals belong to the same choice when no repair contains both — taking
# one of them is what makes the other unnecessary — so the choices are the
# connected components of that relation. The family is exactly their product
# when every repair takes one literal from each and there are as many repairs as
# there are ways of doing that: a repair is determined by its literals, so the
# map from repairs to combinations is injective already, and once the counts
# agree it is onto. Every combination is then a repair, which is what lets the
# menus be offered separately.
#
# The second answer is whether the choices cover the whole family. When they do
# not, what comes back is the largest *rectangle* in it — a set of literals
# every repair in the rectangle shares, and the alternatives that complete it —
# and the family holds smallest repairs the menus do not reach.
function factor_repairs(repairs::Vector{Vector{Int}})
    lits = sort!(unique!(reduce(vcat, repairs; init = Int[])))
    n = length(lits)
    at = Dict{Int,Int}(l => i for (i, l) in enumerate(lits))
    shared = falses(n, n) # do the two literals appear in a repair together?
    for r in repairs, a in r, b in r
        shared[at[a], at[b]] = true
    end
    group = zeros(Int, n)
    ngroups = 0
    for i = 1:n
        group[i] == 0 || continue
        group[i] = (ngroups += 1)
        stack = Int[i]
        while !isempty(stack)
            u = pop!(stack)
            for v = 1:n
                group[v] == 0 && !shared[u, v] || continue
                group[v] = ngroups
                push!(stack, v)
            end
        end
    end
    choices = Vector{Int}[Int[lits[i] for i = 1:n if group[i] == g]
                          for g = 1:ngroups]
    length(repairs) == prod(length, choices; init = 1) &&
        all(r -> all(c -> count(in(c), r) == 1, choices), repairs) &&
        return choices, true
    return largest_rectangle(repairs)
end

# The largest rectangle in a family of same-sized repairs: a core of literals
# every repair in it contains, and the alternatives that complete it. A core one
# literal short of a repair is completed by exactly one literal of every repair
# that contains it, so every such core is a rectangle and there are few enough
# of them to try each one. The cores are tried in an order of their own, so that
# which rectangle wins a tie does not depend on the order the repairs were
# discovered in.
function largest_rectangle(repairs::Vector{Vector{Int}})
    cores = sort!(unique!(Vector{Int}[sort!(Int[k for k in r if k ≠ l])
                                      for r in repairs for l in r]))
    core = alts = Int[]
    for c in cores
        a = sort!(Int[k for s in repairs if c ⊆ s for k in s if k ∉ c])
        length(a) > length(alts) || continue
        core, alts = c, a
    end
    choices = push!(Vector{Int}[Int[k] for k in core], alts)
    return choices, length(alts) == length(repairs)
end

# Presentation order for the choices and the alternatives within them, given the
# order the repairs were asked for in. Within a choice, relaxing a constraint
# comes before dropping a requirement, which is the same preference the
# candidate order states by putting the requirements first; the choices
# themselves go in the order their earliest alternative was offered in.
function order_choices!(choices::Vector{Vector{Int}}, cands::Vector{Int},
                        nreqs::Int)
    at = Dict{Int,Int}(l => i for (i, l) in enumerate(cands))
    for c in choices
        sort!(c; by = l -> (at[l] ≤ nreqs, at[l]))
    end
    sort!(choices; by = c -> minimum(at[l] for l in c))
    return choices
end

# The withdrawal a set of actions asks for, in the two parts `relax` takes it:
# the requirements to stop requiring, and per constraint kind, the packages to
# lift that kind for. An action is already one or the other, so nothing here
# translates anything — which is the whole reason `Action` speaks in the
# problem's own kinds rather than a vocabulary of its own.
function withdrawal(actions::Vector{Action{P}}) where {P}
    drop_reqs = P[]
    drop_constraints = Dict{Symbol,Set{P}}()
    for a in actions
        a.kind === :drop ?
            push!(drop_reqs, a.pkg) :
            push!(get!(Set{P}, drop_constraints, a.kind), a.pkg)
    end
    return drop_reqs, drop_constraints
end

"""
    Resolver.Diagnostics.diagnose(sat, prob, univ; by, order) :: Diagnosis

Why the resolve of `prob` against the instance `sat` failed, and what would
make it succeed. `univ` is the universe `sat` was built from — the one
`prepare_pkg_info(info, prob; order)` produced, filtered for `prob` — and every
fix comes with a witness: the relaxation of `prob` that carrying it out (along
with the first fix of every other conflict) gives, resolved on `univ` with the
same `by` and `order` the failed resolve used. Filtering for `prob` deletes
nothing such a relaxation needs, which is the manual's Theorem C, so no fix
needs a universe of its own.

`sat` must be the instance exactly as the failed query posed it — the
deactivation frame in place and no clause added since — and is left that way.
"""
function diagnose(
    sat  :: SAT{P,V},
    prob :: Problem{P},
    univ :: Universe{P,V};
    by   :: Function = identity,
    order = nothing,
) where {P,V}
    vm = VarMap(sat)
    # the same package required twice is one thing to assume and one to drop
    reqs = unique(prob.reqs)
    # a requirement the universe holds no installable version of is not
    # something the instance can be asked about: it has no variable, and no
    # constraint of this query is what took it away, since a constraint empties
    # classes rather than deleting them. It is a conflict of its own, and
    # dropping it is the only thing that could settle that conflict
    gone = P[p for p in reqs if !haskey(vm, p)]
    req_lits = Int[installed_lit(sat, p) for p in reqs if haskey(vm, p)]

    # per conflict, in report order: the story it tells and its menu, each
    # alternative on it a set of actions
    stories = Tuple{Vector{P},Vector{Fact}}[
        (P[p], Fact[Requirement(p),
            Availability{P,V}(p, V[], Vector{Symbol}[])]) for p in gone]
    menus = Vector{Vector{Action{P}}}[
        Vector{Action{P}}[Action{P}[Action{P}(:drop, p)]] for p in gone]

    others = :none
    with_classes_relaxed(sat) do
        with_emptied_packages(sat, vm) do pkg_lits, pkgs
            emptied = Dict{Int,P}(zip(pkg_lits, pkgs))
            # requirements first, so that a repair keeps them and prefers
            # relaxing a constraint to dropping a dependency, and so that a
            # story drops the requirements it does not need
            cands = [req_lits; pkg_lits]
            repairs = sat_mcses(sat, cands; limit = MAX_REPAIRS)
            # the instance is satisfiable with nothing assumed — install
            # nothing — so there is always a repair to be found
            @assert !isempty(repairs)
            smallest = minimum(length, repairs)
            cheapest = filter(r -> length(r) == smallest, repairs)
            choices, whole = factor_repairs(cheapest)
            order_choices!(choices, cands, length(req_lits))
            # what the menus leave out: nothing when they reach every cheapest
            # repair and no repair costs more than the cheapest. A truncated
            # enumeration has both left to prove and repairs it never saw
            others =
                length(repairs) ≥ MAX_REPAIRS || !whole ? :some :
                length(cheapest) < length(repairs) ? :larger : :none
            # A story explains one choice by settling all the others: drop the
            # rest of one repair from what may be assumed, and what is left is
            # still unsatisfiable, minimally so once this choice is made. The
            # repair's own minimality is what puts this choice's literal in the
            # conflict and keeps the other choices' literals out of it.
            settled = Int[first(c) for c in choices]
            for (i, choice) in enumerate(choices)
                rest = Set{Int}(settled[j] for j in eachindex(settled) if j ≠ i)
                push!(stories, conflict_story(sat, prob, vm, emptied,
                    Int[l for l in cands if l ∉ rest]))
                push!(menus, Vector{Action{P}}[
                    correction_actions(sat, prob, vm, emptied, l)
                    for l in choice])
            end
        end
    end

    # the instance is back as the query posed it, which is what answering a
    # relaxation on it needs — the deactivation frame in place for it to lift.
    # A fix repairs one conflict, so what it is resolved with is the first fix
    # of every other one; those combinations repeat, and a solution is a
    # resolve, so they are answered once each
    defaults = Vector{Action{P}}[first(menu) for menu in menus]
    solutions = Dict{Vector{Action{P}},Dict{P,V}}()
    conflicts = Conflict{P,V}[]
    for (i, (creqs, chain)) in enumerate(stories)
        fixes = Fix{P,V}[]
        for actions in menus[i]
            together = combined_actions(Vector{Action{P}}[
                j == i ? actions : defaults[j] for j in eachindex(menus)])
            sol = get!(solutions, together) do
                s = resolve(sat,
                    relax(univ, prob, withdrawal(together)...; order); by)
                # relaxing a kind admits at least what the repair asked for and
                # never less, so a repair the instance found resolves here too
                @assert s !== nothing
                s
            end
            push!(fixes, Fix{P,V}(actions, sol))
        end
        push!(conflicts, Conflict{P,V}(creqs, chain, fixes))
    end
    return Diagnosis{P,V}(conflicts, others)
end

end # module
