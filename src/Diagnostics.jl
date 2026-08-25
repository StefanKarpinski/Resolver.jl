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
    Dependency, Incompatibility, unavailable, action_phrase, diagnose

using ..Resolver: SAT, PkgInfo, Problem, Universe, relax, resolve, nclasses,
    installed_lit, with_temp_clauses, with_classes_relaxed, sat_new_variable,
    sat_add_var, sat_add, sat_solve, exclusion_kinds, class_exclusions
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
explain why the requirements cannot hold together. The kinds are
[`Requirement`](@ref Resolver.Diagnostics.Requirement) — "you require this
package"; [`Availability`](@ref Resolver.Diagnostics.Availability) — "this is
what your constraints leave of the package's versions";
[`Dependency`](@ref Resolver.Diagnostics.Dependency) — "these versions of this
package need that one, at these versions"; and
[`Incompatibility`](@ref Resolver.Diagnostics.Incompatibility) — "these
versions of the two do not go together".

The first two are about what the query asked for, and the last two about what
the registry says, which is the step from one to the next that a reader could
not have worked out for themselves.

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
    Resolver.Diagnostics.Dependency{P,V}(pkg, versions, dep, offering, allowed, newest, oldest)

The fact that some versions of `pkg` need `dep`, and which versions of it they
will take:

  * `versions` — the versions of `pkg` this fact speaks for: a run of what the
    package offers, best first. Every one of them depends on `dep` and takes
    the same versions of it, and they are versions this query admits, so the
    fact is about the choice the resolve actually had.
  * `offering` — every version of `dep` the resolve could have chosen, in the
    order `dep` lists them.
  * `allowed` — per version of the offering, whether `versions` admit it.
    Nothing of the query is in this: it is what `pkg` says about `dep`. It is
    `nothing` where no bound of `pkg`'s is what the conflict turns on — a step
    on the way to the package that has nothing left, or a `dep` the query has
    left nothing of, which every bound would starve alike — and then the fact
    says only that `dep` is needed, which is all it is there to say.

The direction is a claim, not a convention. A registry records compatibility
symmetrically, so which side declared a bound is not generally recoverable, and
this fact is stated only where the registry does record `pkg` as depending on
`dep`. Where it does not, all that can be said is that the two do not go
together, which is an
[`Incompatibility`](@ref Resolver.Diagnostics.Incompatibility).
"""
struct Dependency{P,V} <: Fact
    pkg      :: P
    versions :: Vector{V}
    dep      :: P
    offering :: Vector{V}
    allowed  :: Union{Nothing,Vector{Bool}}
    # whether `versions` reaches the newest / oldest version of `pkg` this
    # query admits, so the range can be stated as an open bound
    newest   :: Bool
    oldest   :: Bool
end

"""
    Resolver.Diagnostics.Incompatibility{P,V}(pkg, versions, other, offering, allowed, newest, oldest)

The fact that some versions of `pkg` go with only some versions of `other`,
with nothing said about which of the two declared it. The fields are
[`Dependency`](@ref Resolver.Diagnostics.Dependency)'s: `versions` is a run of
`pkg`'s versions this query admits that all agree about `other`, and `allowed`
says, per version of `other`'s `offering`, whether they admit it.

A registry records compatibility symmetrically — that these two versions
cannot be installed together, not that one of them said so — so this is what
is left to say when the versions it speaks for do not depend on the package it
speaks about: there is no side to put the bound on, and it puts it on neither.
Reading it the other way round is just as true.
"""
struct Incompatibility{P,V} <: Fact
    pkg      :: P
    versions :: Vector{V}
    other    :: P
    offering :: Vector{V}
    allowed  :: Vector{Bool}
    newest   :: Bool
    oldest   :: Bool
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

Base.:(==)(a::Dependency, b::Dependency) =
    a.pkg == b.pkg && a.versions == b.versions && a.dep == b.dep &&
    a.offering == b.offering && a.allowed == b.allowed
Base.hash(a::Dependency, h::UInt) =
    hash(a.allowed, hash(a.offering, hash(a.dep,
        hash(a.versions, hash(a.pkg, hash(:Dependency, h))))))

Base.:(==)(a::Incompatibility, b::Incompatibility) =
    a.pkg == b.pkg && a.versions == b.versions && a.other == b.other &&
    a.offering == b.offering && a.allowed == b.allowed
Base.hash(a::Incompatibility, h::UInt) =
    hash(a.allowed, hash(a.offering, hash(a.other,
        hash(a.versions, hash(a.pkg, hash(:Incompatibility, h))))))

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
which cannot all hold, `chain`, a short sequence of
[`Fact`](@ref Resolver.Diagnostics.Fact)s about the packages that says why, and
`fixes`, the menu of [`Fix`](@ref Resolver.Diagnostics.Fix)es that settle it.

`reqs` is every requirement this conflict's fixes rescue, which need not be one
thing that fails: requirements can fail together, or each on its own for a
reason they have in common, and a conflict gathers both. So the chain is
unsatisfiable and nothing in it is a passenger — every fact belongs to some
minimal reason the requirements cannot hold — but it is a *union* of those
reasons rather than a single one, and taking a fact away need not dissolve all
of them.

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
# When a run reaches the end of what is on offer, an open bound says so: "≥
# 1.2.0" tells the reader nothing newer is at issue, where "1.2.0–1.9.3" leaves
# them wondering what is wrong with 1.9.4. `high`/`low` say whether the run
# reaches the newest/oldest version there is to reach.
function range_phrase(members, i::Int, j::Int; high::Bool = false, low::Bool = false)
    i == j && !(high && low) && return string(members[i])
    high && low && return "any version"
    high && return "≥ $(members[j])"
    low && return "≤ $(members[i])"
    return "$(members[j])–$(members[i])"
end

# join a list of phrases into an English list
function list_phrase(parts::Vector{String})
    length(parts) ≤ 1 && return join(parts)
    length(parts) == 2 && return "$(parts[1]) and $(parts[2])"
    return join(parts[1:end-1], ", ") * " and " * parts[end]
end

# An availability fact as one line, stated positively: what the query leaves,
# not what it took. The negative form made the reader compute a complement
# against a version list they do not have -- "≤ 0.6.77 excluded" gives the next
# line's "0.7.5" only to someone who knows what exists in between -- and the
# positive form is the same set the next line is about.
#
# The kinds are still named, since "your compat" is the part the user can act
# on, and a run of versions the query admits is what survives them.
function availability_phrase(f::Availability)
    members, excluded = f.members, f.excluded
    unavailable(f) && return isempty(members) ?
        "no version of $(f.pkg) is available" :
        "no version of $(f.pkg) is left: " *
            join(exclusion_clauses(members, excluded), ", ")
    # every maximal run the query admits, as a range
    parts = String[]
    i = 1
    while i ≤ length(members)
        if isempty(excluded[i])
            j = i
            while j < length(members) && isempty(excluded[j+1])
                j += 1
            end
            push!(parts, range_phrase(members, i, j;
                high = i == 1, low = j == length(members)))
            i = j
        end
        i += 1
    end
    kinds = unique!(reduce(vcat, excluded; init = Symbol[]))
    left = list_phrase(parts)
    # What is left is credited to the query's constraints only when they are
    # the whole reason it is what it is. Where redundancy elimination took
    # versions too, the range is still what the resolve had to choose between,
    # but saying the compat left it there would be a claim about the compat
    # that is not true of it alone
    isempty(kinds) && return "$(f.pkg) can be $left"
    return "$(kinds_phrase(kinds)) leaves $(f.pkg) at $left"
end

# the excluded runs, as the negative clauses the unavailable case still needs:
# there is no surviving range to name, so what took each version away is all
# there is to say
function exclusion_clauses(members, excluded)
    clauses = String[]
    i = 1
    while i ≤ length(members)
        j = i
        while j < length(members) && excluded[j+1] == excluded[i]
            j += 1
        end
        isempty(excluded[i]) || push!(clauses,
            "$(range_phrase(members, i, j; high = open_high(i, j, members),
                low = open_low(i, j, members))) by $(kinds_phrase(excluded[i]))")
        i = j + 1
    end
    return clauses
end

# The versions a mask picks out of a package's offering, as an English list of
# ranges. Packages list their versions best first, so a run reads low to high;
# a run is a run of the offering, so a range never spans a version the mask
# leaves out.
function offering_phrase(offering::Vector, allowed::Vector{Bool})
    parts = String[]
    i = 1
    while i ≤ length(allowed)
        if allowed[i]
            j = i
            while j < length(allowed) && allowed[j+1]
                j += 1
            end
            push!(parts, range_phrase(offering, i, j;
                high = i == 1, low = j == length(allowed)))
            i = j
        end
        i += 1
    end
    return join(parts, ", ")
end

# What a relation is about. The versions are named whatever they are: they are
# the versions this query left standing, so leaving them out would let the
# sentence be read as one about the package rather than about the choice the
# resolve had — and a version range is shorter than saying so.
# A run reaching one end of the offering is stated as an open bound; a run
# reaching *both* is not "any version" — it is every version there is, and the
# sentence around it already says so — so it is named in full.
open_high(i::Int, j::Int, members) = i == 1 && j != length(members)
open_low(i::Int, j::Int, members) = j == length(members) && i != 1

# An open bound is stated only at one end: a run reaching *both* ends of what
# the query admits is still not "any version" — the availability fact just said
# which versions were taken away — so it is named in full.
relation_subject(f) =
    "$(f.pkg) $(range_phrase(f.versions, 1, length(f.versions);
        high = f.newest & !f.oldest, low = f.oldest & !f.newest))"

# a dependency as one line, with the bound where there is one to state
function dependency_phrase(f::Dependency)
    subject = relation_subject(f)
    (f.allowed === nothing || all(f.allowed)) &&
        return "$subject requires $(f.dep)"
    any(f.allowed) || return "$subject requires $(f.dep) at no version of it"
    return "$subject requires $(f.dep) at " *
        offering_phrase(f.offering, f.allowed)
end

# An incompatibility as one line. Nothing here says which package's compat
# entry it came from, because nothing in the universe does
function incompatibility_phrase(f::Incompatibility)
    subject = relation_subject(f)
    any(f.allowed) || return "$subject works with no version of $(f.other)"
    return "$subject works with $(f.other) only at " *
        offering_phrase(f.offering, f.allowed)
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

# a fact as the line the report prints for it; the requirement is the one kind
# with no versions in it, so it is the caller's to print
fact_phrase(f::Availability)    = availability_phrase(f)
fact_phrase(f::Dependency)      = dependency_phrase(f)
fact_phrase(f::Incompatibility) = incompatibility_phrase(f)

# the packages a fact speaks about, in the order it names them: what a report
# shows versions of is what its facts are about
fact_packages(f::Requirement)     = (f.pkg,)
fact_packages(f::Availability)    = (f.pkg,)
fact_packages(f::Dependency)      = (f.pkg, f.dep)
fact_packages(f::Incompatibility) = (f.pkg, f.other)

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
    others === :larger ? "Larger solutions also exist." :
                         "Other solutions also exist."

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
        length(d.conflicts) > 1 ? ", each of which must be fixed:" : ":")
    for (k, c) in enumerate(d.conflicts)
        println(io)
        subject = list_phrase(String[string(p) for p in c.reqs])
        # a conflict can gather requirements that fail together and
        # requirements that each fail on their own, so what is said of them is
        # what is true of both: the conjunction is what cannot hold
        println(io, "Conflict $k: $subject cannot ",
            length(c.reqs) == 1 ? "be satisfied." :
            length(c.reqs) == 2 ? "both be satisfied." : "all be satisfied.")
        for f in c.chain
            f isa Requirement ?
                println(io, "  • you require $(f.pkg)") :
                print_wrapped(io, fact_phrase(f), "  • ", "    ")
        end
        # the versions worth showing are the ones this conflict is about
        named = P[]
        for f in c.chain, p in fact_packages(f)
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

# How many of the smallest repairs to enumerate before settling for the ones
# already found. What is enumerated over is the requirements and the packages
# this query left short, of which there are few, and only the smallest repairs
# among them — a median of 2 of them on the registry queries this was measured
# on — so it is a ceiling rather than the usual stopping point.
# Independent conflicts multiply, though, so there can be exponentially many,
# and what a truncated enumeration costs is not a fix but a claim: the repairs
# found still factor into conflicts and their menus still repair the query, but
# the report can no longer say that nothing else would, nor ask whether
# anything larger exists — that question is put by ruling every one of the
# smallest repairs out, so it needs all of them.
#
# Truncation is deliberately reported as `:some` — the same thing a family that
# is not a product reports — because the sentence is the same either way: there
# are repairs the menus do not offer. It says less than it could, in a case that
# is rare enough (the ceiling is an order of magnitude above what was measured)
# that a third sentence would be a shape of report nobody has read.
const MAX_REPAIRS = 256

# "At most `k` of `lits` hold", as clauses over counting variables: Sinz's
# sequential counter, where the variable for `(i, j)` says that at least `j` of
# the first `i` literals hold. Nothing forces one of those variables down, only
# up, so the last clause of each row — which refuses a literal once `k` are
# already counted — is the whole of what is being asked, and any assignment with
# `k` or fewer literals true extends to a model by counting honestly.
#
# The counter is built for the bound it is asked about rather than for the whole
# range, so a small bound costs a small counter, and nothing is asked at all
# when the bound leaves the literals alone. `counts` is the caller's own store
# of counting variables, extended here when a bound needs more than it holds: a
# variable outlives the frame that gave it meaning and the solver goes on
# deciding it afterwards, so a search that raises the bound step by step hands
# the same store back and leaves one counter behind rather than one per step.
function add_at_most!(sat::SAT, lits::Vector{Int}, k::Int, counts::Vector{Int})
    n = length(lits)
    @assert k ≥ 1
    k < n || return sat
    while length(counts) < (n-1)*k
        push!(counts, sat_new_variable(sat))
    end
    at(i, j) = counts[(i-1)*k + j]
    function clause(ls::Int...)
        for l in ls
            sat_add_var(sat, l)
        end
        sat_add(sat)
    end
    clause(-lits[1], at(1, 1))
    for j = 2:k
        clause(-at(1, j))
    end
    for i = 2:n-1
        clause(-lits[i], at(i, 1))
        clause(-at(i-1, 1), at(i, 1))
        for j = 2:k
            clause(-lits[i], -at(i-1, j-1), at(i, j))
            clause(-at(i-1, j), at(i, j))
        end
        clause(-lits[i], -at(i-1, k))
    end
    clause(-lits[n], -at(n-1, k))
    return sat
end

# Every repair of the smallest size there is, and nothing larger paid for along
# the way.
#
# A repair is the candidates a model leaves unsatisfied, so "no repair larger
# than `k`" is a bound on how many of them a model may leave — one at-most-`k`
# constraint over their negations. Under that bound the instance is satisfiable
# exactly when a repair that small exists, and its minimal correction sets are
# exactly the query's own repairs of that size or less: a model the bound
# admits is one the query admits, and anything smaller than a repair the bound
# admits fits under the bound too, so minimal here and minimal at large are the
# same thing.
#
# So the size of the smallest repair is found by raising the bound until the
# instance gives way, and enumerating at that bound enumerates exactly the
# smallest repairs. The sizes seen in practice are 2 to 6, so counting up costs
# a handful of solves and keeps every counter as small as the answer; a bound of
# `k-1` that failed is what makes `k` the smallest size, so what comes back at
# `k` is all of the smallest repairs, not merely some.
#
# Dropping every candidate is a repair, so the last bound worth trying is the
# one that forbids nothing — which is the plain enumeration, and also what a
# query with no candidate to give up at all wants.
function smallest_repairs(sat::SAT, candidates::Vector{Int}, limit::Int)
    relaxed = Int[-l for l in candidates]
    counts = Int[]
    for k = 1:length(candidates)-1
        repairs = with_temp_clauses(sat) do
            add_at_most!(sat, relaxed, k, counts)
            sat_solve(sat) ? sat_mcses(sat, candidates; limit) : nothing
        end
        repairs === nothing || return repairs
    end
    return sat_mcses(sat, candidates; limit)
end

# Whether the query has a repair that costs more than the smallest ones do.
#
# `smallest` has to be all of the smallest repairs. Ask for a model that leaves
# a candidate of each of them satisfied: one clause per repair, the plain
# disjunction of its candidates, since the repair a model witnesses is the
# candidates it does *not* satisfy. What that model witnesses is then a repair
# holding none of the smallest ones, so the minimal repair inside it is none of
# them either, and — all of the smallest being here — is a larger one.
#
# There is such a model whenever a larger repair exists: the model witnessing
# that repair leaves exactly it unsatisfied, and no smallest repair sits inside
# it, since one minimal repair never contains another. So the answer is exact.
#
# Each candidate says something the query asked for — this package installed,
# this package's emptied classes left empty — so the candidates a model
# satisfies are the ones the repair it witnesses does not touch, and ruling a
# repair out is a disjunction of that repair's own literals rather than of their
# negations.
function larger_repair_exists(sat::SAT, smallest::Vector{Vector{Int}})
    with_temp_clauses(sat) do
        for repair in smallest
            for l in repair
                sat_add_var(sat, l)
            end
            sat_add(sat)
        end
        sat_solve(sat)
    end
end

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
function availability_fact(sat::SAT{P,V}, prob::Problem{P}, p::P) where {P,V}
    info_p = sat.info[p]
    # Shadows belong in this fact and nowhere else. They are versions
    # redundancy elimination took of its own accord, so a fact built from what
    # survives would say the query leaves less than it does -- and the whole
    # point of stating what is left is that the user can act on it. Whatever
    # excludes the class excludes everything it shadows, so a chain that rules
    # the class out rules them out too, and each one's own kinds are a question
    # for the problem rather than for the class hosting it.
    #
    # The order is the version type's, not the provider's. A package lists its
    # versions best first, which is a claim about preference and is the
    # provider's to make; a range is read low to high, which is a claim about
    # order and is the type's. They coincide for version numbers and need not
    # in general, and it is the second one a report wants.
    members = copy(info_p.versions)
    for sh in info_p.shadows
        append!(members, sh)
    end
    sort!(members; rev = true)
    return Availability{P,V}(p, members,
        Vector{Symbol}[exclusion_kinds(prob, p, v) for v in members])
end

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

## why one package matters to another
#
# A chain says what the query asked for and what its constraints left of the
# packages, and between the two there is a step it cannot take on its own: that
# some versions of one package need another, at some versions of it. The user
# wrote the requirement and wrote the bound; the dependency is the one link in
# the story they could not have written, and until it is said the report is a
# pair of facts with nothing between them.
#
# It is read off the universe, which records per class which packages that
# class depends on and which classes of them it admits. What has to be found is
# *which* steps to say, and that is a walk: start at the conflict's
# requirements with the versions this query admits of them, follow the
# dependencies every one of those versions has — the ones no choice of version
# can avoid — and narrow each package reached by what the packages already
# reached admit of it. A package the walk leaves with no version is where the
# story ends, and the steps that got there are the story.
#
# The walk answers with a set of versions per package rather than pairing them
# up, so it says less than the solver does: a package it leaves with nothing
# really has nothing, and a package it leaves something with may still be
# unreachable. That is the safe direction. Nothing it emits is a claim about
# satisfiability — the chain's requirements and availabilities carry that, and
# the solver checked them — and every fact it does emit is a statement about
# the registry and this query's constraints that can be read straight back off
# the universe.

# the classes of `p` this query admits — the ones it left a version in
admitted_classes(sat::SAT{P}, p::P) where {P} =
    Bool[!iszero(r) for r in sat.reps[p]]

# where `q` sits in `p`'s dependency columns, or nothing if `p` never needs it
dep_column(info_p::PkgInfo{P}, q::P) where {P} =
    (k = searchsortedfirst(info_p.depends, q);
     k ≤ length(info_p.depends) && info_p.depends[k] == q ? k : nothing)

# does every class in `live` of `p` depend on `q`? Only then is `q` forced by
# `p`: a class that does without it is a choice the story would have to rule
# out before it could claim `q` has to be installed
function forces(info_p::PkgInfo{P}, q::P, live::Vector{Bool}) where {P}
    k = dep_column(info_p, q)
    k === nothing && return false
    return all(!live[c] || info_p.conflicts[c, k] for c in eachindex(live))
end

# the classes of `q` that some class in `live` of `p` admits. A package the
# two never interact over is one `p` says nothing about, so all of it is
# admitted
function admissible(sat::SAT{P}, p::P, q::P, live::Vector{Bool}) where {P}
    info_p = sat.info[p]
    n_q = nclasses(sat.info[q])
    # the first partner's block of columns starts at offset 0, so whether
    # there is a block at all is a question about the keys
    haskey(info_p.interacts, q) || return fill(true, n_q)
    b = info_p.interacts[q]
    return Bool[any(live[c] && !info_p.conflicts[c, b+j] for c in eachindex(live))
                for j = 1:n_q]
end

# The walk. Every pass reads one snapshot and writes the next, so what a
# package's set was when it starved is a set the report can quote; taking each
# package's answer as it comes would leave the story depending on the order the
# packages were visited in.
#
# A pass first extends the forced set — a package every live class of a forced
# package depends on has to be installed too — and then narrows every forced
# package by every other forced package's bounds on it, dependency or not,
# since two packages that are both installed have to agree however the bound
# got there. Both directions are monotone, so the walk settles; it stops at the
# first pass that leaves a package with nothing, which is a complete story and
# the shortest one this walk has.
#
# Returns the snapshot the last pass read, the package each was first forced
# by, and, per starved package, the classes it had before and the packages
# whose bounds took them.
function forced_descent(sat::SAT{P}, reqs::Vector{P}) where {P}
    live = Dict{P,Vector{Bool}}(p => admitted_classes(sat, p) for p in reqs)
    parent = Dict{P,P}()
    forced = copy(reqs)
    starved = Dict{P,Tuple{Vector{Bool},Vector{P}}}()
    for _ = 1:length(sat.info)
        snap = live
        # to closure, so that a package forced by one forced this pass is
        # reached this pass too — the narrowing below is what a pass is for,
        # and it has nothing to say about a package it has not met
        i = 1
        while i ≤ length(forced)
            p = forced[i]
            i += 1
            any(snap[p]) || continue
            for q in sat.info[p].depends
                haskey(snap, q) && continue
                forces(sat.info[p], q, snap[p]) || continue
                push!(forced, q)
                snap[q] = admitted_classes(sat, q)
                parent[q] = p
            end
        end
        live = Dict{P,Vector{Bool}}()
        changed = false
        for q in forced
            admits = admitted_classes(sat, q)
            new = copy(admits)
            blame = P[]
            for p in forced
                (p == q || !any(snap[p])) && continue
                haskey(sat.info[p].interacts, q) || continue
                a = admissible(sat, p, q, snap[p])
                # everything that takes a version this query still admits,
                # whether or not one is left by the time it is asked. Which of
                # them the story needs is a question for later, and stopping
                # at the first would name whichever came up first
                any(admits[c] && !a[c] for c in eachindex(a)) && push!(blame, p)
                new .&= a
            end
            live[q] = new
            !any(new) && (starved[q] = (snap[q], blame))
            new == snap[q] || (changed = true)
        end
        isempty(starved) || return snap, parent, starved
        changed || break
    end
    return live, parent, starved
end

# Do the bounds of `blamed` between them leave `q` nothing? What a story about
# `q` claims, and what makes each of the packages it names load-bearing
function starves(sat::SAT{P}, snap::Dict{P,Vector{Bool}}, q::P,
                 blamed::Vector{P}) where {P}
    left = admitted_classes(sat, q)
    for r in blamed
        left .&= admissible(sat, r, q, snap[r])
    end
    return !any(left)
end

# The fewest of `blame` that still leave `q` nothing, so that a story names no
# package it could have done without. Dropping is attempted on the packages the
# query never asked for first, so that what survives is what the user can act
# on rather than whichever package the walk happened to meet first
function fewest_blamed(sat::SAT{P}, snap::Dict{P,Vector{Bool}}, q::P,
                       blame::Vector{P}, reqs::Vector{P}) where {P}
    keep = copy(blame)
    order = [P[r for r in blame if r ∉ reqs]; P[r for r in blame if r ∈ reqs]]
    for r in order
        rest = P[x for x in keep if x != r]
        starves(sat, snap, q, rest) && (keep = rest)
    end
    return keep
end

# the forcing steps from a requirement down to `p`, in that order
function forcing_path(parent::Dict{P,P}, p::P) where {P}
    path = Tuple{P,P}[]
    while haskey(parent, p)
        pushfirst!(path, (parent[p], p))
        p = parent[p]
    end
    return path
end

# What `pkg` says about `q`, as facts. One fact per set of `pkg`'s live classes
# that agree about `q`, split again wherever the versions they stand for are
# not a run of what `pkg` offers — so a fact's versions read as one range and
# every version in that range says what the fact says. Versions this query
# forbids are left out: the fact is about the choice the resolve had.
#
# `bound` is whether the bound is worth stating at all. On the last step of a
# story it is the whole point — it is what the query's own constraints have
# left nothing of. On the steps that lead there it is decoration, and worse
# than that: versions that agree about needing a package often disagree about
# which of it they will take, so a step that carried its bound would be a
# paragraph where a sentence was wanted, and the versions it names are named
# again as the subject of the next step anyway.
#
# Whether the versions a fact speaks for depend on `q` is part of what they
# have to agree about, because it is what decides which way round the fact may
# be said. A registry records compatibility symmetrically, so a bound can be
# attributed to the package that declared the dependency and to no other; the
# versions that merely have to get along with `q` get the symmetric sentence,
# which claims nothing about who said so.
function relation_facts(
    sat  :: SAT{P,V},
    prob :: Problem{P},
    pkg  :: P,
    live :: Vector{Bool},
    q    :: P,
    bound :: Bool,
) where {P,V}
    info_p, info_q = sat.info[pkg], sat.info[q]
    k = dep_column(info_p, q)
    # a package this query has left nothing of is starved by every bound
    # alike, so which bound `pkg` declares is not what the story turns on
    bounded = bound && any(admitted_classes(sat, q))
    # what a group of classes says about `q`, and the versions it stands for
    Says = Tuple{Bool,Union{Nothing,Vector{Bool}}}
    groups = Vector{Pair{Says,Vector{Bool}}}()
    for c in eachindex(live)
        live[c] || continue
        needs_it = k !== nothing && info_p.conflicts[c, k]
        a = nothing
        if bounded
            cs = admissible(sat, pkg, q, Bool[i == c for i in eachindex(live)])
            a = fill(false, length(info_q.versions))
            for (j, ok) in enumerate(cs), i in info_q.members[j]
                a[i] = ok
            end
        end
        key = (needs_it, a)
        g = findfirst(x -> first(x) == key, groups)
        g === nothing && (push!(groups, key => fill(false, length(info_p.versions)));
                          g = length(groups))
        for i in info_p.members[c]
            isempty(exclusion_kinds(prob, pkg, info_p.versions[i])) &&
                (last(groups[g])[i] = true)
        end
    end
    facts = Fact[]
    offering = copy(info_q.versions)
    # a run reaching the newest (or oldest) version this query admits can be
    # stated as an open bound: there is nothing above (or below) it to wonder
    # about. What the query excludes is said by the availability fact instead.
    admitted = Bool[isempty(exclusion_kinds(prob, pkg, v)) for v in info_p.versions]
    hi = something(findfirst(admitted), 1)
    lo = something(findlast(admitted), length(admitted))
    for ((needs_it, allowed), mask) in groups
        # versions that do not need `q`, and have nothing to say about it
        # either, contribute nothing to the story
        needs_it || (allowed !== nothing && !all(allowed)) || continue
        i = 1
        while i ≤ length(mask)
            if mask[i]
                j = i
                while j < length(mask) && mask[j+1]
                    j += 1
                end
                vers = info_p.versions[i:j]
                push!(facts, needs_it ?
                    Dependency{P,V}(pkg, vers, q, offering, allowed,
                        i == hi, j == lo) :
                    Incompatibility{P,V}(pkg, vers, q, offering, allowed,
                        i == hi, j == lo))
                i = j
            end
            i += 1
        end
    end
    return facts
end

# Why the packages a conflict names have anything to do with each other, as
# facts about the registry. The walk from the conflict's requirements ends at a
# package with nothing left; the story is the packages whose bounds left it
# nothing, and the forcing steps that put each of them in the query's way.
#
# A step is said in the direction of the dependency where there is one, since
# that is the only direction a bound can honestly be attributed in. Where a
# package is starved by one other and depends on it, the same relation reads
# better the other way about — "this is what these versions need" rather than
# "this is what those versions leave" — and the walk's own narrowing is what
# makes the two the same fact.
#
# The chain already says which packages the query emptied, so a walk that ends
# at one of them is telling this conflict's story rather than some other one;
# that is what `emptied` picks between when several packages starve at once.
function chain_relations(
    sat     :: SAT{P,V},
    prob    :: Problem{P},
    reqs    :: Vector{P},
    emptied :: Vector{P},
) where {P,V}
    snap, parent, starved = forced_descent(sat, reqs)
    isempty(starved) && return Fact[]

    # Which story to tell, when the walk leaves several packages with nothing.
    # Bounds are symmetric, so a package and the one that starves it usually
    # starve each other, and the choice is which way round to say it: a
    # requirement is where the story starts rather than where it ends, and a
    # story about a package this query emptied is a story about this conflict
    # rather than some other one. Then the shortest, and then by name, so that
    # a tie does not depend on iteration order.
    best = nothing
    for (q, (before, blame)) in starved
        blamed = fewest_blamed(sat, snap, q, blame, reqs)
        steps = length(blamed) + length(forcing_path(parent, q)) +
            sum(r -> length(forcing_path(parent, r)), blamed; init = 0)
        touches = q ∈ emptied || any(∈(emptied), blamed)
        key = (q ∈ reqs, !touches, steps, string(q))
        (best === nothing || key < best[1]) && (best = (key, q, before, blamed))
    end
    _, q, before, blamed = best
    # ... and, of the requirements the fewest left out, the ones that starve
    # `q` all by themselves. Any one of them is enough to make it a conflict,
    # so minimality names one and drops the rest — and a reader would fix that
    # one bound and never learn their own other requirement had been broken
    # too. This is the same union of minimal reasons the requirements of a
    # conflict are gathered by. A package the query never asked for earns no
    # such sentence: it is not something the reader is owed an account of, and
    # a story that named every one of them would be a list, not a story
    append!(blamed, P[r for r in starved[q][2] if r ∈ reqs && r ∉ blamed &&
                      starves(sat, snap, q, P[r])])
    # ... back in the order the walk met them, which starts with the
    # requirements in the order the query gave them
    blamed = P[r for r in starved[q][2] if r ∈ blamed]

    # one relation per package blamed, said the better way round when there is
    # only one of them: with several, what matters is that their demands do not
    # overlap, and that is only visible with each stated as a demand on `q`
    pairs = Tuple{P,Vector{Bool},P}[]
    for r in blamed
        length(blamed) == 1 && forces(sat.info[q], r, before) ?
            push!(pairs, (q, before, r)) :
            push!(pairs, (r, snap[r], q))
    end
    # A package the story speaks of has to be in the query's way to begin with,
    # and the forcing steps are what put it there. The package a relation is
    # *about* needs no path of its own when the relation is a dependency on it,
    # since that is what puts it there; where none of them is, `q` gets a
    # walk of its own, and nothing is said about it beyond how it got here
    any(b == q && forces(sat.info[a], q, live) for (a, live, b) in pairs) ||
        pushfirst!(pairs, (q, before, q))
    facts = Fact[]
    said = Set{Tuple{P,P}}((a, b) for (a, _, b) in pairs)
    for (a, live, b) in pairs
        for (x, y) in forcing_path(parent, a)
            (x, y) in said && continue
            push!(said, (x, y))
            # with the bound: the next step's subject is the run this one
            # leaves, so a step that keeps its bound to itself makes that run
            # appear from nowhere. Where there is no bound to state -- either
            # because `pkg` declares none, or because the query has left `y`
            # nothing and every bound would starve it alike -- `relation_facts`
            # says only that `y` is needed, which is what it said before
            append!(facts, relation_facts(sat, prob, x, snap[x], y, true))
        end
        a == b || append!(facts, relation_facts(sat, prob, a, live, b, true))
    end
    return facts
end

# The story of one conflict: the query's own facts that this conflict's fix is
# what rescues, and why each of them needs rescuing.
#
# The smallest unsatisfiable set of them is where it starts. Deletion is
# attempted in the order given, so putting the requirements first drops the ones
# that set does not need and keeps the facts the user can act on.
#
# But "does not need" is the trap. Several requirements can fail for one and the
# same reason, and one of them is enough to make that reason a conflict — so the
# smallest set names one of them, arbitrarily, and the others go unmentioned
# though this conflict's fix is what rescues them too. A user would fix the
# bound and never learn the rest had been broken. So every requirement the set
# left out is asked whether it can be satisfied on its own against the packages
# as this conflict finds them; the ones that cannot are exactly the ones this
# fix rescues, since the fixes of the other conflicts are already in force here
# and carrying this one out satisfies the query. Each of those brings its own
# smallest explanation, and the chain is all of them at once.
#
# One `sat_mus` call per requirement, which answers both questions: empty when
# the requirement is satisfiable, and its explanation when it is not. What comes
# back is a subsequence of what was asked, so the requirements are in the
# problem's requirement order and the packages in package order.
function conflict_story(
    sat        :: SAT{P,V},
    prob       :: Problem{P},
    vm         :: VarMap{P},
    emptied    :: Dict{Int,P}, # per package literal, the package
    candidates :: Vector{Int}, # the literals a repair may withdraw
) where {P,V}
    named = Set{Int}(sat_mus(sat, candidates))
    pkgs = Int[l for l in candidates if haskey(emptied, l)]
    for l in candidates
        (haskey(emptied, l) || l in named) && continue
        union!(named, sat_mus(sat, Int[l; pkgs]))
    end
    reqs = P[]
    chain = Fact[]
    avails = Fact[]
    for l in candidates
        l in named || continue
        p = get(emptied, l, nothing)
        if p === nothing
            q, _ = decode(vm, l)
            push!(reqs, q)
            push!(chain, Requirement(q))
        else
            push!(avails, availability_fact(sat, prob, p))
        end
    end
    # the registry's own part of the story goes between the packages it speaks
    # for and the ones it speaks about, so that the chain reads as a chain:
    # you asked for this, this is what it needs, this is what you left of it
    rels = chain_relations(sat, prob, reqs, P[f.pkg for f in avails])
    subjects = Set{P}(f.pkg for f in rels)
    for f in avails
        (f::Availability).pkg in subjects && push!(chain, f)
    end
    append!(chain, rels)
    for f in avails
        (f::Availability).pkg in subjects || push!(chain, f)
    end
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
function order_choices!(
    choices    :: Vector{Vector{Int}},
    candidates :: Vector{Int},
    nreqs      :: Int,
)
    at = Dict{Int,Int}(l => i for (i, l) in enumerate(candidates))
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
            candidates = [req_lits; pkg_lits]
            cheapest = smallest_repairs(sat, candidates, MAX_REPAIRS)
            choices, whole = factor_repairs(cheapest)
            order_choices!(choices, candidates, length(req_lits))
            # what the menus leave out: nothing when they reach every cheapest
            # repair and no repair costs more than the cheapest. The second is
            # a question to put to the instance, and putting it takes every
            # cheapest repair, so a truncated enumeration has to leave both
            # unanswered
            others =
                length(cheapest) ≥ MAX_REPAIRS || !whole ? :some :
                larger_repair_exists(sat, cheapest) ? :larger : :none
            # A story explains one choice by settling all the others: drop the
            # rest of one repair from what may be assumed, and what is left is
            # still unsatisfiable, minimally so once this choice is made. The
            # repair's own minimality is what puts this choice's literal in the
            # conflict and keeps the other choices' literals out of it.
            settled = Int[first(c) for c in choices]
            for (i, choice) in enumerate(choices)
                rest = Set{Int}(settled[j] for j in eachindex(settled) if j ≠ i)
                push!(stories, conflict_story(sat, prob, vm, emptied,
                    Int[l for l in candidates if l ∉ rest]))
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
