"""
    Resolver.Diagnostics

Why an unsatisfiable query fails, and what its user could change.

A report owes three things: the **fixes** — how many independent things are
wrong and, for each, a menu of alternatives such that any combination of
choices, one per menu, repairs the query while giving up as little as any
repair can; the **explanations** — for each conflict, which of the user's own
facts collide, at which package, and through which of the registry's
statements; and the **witnesses** — what taking each fix yields.

The facts a query is made of are per package: *requiring* it, and the query's
own limits on *which versions of it* are admissible. Package granularity is the
finest grain that stays a user action — no edit gives back one excluded version
while keeping its siblings excluded. A **reason** is a minimal unsatisfiable set
of facts, a **repair** a minimal set whose withdrawal is satisfiable, and the
two are hitting-set duals; the menus are the factors the cheapest repairs
decompose into, and each conflict explains the reasons its whole menu is
contained in.

See the manual's *Explaining an unsatisfiable resolve* for the theory this
implements, including the proof of every claim above.
"""
module Diagnostics

using ..Resolver: Resolver, SAT, Problem, PkgInfo, Universe, PicoSAT, Relation,
    nclasses, installed_lit, forbidden_lit, sat_assume_var, sat_solve,
    sat_new_variable, sat_add_var, sat_add, with_classes_relaxed,
    with_temp_clauses, exclusion_kinds, relax, resolve
using ..Resolver.Clauses: Clauses, Clause, Lit, literal, clause, packages,
    isbottom, subsumes, absent, present, resolve_raw, resolve_on, clause_phrase,
    range_phrase, selected, unselected, nversions
using ..Resolver.UnsatCores: sat_mus

export Diagnosis, Conflict, Fix, Action, Line, action_phrase

## what a report is made of

"""
    Action(kind, pkg)

One thing a user could do: `:drop` a requirement on `pkg`, or lift the
constraint of kind `kind` for `pkg`. The kinds are the query's own
([`Problem`](@ref Resolver.Problem)), so an action is always something the
reader can carry out by editing what they wrote.
"""
struct Action{P}
    kind :: Symbol
    pkg  :: P
end

Base.:(==)(a::Action, b::Action) = a.kind === b.kind && a.pkg == b.pkg
Base.hash(a::Action, h::UInt) = hash(a.pkg, hash(a.kind, hash(:Action, h)))
Base.show(io::IO, a::Action) = print(io, "Action(", repr(a.kind), ", ",
                                     repr(a.pkg), ")")

"""
    Line(clause, through, given, proof = 0, pivot = nothing)

One printed statement. `clause` is the whole of what it says; `through` names
the packages an elimination reached it by — a courtesy pointer, carrying no
claim of its own, since a flat page cannot show how two packages reach each
other. `given` marks the query's own facts, which every proof on the page
shares; `proof` numbers the reason a derived line argues for, and `pivot` the
package its meet is taken at.
"""
struct Line{P}
    clause  :: Clause{P}
    through :: Vector{P}
    given   :: Bool
    proof   :: Int
    pivot   :: Union{Nothing,P}
end

Line{P}(c::Clause{P}, t::Vector{P}, g::Bool) where {P} =
    Line{P}(c, t, g, 0, nothing)
Line{P}(c::Clause{P}, t::Vector{P}, g::Bool, n::Integer) where {P} =
    Line{P}(c, t, g, Int(n), nothing)

"""
    Fix(actions, solution)

One entry of a menu: the actions to carry out, and what the resolver answers
once they are — with every other conflict settled the first way its own menu
offers. The versions are the resolver's optimising answer for the withdrawn
query, never a model found during diagnosis: diagnosis decides what is true,
the resolver decides what is chosen.
"""
struct Fix{P,V}
    actions  :: Vector{Action{P}}
    solution :: Dict{P,V}
end

"""
    Conflict(reqs, lines, versions, excluded, fixes)

One independent thing that is wrong: the requirements it answers for, the lines
that prove it, the version list each named package is spoken of in, which of
the query's constraint kinds exclude which of those versions, and the menu of
alternatives that settle it.
"""
struct Conflict{P,V}
    reqs     :: Vector{P}
    lines    :: Vector{Line{P}}
    versions :: Dict{P,Vector{V}}
    excluded :: Dict{P,Vector{Vector{Symbol}}}
    fixes    :: Vector{Fix{P,V}}
end

"""
    Diagnosis

What [`resolve`](@ref Resolver.resolve) answers when the query cannot be
satisfied: the conflicts, each of which must be fixed, and what the menus leave
out — `:none` when every repair is a combination of them, `:larger` when the
ones outside all give up more, and `:some` when equally cheap repairs lie
outside. `truncated` records that the search for reasons was cut short, so the
account of some conflict may be incomplete.

`show`ing one prints the report.
"""
struct Diagnosis{P,V}
    conflicts :: Vector{Conflict{P,V}}
    others    :: Symbol # :none, :larger, :some
    truncated :: Bool
end

# a diagnosis rebuilt by a caller — renamed, filtered, whatever — is not one
# whose search was cut short, so the disclosure defaults off
Diagnosis(conflicts::Vector{Conflict{P,V}}, others::Symbol) where {P,V} =
    Diagnosis{P,V}(conflicts, others, false)

## the universe as the clause logic sees it

"""
    clause_versions(sat, p) :: Vector{V}

The versions of `p` a clause's literal indexes — the ones the universe this
query was run against still holds. A literal has one more slot than this, for
⊥.
"""
clause_versions(sat::SAT{P,V}, p) where {P,V} =
    haskey(sat.info, p) ? sat.info[p].versions : V[]

"""
    clause_of(sat, r::Relation) :: Union{Clause{P}, Nothing}

One registry statement as the clause it is: the versions of `r.pkg` it speaks
for, complemented and admitting ⊥ — the statement binds only where that package
is at one of them — together with what it asks of `r.other`, which is presence
(a dependency) or everything but a run of versions (a compatibility bound).
"""
function clause_of(sat::SAT{P,V}, r::Relation{P}) where {P,V}
    ip, iq = sat.info[r.pkg], sat.info[r.other]
    np, nq = length(ip.versions), length(iq.versions)
    a = falses(np + 1)
    for c in r.classes, j in ip.members[c]
        a[j] = true
    end
    for j = 1:np
        a[j] = !a[j]
    end
    a[np+1] = true
    b = falses(nq + 1)
    if r.dep
        b[1:nq] .= true
    else
        for c in r.others, j in iq.members[c]
            b[j] = true
        end
        for j = 1:nq
            b[j] = !b[j]
        end
        b[nq+1] = true
    end
    return clause([r.pkg => Lit(a), r.other => Lit(b)])
end

"""
    clauses_satisfiable(sat, clauses) :: Bool

Can these clauses hold together, over the versions the universe leaves each
package they name? Asked of the clauses alone — no registry, no query — because
whether a set of statements contradicts is a question about the statements.
"""
function clauses_satisfiable(sat::SAT{P,V}, cs) where {P,V}
    ps = P[]
    for c in cs, p in packages(c)
        p in ps || push!(ps, p)
    end
    isempty(ps) && return !any(isbottom, cs)
    pico = PicoSAT.init()
    try
        var = Dict{P,Int}()
        n = 0
        for p in ps
            var[p] = n
            n += length(clause_versions(sat, p)) + 1
        end
        PicoSAT.adjust(pico, max(n, 1))
        function add(lits)
            for l in lits
                PicoSAT.add(pico, l)
            end
            PicoSAT.add(pico, 0)
        end
        for p in ps
            k = length(clause_versions(sat, p)) + 1
            add(Int[var[p] + i for i = 1:k])
            for i = 1:k, j = i+1:k
                add(Int[-(var[p] + i), -(var[p] + j)])
            end
        end
        for c in cs
            add(Int[var[p] + i for (p, m) in c.lits for i in eachindex(m) if m[i]])
        end
        return PicoSAT.sat(pico) != PicoSAT.UNSATISFIABLE
    finally
        PicoSAT.reset(pico)
    end
end

# Does `base` entail `c`? Deny the clause — every package it names confined to
# what its literal does not allow — and ask whether the rest can still hold.
function entails(sat::SAT{P,V}, base::Vector{Clause{P}}, c::Clause{P}) where {P,V}
    deny = Clause{P}[]
    for (p, m) in c.lits
        d = clause([p => Lit(.~m.bits)])
        d === nothing && return true
        push!(deny, d)
    end
    return !clauses_satisfiable(sat, Clause{P}[base; deny])
end

# The smallest vocabulary outside a line that still proves it, in the order the
# route reaches it. The packages an elimination passed through are not all
# packages the line needs, and the list is what the report offers in place of
# the argument, so it should be the argument's own vocabulary and no more. And
# order matters as much as membership: named alphabetically the list is a set,
# and a set is not what the reader wants — they want the way from what the
# line argues from to what it concludes, so the route is walked out from the
# line's own packages over the statements the proof may use, nearest first.
function trim_route(sat::SAT{P,V}, base::Vector{Clause{P}}, c::Clause{P},
                    through::Vector{P}) where {P,V}
    named = Set{P}(packages(c))
    ok(v) = entails(sat, Clause{P}[b for b in base
        if all(q -> q in named || q in v, packages(b))], c)
    v = Set{P}(q for q in through if q ∉ named)
    isempty(v) && return P[]
    if ok(v)
        for q in sort!(collect(v))
            w = setdiff(v, [q])
            ok(w) && (v = w)
        end
    end
    isempty(v) && return P[]
    usable = Clause{P}[b for b in base
                       if all(q -> q in named || q in v, packages(b))]
    seen = copy(named); out = P[]
    front = sort!(collect(named))
    while !isempty(front)
        step = P[]
        for b in usable
            ps = packages(b)
            any(q -> q in front, ps) || continue
            for q in ps
                q in v && q ∉ seen && q ∉ step && push!(step, q)
            end
        end
        sort!(step)
        append!(out, step); union!(seen, step)
        front = step
    end
    for q in sort!(collect(v))
        q in seen || push!(out, q)
    end
    return out
end

## elimination

# a clause on its way through a projection: what it says, the packages whose
# elimination derived it, and which of the reason's facts its derivation used.
# The mask over-attributes — a derivation that happened to use a fact is a
# witness that the fact may be needed, not a proof of it — which is sound for
# everything printed, since a side stated from a larger support is still true.
struct Item{P}
    clause :: Clause{P}
    route  :: Vector{P}
    mask   :: UInt64
end

# how much subset enumeration one elimination may do before the projection is
# abandoned — on overflow the caller falls back to printing the core, a true
# but unshaped account being better than a shaped partial one
const ELIM_NODES = 30_000
const ELIM_CLAUSES = 2_000
const ELIM_ARITY = 40
const ELIM_SET = 5

# Every minimal subset of `bits` whose intersection is empty, up to
# `ELIM_SET` members. Minimal subsets and not merely emptying ones: the
# resolvent of a larger emptying set is subsumed by that of a minimal one
# inside it, so restricting to these loses nothing — and pairs do not suffice,
# since three set-valued literals can have empty intersection while every two
# overlap.
#
# Three prunes keep this from walking the powerset, two of them lossless. A
# clause that does not shrink the running intersection is redundant in every
# minimal set the branch could reach, so it is skipped where it stands. A
# branch whose intersection survives even the conjunction of everything still
# ahead of it can never empty, so it is cut — this is what kills the long
# chains of statements that all admit ⊥. The size cap alone loses something:
# a minimal emptying set wider than `ELIM_SET` is not found, and the
# projection may then fail to close at this pivot — sound (nothing false is
# ever derived, only less), and the search moves to the next pivot or falls
# back to the core. Convex bounds cannot need more than three (the manual's
# proposition on the shape of convex meets), so the cap is generous.
function minimal_emptying(bits::Vector{BitVector}, budget::Ref{Int})
    n = length(bits)
    sols = Vector{Vector{Int}}()
    chosen = Int[]
    suffix = [trues(length(bits[1])) for _ = 1:n+1]
    for i = n:-1:1
        suffix[i] = suffix[i+1] .& bits[i]
    end
    function is_minimal()
        for j in eachindex(chosen)
            cur = trues(length(bits[1]))
            for k in eachindex(chosen)
                k == j || (cur .&= bits[chosen[k]])
            end
            any(cur) || return false
        end
        return true
    end
    function rec(start::Int, cur::BitVector)
        any(cur .& suffix[start]) && return true
        length(chosen) ≥ ELIM_SET && return true
        for i = start:n
            budget[] -= 1
            budget[] ≤ 0 && return false
            nxt = cur .& bits[i]
            nxt == cur && continue
            push!(chosen, i)
            if !any(nxt)
                is_minimal() && push!(sols, copy(chosen))
            elseif i < n
                rec(i + 1, nxt) || (pop!(chosen); return false)
            end
            pop!(chosen)
        end
        return true
    end
    isempty(bits) && return sols
    return rec(1, trues(length(bits[1]))) ? sols : nothing
end

# drop what another item already says: an item saying at least as much on no
# more of the reason's facts makes the other one redundant on every page
function reduce_items(items::Vector{Item{P}}) where {P}
    keep = trues(length(items))
    for i in eachindex(items)
        keep[i] || continue
        for j in eachindex(items)
            (i == j || !keep[j]) && continue
            subsumes(items[i].clause, items[j].clause) || continue
            iszero(items[i].mask & ~items[j].mask) || continue
            if items[i].clause == items[j].clause && items[i].mask == items[j].mask
                (length(items[i].route), i) ≤ (length(items[j].route), j) || continue
            end
            keep[j] = false
        end
    end
    return Item{P}[items[i] for i in eachindex(items) if keep[i]]
end

# Two statements that agree everywhere but one package are one statement about
# the union: `A 1.0.0 requires C S` beside `A 1.1.0–1.3.0 requires C S` is
# `A 1.0.0–1.3.0 requires C S`, by resolution on `A` like anything else.
# Carried separately into an elimination each becomes a case of its own, the
# subset enumeration pays for every combination, and the reader is handed a
# dozen lines with the same right-hand side; merged first, the split is over
# what actually differs. A merged support is the union of the parents' — wider
# than either may have needed, which the mask is allowed to be.
function factor_items(items::Vector{Item{P}}, keep = ()) where {P}
    items = copy(items)
    while true
        done = true
        for i in eachindex(items), j in eachindex(items)
            i < j || continue
            c, e = items[i].clause, items[j].clause
            length(c.lits) == length(e.lits) || continue
            diff = P[q for (q, m) in c.lits
                     if (n = e[q]; n === nothing || n != m)]
            length(diff) == 1 || continue
            # merge only across a package on its way out: resolving two of a
            # meet's own sides on the kept package would fold the meet itself
            # into the empty clause, which is a derivation, not a merge
            q = only(diff)
            q in keep && continue
            r = resolve_on(c, e, q)
            (r === nothing || isbottom(r)) && continue
            route = copy(items[i].route)
            for q in items[j].route
                q in route || push!(route, q)
            end
            items[i] = Item{P}(r, sort!(route), items[i].mask | items[j].mask)
            deleteat!(items, j)
            done = false; break
        end
        done && return reduce_items(items)
    end
end

# Eliminate `q`: keep the clauses that do not mention it, and resolve every
# minimal set of the ones that do that leaves it nothing to be. That is the
# coverage condition stated as the rule itself — a version of `q` no resolvent
# accounts for is a gap in the argument, and this cannot leave one.
function eliminate_one(items::Vector{Item{P}}, q::P, budget::Ref{Int}) where {P}
    hit = Item{P}[]
    out = Item{P}[]
    for it in items
        push!(it.clause[q] === nothing ? out : hit, it)
    end
    isempty(hit) && return out
    length(hit) > ELIM_ARITY && return nothing
    sets = minimal_emptying(BitVector[hit[i].clause[q].bits for i in eachindex(hit)],
                            budget)
    sets === nothing && return nothing
    for S in sets
        r = resolve_raw(Clause{P}[hit[i].clause for i in S], q)
        r === nothing && continue
        route = P[q]
        mask = UInt64(0)
        for i in S
            append!(route, hit[i].route)
            mask |= hit[i].mask
        end
        push!(out, Item{P}(r, sort!(unique!(route)), mask))
    end
    return reduce_items(out)
end

# eliminate everything outside `keep`, cheapest package first — any order
# preserves satisfiability, and taking the least-mentioned one keeps the
# intermediate sets small
function eliminate_to(items::Vector{Item{P}}, keep) where {P}
    work = factor_items(items, keep)
    while true
        counts = Dict{P,Int}()
        for it in work, p in packages(it.clause)
            p in keep && continue
            counts[p] = get(counts, p, 0) + 1
        end
        isempty(counts) && return work
        q = first(sort!(collect(counts); by = x -> (last(x), first(x))))[1]
        work = eliminate_one(work, q, Ref(ELIM_NODES))
        work === nothing && return nothing
        work = factor_items(work, keep)
        length(work) > ELIM_CLAUSES && return nothing
    end
end

"""
    project(lines, keep) :: Union{Vector{Line{P}}, Nothing}

The statements `lines` make about the packages in `keep` alone: every package
outside `keep` eliminated, one at a time, by resolving every minimal set of
clauses that leaves it nothing to be. Satisfiability is preserved, so what comes
out contradicts exactly when what went in did.

`nothing` when the elimination would cost more work than it is worth; the
caller then has the clauses it started with, which are true whether or not they
have been shaped.
"""
function project(lines::Vector{Line{P}}, keep) where {P}
    out = eliminate_to(
        Item{P}[Item{P}(l.clause, copy(l.through), UInt64(0)) for l in lines], keep)
    out === nothing && return nothing
    return Line{P}[Line{P}(it.clause, it.route, false) for it in out]
end

## finding the conflicts
#
# The facts a query is made of, the cheapest repairs, the menus those factor
# into, and one explanation per reason each menu owns.
#
# Every question here is a solve under a subset of the *assumptions* standing
# for the query's facts, on an instance whose own deactivations have been
# lifted for the duration — so the registry is what the solver holds and the
# query is what it is asked about. Blame is settled on that instance first,
# with the registry inviolable; only then is the registry's share asked for, on
# a second instance where its statements are individually switchable. Asking
# the other way round lets a minimal answer drop a registry statement and blame
# a requirement in its place, which is minimal and a lie about what the user
# can do.

# a pivot with more classes than this has its sides taken as the elimination
# derived them; below it, each class is put to the instance, so that a side
# states everything the registry entails from its support rather than only what
# this reason's core happened to need
const TIGHTEN_CLASSES = 24

# how many of each enumeration is worth having. The reason walk is exponential
# and the repair enumeration can be too; both truncate, and truncation is
# disclosed rather than hidden (a conflict simply explains fewer of the reasons
# it owns).
const REPAIR_CAP = 64
const REASON_CAP = 24    # reasons the walk records before it stops
const REASON_NODES = 400 # calls the walk makes before it stops
const CONFLICT_REASONS = 4 # reasons one conflict sets out as proofs
# a reason wider than this has more facts than a bitmask carries; its
# explanation falls back to printing the core
const MASK_WIDTH = 62

"""
    VarMap(sat)

Which packages a diagnosis can ask about, and how. `haskey(vm, p)` is whether
the universe holds `p` at all — a requirement it does not is one nothing can
give back. `pkgs` are the packages this query has emptied classes of, in a
fixed order, with the literals that forbid them: those are the only packages a
constraint can be blamed for, since a constraint that spares a member of every
class it touches has changed nothing the registry can tell.
"""
struct VarMap{P}
    known :: Set{P}
    pkgs  :: Vector{P}
    lits  :: Dict{P,Vector{Int}}
end

function VarMap(sat::SAT{P,V}) where {P,V}
    pkgs = P[]
    lits = Dict{P,Vector{Int}}()
    for p in sort!(collect(keys(sat.reps)))
        reps = sat.reps[p]
        ls = Int[forbidden_lit(sat, p, c) for c in eachindex(reps)
                 if iszero(reps[c])]
        isempty(ls) && continue
        push!(pkgs, p)
        lits[p] = ls
    end
    return VarMap{P}(Set{P}(keys(sat.info)), pkgs, lits)
end

Base.haskey(vm::VarMap{P}, p) where {P} = p in vm.known

"""
    with_emptied_packages(body, sat, vm) -> body(pkg_lits, pkg_of)

Run `body` with one fresh literal per package this query emptied classes of,
each of which says "this query's limits on that package are in force". The
definitions live in a frame of their own and are retracted afterwards, so the
instance is as reusable as it was found.

Must be run inside `with_classes_relaxed`:
the point of the literal is to impose by assumption what the query's own frame
imposes unconditionally, and it can only do that once that frame is off.
"""
function with_emptied_packages(body::Function, sat::SAT{P,V},
                               vm::VarMap{P}) where {P,V}
    isempty(vm.pkgs) && return body(Int[], Dict{Int,P}())
    return with_temp_clauses(sat) do
        pkg_lits = Int[]
        pkg_of = Dict{Int,P}()
        for p in vm.pkgs
            v = sat_new_variable(sat)
            for l in vm.lits[p]
                sat_add_var(sat, -v)
                sat_add_var(sat, l)
                sat_add(sat)
            end
            push!(pkg_lits, v)
            pkg_of[v] = p
        end
        body(pkg_lits, pkg_of)
    end
end

# One of the query's facts: requiring a package, or the query's own limits on
# which of its versions are admissible. `lit` assumes it on the blame instance,
# `litx` on the one that can switch registry statements off.
#
# The canonical order is constraints before requirements, each by package: a
# menu prints in it, so a fix that gives something back is offered before one
# that gives a requirement up.
struct Fact{P}
    req  :: Bool
    pkg  :: P
    lit  :: Int
    litx :: Int
end

# solve with exactly `lits` assumed
function assuming(sat::SAT, lits)
    for l in lits
        sat_assume_var(sat, l)
    end
    return sat_solve(sat)
end

# What the registry entails about `pivot` from `support` alone — the strongest
# bound, not merely the one this reason's core happened to give. A core is
# minimal for the whole reason, so a statement the contradiction can do without
# is left out of it, and a side derived from what is left can be weaker than
# the truth: it may permit the package's absence where the registry demands it,
# which is the difference between *constrains* and *requires*. Asking closes
# that gap, and asking costs one solve per value.
function tighten!(sat::SAT{P,V}, bits::BitVector, pivot::P,
                  support::Vector{Int}) where {P,V}
    if bits[end] && !assuming(sat, Int[support; -installed_lit(sat, pivot)])
        bits[end] = false
    end
    info = sat.info[pivot]
    nclasses(info) ≤ TIGHTEN_CLASSES || return bits
    for c = 1:nclasses(info)
        mem = info.members[c]
        any(j -> bits[j], mem) || continue
        assuming(sat, Int[support; -forbidden_lit(sat, pivot, c)]) && continue
        for j in mem
            bits[j] = false
        end
    end
    return bits
end

## the bounded oracle

# "at most `k` of these facts are violated", as clauses over the fact literals
# (Sinz's sequential counter). A model of the instance under this constraint
# violates a correction set of at most `k` facts, and the least `k` at which one
# exists is the size of the cheapest repairs.
function add_at_most!(sat::SAT, lits::Vector{Int}, k::Int)
    n = length(lits)
    (n == 0 || k ≥ n) && return
    unit(a) = (sat_add_var(sat, a); sat_add(sat))
    two(a, b) = (sat_add_var(sat, a); sat_add_var(sat, b); sat_add(sat))
    three(a, b, c) =
        (sat_add_var(sat, a); sat_add_var(sat, b); sat_add_var(sat, c); sat_add(sat))
    if k == 0
        for l in lits
            unit(l)
        end
        return
    end
    s = [Int[sat_new_variable(sat) for _ = 1:k] for _ = 1:n-1]
    two(lits[1], s[1][1])
    for j = 2:k
        unit(-s[1][j])
    end
    for i = 2:n-1
        two(lits[i], s[i][1])
        two(-s[i-1][1], s[i][1])
        for j = 2:k
            three(lits[i], -s[i-1][j-1], s[i][j])
            two(-s[i-1][j], s[i][j])
        end
        two(lits[i], -s[i-1][k])
    end
    two(lits[n], -s[n-1][k])
    return
end

# The cheapest repairs, and how much they cost. Raising the bound one at a time
# stops at the first `k` a model exists at, which is the least size a correction
# set has; at that `k` every model's violated set *is* a repair, and blocking
# each one found enumerates them all.
function min_repairs(sat::SAT, lits::Vector{Int})
    n = length(lits)
    for k = 0:n
        found = with_temp_clauses(sat) do
            add_at_most!(sat, lits, k)
            out = Vector{Vector{Int}}()
            while sat_solve(sat)
                viol = Int[i for i = 1:n if PicoSAT.deref(sat.pico, lits[i]) < 0]
                push!(out, viol)
                length(out) ≥ REPAIR_CAP && break
                for i in viol
                    sat_add_var(sat, lits[i])
                end
                sat_add(sat)
            end
            return out
        end
        isempty(found) || return k, found
    end
    return n, Vector{Int}[]
end

# Is there a repair holding none of the cheapest ones? One question, and its
# answer is the whole of what a report may say about what it leaves out.
function larger_repairs(sat::SAT, lits::Vector{Int}, fmin::Vector{Vector{Int}})
    return with_temp_clauses(sat) do
        for m in fmin
            for i in m
                sat_add_var(sat, lits[i])
            end
            sat_add(sat)
        end
        return sat_solve(sat)
    end
end

## menus

# The factors the cheapest repairs decompose into, or `nothing` where they do
# not decompose. Two facts lie in one factor exactly when they never share a
# repair, so the candidate partition is the components of that relation and
# needs no search; a family that counts right — one fact of each factor in every
# member, and as many members as the product has tuples — *is* the product.
function product_menus(fmin::Vector{Vector{Int}}, used::Vector{Int})
    n = length(used)
    at = Dict{Int,Int}(f => i for (i, f) in enumerate(used))
    shares = falses(n, n)
    for m in fmin, a in m, b in m
        a == b || (shares[at[a], at[b]] = true)
    end
    parent = collect(1:n)
    function root(x::Int)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        return x
    end
    for i = 1:n, j = i+1:n
        shares[i, j] && continue
        ri, rj = root(i), root(j)
        ri == rj || (parent[ri] = rj)
    end
    groups = Dict{Int,Vector{Int}}()
    for i = 1:n
        push!(get!(Vector{Int}, groups, root(i)), used[i])
    end
    factors = Vector{Int}[sort!(g) for g in values(groups)]
    k = length(factors)
    for m in fmin
        length(m) == k || return nothing
        all(f -> count(x -> x in f, m) == 1, factors) || return nothing
    end
    length(fmin) == prod(length, factors; init = 1) || return nothing
    sort!(factors; by = first)
    return factors
end

# The largest part of the family that *is* a product: `k-1` singleton menus and
# one free menu, every selection of which is a member. Any such rectangle is a
# legitimate offer; this takes one of maximal coverage among those searched and
# breaks ties on the fact order, so the same family always yields the same page.
function rectangle_menus(fmin::Vector{Vector{Int}}, used::Vector{Int})
    members = Set{Vector{Int}}(sort(m) for m in fmin)
    best = nothing
    bestkey = nothing
    for m in fmin, c in m
        core = sort!(Int[x for x in m if x != c])
        free = Int[d for d in used if d ∉ core && sort!([core; d]) in members]
        key = (-length(free), core, free)
        if bestkey === nothing || key < bestkey
            best, bestkey = (core, free), key
        end
    end
    best === nothing && return Vector{Int}[], 0
    core, free = best
    menus = Vector{Int}[Int[x] for x in core]
    push!(menus, free)
    sort!(menus; by = first)
    return menus, length(free)
end

## reasons

# Reasons are enumerated by removal, over the whole pool and never a restricted
# one: a reason can hold this conflict's menu *and* facts another menu offers,
# and a pool with those held out cannot contain it. Every reason in the pool is
# recorded — take a minimal unsatisfiable subset, then recurse with each of its
# members removed — until the budget stops the walk, which is disclosed.
function reason_walk(sat::SAT, pool::Vector{Int}, lits::Vector{Int})
    found = Vector{Vector{Int}}()
    seen = Set{Vector{Int}}()
    budget = Ref(REASON_NODES)
    complete = Ref(true)
    index = Dict{Int,Int}(l => i for (i, l) in enumerate(lits))
    function walk(p::Vector{Int})
        if length(found) ≥ REASON_CAP || budget[] ≤ 0
            complete[] = false
            return
        end
        budget[] -= 1
        m = sat_mus(sat, Int[lits[i] for i in p])
        isempty(m) && return
        r = sort!(Int[index[l] for l in m])
        if !(r in seen)
            push!(seen, r)
            push!(found, r)
        end
        for x in r
            walk(Int[y for y in p if y != x])
        end
    end
    walk(pool)
    return found, complete[]
end

## explanation

# The query's own fact, as the clause it is: a requirement demands one of the
# package's versions, a constraint admits the versions it leaves and absence.
# Read straight off the query — which versions each kind rules out is data, not
# a solver question — so it is finer than the instance, which can only empty a
# whole class. Finer is sound: a stronger premise cannot make a contradiction
# less of one.
function fact_clause(sat::SAT{P,V}, prob::Problem{P}, f::Fact{P}) where {P,V}
    vs = clause_versions(sat, f.pkg)
    n = length(vs)
    b = falses(n + 1)
    if f.req
        b[1:n] .= true
    else
        for j = 1:n
            b[j] = isempty(exclusion_kinds(prob, f.pkg, vs[j]))
        end
        b[n+1] = true
    end
    return clause([f.pkg => Lit(b)])
end

# One meet: the pivot, the sides that close at it, and what the page costs
struct Meet{P}
    pivot :: P
    sides :: Vector{Tuple{Clause{P},Vector{P}}}
    key   :: Tuple{Bool,Int,Int,Int,Int,Int,String}
end

# The strongest bound on `pivot` each support entails, and the smallest
# subfamily of them that closes.
#
# The facts on the pivot itself are held out of the derivation and put back as
# sides: they are the query's own lines, already on the page, and leaving them
# in would let a side be stated from a support the printed implication cannot
# name.
function meet_at(
    sat   :: SAT{P,V},
    pivot :: P,
    reason:: Vector{Int},
    facts :: Vector{Fact{P}},
    fcl   :: Dict{Int,Clause{P}},
    core  :: Vector{Clause{P}},
    menu  :: Vector{Int},
) where {P,V}
    bit = Dict{Int,Int}(f => i - 1 for (i, f) in enumerate(reason))
    items = Item{P}[]
    on_pivot = Int[]
    for i in reason
        if facts[i].pkg == pivot
            push!(on_pivot, i)
        else
            push!(items, Item{P}(fcl[i], P[], UInt64(1) << bit[i]))
        end
    end
    for c in core
        push!(items, Item{P}(c, P[], UInt64(0)))
    end
    proj = eliminate_to(items, Set{P}([pivot]))
    proj === nothing && return nothing
    nv = length(clause_versions(sat, pivot))
    sigma = Dict{UInt64,BitVector}()
    routes = Dict{UInt64,Vector{P}}()
    for it in proj
        isbottom(it.clause) && return nothing
        m = it.clause[pivot]
        m === nothing && continue
        if haskey(sigma, it.mask)
            sigma[it.mask] .&= m.bits
            append!(routes[it.mask], it.route)
        else
            sigma[it.mask] = copy(m.bits)
            routes[it.mask] = copy(it.route)
        end
    end
    for (m, b) in sigma
        tighten!(sat, b, pivot, Int[facts[i].lit for i in reason
                                    if (m >> bit[i]) & 1 == 1])
    end
    base = trues(nv + 1)
    for i in on_pivot
        base .&= fcl[i][pivot].bits
    end
    # trim to an irredundant family, dropping the sides that cost most to read
    # first; nothing that matters can be lost this way, since every fact of the
    # reason roots a side of *whatever* family closes
    keep = Set{UInt64}(keys(sigma))
    for m in sort!(collect(keys(sigma)); by = x -> (-count_ones(x), -Int(x)))
        trial = copy(base)
        for x in keep
            x == m || (trial .&= sigma[x])
        end
        any(trial) || delete!(keep, m)
    end
    closed = copy(base)
    for x in keep
        closed .&= sigma[x]
    end
    any(closed) && return nothing
    # each side, stated from its support: the packages the support names, held
    # where the support holds them, imply the bound
    sides = Tuple{Clause{P},Vector{P}}[]
    menubits = UInt64(0)
    for i in menu
        haskey(bit, i) && (menubits |= UInt64(1) << bit[i])
    end
    forcing = false
    cond = cond_menu = fracture = routelen = 0
    for m in sort!(collect(keep))
        ant = Dict{P,BitVector}()
        for i in reason
            (m >> bit[i]) & 1 == 1 || continue
            p = facts[i].pkg
            b = fcl[i][p].bits
            haskey(ant, p) ? (ant[p] .&= b) : (ant[p] = copy(b))
        end
        pairs = Pair{P,Lit}[p => Lit(.~b) for (p, b) in ant]
        push!(pairs, pivot => Lit(copy(sigma[m])))
        cl = clause(pairs)
        cl === nothing && continue
        route = P[q for q in routes[m] if q ∉ packages(cl)]
        push!(sides, (cl, sort!(unique!(route))))
        sigma[m][end] || (forcing = true)
        wide = count_ones(m) > 1
        wide && (cond += 1)
        wide && !iszero(m & menubits) && (cond_menu += 1)
        lit = cl[pivot]
        # a support the registry cannot satisfy at all leaves the pivot out of
        # the statement: what the side says is then about its own packages
        lit === nothing ||
            (fracture += length(Clauses.selected_runs(selected(lit))))
        routelen += length(route)
    end
    key = (!forcing, cond_menu, cond, length(sides), routelen, fracture,
           string(pivot))
    return Meet{P}(pivot, sides, key)
end

# Every package the argument touches is a lawful pivot, so the choice is made
# for the reader by ranking what each would print: a meet with a forcing side
# first (a meet all of whose registry sides admit absence has elided the
# registry's half of the argument), then separated sides over conditional ones,
# then fewer sides, shorter routes and less version fracture.
function best_meet(
    sat   :: SAT{P,V},
    reason:: Vector{Int},
    facts :: Vector{Fact{P}},
    fcl   :: Dict{Int,Clause{P}},
    core  :: Vector{Clause{P}},
    menu  :: Vector{Int},
) where {P,V}
    cands = P[]
    for c in core, p in packages(c)
        p in cands || push!(cands, p)
    end
    for i in reason
        facts[i].pkg in cands || push!(cands, facts[i].pkg)
    end
    sort!(cands)
    best = nothing
    for p in cands
        m = meet_at(sat, p, reason, facts, fcl, core, menu)
        m === nothing && continue
        (best === nothing || m.key < best.key) && (best = m)
    end
    return best
end

# The registry's share of one reason: a minimal set of its statements that the
# reason cannot live with. Empty is a value here — a requirement whose package
# the query leaves nothing of contradicts it with no help from the registry.
function reason_core(satx::SAT{P,V}, reason::Vector{Int},
                     facts::Vector{Fact{P}}, selectors::Vector{Int}) where {P,V}
    isempty(selectors) && return Clause{P}[]
    fixed = Int[facts[i].litx for i in reason]
    out = Clause{P}[]
    for v in sat_mus(satx, fixed, selectors)
        c = clause_of(satx, satx.why[v])
        c === nothing || push!(out, c)
    end
    return out
end

# no line may say less than another beside it in the same proof: a covered line
# is one the reader has already been told
function drop_covered(lines::Vector{Line{P}}) where {P}
    keep = trues(length(lines))
    for i in eachindex(lines)
        keep[i] || continue
        for j in eachindex(lines)
            (i == j || !keep[j]) && continue
            subsumes(lines[i].clause, lines[j].clause) || continue
            lines[i].clause == lines[j].clause && i > j && continue
            keep[j] = false
        end
    end
    return Line{P}[lines[i] for i in eachindex(lines) if keep[i]]
end

## putting it together

# what the analysis settles, on the instance, before anything is resolved
struct Plan{P}
    menus     :: Vector{Vector{Int}}
    reqs      :: Vector{Vector{P}}
    lines     :: Vector{Vector{Line{P}}}
    others    :: Symbol
    truncated :: Bool
end

function analyse(
    sat   :: SAT{P,V},
    satx  :: SAT{P,V},
    prob  :: Problem{P},
    facts :: Vector{Fact{P}},
) where {P,V}
    lits = Int[f.lit for f in facts]
    k, fmin = min_repairs(sat, lits)
    used = sort!(unique!(reduce(vcat, fmin; init = Int[])))
    menus = Vector{Int}[]
    covered = length(fmin)
    if !isempty(used)
        factors = product_menus(fmin, used)
        if factors === nothing
            menus, covered = rectangle_menus(fmin, used)
        else
            menus = factors
        end
    end
    others = (covered < length(fmin) || length(fmin) ≥ REPAIR_CAP) ? :some :
        (larger_repairs(sat, lits, fmin) ? :larger : :none)

    fcl = Dict{Int,Clause{P}}()
    for i in eachindex(facts)
        c = fact_clause(sat, prob, facts[i])
        c === nothing || (fcl[i] = c)
    end
    selectors = sort!(collect(keys(satx.why)))
    reqs = Vector{P}[]
    lines = Vector{Line{P}}[]
    truncated = false
    # one walk, over the whole fact set: a reason can hold this conflict's menu
    # *and* facts another menu offers, and a pool with those held out cannot
    # contain it. Ownership is then a filter, and a reason two conflicts own is
    # set out under both — redundant, never wrong.
    found = Vector{Vector{Int}}()
    if !isempty(menus)
        found, ok = reason_walk(sat, collect(eachindex(facts)), lits)
        ok || (truncated = true)
    end
    for (i, menu) in enumerate(menus)
        reasons = Vector{Int}[r for r in found if menu ⊆ r]
        sort!(reasons; by = r -> (length(r), r))
        if isempty(reasons)
            # every conflict owns a reason of its very own, whatever the walk
            # got to: with the other menus settled one way, what is left is
            # still unsatisfiable and every reason in it is this one's
            held = Set{Int}(first(menus[j]) for j in eachindex(menus) if j != i)
            pool = Int[lits[x] for x in eachindex(facts) if x ∉ held]
            at = Dict{Int,Int}(l => j for (j, l) in enumerate(lits))
            r = sort!(Int[at[l] for l in sat_mus(sat, pool)])
            isempty(r) || push!(reasons, r)
            isempty(reasons) && push!(reasons, sort!(unique!(copy(menu))))
            truncated = true
        end
        if length(reasons) > CONFLICT_REASONS
            resize!(reasons, CONFLICT_REASONS)
            truncated = true
        end
        given = sort!(unique!(reduce(vcat, reasons; init = Int[])))
        ls = Line{P}[Line{P}(fcl[j], P[], true) for j in given if haskey(fcl, j)]
        for (n, r) in enumerate(reasons)
            core = reason_core(satx, r, facts, selectors)
            meet = length(r) ≤ MASK_WIDTH ?
                best_meet(sat, r, facts, fcl, core, menu) : nothing
            derived = meet === nothing ?
                Line{P}[Line{P}(it.clause, P[], false, n, nothing)
                        for it in factor_items(
                            Item{P}[Item{P}(c, P[], UInt64(0)) for c in core])] :
                Line{P}[Line{P}(cl, trim_route(sat, core, cl, route),
                                false, n, meet.pivot)
                        for (cl, route) in meet.sides]
            append!(ls, drop_covered(derived))
        end
        push!(lines, ls)
        push!(reqs, sort!(unique!(P[facts[j].pkg for j in given if facts[j].req])))
    end
    return Plan{P}(menus, reqs, lines, others, truncated)
end

# The kinds to lift so that a package this query emptied is choosable again.
# Its classes are what a lift has to give back — a version is not something the
# universe can be asked about — so the cheapest lift is the smallest set of
# kinds restoring every class the full lift would, and among those the one
# whose best restored version is best.
function lift_actions(prob::Problem{P}, sat::SAT{P,V}, univ::Universe{P,V},
                      p::P) where {P,V}
    info = sat.info[p]
    reps = univ.reps[p]
    vs = info.versions
    excl = Vector{Symbol}[exclusion_kinds(prob, p, v) for v in vs]
    kinds = Symbol[]
    for ks in excl, k in ks
        k in kinds || push!(kinds, k)
    end
    sort!(kinds)
    dead = Int[c for c in eachindex(reps) if iszero(reps[c])]
    admits(S, j) = isempty(setdiff(excl[j], S))
    back(S) = Set{Int}(c for c in dead if any(j -> admits(S, j), info.members[c]))
    full = back(kinds)
    best = Symbol[]
    for size = 1:length(kinds)
        bestkey = nothing
        for S in subsets(kinds, size)
            back(S) == full || continue
            j = findfirst(j -> admits(S, j) && !isempty(excl[j]), eachindex(vs))
            key = (something(j, length(vs) + 1), S)
            bestkey === nothing || key < bestkey || continue
            best, bestkey = S, key
        end
        isempty(best) || break
    end
    isempty(best) && (best = kinds)
    return Action{P}[Action(k, p) for k in best]
end

# the size-`k` subsets of `v`, in order
function subsets(v::Vector{T}, k::Int) where {T}
    out = Vector{T}[]
    n = length(v)
    k > n && return out
    idx = collect(1:k)
    while true
        push!(out, T[v[i] for i in idx])
        i = k
        while i ≥ 1 && idx[i] == n - k + i
            i -= 1
        end
        i == 0 && break
        idx[i] += 1
        for j = i+1:k
            idx[j] = idx[j-1] + 1
        end
    end
    return out
end

fix_actions(prob::Problem{P}, sat::SAT{P,V}, univ::Universe{P,V},
            f::Fact{P}) where {P,V} =
    f.req ? Action{P}[Action(:drop, f.pkg)] : lift_actions(prob, sat, univ, f.pkg)

# What carrying `actions` out gets you: the resolver's own answer for the
# withdrawn query, on the universe the failed resolve was run against. A model
# found during diagnosis witnesses satisfiability; which versions the user
# would *get* is the resolver's optimising answer, and nothing else may stand
# in for it.
function witness(sat::SAT{P,V}, univ::Universe{P,V}, prob::Problem{P},
                 actions::Vector{Action{P}}; by, order) where {P,V}
    drop_reqs = P[a.pkg for a in actions if a.kind === :drop]
    drop_constraints = Dict{Symbol,Set{P}}()
    for a in actions
        a.kind === :drop && continue
        push!(get!(Set{P}, drop_constraints, a.kind), a.pkg)
    end
    sol = resolve(sat, relax(univ, prob, drop_reqs, drop_constraints; order); by)
    return sol === nothing ? Dict{P,V}() : sol
end

"""
    diagnose(sat, prob, univ; by, order) :: Diagnosis

Why `prob` cannot be satisfied on `univ`, and what its user could change. The
instance is left exactly as it was found, so a caller may go on using it.

Requirements the universe holds nothing of are their own conflicts: nothing
this query did took them away, so dropping them is all that could help. What is
left is settled on the facts — the requirements and the query's own version
limits — and one relaxation is resolved per fix on the menus offered.
"""
function diagnose(
    sat  :: SAT{P,V},
    prob :: Problem{P},
    univ :: Universe{P,V};
    by   :: Function = identity,
    order = nothing,
) where {P,V}
    vm = VarMap(sat)
    seen = P[]
    for p in prob.reqs
        p in seen || push!(seen, p)
    end
    sort!(seen)
    gone = P[p for p in seen if !haskey(vm, p)]
    live = P[p for p in seen if haskey(vm, p)]

    facts = Fact{P}[]
    satx = SAT(univ; explain = true)
    plan = try
        vmx = VarMap(satx)
        with_classes_relaxed(sat) do
            with_emptied_packages(sat, vm) do lits, _
                with_classes_relaxed(satx) do
                    with_emptied_packages(satx, vmx) do litsx, _
                        for i in eachindex(vm.pkgs)
                            push!(facts, Fact{P}(false, vm.pkgs[i], lits[i], litsx[i]))
                        end
                        for p in live
                            push!(facts, Fact{P}(true, p, installed_lit(sat, p),
                                                 installed_lit(satx, p)))
                        end
                        analyse(sat, satx, prob, facts)
                    end
                end
            end
        end
    finally
        Resolver.finalize(satx)
    end

    # the menus, as the actions they ask for: the missing requirements first,
    # each its own forced choice, then the factors of the cheapest repairs
    menus = Vector{Vector{Action{P}}}[]
    for p in gone
        push!(menus, Vector{Action{P}}[Action{P}[Action(:drop, p)]])
    end
    for menu in plan.menus
        push!(menus, Vector{Action{P}}[fix_actions(prob, sat, univ, facts[i])
                                       for i in menu])
    end

    conflicts = Conflict{P,V}[]
    for (i, entries) in enumerate(menus)
        fixes = Fix{P,V}[]
        for e in entries
            acts = copy(e)
            for j in eachindex(menus)
                j == i || append!(acts, first(menus[j]))
            end
            unique!(acts)
            push!(fixes, Fix{P,V}(e, witness(sat, univ, prob, acts; by, order)))
        end
        if i ≤ length(gone)
            push!(conflicts, Conflict{P,V}(P[gone[i]], Line{P}[],
                Dict{P,Vector{V}}(), Dict{P,Vector{Vector{Symbol}}}(), fixes))
            continue
        end
        n = i - length(gone)
        lines = plan.lines[n]
        pkgs = P[]
        for l in lines, p in packages(l.clause)
            p in pkgs || push!(pkgs, p)
        end
        versions = Dict{P,Vector{V}}(p => clause_versions(sat, p) for p in pkgs)
        excluded = Dict{P,Vector{Vector{Symbol}}}()
        for p in pkgs
            ks = Vector{Symbol}[exclusion_kinds(prob, p, v) for v in versions[p]]
            any(!isempty, ks) && (excluded[p] = ks)
        end
        push!(conflicts, Conflict{P,V}(plan.reqs[n], lines, versions, excluded, fixes))
    end
    return Diagnosis{P,V}(conflicts, plan.others, plan.truncated)
end

## the report
#
# What a diagnosis says, and the rules that keep every sentence of it true.
#
# The page is flat. The query's own facts come first, said as the user's ("you
# require A"; "your compat leaves A 1.2"), each package once; then the
# registry's statements, each said as the implication from the facts it rests on
# to the bound it puts on the package the argument meets at, with the packages
# an elimination reached it through in parentheses. Nothing on the page is
# derived from anything else on the page, which is why a line can carry a route
# and why the meet — the sides whose intersection is empty — prints last and
# together.
#
# Only a query line may say "your". Everything else is the registry's, and a
# registry statement never attributes a bound to a `Project.toml` it cannot see.

# print `text` filled to `width` columns, `lead` before the first line and
# `rest` before every later one — a report is read in a terminal, and a line
# that runs off the edge of one is a line the reader scrolls sideways for
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

"""
    action_phrase(a::Action) :: String

One action, said as something the reader could carry out. Whatever a constraint
kind is called inside the resolver, what it reads as here is an edit.
"""
function action_phrase(a::Action)
    a.kind === :drop && return "drop requirement $(a.pkg)"
    a.kind === :compat && return "relax your compat on $(a.pkg)"
    a.kind === :pin && return "unpin $(a.pkg)"
    return "allow $(a.kind) versions of $(a.pkg)"
end

join_and(xs) = join(xs, ", ", " and ")

fix_phrase(f::Fix) = join_and(String[action_phrase(a) for a in f.actions])

# What a conflict is about, in one sentence. A requirement whose package the
# universe holds nothing of has no argument to make, so what it says is the
# whole of what happened to it.
function conflict_heading(c::Conflict)
    rs = c.reqs
    isempty(c.lines) && length(rs) == 1 &&
        return "no version of $(only(rs)) is available."
    isempty(rs) && return "the requirements cannot all be satisfied."
    length(rs) == 1 && return "$(only(rs)) cannot be satisfied."
    length(rs) == 2 && return "$(rs[1]) and $(rs[2]) cannot both be satisfied."
    return join_and(String[string(r) for r in rs]) * " cannot all be satisfied."
end

# Which of the query's kinds took versions of `p` away, and what they left.
# Read straight off the query, so no line here needs a solver's licence; and
# named as the user's, which nothing else on the page may be.
function constraint_phrase(c::Conflict{P,V}, p::P, l::Line{P},
                           named::Bool) where {P,V}
    kinds = Symbol[]
    for ks in c.excluded[p], k in ks
        k in kinds || push!(kinds, k)
    end
    sort!(kinds)
    lead = join(String["your $k" for k in kinds], " and ")
    m = l.clause[p]
    sel = selected(m)
    any(sel) || return "$lead leaves no version of $p"
    r = range_phrase(c.versions[p], sel)
    isempty(r) && return "$lead leaves $p"
    return named ? "$lead leaves $p $r" : "$lead leaves $r"
end

# is this given line the requirement itself, rather than a limit on it?
is_requirement(l::Line{P}, p::P) where {P} =
    (m = l.clause[p]; m !== nothing && !absent(m) &&
     all(m[i] for i = 1:nversions(m)))

function line_phrase(l::Line{P}, vers, names) where {P}
    s = clause_phrase(l.clause, vers, names)
    isempty(l.through) && return s
    return s * " (through " * join_and(String[names(q) for q in l.through]) * ")"
end

# do these lines leave the package they meet at nothing to be? Only then may
# the page say so: two lines that agree about a package leave it something, and
# claiming otherwise would claim more than the page shows.
function meet_is_empty(group::Vector{Line{P}}, pivot) where {P}
    acc = nothing
    for l in group
        m = l.clause[pivot]
        m === nothing && return false
        acc = acc === nothing ? copy(m.bits) : (acc .& m.bits)
    end
    return acc !== nothing && !any(acc)
end

function print_given(io::IO, c::Conflict{P,V}, given::Vector{Line{P}},
                     vers, names) where {P,V}
    order = P[]
    single = Dict{P,Vector{Line{P}}}()
    other = Line{P}[]
    for l in given
        ps = packages(l.clause)
        if length(ps) == 1
            push!(get!(Vector{Line{P}}, single, ps[1]), l)
        else
            push!(other, l)
        end
    end
    for p in c.reqs
        haskey(single, p) && p ∉ order && push!(order, p)
    end
    for l in given
        ps = packages(l.clause)
        length(ps) == 1 && ps[1] ∉ order && push!(order, ps[1])
    end
    for p in order
        req = nothing
        con = nothing
        extra = Line{P}[]
        for l in single[p]
            if is_requirement(l, p) && p in c.reqs
                req = l
            elseif absent(l.clause[p]) && haskey(c.excluded, p)
                con = l
            else
                push!(extra, l)
            end
        end
        if req !== nothing
            s = "you require $p"
            con === nothing || (s *= ", and " * constraint_phrase(c, p, con, false))
            print_wrapped(io, s, "  • ", "    ")
        elseif con !== nothing
            print_wrapped(io, constraint_phrase(c, p, con, true), "  • ", "    ")
        end
        for l in extra
            print_wrapped(io, line_phrase(l, vers, names), "  • ", "    ")
        end
    end
    for l in other
        print_wrapped(io, line_phrase(l, vers, names), "  • ", "    ")
    end
end

function print_derived(io::IO, derived::Vector{Line{P}}, vers, names) where {P}
    i = 1
    while i ≤ length(derived)
        l = derived[i]
        j = i
        if l.pivot !== nothing
            while j < length(derived) && derived[j+1].pivot == l.pivot &&
                  derived[j+1].proof == l.proof
                j += 1
            end
        end
        group = derived[i:j]
        if length(group) ≥ 2 && meet_is_empty(group, l.pivot)
            println(io, "  • no version of ", names(l.pivot), " is all of these:")
            for g in group
                print_wrapped(io, line_phrase(g, vers, names), "      — ", "        ")
            end
        else
            for g in group
                print_wrapped(io, line_phrase(g, vers, names), "  • ", "    ")
            end
        end
        i = j + 1
    end
end

# What the fix gets you, of the packages this conflict is about: the reader sees
# the witness land where the opened meet says it can, which is what the versions
# are on the page for.
function print_allows(io::IO, c::Conflict{P,V}, f::Fix{P,V},
                      indent::String) where {P,V}
    ps = sort!(P[p for p in keys(c.versions) if haskey(f.solution, p)])
    isempty(ps) && return
    print_wrapped(io, join(String["$p $(f.solution[p])" for p in ps], ", "),
                  indent * "→ allows: ", indent * "  ")
end

# A menu of one has exactly three honest wordings, and which one is the whole of
# what the reader learns about the gap. Never derived from the length of a
# vector; derived from the two decided questions — whether anything larger
# exists, and whether the menus reach every repair as cheap as theirs.
function print_menu(io::IO, c::Conflict{P,V}, others::Symbol) where {P,V}
    isempty(c.fixes) && return
    if length(c.fixes) == 1
        word = others === :none ? "The only fix" :
               others === :larger ? "The only minimal fix" : "One fix"
        println(io, "  ", word, ": ", fix_phrase(c.fixes[1]))
        print_allows(io, c, c.fixes[1], "    ")
    else
        println(io, "  Fix it by any one of:")
        for (i, f) in enumerate(c.fixes)
            println(io, "    ", i, ". ", fix_phrase(f))
            print_allows(io, c, f, "       ")
        end
    end
end

"""
    print_conflict(io, c, index = nothing; others = :some)

One conflict's page: its heading (where it is numbered), the lines that prove
it, and the menu that settles it. `others` is what the whole diagnosis knows
about the repairs its menus do not reach, which is what a menu of one is
entitled to say about itself.
"""
function print_conflict(io::IO, c::Conflict{P,V}, index = nothing;
                        others::Symbol = :some) where {P,V}
    index === nothing ||
        println(io, "Conflict ", index, ": ", conflict_heading(c))
    vers(p) = c.versions[p]
    names(p) = string(p)
    print_given(io, c, Line{P}[l for l in c.lines if l.given], vers, names)
    print_derived(io, Line{P}[l for l in c.lines if !l.given], vers, names)
    print_menu(io, c, others)
end

function Base.show(io::IO, d::Diagnosis)
    n = length(d.conflicts)
    f = prod(length(c.fixes) for c in d.conflicts; init = 1)
    print(io, "Diagnosis: ", n, n == 1 ? " conflict, " : " conflicts, ",
          f, f == 1 ? " fix" : " fixes")
end

function Base.show(io::IO, ::MIME"text/plain", d::Diagnosis)
    n = length(d.conflicts)
    print(io, "Unsatisfiable — ", n, n == 1 ? " conflict" : " conflicts")
    n > 1 && print(io, ", each of which must be fixed")
    println(io, ":")
    for (i, c) in enumerate(d.conflicts)
        println(io)
        print_conflict(io, c, i; others = d.others)
    end
    if d.others === :larger
        println(io)
        println(io, "Larger solutions also exist.")
    elseif d.others === :some
        println(io)
        println(io, "Other solutions also exist.")
    end
    if d.truncated
        println(io)
        println(io, "There may be more to say about some of these.")
    end
end

## verification

# does `sol` satisfy `c`, reading a package it does not name as absent?
function models(c::Clause{P}, sol::Dict{P,V},
                versions::Dict{P,Vector{V}}) where {P,V}
    for (p, m) in c.lits
        v = get(sol, p, nothing)
        if v === nothing
            absent(m) && return true
        else
            i = findfirst(==(v), get(versions, p, V[]))
            i === nothing || m[i] && return true
        end
    end
    return false
end

# the packages a fix's withdrawal takes the query's word about
touched(f::Fix{P,V}) where {P,V} = Set{P}(a.pkg for a in f.actions)

"""
    report_problems(d) :: Vector{String}

Everything Section 8's checker can decide without asking the solver:

  * **(V2) visible closure** — each meet's sides really do intersect emptily,
    so the contradiction is on the page rather than behind it;
  * **(V3) source coverage** — every requirement the page answers for, and
    every package its menu asks the reader to act on, is named by a line;
  * **(V5) witness coherence** — each fix's witness lands inside every line the
    fix's own withdrawal leaves standing. Silent breakage here is invisible to
    every other check, which is exactly why this one exists.

Empty when the report is sound. The remaining obligation — that each printed
line is true of the universe this query left (V1) — is one entailment query per
line and belongs to whoever holds the instance; the menu wordings (V6) are read
off the two decided questions when the page is printed, so there is no second
place for them to disagree.

Every check is per explanation, never against a union: where two explanations'
lines are `S₁ ∪ S₂` and `S₂` alone contradicts, the union stays contradictory
whatever is deleted from `S₁`, so a union-level check would pass a page that
silently destroyed the whole account of `S₁`.
"""
function report_problems(d::Diagnosis{P,V}) where {P,V}
    bad = String[]
    # V5 is stated against the *full* withdrawal, never the single entry: an
    # owned reason can hold other conflicts' facts, and the witness respects
    # only the sides whose supports survive everything withdrawn. So each fix
    # is judged with every other conflict settled the way its own witness was
    # taken — the first entry of its menu.
    firsts = Vector{Action{P}}[isempty(c.fixes) ? Action{P}[] :
                               first(c.fixes).actions for c in d.conflicts]
    for (n, c) in enumerate(d.conflicts)
        rest = Set{P}(a.pkg for j in eachindex(firsts) if j != n
                             for a in firsts[j])
        for s in conflict_problems(c, rest)
            push!(bad, "conflict $n: $s")
        end
    end
    return bad
end

function conflict_problems(c::Conflict{P,V}, rest::Set{P} = Set{P}()) where {P,V}
    bad = String[]
    given = Line{P}[l for l in c.lines if l.given]
    derived = Line{P}[l for l in c.lines if !l.given]

    # (V2) every meet closes, on its own lines and the query's, never on
    # another proof's
    for n in unique!(Int[l.proof for l in derived])
        mine = Line{P}[l for l in derived if l.proof == n]
        for pivot in unique(P[l.pivot for l in mine if l.pivot !== nothing])
            haskey(c.versions, pivot) || continue
            acc = trues(length(c.versions[pivot]) + 1)
            closed = false
            for l in Line{P}[mine; given]
                m = l.clause[pivot]
                if m === nothing
                    # a side stated from a support the registry cannot satisfy
                    # names no bound at the pivot: it closes on its own
                    l.given || l.pivot != pivot || (closed = true)
                else
                    acc .&= m.bits
                end
            end
            (closed || !any(acc)) && continue
            push!(bad,
                "proof $n leaves $pivot something every one of its lines admits")
        end
    end

    # (V3) the page names what it answers for and what it asks to be changed
    named = Set{P}(p for l in c.lines for p in packages(l.clause))
    if !isempty(c.lines)
        for r in c.reqs
            r in named || push!(bad, "names $r and no line mentions it")
        end
        for f in c.fixes, a in f.actions
            a.pkg in named || push!(bad, "offers $(a.pkg) and no line mentions it")
        end
    end

    # (V5) the witness lands where the opened meet says it can
    for f in c.fixes
        isempty(f.solution) && continue
        gone = union(touched(f), rest)
        for l in c.lines
            any(p -> p in gone, packages(l.clause)) && continue
            models(l.clause, f.solution, c.versions) && continue
            push!(bad, "the witness for " * fix_phrase(f) *
                  " does not satisfy: " * line_phrase(l, p -> c.versions[p], string))
        end
    end
    return bad
end

end # module Diagnostics
