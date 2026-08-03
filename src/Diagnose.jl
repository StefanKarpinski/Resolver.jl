# UNSAT diagnostics: explain why a resolution problem is unsatisfiable and list
# verified fixes.
#
# The vocabulary below is the whole public surface of a diagnosis: a `Diagnosis`
# is a list of independent `Conflict`s — each a requirement-level cluster, a
# chain of `Fact`s telling that cluster's story, and any upstream suggestions —
# plus a global list of verified `Fix`es. Everything is structured data;
# rendering is a separate layer at the bottom of this file, so a caller (Pkg,
# say) can render its own.

# facts — the vocabulary of explanations and fixes

abstract type Fact end

"""
    Requirement(pkg)

The user asked for `pkg`. Relaxing it means dropping the requirement.
"""
struct Requirement{P} <: Fact
    pkg :: P
end

"""
    Uninstallable(pkg)

`pkg` has no installable version at all — every version of it depends,
transitively, on a package with no versions. Nothing the user can relax fixes
this; the only correction is to stop requiring `pkg`.
"""
struct Uninstallable{P} <: Fact
    pkg :: P
end

"""
    UserCompat(pkg, allowed, [label])

The user's compat bound on `pkg`, as the list of that package's versions it
admits. Relaxing it means loosening (or deleting) the bound.

`label` says what kind of restriction it was, for reports: `:compat` (the
default) reads as "your compat restricts …", `:requested` as "you requested …".
A caller sets it through [`Problem`](@ref)'s `labels`, because "your compat" is
wrong when the bound came from `Pkg.add(name, version)` rather than a compat
section. Unknown labels get a neutral phrasing rather than an error.
"""
struct UserCompat{P,V} <: Fact
    pkg     :: P
    allowed :: Vector{V}
    label   :: Symbol
end

UserCompat{P,V}(pkg, allowed) where {P,V} =
    UserCompat{P,V}(pkg, allowed, :compat)
UserCompat(pkg::P, allowed::Vector{V}) where {P,V} =
    UserCompat{P,V}(pkg, allowed, :compat)

"""
    Pin(pkg, version)

`pkg` is pinned at `version`. Relaxing it means unpinning.
"""
struct Pin{P,V} <: Fact
    pkg     :: P
    version :: V
end

"""
    Admission(kind, pkg, forbidden)

An *admission knob* ruling versions of `pkg` out: `kind` names the source
(`:prerelease`, `:yanked`, …) and `forbidden` lists the versions it excludes.
Unlike a compat bound or a pin, a knob is stated about versions rather than
about a package, so the same one appears once per package it touches.

Relaxing it means admitting that class of version — passing `--allow-pre`, say.
"""
struct Admission{P,V} <: Fact
    kind      :: Symbol
    pkg       :: P
    forbidden :: Vector{V}
end

"""
    Bound(pkg, versions, dep, allowed)

A third-party compatibility bound: the listed `versions` of `pkg` are
incompatible with every version of `dep` outside `allowed`. Only an upstream
release can relax one.

The resolver's package info records compatibility symmetrically — a conflict
between `pkg@v` and `dep@w` is one fact, not two directed ones — so a `Bound`
names an incompatible *pair* of version groups rather than saying which side
declared the bound. See the manual's diagnostics page.
"""
struct Bound{P,V} <: Fact
    pkg      :: P
    versions :: Vector{V}
    dep      :: P
    allowed  :: Vector{V}
    # true when every version of `pkg` this fact is about is one the query's own
    # constraints already rule out. Such a fact is needed for the conflict to be
    # *minimal* -- something has to close off the old versions -- but it is not
    # what a reader wants at the same level as the lines that name versions they
    # could actually get. Reports demote them; see `pair_facts`.
    incidental :: Bool
end

Bound{P,V}(pkg, versions, dep, allowed) where {P,V} =
    Bound{P,V}(pkg, versions, dep, allowed, false)
Bound(pkg::P, versions::Vector{V}, dep::P, allowed::Vector{V}) where {P,V} =
    Bound{P,V}(pkg, versions, dep, allowed, false)

fact_pkg(f::Fact) = f.pkg

# value equality for facts (immutable, compared field-wise)
Base.:(==)(a::F, b::F) where {F<:Fact} =
    all(getfield(a, i) == getfield(b, i) for i = 1:nfields(a))
function Base.hash(a::Fact, h::UInt)
    for i = 1:nfields(a)
        h = hash(getfield(a, i), h)
    end
    hash(typeof(a), h)
end

# deterministic ordering: actionability first, then by package. requirements
# lead (dropping one is always available), uninstallability follows its
# requirement, then the user's own knobs, then the bounds only upstream can move
kind_rank(::Requirement)   = 1
kind_rank(::Uninstallable) = 2
kind_rank(::Pin)           = 3
kind_rank(::UserCompat)    = 4
kind_rank(::Admission)     = 5
kind_rank(::Bound)         = 6

order_facts(fs::Vector{Fact}) =
    sort(fs; by = f -> (kind_rank(f), fact_pkg(f)))

"""
    Fix(actions, solution)

A verified minimal correction: relaxing every fact in `actions` — dropping a
`Requirement`, loosening a `UserCompat`, lifting a `Pin` — makes the problem
solvable, and `solution` is the resolution it then produces.
"""
struct Fix{P,V}
    actions  :: Vector{Fact}
    solution :: Dict{P,V}
end

"""
    UpstreamFix(bound, solution)

A conflict that a new upstream release could fix on its own: were `bound`
relaxed, the conflict would go away and `solution` is what would resolve.
"""
struct UpstreamFix{P,V}
    bound    :: Bound{P,V}
    solution :: Dict{P,V}
end

"""
    Conflict(reqs, chain, upstream)

One independent reason the problem is unsatisfiable: `reqs` is a minimal set of
requirements that cannot be satisfied together, `chain` tells the story as an
ordered list of facts, and `upstream` lists bounds whose relaxation would
resolve this conflict.
"""
struct Conflict{P,V}
    reqs     :: Vector{P}
    chain    :: Vector{Fact}
    upstream :: Vector{UpstreamFix{P,V}}
end

"""
    Diagnosis(reqs, conflicts, fixes, versions)

Why a resolution problem is unsatisfiable: the `conflicts` are independent, the
`fixes` are global (each one repairs every conflict at once), and `versions`
records the version list of every package the report mentions, so a renderer
can compress version sets into ranges.

An empty `fixes` means there is nothing to propose: every correction would have
to drop a requirement that dropping does not help, which happens when only an
upstream release could resolve the conflict.
"""
struct Diagnosis{P,V}
    reqs      :: Vector{P}
    conflicts :: Vector{Conflict{P,V}}
    fixes     :: Vector{Fix{P,V}}
    versions  :: Dict{P,Vector{V}}
end

# a fix dominates another if it relaxes a proper subset of what the other
# relaxes; dominated fixes carry no information and are suppressed
function filter_dominated(fixes::Vector{Fix{P,V}}) where {P,V}
    sets = [Set{Fact}(fix.actions) for fix in fixes]
    keep = [!any(sets[j] ⊊ sets[i] for j in eachindex(sets))
            for i in eachindex(sets)]
    return fixes[keep]
end

# order a cluster's facts into story order: requirements first (sorted), then a
# BFS from the required packages along package links — introducing a package
# emits its uninstallability, pin and compat facts, then its bounds (each bound
# introduces the package on its other side). Deterministic: the queue is
# processed in sorted order.
function order_chain(facts::Vector{Fact}, ::Type{P}, ::Type{V}) where {P,V}
    reqfacts = sort!(Fact[f for f in facts if f isa Requirement]; by = fact_pkg)
    Uby = Dict{P,Uninstallable{P}}()
    Cby = Dict{P,UserCompat{P,V}}()
    Hby = Dict{P,Pin{P,V}}()
    Aby = Dict{P,Vector{Admission{P,V}}}()
    Bby = Dict{P,Vector{Bound{P,V}}}()
    for f in facts
        f isa Uninstallable && (Uby[f.pkg] = f)
        f isa UserCompat    && (Cby[f.pkg] = f)
        f isa Pin           && (Hby[f.pkg] = f)
        f isa Admission     && push!(get!(() -> Admission{P,V}[], Aby, f.pkg), f)
        f isa Bound         && push!(get!(() -> Bound{P,V}[], Bby, f.pkg), f)
    end
    foreach(v -> sort!(v; by = a -> a.kind), values(Aby))
    foreach(v -> sort!(v; by = b -> (b.dep, b.versions)), values(Bby))
    chain = copy(reqfacts)
    visited = Set{P}()
    queue = P[fact_pkg(f) for f in reqfacts]
    while !isempty(queue)
        sort!(queue)
        p = popfirst!(queue)
        p in visited && continue
        push!(visited, p)
        haskey(Uby, p) && push!(chain, Uby[p])
        haskey(Hby, p) && push!(chain, Hby[p])
        haskey(Cby, p) && push!(chain, Cby[p])
        append!(chain, get(Aby, p, Admission{P,V}[]))
        for b in get(Bby, p, Bound{P,V}[])
            push!(chain, b)
            push!(queue, b.dep)
        end
    end
    # emit any facts the BFS didn't reach (deterministically)
    for f in order_facts(Fact[f for f in facts if f ∉ chain])
        push!(chain, f)
    end
    return chain
end

# diagnosis on the failed production instance
#
# The instance `resolve` just failed on is convertible in place: one `sat_pop`
# frees every user-constraint selector, and from then on a diagnostic query is
# an assumption set -- some requirement package variables, some selectors --
# and nothing else. That matters: the theory licensing all of this (see the
# manual's diagnostics page) covers *relaxations* only, so no query may add a
# clause. Fix enumeration force-keeps assumptions instead of blocking.
#
# It also fixes the granularity. Theorem D1 holds for column-closed groups, and
# a package's user constraints are column-closed only taken *together* -- its
# compat bound and its pin share one exclusion column, and relaxing half of one
# is exactly the counterexample of Proposition D1'. So the relaxable unit here
# is a package: a fix touching a package with both sources relaxes both.

# one relaxable group: the assumption literals that switch it on
struct Relaxable{P}
    kind :: Symbol       # :req or :user
    pkg  :: P
    lits :: Vector{Int}
end

"""
    Relaxation

The knobs a diagnosis may turn on the failed instance: `rel` is every group,
`reqs` the requirement groups a fix may drop, and `user` the constraint groups
a fix may relax.

`user` is one group per constrained package, holding all of its selectors --
its whole exclusion column. Switching a whole column off is a column-closed
relaxation whatever the sources in it, which is the granularity Theorem D1
licenses and the only one anything here asks for.

No group is assumed implicitly: a group left out of an assumption set is
*free*, and the solver will happily switch it off, so each caller names exactly
the context its query means to hold.
"""
struct Relaxation{P}
    rel  :: Vector{Relaxable{P}}
    reqs :: Vector{Int}
    user :: Vector{Int}
end

# switch exactly the given groups on and test satisfiability; everything else is
# free, and the solver may set it false
function relax_sat(sat::SAT, R::Relaxation, on)
    for i in on, l in R.rel[i].lits
        sat_assume_var(sat, l)
    end
    is_satisfiable(sat)
end

# group MUS with an actionability-biased deletion order.
#
# `on` are the groups switched on; `shrink` is the (ordered) subset we may
# delete, least actionable first, so that when the same conflict can be told
# two ways the surviving MUS keeps the explanation the user can act on. Groups
# in `on ∖ shrink` are context and stay on.
#
# A single ordered pass suffices -- no restart after deletions. Satisfiability
# is monotone under removing assumptions: if deleting `g` once made the
# instance satisfiable (so `g` was kept), any later deletion of `g` would test
# a subset of that satisfiable assumption set, which is still satisfiable. Kept
# groups provably stay needed, and the result is minimal for the same reason.
# We start from the full assumption set rather than the failed-assumption core,
# which is nondeterministic and can drop the actionable explanation before the
# biased order gets to prefer it.
function group_mus(
    sat    :: SAT,
    R      :: Relaxation,
    on     :: AbstractVector{Int},
    shrink :: AbstractVector{Int},
)
    relax_sat(sat, R, on) && return Int[]
    mus = collect(on)
    for g in shrink
        k = findfirst(==(g), mus)
        k === nothing && continue
        deleteat!(mus, k)
        relax_sat(sat, R, mus) && insert!(mus, k, g) # needed: keep it
    end
    return sort!(mus)
end

# partition the requirements into independent conflicts (requirement-level
# MICE): repeatedly extract a requirement-level MUS with every user group held
# on as context, remove its members, repeat until what is left is satisfiable.
# Each disjoint MUS is one cluster, and one story.
function cluster_reqs(sat::SAT, R::Relaxation, reqs::Vector{Int})
    remaining = copy(reqs)
    reqset = Set(reqs)
    clusters = Vector{Int}[]
    while !isempty(remaining)
        mus = group_mus(sat, R, Int[remaining; R.user], remaining)
        rmus = filter(in(reqset), mus)
        isempty(rmus) && break
        push!(clusters, rmus)
        setdiff!(remaining, rmus)
    end
    return clusters
end

# the resolution a fix produces: re-resolve this same instance with the kept
# groups asserted and the dropped ones denied, inside a frame that unwinds
# afterwards. Theorem D1 preserves the layered *answer*, not just the verdict,
# so these are the versions the raw problem would give under the fix.
function fix_solution(
    sat  :: SAT{P,V},
    R    :: Relaxation{P},
    kept :: AbstractVector{Int},
    by   :: Function,
) where {P,V}
    keptset = Set(kept)
    reqs = P[R.rel[i].pkg for i in kept if R.rel[i].kind === :req]
    with_temp_clauses(sat) do
        for (i, r) in enumerate(R.rel)
            r.kind === :user || continue
            on = i in keptset
            for l in r.lits
                sat_add_var(sat, on ? l : -l)
                sat_add(sat)
            end
        end
        sol = resolve_core(sat, sort!(reqs); by, restore = false)
        sol === nothing ? Dict{P,V}() :
            Dict{P,V}(p => sat.info[p].versions[i] for (p, i) in sol)
    end
end

# enumerate verified fixes as minimal correction sets, MARCO-style.
#
# A greedy pass over the relaxables in keep-order (requirements first, then the
# user's own knobs) yields a maximal satisfiable subset; its complement is a
# minimal correction set. To get the *next* fix we force-keep one of the groups
# a previous fix dropped -- an assumption, not a clause -- and run the greedy
# pass again. Working through the force-keep sets breadth-first enumerates
# undominated alternatives without ever adding a constraint to the instance,
# which is what the licensing theorem requires.
function enumerate_fixes(
    sat   :: SAT{P,V},
    R     :: Relaxation{P},
    order :: Vector{Int}, # keep-order over the relaxable group ids
    extra :: Vector{Fact}; # actions every fix must take (uninstallable reqs)
    max_fixes :: Integer = 8,
    by :: Function = identity,
    acceptable :: Function = (actions, sol) -> true,
) where {P,V}
    seed = Int[]
    fixes = Fix{P,V}[]
    found = Set{Vector{Int}}()
    queued = Set{Vector{Int}}([seed])
    queue = Vector{Int}[seed]
    # most force-keep sets yield a fix, but an infeasible one yields none and
    # costs a solve, so bound the passes rather than the fixes alone
    passes = 0
    while !isempty(queue) && length(fixes) < max_fixes && passes < 16max_fixes
        passes += 1
        force = popfirst!(queue)
        # nothing can be fixed while keeping all of `force`
        relax_sat(sat, R, force) || continue
        # greedy maximal satisfiable subset in keep-order
        kept = copy(force)
        for g in order
            g in force && continue
            push!(kept, g)
            relax_sat(sat, R, kept) || pop!(kept)
        end
        mcs = sort!(setdiff(order, kept))
        # an empty correction set is still a fix when there are forced actions
        # (uninstallable requirements); with none it means nothing was broken
        isempty(mcs) && isempty(extra) && continue
        if mcs ∉ found
            push!(found, mcs)
            actions = copy(extra)
            for i in mcs
                append!(actions, group_facts(sat, R.rel[i]))
            end
            actions = order_facts(actions)
            sol = fix_solution(sat, R, kept, by)
            # a correction the caller will not stand behind is not a fix: the
            # BFS still expands past it, so the alternatives get found
            acceptable(actions, sol) && push!(fixes, Fix{P,V}(actions, sol))
        end
        # the next fixes keep something this one dropped
        for g in mcs
            next = sort!(push!(copy(force), g))
            next ∈ queued && continue
            push!(queued, next)
            push!(queue, next)
        end
    end
    return filter_dominated(fixes)
end

# the facts a relaxable group stands for.
#
# A group is a whole package's constraints -- compat bound, pin, and every
# admission knob that touches it -- because they share one exclusion column and
# Theorem D1 only licenses moving a column as a whole. What the facts name is
# the sources that put a cell in that column: a source forbidding no version
# this instance still has contributes nothing to it, so relaxing the group
# neither needs it nor is affected by it, and naming it would ask the user to
# make a change that had no part in the verification. (`SAT` decides which
# sources get a selector by the same test, so facts and selectors agree.)
function group_facts(sat::SAT{P,V}, r::Relaxable{P}) where {P,V}
    r.kind === :req && return Fact[Requirement(r.pkg)]
    facts = Fact[]
    prob = sat.prob
    prob === nothing && return facts
    p = r.pkg
    vers = sat.info[p].versions
    # driven by the same enumeration `SAT` builds its selectors from, so the
    # facts and the things they stand for cannot drift apart
    for (kind, forbids) in exclusion_sources(prob, p)
        forbidden = V[v for v in vers if forbids(v)]
        isempty(forbidden) && continue
        push!(facts,
            kind === :pin    ? Pin{P,V}(p, prob.pins[p]) :
            kind === :compat ? UserCompat{P,V}(p,
                                   V[v for v in vers if v ∈ prob.compat[p]],
                                   get(prob.labels, p, :compat)) :
                               Admission{P,V}(kind, p, forbidden))
    end
    return facts
end

# the groups of a failed instance, in keep-order: requirements first (dropping
# one is the coarsest change a user can make), then the packages whose
# constraints the selectors guard.
function relaxation(sat::SAT{P,V}, reqs::AbstractVector{P}) where {P,V}
    rel = Relaxable{P}[]
    rids = Int[]
    for p in reqs
        push!(rel, Relaxable{P}(:req, p, Int[sat.vars[p]]))
        push!(rids, length(rel))
    end
    bypkg = Dict{P,Vector{Int}}()
    for (sel, (_, p)) in sat.sels
        push!(get!(() -> Int[], bypkg, p), sel)
    end
    uids = Int[]
    for p in sort!(collect(keys(bypkg)))
        push!(rel, Relaxable{P}(:user, p, sort!(bypkg[p])))
        push!(uids, length(rel))
    end
    return Relaxation{P}(rel, rids, uids)
end

"""
    diagnose(sat, reqs; by = identity, max_fixes = 8) :: Diagnosis

Explain why `reqs` cannot be satisfied on `sat`, and list verified fixes.

Must be called on an instance whose satisfiability check has just failed, with
no clauses added since: the frame holding the user-constraint selectors is
popped for the duration and re-asserted afterwards, so the instance is left
exactly as it was found and can be reused.
"""
function diagnose(
    sat  :: SAT{P,V},
    reqs :: AbstractVector{P} = collect(keys(sat.info));
    by   :: Function = identity,
    max_fixes :: Integer = 8,
) where {P,V}
    # a requirement the filter dropped has no installable version at all, so it
    # is not in the instance and cannot be reasoned about there: it is its own
    # conflict, and every fix has to drop it
    absent = sort!(P[p for p in reqs if !haskey(sat.info, p)])
    present = sort!(P[p for p in reqs if haskey(sat.info, p)])
    conflicts = Conflict{P,V}[
        Conflict{P,V}([p], Fact[Requirement(p), Uninstallable(p)],
                      UpstreamFix{P,V}[])
        for p in absent]
    forced = Fact[Requirement(p) for p in absent]

    fixes, versions = with_relaxed_selectors(sat) do
        R = relaxation(sat, present)
        rel = R.rel
        prob = something(sat.prob, Problem(P[]))
        # the whole-column group of each constrained package, for story facts
        bygroup = Dict{P,Int}(rel[i].pkg => i for i in R.user)
        for cluster in cluster_reqs(sat, R, R.reqs)
            creqs = sort!(P[rel[i].pkg for i in cluster])
            # Phase B: which package-pair incompatibilities tell this story,
            # and which single upstream release would end it. Third-party
            # bounds are hard clauses here -- no selector overhead on the hot
            # path -- so it runs on a closure sub-instance (see Closure.jl),
            # which also lets it shrink the user groups in the right order:
            # bounds first, so the actionable explanation survives
            bfacts, ups, kept = bound_story(
                sat.info, prob, creqs, P[rel[i].pkg for i in R.user], by)
            if kept === nothing
                # the closure was too big to explore: shrink the columns
                # against the production instance instead, without the bias
                mus = group_mus(sat, R, Int[cluster; R.user], R.user)
                kept = P[rel[i].pkg for i in mus if rel[i].kind === :user]
            end
            # the packages this story is about: the cluster's requirements, the
            # ones whose constraints the MUS needs, and the ones the
            # incompatibilities name
            told = Set{P}(creqs)
            union!(told, kept)
            for f in bfacts
                push!(told, fact_pkg(f))
                f isa Bound && push!(told, f.dep)
            end
            facts = Fact[Requirement(p) for p in creqs]
            # A package the story already blames shows every source that rules a
            # version of it out, the firm ones included: "the version that would
            # have worked is yanked" explains something, it just cannot be
            # offered as a fix. A package the story does *not* mention shows
            # nothing -- there are plenty of yanked versions in a registry and
            # none of them are news. Emitted once per package, since
            # `group_facts` enumerates all of a package's sources at once.
            for p in sort!(collect(told))
                i = get(bygroup, p, 0)
                i == 0 && continue
                append!(facts, group_facts(sat, rel[i]))
            end
            append!(facts, bfacts)
            push!(conflicts,
                  Conflict{P,V}(creqs, order_chain(facts, P, V), ups))
        end
        order = Int[R.reqs; R.user]
        fixes = enumerate_fixes(sat, R, order, forced; max_fixes, by,
                                acceptable = fix_acceptable(sat))
        # the version lists a renderer needs to compress fact version sets
        # against: only the packages the report actually names, so a
        # registry-scale diagnosis doesn't carry the whole universe around
        vers = Dict{P,Vector{V}}()
        record(p) = haskey(sat.info, p) &&
            (vers[p] = copy(sat.info[p].versions))
        for c in conflicts, f in c.chain
            record(fact_pkg(f))
            f isa Bound && record(f.dep)
        end
        fixes, vers
    end
    sort!(conflicts; by = c -> c.reqs)
    return Diagnosis{P,V}(sort!(collect(reqs)), conflicts, fixes, versions)
end

"""
    superseded(v, versions) :: Bool

Is `v` a prerelease that a release has already replaced — a version `x.y.z-…`
for which `x.y.z` itself is in `versions`?

Telling someone to accept `1.2.3-alpha1` when `1.2.3` exists is not advice, so
admitting prereleases is only ever suggested for versions this returns `false`
for. The test is deliberately restricted to the *same base version*: the
generalization that a newer release supersedes an older prerelease across base
versions was considered and left out, pending evidence that it is what people
want (an unreleased `2.0.0-rc1` is a real thing to want even when `1.9.0` is
out).

Version types with no notion of a prerelease are never superseded, so this costs
nothing outside the `VersionNumber` world.
"""
superseded(v, versions) = false

function superseded(v::VersionNumber, versions)
    isempty(v.prerelease) && return false
    any(versions) do w
        isempty(w.prerelease) && w.major == v.major &&
            w.minor == v.minor && w.patch == v.patch
    end
end

# Would this correction be advice worth giving? A prerelease admission whose
# newly-chosen version has since been released should not be offered. Checked
# against the witness, since that is the version the user would actually end up
# with.
fix_acceptable(sat::SAT) = function (actions, sol)
    for a in actions
        a isa Admission && a.kind === :prerelease || continue
        haskey(sol, a.pkg) || continue
        superseded(sol[a.pkg], sat.info[a.pkg].versions) && return false
    end
    return true
end

# rendering

# compress a subset of a package's versions into runs against the full list,
# e.g. [1,2,3,5] against 1:5 → "1–3, 5"; the whole list → "all versions"
function format_versions(whole::AbstractVector, subset::AbstractVector)
    isempty(subset) && return "no versions"
    idx = sort!(Int[i for i in (findfirst(==(v), whole) for v in subset)
                    if i !== nothing])
    isempty(idx) && return "no versions"
    length(idx) == length(whole) && return "all versions"
    # compress consecutive-in-`whole` versions into runs, then order the runs
    # and each run's endpoints by version value -- independent of whether
    # `whole` is sorted ascending or descending (the registry provider sorts
    # versions descending, which otherwise prints ranges backwards, e.g.
    # "1.11.0–1.0.0")
    runs = Tuple{eltype(whole),eltype(whole)}[]
    i = 1
    while i ≤ length(idx)
        j = i
        while j < length(idx) && idx[j+1] == idx[j] + 1
            j += 1
        end
        a, b = whole[idx[i]], whole[idx[j]]
        push!(runs, (min(a, b), max(a, b)))
        i = j + 1
    end
    sort!(runs; by = first)
    return join((lo == hi ? string(lo) : string(lo, "–", hi)
                 for (lo, hi) in runs), ", ")
end

# the version list to compress against; a package the diagnosis never recorded
# (possible only for hand-built facts) falls back to the subset itself
versions_of(d::Diagnosis{P,V}, p::P) where {P,V} = get(d.versions, p, V[])

render_fact(io::IO, f::Requirement, d::Diagnosis) =
    print(io, "you require ", f.pkg)
render_fact(io::IO, f::Uninstallable, d::Diagnosis) =
    print(io, "no version of ", f.pkg, " is installable")
function render_fact(io::IO, f::UserCompat, d::Diagnosis)
    vers = format_versions(versions_of(d, f.pkg), f.allowed)
    f.label === :requested ?
        print(io, "you requested ", f.pkg, " at ", vers) :
    f.label === :compat ?
        print(io, "your compat restricts ", f.pkg, " to ", vers) :
        print(io, f.label, " restricts ", f.pkg, " to ", vers)
end
render_fact(io::IO, f::Pin, d::Diagnosis) =
    print(io, f.pkg, " is pinned at ", f.version)
# admission knobs get a sentence apiece where one reads better than the
# generic phrasing; an unknown kind still renders sensibly
function render_fact(io::IO, f::Admission, d::Diagnosis)
    whole = versions_of(d, f.pkg)
    n = length(f.forbidden)
    # "all versions" is what `format_versions` says when the set is everything,
    # which does not fit inside these sentences -- give those their own
    if !isempty(whole) && n == length(whole)
        f.kind === :prerelease ?
            print(io, "every version of ", f.pkg,
                  " is a prerelease, which is not allowed") :
        f.kind === :yanked ?
            print(io, "every version of ", f.pkg, " is yanked") :
            print(io, "every version of ", f.pkg, " is ", f.kind,
                  ", which is not allowed")
        return
    end
    vers = format_versions(whole, f.forbidden)
    if f.kind === :prerelease
        print(io, "prereleases of ", f.pkg, " are not allowed (", vers, ")")
    elseif f.kind === :yanked
        print(io, n == 1 ? "version " : "versions ", vers, " of ", f.pkg,
              n == 1 ? " is yanked" : " are yanked")
    else
        print(io, f.kind, " versions of ", f.pkg,
              " are not allowed (", vers, ")")
    end
end
function render_fact(io::IO, f::Bound, d::Diagnosis)
    all_p = versions_of(d, f.pkg)
    src = length(f.versions) == length(all_p) ? "(all versions)" :
          format_versions(all_p, f.versions)
    print(io, f.pkg, " ", src, " works with ", f.dep, " only at ",
          format_versions(versions_of(d, f.dep), f.allowed))
end

render_action(f::Requirement) = string("drop requirement ", f.pkg)
render_action(f::UserCompat)  =
    f.label === :requested ? string("relax the version you requested for ", f.pkg) :
    f.label === :compat    ? string("relax your compat on ", f.pkg) :
                             string("relax the ", f.label, " on ", f.pkg)
render_action(f::Pin)         = string("unpin ", f.pkg)
render_action(f::Admission)   =
    f.kind === :prerelease ? string("allow prereleases of ", f.pkg) :
    f.kind === :yanked ?
        string("allow the yanked version", length(f.forbidden) == 1 ? " " : "s ",
               join(f.forbidden, ", "), " of ", f.pkg) :
        string("allow ", f.kind, " versions of ", f.pkg)

# packages worth naming in a fix's resulting solution: the requirements plus the
# dependencies actually implicated in some conflict. A verified solution lists
# the entire transitive closure at whatever versions it happened to pick, which
# is mostly noise -- the contested packages are what a fix meaningfully changes.
# The complete assignment is still available on `Fix.solution` /
# `UpstreamFix.solution` for programmatic use (e.g. emitting a manifest).
function relevant_pkgs(d::Diagnosis{P,V}) where {P,V}
    rel = Set{P}(d.reqs)
    union!(rel, story_pkgs(d))
    return rel
end

# just the packages some conflict's story names
function story_pkgs(d::Diagnosis{P,V}) where {P,V}
    rel = Set{P}()
    for c in d.conflicts, f in c.chain
        push!(rel, fact_pkg(f))
        f isa Bound && push!(rel, f.dep)
    end
    return rel
end

# Packages worth naming in a witness: the requirements, and the ones some
# conflict's story is about. A platform package the user never chose -- julia --
# stays out of this by being modelled as a dependency edge rather than a
# requirement, which is the convention integrators are asked to follow; it then
# appears only when a story actually blames it.
function witness_pkgs(d::Diagnosis{P,V}, sol::AbstractDict{P,V}) where {P,V}
    rel = relevant_pkgs(d)
    P[p for p in keys(sol) if p ∈ rel]
end

function render_solution(d::Diagnosis{P,V}, sol::AbstractDict{P,V}) where {P,V}
    pkgs = witness_pkgs(d, sol)
    isempty(pkgs) && return "nothing"
    reqset = Set(d.reqs)
    sort!(pkgs)
    sort!(pkgs; by = p -> p ∉ reqset) # requirements first (stable)
    join((string(p, " ", sol[p]) for p in pkgs), ", ")
end

# the demoted bound facts of a conflict, as package => versions, in a
# deterministic order
function incidental_versions(c::Conflict{P,V}) where {P,V}
    byp = Dict{P,Vector{V}}()
    for f in c.chain
        f isa Bound && f.incidental || continue
        append!(get!(() -> V[], byp, f.pkg), f.versions)
    end
    return [p => sort!(unique!(byp[p])) for p in sort!(collect(keys(byp)))]
end

# compact summary
function Base.show(io::IO, d::Diagnosis)
    print(io, "Diagnosis(", summary_counts(d), ")")
end

# "1 conflict, 2 fixes" -- the shape of the answer, for either show
function summary_counts(d::Diagnosis)
    nc = length(d.conflicts)
    nf = length(d.fixes)
    string(nc, nc == 1 ? " conflict, " : " conflicts, ",
           nf, nf == 1 ? " fix" : " fixes")
end

# full report. It opens by saying what it is: a `Diagnosis` *is* the answer
# `resolve` gives to requirements it cannot satisfy, so one showing up at the
# REPL has to read as that answer and not as a stray object.
function Base.show(io::IO, ::MIME"text/plain", d::Diagnosis)
    println(io, "Unsatisfiable — ", summary_counts(d),
            isempty(d.conflicts) ? "" : ":")
    for (n, c) in enumerate(d.conflicts)
        n > 1 && println(io)
        rs = join(c.reqs, ", ")
        println(io, "Conflict ", n, ": ", rs,
                length(c.reqs) == 1 ? " cannot be satisfied." :
                " cannot be satisfied together.")
        for f in c.chain
            f isa Bound && f.incidental && continue
            print(io, "  • ")
            render_fact(io, f, d)
            println(io)
        end
        # the facts that only close off versions the query already excludes:
        # one summary line per package, rather than a bullet apiece competing
        # with the lines that name versions the reader could actually get
        for (p, vs) in incidental_versions(c)
            println(io, "    (likewise for ", p, " ",
                    format_versions(versions_of(d, p), vs),
                    ", which your constraints already rule out)")
        end
    end
    ups = [u for c in d.conflicts for u in c.upstream]
    if !isempty(d.fixes)
        println(io)
        println(io, "Verified fixes:")
        for (n, fix) in enumerate(d.fixes)
            println(io, "  ", n, ". ",
                    join(map(render_action, fix.actions), " and "))
            println(io, "     → allows: ", render_solution(d, fix.solution))
        end
    elseif !isempty(d.conflicts)
        println(io)
        println(io, isempty(ups) ?
            "There is nothing you can change that would fix this." :
            "There is nothing you can change that would fix this — only an upstream release could.")
    end
    if !isempty(ups)
        println(io)
        println(io, "Upstream fixes:")
        for u in Iterators.take(ups, MAX_UPSTREAM_SHOWN)
            println(io, "  • a release of ", u.bound.pkg, " or ", u.bound.dep,
                    " relaxing their compat on each other")
            println(io, "    → allows: ", render_solution(d, u.solution))
        end
        extra = length(ups) - MAX_UPSTREAM_SHOWN
        extra > 0 && println(io, "  (", extra,
            " more possible upstream fix", extra == 1 ? "" : "es",
            " not shown)")
    end
end

# How many upstream suggestions a report prints. A tangled conflict can produce
# a dozen, all true and most useless -- ten for one real project -- and
# they are ranked (see `rank_upstream!`), so the tail is the least plausible.
# The structured `Diagnosis` keeps every one of them.
const MAX_UPSTREAM_SHOWN = 3
