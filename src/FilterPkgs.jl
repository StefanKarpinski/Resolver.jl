"""
    Universe{P,V}

The universe a single resolve is run against: one
[`PkgInfo`](@ref Resolver.PkgInfo) per package, plus `reps` — per package, per
class, the index into that package's `versions` of the member the class stands
for under this query, or `0` when the query admits none of them.

`info` carries the package universe ([`pkg_info`](@ref Resolver.pkg_info)) and
`reps` carries everything the query determines about it. A user constraint acts
by forbidding *members* of classes. Class members are indistinguishable to
everything in the registry, so forbidding some of them only changes which one
stands for the class; forbidding all of them leaves the class **empty**, hence
**deactivated** — it cannot be selected.

Deactivated classes are **not** deleted. They stay in the matrices as the
columns a later relaxation of whatever emptied them would need, and the passes
read the deactivation the way they read a conflict with an always-present
version: reachability continues its prefix past a deactivated class while still
following its dependencies, redundancy elimination neither deletes one nor lets
one license a deletion, and the SAT instance forbids each of them with one
relaxable unit clause. BitKernels.jl's header spells out how a class's
*in-universe* flag in the matrix and its *activated* state here divide the work.

## Why deactivation is read off `reps`

`reps` has to exist regardless: a class's representative is what names the
version in the answer, and it is what the class competes at when the classes
are ranked. A class the query admits no member of simply has no representative,
so "unavailable for this query" costs nothing extra to represent. A second flag
bit in the matrix would work equally well; it would just encode a fact `reps`
already carries.

## Why the query's state is separate from the artifact

The invariant worth protecting is that the artifact `info` is
*query-independent*: it is a function of the registry alone, `PkgInfoFiles`
writes and reads it, and one of them answers every query, whatever that query
constrains or prefers. `reps` is the opposite — one query's worth of state,
recomputed per resolve.

Keeping the two in separate values is what makes that invariant checkable
rather than merely intended. Were `reps` a field of `PkgInfo`, the persisted
type would contain query state: serialization would have to know to skip it,
and a freshly loaded artifact would carry a field whose only honest value is
"not computed yet" — with `0`, which means *deactivated*, as its natural
uninitialized value, so forgetting to fill it in yields a universe where
nothing can be selected rather than an error. Whatever shape this ends up
taking, that is the property it has to preserve: nothing a query decides may
be reachable from the thing that gets cached and shared.
"""
struct Universe{P,V}
    info :: Dict{P, PkgInfo{P,V}}
    reps :: Dict{P, Vector{Int}}
end

# the universe an unconstrained query makes of an artifact: every class stands
# for its first member (the canonical order's best), nothing is deactivated
default_reps(info::AbstractDict{P, PkgInfo{P,V}}) where {P,V} =
    Dict{P,Vector{Int}}(
        p => Int[first(mem) for mem in info_p.members]
        for (p, info_p) in info)

Universe(info::Dict{P, PkgInfo{P,V}}) where {P,V} =
    Universe{P,V}(info, default_reps(info))
Universe(info::AbstractDict{P, PkgInfo{P,V}}) where {P,V} =
    Universe{P,V}(Dict{P, PkgInfo{P,V}}(info), default_reps(info))

# per package, a mask of its deactivated (empty) classes — only for the
# packages that have one, so that an unconstrained resolve, where there are
# none at all, does no per-class work downstream
function deactivations(univ::Universe{P,V}) where {P,V}
    deact = Dict{P,BitVector}()
    for (p, reps_p) in univ.reps
        any(iszero, reps_p) || continue
        deact[p] = BitVector(iszero(r) for r in reps_p)
    end
    return isempty(deact) ? EmptyDict{P,BitVector}() : deact
end

"""
    version_permutations(info, order) :: Union{Nothing, Dict{P, Vector{Int}}}

Per package, the permutation that puts its versions into the ordering `order`
asks for: `perm[r]` is the index, in the T1 version list, of the `r`-th best
version. Only the packages that actually need reordering get an entry, and the
whole result is `nothing` when none of them do — so the *canonical* order, which
is the one the T1 artifact carries, costs nothing.

`order` is `nothing` for the canonical order, or a callable mapping a package to
a `lt` comparator over its versions (`lt(u, v)` = "u is preferred to v"). The
ordering is the one query knob that is not a *constraint*: it does not change
which solutions are valid, only which valid solution is optimal. That is why it
is a `resolve` parameter rather than part of the [`Problem`](@ref), and why the
T1 artifact can be built and cached without knowing it (see the manual's Theory
section, on order-free deletion and on interchangeability classes being
ordering-independent while the representative choice is not).

This is the ordering of *versions*; what the universe is laid out in is the
ordering of *classes* it induces, which is
[`class_ranking`](@ref Resolver.class_ranking)'s business.
"""
version_permutations(
    info  :: AbstractDict{P, PkgInfo{P,V}},
    order :: Nothing,
) where {P,V} = nothing

function version_permutations(
    info  :: AbstractDict{P, PkgInfo{P,V}},
    order,
) where {P,V}
    perms = Dict{P, Vector{Int}}()
    for (p, info_p) in info
        vers = info_p.versions
        lt = order(p)
        # the common case by far — a comparator that agrees with the canonical
        # order, e.g. "newest first" on a registry that already lists versions
        # that way — is detected in one linear pass and costs nothing more
        issorted(vers; lt) && continue
        perms[p] = sortperm(vers; lt, alg = MergeSort)
    end
    return isempty(perms) ? nothing : perms
end

"""
    class_ranking(info, prob, perms) :: (reps, perms′)

What a query makes of an artifact's classes, before anything is materialized:
per package, which member each class stands for (`reps`, a version index per
class, `0` for a class with no admissible member) and how the classes have to be
reordered to be in rank order (`perms′`, an entry per package that needs it,
`nothing` when none does).

**Representative.** A class's members are indistinguishable to the registry, so
the only thing a user constraint can do to a class is forbid some of its
members. What is left, the class *admits*; the best of those in the active
version ordering — `perms`, see
[`version_permutations`](@ref Resolver.version_permutations) — represents it.
A class admitting nothing is empty, hence deactivated, and has no
representative.

**Rank.** A class ranks where the member it stands for ranks, because that is
the version choosing it yields: a class whose best member is forbidden competes
at its best *admissible* member, and can therefore fall behind classes it
outranks in the artifact. An empty class ranks where its best member would, the
position it would compete from if the constraint that emptied it were lifted.
Every class's key is a distinct version rank, so this is a total order.

Both outputs are per-resolve state: the partition they index is the registry's,
but which member stands for a class, and hence how the classes are ordered,
belongs to the query.
"""
function class_ranking(
    info  :: AbstractDict{P, PkgInfo{P,V}},
    prob  :: Problem{P} = Problem(P[]),
    perms :: Union{Nothing, AbstractDict{P, Vector{Int}}} = nothing,
) where {P,V}
    excl = exclusion_masks(info, prob)
    reps = Dict{P, Vector{Int}}()
    perms′ = Dict{P, Vector{Int}}()
    key = Int[] # scratch: per class, the rank of the member it competes at
    for (p, info_p) in info
        m = nclasses(info_p)          # one entry of `reps` and `key` per class
        nv = length(info_p.versions)  # ... walking this many versions to fill them
        cls = info_p.classes
        e = get(excl, p, nothing)
        perm = perms === nothing ? nothing : get(perms, p, nothing)
        reps_p = zeros(Int, m)
        resize!(key, m)
        fill!(key, 0)
        # walk the versions best first, so the first admissible member of a
        # class found is the one the active order ranks best. i indexes
        # versions, j the class each one is in
        for r = 1:nv
            i = perm === nothing ? r : perm[r]
            j = cls[i]
            reps_p[j] == 0 || continue
            e === nothing || !e[i] || continue
            reps_p[j] = i
            key[j] = r
        end
        if any(iszero, reps_p)
            # the empty classes rank where their best member would
            for r = 1:nv
                j = cls[perm === nothing ? r : perm[r]]
                key[j] == 0 || continue
                key[j] = r
            end
        end
        reps[p] = reps_p
        issorted(key) || (perms′[p] = sortperm(key))
    end
    return reps, isempty(perms′) ? nothing : perms′
end

"""
    prepare_pkg_info(info, prob, [info′]; order = nothing) :: Universe

The per-resolve half of preprocessing, run on a T1 artifact (see
[`pkg_info`](@ref Resolver.pkg_info)) to produce the universe the SAT instance
is built over. Two steps:

1. **Rank** each package's classes for this query (`class_ranking`): pick the
   member each class stands for and lay the classes out in the order those
   members put them in. Under the canonical ordering and no constraints this is
   the layout the artifact already has, and costs nothing.
2. **Filter** the result for reachability and redundancy against the actual
   requirements and constraints (`filter_pkg_info!`).

Filtering has to come second: both of its passes are stated in terms of the rank
order, so the layout has to exist first.

The user's constraints reach the matrices through nothing but step 1's
representatives. A constrained query and an unconstrained one produce the same
matrices, up to the class reordering their representatives induce; the whole of
what a constraint does is empty classes, and an empty class is deactivated
rather than deleted, split or rewritten. That is what makes a class's
deactivation the only thing there is to relax.

The result is a fresh universe — `info` is left untouched, so a T1 artifact
stays reusable across resolves — unless the caller passes its own scratch dict
as `info′`, or `info` itself when it owns it (which reordering rules out, since
it rebuilds the matrices).
"""
prepare_pkg_info(
    info  :: AbstractDict{P, PkgInfo{P,V}},
    prob  :: Problem{P},
    info′ :: Dict{P, PkgInfo{P,V}} = Dict{P, PkgInfo{P,V}}();
    order = nothing, # nothing = canonical, else package -> `lt` comparator
) where {P,V} =
    filter_pkg_info!(rank_pkg_info(info, prob, info′; order), prob)

# step 1 on its own: the query's universe before anything is filtered out of
# it. This is where a constraint's entire effect lands — the representatives
# and the class order they induce — so it is also the shape in which "no user
# constraint reaches a matrix" is a statement one can check.
function rank_pkg_info(
    info  :: AbstractDict{P, PkgInfo{P,V}},
    prob  :: Problem{P} = Problem(P[]),
    info′ :: Dict{P, PkgInfo{P,V}} = Dict{P, PkgInfo{P,V}}();
    order = nothing,
) where {P,V}
    perms = version_permutations(info, order)
    reps, cperms = class_ranking(info, prob, perms)
    # reordering reads every matrix while writing a differently laid out one,
    # so the in-place shortcut is only available when nothing moves
    cperms === nothing || info′ !== info ||
        (info′ = Dict{P, PkgInfo{P,V}}())
    univ = Universe{P,V}(info′, Dict{P, Vector{Int}}())
    return copy_ranked!(univ, info, reps, cperms)
end

# Copy `info` into `univ`, laid out in the class order `perms` gives (the order
# it already has when that is `nothing`), with every class and column marked
# active — deactivation is `univ.reps`, not a flag, precisely so that the
# filtering passes below can tell "nothing can choose this class" apart from
# "this class is not in the universe". When the destination dict is the source,
# only the flags are normalized.
function copy_ranked!(
    univ  :: Universe{P,V},
    info  :: AbstractDict{P, PkgInfo{P,V}},
    reps  :: Dict{P, Vector{Int}},
    perms :: Union{Nothing, AbstractDict{P, Vector{Int}}} = nothing,
) where {P,V}
    info′ = univ.info
    reps′ = univ.reps
    src = Int[] # scratch: output column => source column
    for (p, info_p) in info
        X = info_p.conflicts
        m = nclasses(info_p)
        n = size(X, 2) - 1
        if info′ === info
            @assert perms === nothing
            X[1:m, end] .= true
            X[m+1, 1:n] .= true
            reps′[p] = reps[p]
            continue
        end
        perm = perms === nothing ? nothing : get(perms, p, nothing)
        # a reordered partner's interaction block is indexed by that partner's
        # classes, so its columns move with them; every other column stays put
        cols = nothing
        if perms !== nothing
            for (q, off) in info_p.interacts
                perm_q = get(perms, q, nothing)
                perm_q === nothing && continue
                if cols === nothing
                    cols = resize!(src, n)
                    for j = 1:n
                        cols[j] = j
                    end
                end
                for (r, j) in enumerate(perm_q)
                    cols[off + r] = off + j
                end
            end
        end
        X′ = perm === nothing && cols === nothing ? copy(X) :
             relaid_conflicts(X, m, n, perm, cols)
        X′[1:m, end] .= true
        X′[m+1, 1:n] .= true
        if perm === nothing
            info′[p] = PkgInfo(
                copy(info_p.versions), copy(info_p.classes),
                copy(info_p.members), copy(info_p.depends),
                copy(info_p.interacts), X′)
            reps′[p] = copy(reps[p])
        else
            # the partition follows the classes: class `r` of the new layout is
            # class `perm[r]` of the old one
            inv = invperm(perm)
            info′[p] = PkgInfo(
                copy(info_p.versions), Int[inv[i] for i in info_p.classes],
                info_p.members[perm], copy(info_p.depends),
                copy(info_p.interacts), X′)
            reps′[p] = reps[p][perm]
        end
    end
    return univ
end

# A conflicts matrix laid out under a reordering: row `r` reads source row
# `perm[r]` (this package's classes) and column `j` reads source column
# `cols[j]` (a reordered partner's block); `nothing` for either means unchanged.
# Only the class rows of the `n` conflict columns are written — the caller
# sets the flag row and the class-flag column. Columns whose rows do not move
# are blitted word by word, which is the whole matrix whenever only partners
# were reordered.
function relaid_conflicts(
    X    :: BitMatrix,
    m    :: Int,
    n    :: Int,
    perm :: Union{Nothing, Vector{Int}},
    cols :: Union{Nothing, Vector{Int}},
)
    X′ = falses(size(X, 1), n + 1)
    W = col_words(X)
    ch = X.chunks
    ch′ = X′.chunks
    @inbounds for r = 1:n
        j = cols === nothing ? r : cols[r]
        base = (j - 1) * W
        base′ = (r - 1) * W
        if perm === nothing
            for w = 1:W
                ch′[base′ + w] = ch[base + w]
            end
        else
            for i′ = 1:m
                i = perm[i′]
                b = (ch[base + ((i - 1) >> 6) + 1] >> ((i - 1) & 63)) & 1
                ch′[base′ + ((i′ - 1) >> 6) + 1] |= b << ((i′ - 1) & 63)
            end
        end
    end
    return X′
end

filter_pkg_info!(
    info :: Dict{P, PkgInfo{P,V}},
    reqs :: SetOrVec{P} = keys(info),
) where {P,V} = filter_pkg_info!(Universe(info), Problem(reqs))

function filter_pkg_info!(
    univ :: Universe{P,V},
    prob :: Problem{P},
) where {P,V}
    info = univ.info
    reqs = prob.reqs
    # the query's constraints reach the passes below as the *deactivated*
    # classes and nothing else — no exclusion column, no virtual package, no
    # row of their own — because a class its members all fail is the only mark
    # a constraint can leave on a universe whose rows are classes. the masks
    # are index-based, so they are rebuilt after every deletion round; only the
    # packages with an empty class cost anything
    #
    # the installability prune first: it needs no reachability marks and it is
    # what makes the deletions safe — dropping a package while a kept class
    # still depends on it would leave a dangling name. then redundancy
    # elimination — it needs no reachability marks and does most of the
    # shrinking — then alternate reachability and redundancy until neither
    # deletes anything: each deleted class can expose more of both kinds
    # of pruning, and every round preserves the resolved answer (see the
    # theory docs on iterating the filter). requirement packages survive
    # filtering unless they have no installable class at all — nothing
    # can satisfy those — so the requirements remain valid across rounds.
    # rounds strictly shrink the total class count, so the loop terminates
    mark_installable!(info)
    mark_necessary!(info, deactivations(univ))
    drop_unmarked!(univ)
    while true
        total = sum(nclasses(i) for i in values(info); init = 0)
        deact = deactivations(univ)
        mark_reachable!(info, reqs, deact)
        mark_necessary!(info, deact)
        drop_unmarked!(univ)
        sum(nclasses(i) for i in values(info); init = 0) < total ||
            break
    end
    return univ
end

"""
    find_reachable(info, reqs) :: Dict{P, Int}

This function finds a minimal "reachable" subset of packages and classes that
could appear in pareto-optimal solutions to version resolution for the given set
of required "root" packages, using the following recursive logic:

- P in reqs => P[1] reachable
- P[i] reachable & P[i] depends on D => D[1] reachable
- P[i] reachable & P[i] conflicts w. reachable => P[i+1] reachable
- D[end] conflicts w. reachable & P[i] depends on D => P[i+1] reachable

The optional `deact` argument gives, per package, a mask of its deactivated
classes — the ones this query admits no member of. Such a class cannot be
selected, so the package cannot be installed at it and the prefix has to
continue past it, exactly as it continues past a class conflicting with
something present:

- P[i] reachable & P[i] deactivated => P[i+1] reachable

A deactivated class's dependencies are followed all the same, which is what
keeps the cone behind it in the universe: it is precisely the cone a relaxation
of whatever emptied the class would need to find still there.

The function returns a dictionary mapping packages to the maximum class index
of that package that could be reached in an optimal solution. If a package
cannot appear in an optimal solution, it will not appear in this dictionary.
"""
function find_reachable(
    info  :: Dict{P, PkgInfo{P,V}},
    reqs  :: SetOrVec{P} = keys(info),
    deact :: AbstractDict{P, BitVector} = EmptyDict{P,BitVector}(),
) where {P,V}
    # intern packages as integer indices so that the fixpoint loop below
    # hashes no package names; all derived tables are indexed by pkg id
    pkgs = sort!(collect(keys(info)))
    N = length(pkgs)
    ix = Dict{P,Int}(p => i for (i, p) in enumerate(pkgs))
    infos = Vector{PkgInfo{P,V}}(undef, N)
    ncls = Vector{Int}(undef, N)
    # interned depends, same order (id 0: dependency absent from info,
    # i.e. a package with no installable class — see below)
    deps = Vector{Vector{Int}}(undef, N)
    partners = Vector{Vector{Tuple{Int,Int}}}(undef, N)
    # interned deactivation masks (nothing: nothing empty here, the norm)
    deacts = Vector{Union{Nothing,BitVector}}(nothing, N)
    for p = 1:N
        info_p = info[pkgs[p]]
        infos[p] = info_p
        ncls[p] = nclasses(info_p)
        deacts[p] = get(deact, pkgs[p], nothing)
        deps[p] = Int[get(ix, q, 0) for q in info_p.depends]
        # (partner id, offset of p's class block in the partner's matrix)
        prt = Tuple{Int,Int}[
            (ix[q], info[q].interacts[pkgs[p]])
            for q in keys(info_p.interacts)
        ]
        partners[p] = sort!(prt)
    end

    # meaning (both map packages to class indices):
    #   - reach tracks fully processed reachable classes
    #   - queue tracks newly reachable classes not yet processed
    reach = zeros(Int, N)
    queue = zeros(Int, N)
    stack = Int[] # packages with queued work
    instack = falses(N)

    function enqueue(p::Int, j::Int)
        queue[p] = j
        if !instack[p]
            instack[p] = true
            push!(stack, p)
        end
    end

    # add next active class of p *after* i to the queue
    # do nothing if there's already class > i in reach/queue
    function next(p::Int, i::Int)
        reach[p] > i && return false
        queue[p] > i && return false
        m = ncls[p]
        X = infos[p].conflicts
        # first active class after i (the flag column is the active set)
        j = col_min_from(X, size(X, 2), i + 1, m)
        # j == 0: we're out of classes (i.e. saturated; see below)
        enqueue(p, j == 0 ? m + 1 : j)
        return true
    end

    rdeps = [Dict{Int,Int}() for _ = 1:N] # reverse dependency map
    # rdeps[p][q] == k means
    #   "k is latest reachable class of q that depends on p"

    for p in reqs
        # a requirement with no installable class reaches nothing
        i = get(ix, p, 0)
        i == 0 || enqueue(i, 1)
    end

    # notation:
    #   - p, q: packages
    #   - info_p = infos[p]
    #   - indices of classes of package p: i, j
    #   - indices of classes of package q: k
    while !isempty(stack)
        # get unprocessed package + class
        p = pop!(stack)
        instack[p] = false
        i = queue[p]
        queue[p] = 0
        info_p = infos[p]
        m = ncls[p]
        # check for saturation
        #   p saturated means: conflicts can force p to be uninstallable
        #   saturation represented by i > ncls[p]
        if i > m
            # p has become saturated
            for (q, k) in rdeps[p]
                # q@k depends on p, therefore
                # q@k conflicts with p being uninstallable
                # p being saturated means that can happen
                next(q, k)
            end
        end
        # process each newly reachable class of p
        deact_p = deacts[p]
        for j = reach[p]+1:min(i, m)
            # p@j is deactivated: it cannot be selected, so it cannot be
            # what p is installed at — the prefix has to continue past it.
            # its dependencies below are followed anyway, which is what keeps
            # the packages behind it in the universe
            deact_p === nothing || !deact_p[j] || next(p, j)
            # dependencies
            for (k, q) in enumerate(deps[p])
                info_p.conflicts[j, k] || continue
                # p@j depends on q
                if q == 0
                    # q has no installable class at all, so neither has
                    # p@j: p can only be installed past it
                    next(p, j)
                    continue
                end
                rdeps_q = rdeps[q]
                rdeps_q[p] = max(get(rdeps_q, p, 0), j)
                next(q, 0) # q can be required
                # check if q is saturated:
                if reach[q] > ncls[q]
                    # p@j depends on q, therefore
                    # p@j conflicts with q being uninstallable
                    # q being saturated means that can happen
                    next(p, j)
                end
            end
            # find p@j's conflicts with reached classes of each partner:
            # conflicts are stored symmetrically, so scan the partner-side
            # column Y[k, c+j] == X[j, b+k], which is contiguous
            for (q, c) in partners[p]
                r = min(reach[q], ncls[q])
                r > 0 || continue
                k = col_max_upto(infos[q].conflicts, c + j, r)
                k == 0 && continue
                # p@j conflicts with q@k (the highest such k; pushing past
                # it subsumes pushing past any lower conflicting class)
                next(p, j)
                next(q, k)
            end
        end
        # update the reach map
        reach[p] = i
    end

    return Dict{P,Int}(pkgs[p] => reach[p] for p = 1:N if reach[p] > 0)
end

function mark_reachable!(
    info  :: Dict{P, PkgInfo{P,V}},
    reqs  :: SetOrVec{P} = keys(info),
    deact :: AbstractDict{P, BitVector} = EmptyDict{P,BitVector}(),
) where {P,V}
    reach = find_reachable(info, reqs, deact)
    for (p, info_p) in info
        r = min(get(reach, p, 0), nclasses(info_p))
        info_p.conflicts[1:r, end] .= true
        info_p.conflicts[r+1:end, end] .= false
    end
    return reach
end

function mark_installable!(
    info :: Dict{P, <: PkgInfo{P}},
) where {P}
    # a package with no in-universe classes cannot be installed, so no class
    # depending on it can be either: take those out too, and repeat until
    # nothing more goes, since each removal can empty another package (the
    # literature calls this arc consistency). It is a sound approximation of
    # deleting the classes no valid solution contains (see the theory docs),
    # and it is also what keeps `depends` from naming a package that
    # `drop_unmarked!` deletes as empty
    #
    # a dependency naming a package that is not in `info` at all counts the
    # same way: it cannot be satisfied either, and a deletion elsewhere (an
    # earlier `drop_unmarked!`, say) is exactly how one comes about. `info`
    # does not change here, so the absentees are found once
    absent = Set{P}()
    for info_p in values(info), q in info_p.depends
        haskey(info, q) || push!(absent, q)
    end
    empties = Set{P}()
    while true
        empty!(empties)
        union!(empties, absent)
        for (p, info_p) in info
            X = info_p.conflicts
            any(X[i, end] for i = 1:nclasses(info_p)) && continue
            push!(empties, p)
        end
        isempty(empties) && return info
        dirty = false
        for info_p in values(info)
            X = info_p.conflicts
            for (k, q) in enumerate(info_p.depends)
                q in empties || continue
                # q is uninstallable, drop the p@i that depend on it
                for i = 1:nclasses(info_p)
                    X[i, end] || continue
                    X[i, k] || continue
                    X[i, end] = false
                    dirty = true
                end
            end
        end
        dirty || return info
    end
end

function mark_necessary!(
    info  :: Dict{P, PkgInfo{P,V}},
    deact :: AbstractDict{P, BitVector} = EmptyDict{P,BitVector}(),
) where {P,V}
    # intern packages as integer indices so that the work loop below
    # hashes no package names (sorted so the processing order — and with
    # it any order-dependent tie — is deterministic)
    pkgs = sort!(collect(keys(info)))
    N = length(pkgs)
    ix = Dict{P,Int}(p => i for (i, p) in enumerate(pkgs))
    infos = Vector{PkgInfo{P,V}}(undef, N)
    ncls = Vector{Int}(undef, N)
    partners = Vector{Vector{NTuple{3,Int}}}(undef, N)
    # interned deactivation masks (nothing: nothing empty here, the norm)
    deacts = Vector{Union{Nothing,BitVector}}(nothing, N)
    for p = 1:N
        info_p = info[pkgs[p]]
        infos[p] = info_p
        ncls[p] = nclasses(info_p)
        deacts[p] = get(deact, pkgs[p], nothing)
        # (partner id, partner's class block offset in p's matrix,
        #  p's class block offset in the partner's matrix)
        prt = NTuple{3,Int}[
            (ix[q], b, info[q].interacts[pkgs[p]])
            for (q, b) in info_p.interacts
        ]
        partners[p] = sort!(prt)
    end
    # some work buffers
    A = UInt64[]        # candidate class mask
    D = UInt64[]        # per-class domination candidate masks
    T = UInt64[]        # mask of classes still tracked by the sweep
    R = Int[]           # redundant indices vector

    # the classes redundancy may reason about: the in-universe ones, less the
    # deactivated ones — the pass is the one place that reads *both* per-class
    # bits (see BitKernels.jl's header). a deactivated class must not be
    # *deleted* — the relaxation that would want it back has to find it still
    # there — and it must not *dominate*, because domination says "this better
    # class will be chosen instead", which is worth nothing when the better
    # class is one nothing can choose. dropping them from the candidate set
    # does both — it takes them off both sides of the test — and it keeps the
    # test itself activation-free: what a row says is compared in full,
    # including conflicts against partner classes this query happens to have
    # emptied
    function candidates!(A::Vector{UInt64}, p::Int)
        info_p = infos[p]
        X = info_p.conflicts
        resize!(A, col_words(X))
        col_copy!(A, X, size(X, 2))
        clear_rows_above!(A, ncls[p])
        d = deacts[p]
        if d !== nothing
            dc = d.chunks
            @inbounds for w = 1:min(length(A), length(dc))
                A[w] &= ~dc[w]
            end
        end
        return A
    end

    # initialize active column flags
    for p = 1:N
        info_p = infos[p]
        X = info_p.conflicts
        m = ncls[p]
        n = size(X, 2) - 1
        candidates!(A, p)
        for j = 1:n
            # we don't have to look at columns that have no
            # conflicts for any candidate class of the package
            X[m+1, j] = col_intersects(X, j, A)
        end
        # conflict columns for deleted partner classes are also
        # irrelevant: no model over the surviving classes can include
        # the partner, so the conflict cannot constrain anything.
        # a partner class that is merely *deactivated* is a different
        # matter — it is still in the universe, so the conflict still
        # counts, and the flag it is read from does not know about
        # deactivation precisely so that it does
        for (q, b, c) in partners[p]
            Y = infos[q].conflicts
            for k = 1:ncls[q]
                X[m+1, b+k] &= Y[k, end]
            end
        end
    end
    # process every package once, cheapest first (stable sort breaks ties
    # by name); packages woken by a partner's deletions are appended to
    # the queue unless already pending
    queue = sort!(collect(1:N), by = p -> length(partners[p]))
    inwork = trues(N)
    head = 1
    while head ≤ length(queue)
        p = queue[head]
        head += 1
        inwork[p] || continue
        inwork[p] = false
        # get conflicts & dimensions
        info_p = infos[p]
        X = info_p.conflicts
        m = ncls[p]
        m > 1 || continue # unique class cannot be redundant
        n = size(X, 2) - 1
        W = col_words(X)
        # candidate class mask
        candidates!(A, p)
        nact = 0
        @inbounds for w = 1:W
            nact += count_ones(A[w])
        end
        nact > 1 || continue # unique class cannot be redundant
        # find redundant classes: i < j and X[i, k] => X[j, k] for all
        # active k means an earlier class is strictly more compatible,
        # therefore i will always be chosen instead of j. the candidates
        # dominated by i are the candidate classes worse than i that have a
        # constraint in every column i does. sweep the active columns
        # once, restricting the candidates of every class with a bit in
        # the column to the column's row set; classes whose candidate
        # set empties drop out of the sweep. a class dominated only by
        # dominated classes is also dominated by their (transitively
        # live) dominators, so batching over the initial candidate set finds
        # exactly the sequentially dominated classes.
        resize!(D, m * W)
        resize!(T, W)
        copyto!(T, A)
        @inbounds for w = 1:W
            c = A[w]
            while !iszero(c)
                i = ((w - 1) << 6) + trailing_zeros(c) + 1
                c &= c - 1
                # candidates of i: candidate classes worse than i
                o = (i - 1) * W
                for w′ = 1:W
                    D[o + w′] = A[w′]
                end
                wt = i >> 6
                for w′ = 1:min(wt, W)
                    D[o + w′] = 0
                end
                r = i & 63
                if r != 0 && wt < W
                    D[o + wt + 1] &= ~((UInt64(1) << r) - 1)
                end
                live = UInt64(0)
                for w′ = 1:W
                    live |= D[o + w′]
                end
                iszero(live) && (T[w] &= ~(UInt64(1) << ((i - 1) & 63)))
            end
        end
        prev = -1 # chunk base of the previous active column
        for k = 1:n
            X[m+1, k] || continue # inactive column
            base = (k - 1) * W
            # compat bounds are ranges, so identical adjacent columns are
            # common — and re-ANDing an identical column changes nothing
            if prev ≥ 0
                same = true
                @inbounds for w = 1:W
                    if X.chunks[base + w] != X.chunks[prev + w]
                        same = false
                        break
                    end
                end
                same && continue
            end
            prev = base
            @inbounds for w = 1:W
                c = X.chunks[base + w] & T[w]
                while !iszero(c)
                    i = ((w - 1) << 6) + trailing_zeros(c) + 1
                    c &= c - 1
                    o = (i - 1) * W
                    live = UInt64(0)
                    for w′ = 1:W
                        live |= (D[o + w′] &= X.chunks[base + w′])
                    end
                    iszero(live) && (T[w] &= ~(UInt64(1) << ((i - 1) & 63)))
                end
            end
        end
        # union of all candidate sets = the dominated classes
        fill!(T, UInt64(0)) # reuse as the union accumulator
        @inbounds for w = 1:W
            c = A[w]
            while !iszero(c)
                i = ((w - 1) << 6) + trailing_zeros(c) + 1
                c &= c - 1
                o = (i - 1) * W
                for w′ = 1:W
                    T[w′] |= D[o + w′]
                end
            end
        end
        empty!(R)
        @inbounds for w = 1:W
            c = T[w]
            while !iszero(c)
                push!(R, ((w - 1) << 6) + trailing_zeros(c) + 1)
                c &= c - 1
            end
        end
        isempty(R) && continue
        # deactivate redundant classes
        X[R, end] .= false
        for (q, b, c) in partners[p]
            infos[q].conflicts[ncls[q] + 1, c .+ R] .= false
            if !inwork[q] # can create new redundancies
                inwork[q] = true
                push!(queue, q)
            end
        end
    end
end

drop_unmarked!(univ::Universe) = drop_unmarked!(univ.info, univ.reps)

function drop_unmarked!(
    info :: Dict{P, <: PkgInfo{P}},
    reps :: Union{Nothing, Dict{P, Vector{Int}}} = nothing,
) where {P}
    # This pass reads the *in-universe* bit and nothing else: it knows about
    # the deleting passes' verdicts and knows nothing about deactivation, which
    # is exactly why deactivation is not stored as that bit — an emptied class
    # must come out of here intact (see BitKernels.jl's header). `reps` is
    # carried along only to be renumbered with the classes it indexes.
    #
    # first pass: per package, compute the kept classes and kept columns
    # from the (definitive) class flags — all of it before any package
    # is rebuilt, since kept columns are read off the partners' matrices
    masks = Dict{P, Tuple{BitVector, BitVector}}()
    A = UInt64[] # active class mask buffer
    for (p, info_p) in info
        X = info_p.conflicts
        m = nclasses(info_p)
        n = size(X, 2) - 1
        W = col_words(X)
        # active classes
        I = X[1:m, end]
        resize!(A, W)
        col_copy!(A, X, n + 1)
        clear_rows_above!(A, m)
        # kept columns: dependency columns some active class depends
        # on, and conflict columns of active partner classes — a
        # partner's column block must stay aligned with its class list,
        # so kept partner classes keep their columns even if empty
        K = falses(n)
        for k = 1:length(info_p.depends)
            K[k] = col_intersects(X, k, A)
        end
        for (q, b) in info_p.interacts
            Y = info[q].conflicts
            copy_col_bits!(K.chunks, b, Y, size(Y, 2), nclasses(info[q]))
        end
        masks[p] = (I, K)
    end
    # save original class counts (the matrices are rebuilt in place)
    N = Dict{P,Int}(p => nclasses(info_p) for (p, info_p) in info)
    # second pass: shrink each PkgInfo
    for p in collect(keys(info))
        info_p = info[p]
        # abbreviate components
        V = info_p.versions
        C = info_p.classes
        M = info_p.members
        D = info_p.depends
        T = info_p.interacts
        X = info_p.conflicts
        m = nclasses(info_p)
        n = size(X, 2) - 1
        W = col_words(X)
        I, K = masks[p]
        m′ = count(I)
        # delete if no active classes
        if m′ == 0
            delete!(info, p)
            reps === nothing || delete!(reps, p)
            continue
        end
        # keep as is if everything is active (with the column flags
        # normalized, as rebuilding would leave them)
        if m′ == m && all(K)
            X[m+1, 1:n] .= true
            continue
        end
        # the surviving classes keep the versions they hold: a deleted class
        # takes its members with it, since nothing is left to name them.
        # dropping classes and columns can only *merge* row-equality classes
        # (fewer rows to distinguish, fewer columns to differ in), so
        # restricting the partition stays sound, if possibly finer than the
        # truth; `pkg_info` is where it is computed exactly
        keep = falses(length(V))
        for i = 1:m
            I[i] || continue
            for j in M[i]
                keep[j] = true
            end
        end
        index = cumsum(keep) # old version index => new one (kept ones only)
        V′ = V[keep]
        C′ = Vector{Int}(undef, length(V′))
        M′ = Vector{Vector{Int}}(undef, m′)
        i′ = 0
        for i = 1:m
            I[i] || continue
            i′ += 1
            M′[i′] = mem = Int[index[j] for j in M[i]]
            for j in mem
                C′[j] = i′
            end
        end
        if reps !== nothing
            r = reps[p]
            reps[p] = Int[iszero(r[i]) ? 0 : index[r[i]] for i = 1:m if I[i]]
        end
        # compute shrunken components
        D′ = D[K[1:length(D)]]
        T′ = Dict{P,Int}()
        b′ = length(D′)
        for (q, b) in sort!(collect(T), by=last)
            n′ = count(K[b .+ (1:N[q])])
            if n′ > 0
                T′[q] = b′
                b′ += n′
            end
        end
        R′ = padded_rows(m′)
        W′ = R′ >> 6
        X′ = falses(R′, b′ + 1)
        src = X.chunks
        dst = X′.chunks
        # kept-class mask words for the gather below
        resize!(A, W)
        col_copy!(A, X, n + 1)
        clear_rows_above!(A, m)
        prefix = findlast(I) == m′ # kept classes are 1:m′
        j′ = 0
        for j = 1:n
            K[j] || continue
            j′ += 1
            sb = (j - 1) * W
            db = (j′ - 1) * W′
            if prefix
                # kept classes are a prefix: straight word copy
                nw = (m′ + 63) >> 6
                @inbounds for w = 1:nw
                    dst[db + w] = src[sb + w]
                end
                r = m′ & 63
                r != 0 && @inbounds (dst[db + nw] &= (UInt64(1) << r) - 1)
            else
                # gather the kept classes' bits
                acc = UInt64(0)
                na = 0
                dw = db + 1
                @inbounds for w = 1:W
                    avail = A[w]
                    iszero(avail) && continue
                    v = src[sb + w]
                    while !iszero(avail)
                        acc |= ((v >> trailing_zeros(avail)) & 1) << na
                        avail &= avail - 1
                        na += 1
                        if na == 64
                            dst[dw] = acc
                            dw += 1
                            acc = UInt64(0)
                            na = 0
                        end
                    end
                end
                na > 0 && @inbounds (dst[dw] = acc)
            end
        end
        @assert j′ == b′
        X′[1:m′, end] .= true  # kept classes are active
        X′[m′+1, 1:b′] .= true # kept columns are active
        # assign new struct into info
        info[p] = PkgInfo(V′, C′, M′, D′, T′, X′)
    end
    return info
end

function check_info_structure(
    info :: Dict{P, PkgInfo{P,V}},
) where {P,V}
    for (p, info_p) in info
        length(info_p.classes) == length(info_p.versions) ||
            return :classes, p, p
        for (i, mem) in enumerate(info_p.members), j in mem
            info_p.classes[j] == i ||
                return :members, p, p
        end
        size(info_p.conflicts, 1) == padded_rows(nclasses(info_p)) ||
            return :rows, p, p
        for q in info_p.depends
            q in keys(info) ||
                return :depends, p, q
        end
        interacts = sort!(collect(info_p.interacts), by=last)
        for (i, (q, b)) in enumerate(interacts)
            if i == 1
                b == length(info_p.depends) ||
                    return :conflicts, p, q
            end
            if i < length(interacts)
                b + nclasses(info[q]) == interacts[i+1][2] ||
                    return :conflicts, p, q
            end
        end
    end
end
