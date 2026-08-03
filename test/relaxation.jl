# Relaxation stability: the empirical check behind Theorem D1 and Proposition
# D1' (see the manual's Theory section, "Diagnosing unsatisfiability").
#
# Diagnosis asks a *family* of questions of one instance: what if this
# requirement were dropped, these constraints relaxed, this incompatibility
# lifted? The instance it asks them of was prepared -- collapsed to class
# representatives, then filtered -- with all of them active. Theorem D1 says
# that is sound, the preparation preserving the layered answer and not merely
# the verdict, provided every relaxable group is *column-closed*: switching it
# off clears whole columns of every package's constraint matrix, never part of
# one. Two granularities qualify, and they are the ones diagnosis uses: all of
# a package's constraints together -- compat bound, pin and every admission
# knob, which share one exclusion column -- and all conflicts between an
# interacting pair.
#
# Below column-closure it is false, and this file pins that down too: the
# one-package counterexample of Proposition D1', the same failure at
# declaration granularity for third-party bounds, and a sweep confirming that
# reachability alone -- the half Lemma R covers unconditionally -- is stable at
# the finer granularity where the rest is not.
#
# Everything here is brute force. For each instance the sweep enumerates EVERY
# relaxation sigma -- every subset of the requirements crossed with every
# subset of the groups -- and compares against a reference resolve of the raw
# data with those groups deleted. Not a sample: every sigma. Both `group`
# settings are swept, so the collapse is under test alongside the filter.

using Resolver: PkgData, PkgInfo, Problem, class_representatives,
    drop_unmarked!, exclusion_masks, is_excluded, mark_installable!,
    mark_reachable!, pkg_info, prepare_pkg_info, resolve_prepared,
    version_permutations

# resolve a universe that is *already* prepared, without preparing it again --
# the whole point being to interrogate the universe the failed resolve used
resolve_asis(info, prob; by = identity) =
    resolve_prepared(info, prob; by, diagnose = false)

## the group granularities

# column-closed: all conflicts between an unordered pair of packages
pair_groups(data::AbstractDict{P}) where {P} =
    sort!(collect(Set{Tuple{P,P}}(
        minmax(p, q)
        for (p, data_p) in data
        for v in data_p.versions if haskey(data_p.compat, v)
        for (q, s) in data_p.compat[v]
        if q ≠ p && haskey(data, q) && any(w ∉ s for w in data[q].versions))))

# NOT column-closed: one group per compat declaration -- the versions of `pkg`
# that share one statement about `dep`. This is the granularity a report would
# like to blame, and Proposition D1' is why it may not.
struct BoundGroup{P,V}
    pkg      :: P
    dep      :: P
    versions :: Vector{V}
end

function decl_groups(data::AbstractDict{P,<:PkgData{P,V}}) where {P,V}
    groups = BoundGroup{P,V}[]
    for p in sort!(collect(keys(data)))
        data_p = data[p]
        byq = Dict{P,Dict{Vector{V},Vector{V}}}() # dep => excluded => versions
        for v in data_p.versions
            haskey(data_p.compat, v) || continue
            for (q, s) in data_p.compat[v]
                (q == p || !haskey(data, q)) && continue
                excl = V[w for w in data[q].versions if w ∉ s]
                isempty(excl) && continue # excludes nothing: no group
                d = get!(() -> Dict{Vector{V},Vector{V}}(), byq, q)
                push!(get!(() -> V[], d, excl), v)
            end
        end
        for q in sort!(collect(keys(byq))), e in sort!(collect(keys(byq[q])))
            push!(groups, BoundGroup{P,V}(p, q, byq[q][e]))
        end
    end
    return groups
end

# `data` with the given conflict groups deleted -- the upstream releases that
# loosened those declarations. Rebuilt with plain containers, since the input
# may use specialized immutable ones (the tiny generators' TinyDict/TinyRange).
function relax_data(
    data :: AbstractDict{P,<:PkgData{P,V}},
    off,
) where {P,V}
    # (package, version) => the deps it no longer says anything about
    drop = Dict{Tuple{P,V},Set{P}}()
    pairs = Set{Tuple{P,P}}()
    for g in off
        if g isa BoundGroup
            for v in g.versions
                push!(get!(() -> Set{P}(), drop, (g.pkg, v)), g.dep)
            end
        else
            a, b = g
            push!(pairs, (a, b))
            push!(pairs, (b, a))
        end
    end
    Deps = Dict{V,Vector{P}}
    Comp = Dict{V,Dict{P,Vector{V}}}
    out = Dict{P,PkgData{P,V,Vector{V},Vector{V},Deps,Comp}}()
    for (p, data_p) in data
        vers = V[v for v in data_p.versions]
        deps = Deps(v => P[q for q in data_p.depends[v]]
                    for v in vers if haskey(data_p.depends, v))
        comp = Comp()
        for v in vers
            haskey(data_p.compat, v) || continue
            d = Dict{P,Vector{V}}()
            gone = get(drop, (p, v), nothing)
            for (q, s) in data_p.compat[v]
                haskey(data, q) || continue # compat on an unknown package: inert
                (p, q) ∈ pairs && continue
                gone === nothing || q ∉ gone || continue
                d[q] = V[w for w in data[q].versions if w ∈ s]
            end
            comp[v] = d
        end
        out[p] = PkgData(vers, deps, comp)
    end
    return out
end

# The layered answer is defined against a version *ranking*, so when the
# prepared universe is laid out in a non-canonical order the reference side has
# to rank versions the same way -- otherwise the two are answering different
# questions. `relax_data` hands back plain vectors, so sorting in place is safe.
function reorder!(data, order)
    order === nothing && return data
    for (p, data_p) in data
        sort!(data_p.versions; lt = order(p), alg = MergeSort)
    end
    return data
end

# do `data`'s compat declarations make p@v and q@w incompatible?
function raw_conflict(data, p, v, q, w)
    (haskey(data[p].compat, v) && haskey(data[p].compat[v], q) &&
        w ∉ data[p].compat[v][q]) ||
    (haskey(data[q].compat, w) && haskey(data[q].compat[w], p) &&
        v ∉ data[q].compat[w][p])
end

# `info` -- prepared with every group active -- with the conflict bits that
# `rdata` no longer implies cleared: the prepared universe *under* the
# relaxation, i.e. the left-hand side of Theorem D1. Dependency columns and the
# version list are untouched; the point is to ask the already-prepared instance
# the relaxed question, which is exactly what diagnosis does.
#
# `classes` is carried over unchanged: clearing conflict bits can only *merge*
# row-equality classes, so the old partition stays sound (possibly finer), which
# is all `PkgInfo` promises of it.
function relax_info(
    info  :: Dict{P,PkgInfo{P,V}},
    rdata :: AbstractDict{P,<:PkgData{P,V}},
) where {P,V}
    out = Dict{P,PkgInfo{P,V}}(
        p => PkgInfo{P,V}(copy(i.versions), copy(i.depends),
                          copy(i.interacts), copy(i.conflicts), copy(i.classes))
        for (p, i) in info)
    for (p, info_p) in out
        X = info_p.conflicts
        for (q, b) in info_p.interacts
            vers_q = out[q].versions
            for (i, v) in enumerate(info_p.versions),
                (j, w) in enumerate(vers_q)
                X[i, b+j] || continue
                raw_conflict(rdata, p, v, q, w) || (X[i, b+j] = false)
            end
        end
    end
    return out
end

## constraint groups, at the two granularities

# the packages `prob` forbids at least one existing version of. Compat bounds
# and pins name their packages; an admission knob reaches every package, so
# this asks the universe rather than the problem.
function constrained_pkgs(data::AbstractDict{P}, prob::Problem{P}) where {P}
    sort!(P[p for p in keys(data)
            if any(is_excluded(prob, p, v) for v in data[p].versions)])
end

# column-closed: all of a package's constraints together -- its compat bound,
# its pin, and every admission knob that touches it. They share one exclusion
# column, which is exactly why they have to move together.
user_pkg_groups(data::AbstractDict{P}, prob::Problem{P}) where {P} =
    constrained_pkgs(data, prob)

# NOT column-closed: one group per constraint *source* on a package
function user_source_groups(data::AbstractDict{P}, prob::Problem{P}) where {P}
    groups = Tuple{Symbol,P}[]
    for p in constrained_pkgs(data, prob)
        vers = data[p].versions
        haskey(prob.compat, p) && any(v ∉ prob.compat[p] for v in vers) &&
            push!(groups, (:compat, p))
        haskey(prob.pins, p) && any(v ≠ prob.pins[p] for v in vers) &&
            push!(groups, (:pin, p))
        for (kind, forbids) in prob.excludes
            any(forbids(p, v)::Bool for v in vers) && push!(groups, (kind, p))
        end
    end
    return groups
end

# `prob` with the given constraint groups switched off. A group named by a bare
# package switches off that package's whole exclusion column, admission knobs
# included -- which needs the knobs' predicates masked, since a knob is stated
# globally and cannot be dropped for one package by deleting an entry.
function relax_problem(prob::Problem{P}, reqs, off) where {P}
    compat = Dict(prob.compat)
    pins   = Dict(prob.pins)
    offpkg = Set{P}()               # whole exclusion columns switched off
    offsrc = Set{Tuple{Symbol,P}}() # single sources switched off
    for g in off
        g isa Tuple ? push!(offsrc, g) : push!(offpkg, g)
    end
    for p in offpkg
        delete!(compat, p)
        delete!(pins, p)
    end
    for (kind, p) in offsrc
        kind === :compat ? delete!(compat, p) :
        kind === :pin    ? delete!(pins, p)   : nothing
    end
    excludes = Pair{Symbol,Any}[
        kind => ((p, v) -> p ∉ offpkg && (kind, p) ∉ offsrc &&
                           forbids(p, v)::Bool)
        for (kind, forbids) in prob.excludes]
    Problem(reqs; compat, pins, excludes, labels = prob.labels)
end

subsets(v::AbstractVector) =
    (v[[isodd(b >> (i-1)) for i = 1:length(v)]] for b = 0:(1 << length(v)) - 1)

## the universes under test

# the production preparation: T1, then collapse, then filter -- everything on
prepared_info(data, prob; group = true, order = nothing) =
    prepare_pkg_info(pkg_info(data, prob.reqs), prob; group, order)

# arc consistency + reachability, no redundancy elimination and no collapse:
# the half of the preparation Lemma R covers for arbitrary group removal
function reach_only_info(data, prob; group = false, order = nothing)
    info = pkg_info(data, prob.reqs; filter = false)
    mark_installable!(info)
    drop_unmarked!(info)
    while true
        total = sum(length(i.versions) for i in values(info); init = 0)
        mark_reachable!(info, prob.reqs, exclusion_masks(info, prob))
        drop_unmarked!(info)
        sum(length(i.versions) for i in values(info); init = 0) < total || break
    end
    return info
end

# One instance, every sigma: the prepared instance under sigma must answer
# exactly as a brute-force resolve of the correspondingly deleted raw data.
# Returns the number of violations, so callers can assert zero or report a rate.
function sweep_relaxations(data, prob::Problem{P}, by::Function;
                           prepare = prepared_info,
                           group = true,
                           order = nothing,
                           ugroups = user_pkg_groups(data, prob),
                           bgroups = pair_groups(data),
                           quiet::Bool = false) where {P}
    info = prepare(data, prob; group, order)
    bad = 0
    for boff in subsets(bgroups)
        rdata = reorder!(relax_data(data, boff), order)
        rinfo = relax_info(info, rdata)
        for uoff in subsets(ugroups), reqs in subsets(prob.reqs)
            σprob = relax_problem(prob, reqs, uoff)
            got  = resolve_asis(rinfo, σprob; by)
            want = ref_resolve(bake(rdata, σprob), reqs; by)
            got == want && continue
            bad += 1
            quiet || @error "relaxation stability violated" reqs uoff boff got want
        end
    end
    return bad
end

sweep_size(reqs, ugroups, bgroups) =
    (1 << length(reqs)) * (1 << length(ugroups)) * (1 << length(bgroups))

# a random admission knob: forbids a random version subset of every package,
# the same shape of constraint as `:prerelease` and `:yanked`
random_kinds(n::Int) = rand(Bool) ? Pair{Symbol,Any}[] :
    let bad = Set{Int}(v for v = 1:n if rand(Bool))
        Pair{Symbol,Any}[:test => (p, v) -> v ∈ bad]
    end

@testset "relaxation stability: brute-force sweep" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    hi = p -> p  # default priority: lower package id first
    lo = p -> -p # reversed priority
    rev = p -> (u, v) -> u > v # a non-canonical version ordering

    # Theorem D1, exhaustive tier: every dependency and compatibility pattern
    # of the smallest grids, a random constraint draw on each, every sigma --
    # with and without the interchangeability collapse
    for m = 1:2, n = 1:2
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for deps_bits = 0:2^d-1
            deps = make_deps(deps_bits)
            for comp_bits = 0:2^c-1
                comp = make_comp(comp_bits)
                fill_data!(m, n, deps, comp, data)
                compat, pins = random_constraints(m, n)
                prob = Problem(collect(1:m);
                               compat, pins, excludes = random_kinds(n))
                for group in (true, false)
                    @test sweep_relaxations(data, prob, rand((hi, lo));
                                            group) == 0
                end
            end
        end
    end

    # Theorem D1, random tier: bigger grids, bounded sigma count per instance,
    # both collapse settings and a non-canonical version ordering
    for (m, n) in ((2, 3), (3, 2), (3, 3), (4, 2))
        (m*n)^2 ≤ 128 || continue
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        tested = 0
        for _ = 1:400
            tested < 20 || break
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            compat, pins = random_constraints(m, n)
            prob = Problem(collect(1:m);
                           compat, pins, excludes = random_kinds(n))
            ug = user_pkg_groups(data, prob)
            bg = pair_groups(data)
            sweep_size(prob.reqs, ug, bg) ≤ 2048 || continue
            tested += 1
            by = rand((hi, lo))
            for (group, order) in ((true, nothing), (false, nothing),
                                   (true, rev))
                @test sweep_relaxations(data, prob, by;
                                        group, order,
                                        ugroups = ug, bgroups = bg) == 0
            end
        end
        @test tested > 0
    end

    # Theorem D1, hand-built shapes the random draws under-sample: a constraint
    # that forbids everything, a pin at a missing version, an incompatibility a
    # relaxation must be able to bring a version back through
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    chain = Dict(
        :A => PkgData([:v2, :v1], Dict(:v2 => [:C], :v1 => [:C]),
                      Dict(:v2 => Dict(:C => [:v1]))),
        :B => PkgData([:v2, :v1], Dict(:v2 => [:C]),
                      Dict(:v2 => Dict(:C => [:v2]))),
        :C => PkgData([:v2, :v1], nodeps, nocomp),
    )
    pre = Pair{Symbol,Any}[:prerelease => (p, v) -> v === :v2]
    cases = [
        (chain, Problem([:A, :B])),
        (chain, Problem([:A, :B]; compat = Dict(:C => [:v1]))),
        (chain, Problem([:A, :B]; compat = Dict(:C => Symbol[]))),
        (chain, Problem([:A, :B]; pins = Dict(:C => :v9))),
        (chain, Problem([:A, :B]; pins = Dict(:A => :v2, :C => :v2))),
        (chain, Problem([:A, :B, :C];
                        compat = Dict(:A => [:v1]), pins = Dict(:C => :v1))),
        (chain, Problem([:A, :B]; excludes = pre)),
        (chain, Problem([:A, :B]; compat = Dict(:C => [:v2]), excludes = pre)),
        (chain, Problem([:A, :B]; pins = Dict(:C => :v2), excludes = pre)),
    ]
    byrev = p -> -Int(first(string(p)))
    for (data, prob) in cases, by in (identity, byrev), group in (true, false)
        @test sweep_relaxations(data, prob, by; group) == 0
    end
end

@testset "relaxation stability: the collapse keeps forbidden versions" begin
    # `class_representatives` refines the T1 interchangeability classes by the
    # exclusion mask *before* picking representatives, so a forbidden
    # sub-class keeps a representative of its own. That is what makes the
    # collapse compatible with diagnosis at all: a version the query forbids is
    # still in the instance, and switching the constraint off can bring it back.
    #
    # Checked for every kind of constraint source, since they all feed the one
    # mask -- which is also why they have to relax together.
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    # :v2 and :v1 are interchangeable: same (empty) dependency set, and nothing
    # anywhere constrains one and not the other
    data = Dict(
        :A => PkgData([:v1], Dict(:v1 => [:B]), nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    t1 = pkg_info(data, [:A])
    @test t1[:B].classes == [1, 1] # one class, as promised

    # forbid the *better* member, so that relaxing has something to show
    sources = [
        ("compat", Problem([:A]; compat = Dict(:B => [:v1]))),
        ("pin",    Problem([:A]; pins = Dict(:B => :v1))),
        ("prerelease",
         Problem([:A]; excludes = [:prerelease => (p, v) -> v === :v2])),
        ("yanked",
         Problem([:A]; excludes = [:yanked => (p, v) -> v === :v2])),
    ]
    for (name, prob) in sources
        keep = class_representatives(pkg_info(data, [:A]), prob)
        # the one T1 class splits in two, and both halves keep a member
        @test count(keep[:B]) == 2
        info = prepared_info(data, prob)
        @test :v2 ∈ info[:B].versions   # forbidden, but present
        @test resolve_asis(info, prob) == Dict(:A => :v1, :B => :v1)
        # ... and switching the source off reactivates it
        σ = relax_problem(prob, [:A], [:B])
        @test resolve_asis(info, σ) == Dict(:A => :v1, :B => :v2)
        @test resolve_asis(info, σ) == bare_resolve(data, [:A])
    end

    # with no constraint at all the class collapses to one representative
    keep = class_representatives(pkg_info(data, [:A]), Problem([:A]))
    @test count(keep[:B]) == 1
end

@testset "relaxation stability: below column-closure it fails" begin
    # Proposition D1': one package, two versions, a compat bound admitting only
    # the worse one and a pin holding the better one. Together they forbid
    # everything, so the two versions have identical constraint rows and
    # redundancy deletes the second -- and relaxing just the *pin* then finds
    # nothing left to satisfy the compat with.
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    data = Dict(:A => PkgData([:v1, :v2], nodeps, nocomp))
    full = Problem([:A]; compat = Dict(:A => [:v2]), pins = Dict(:A => :v1))
    info = prepared_info(data, full)
    @test info[:A].versions == [:v1] # :v2 collapsed into :v1
    σprob = Problem([:A]; compat = Dict(:A => [:v2])) # pin dropped
    @test bare_resolve(data, σprob) == Dict(:A => :v2) # raw: satisfiable
    @test resolve_asis(info, σprob) === nothing        # prepared: not
    # and the whole-package relaxation, which *is* column-closed, is fine
    @test sweep_relaxations(data, full, identity) == 0

    # the same, with an admission knob standing in for the compat bound: the
    # kinds are not special, they are just more sources sharing the one column
    full2 = Problem([:A]; pins = Dict(:A => :v1),
                          excludes = [:prerelease => (p, v) -> v === :v1])
    info2 = prepared_info(data, full2)
    σ2 = Problem([:A]; excludes = [:prerelease => (p, v) -> v === :v1])
    @test bare_resolve(data, σ2) == Dict(:A => :v2)
    @test resolve_asis(info2, σ2) === nothing
    @test sweep_relaxations(data, full2, identity) == 0

    # and at declaration granularity for third-party bounds: two separate
    # compat statements about the same package, only one relaxed
    data = Dict(
        :A => PkgData([:v1, :v2], Dict(:v2 => [:B]), nocomp),
        :B => PkgData([:v1, :v2], nodeps,
                      Dict(:v1 => Dict(:A => [:v1]), :v2 => Dict(:A => Symbol[]))),
    )
    prob = Problem([:A, :B]; compat = Dict(:B => [:v1, :v2]),
                             pins = Dict(:B => :v1, :A => :v2))
    decls = decl_groups(data)
    @test length(decls) == 2 # B@v1's statement and B@v2's, coalesced apart
    # column-closed granularity: sound
    @test sweep_relaxations(data, prob, identity) == 0
    # declaration granularity: not (this is a regression test for the *fact*,
    # so that a future change which makes it sound gets noticed)
    @test sweep_relaxations(data, prob, identity;
                            bgroups = decls, quiet = true) > 0
end

@testset "relaxation stability: reachability alone is unconditionally stable" begin
    # Lemma R holds for arbitrary group removal, so a reach-only universe --
    # no redundancy elimination, no collapse -- is stable even at the finer
    # granularities Proposition D1' breaks, which localizes the whole
    # instability in the answer-preserving-by-domination half of preparation
    Random.seed!(rand(RandomDevice(), UInt64))
    for (m, n) in ((2, 2), (2, 3), (3, 2))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        tested = 0
        for _ = 1:400
            tested < 20 || break
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            compat, pins = random_constraints(m, n)
            prob = Problem(collect(1:m);
                           compat, pins, excludes = random_kinds(n))
            ug = user_source_groups(data, prob)
            bg = decl_groups(data)
            sweep_size(prob.reqs, ug, bg) ≤ 4096 || continue
            tested += 1
            @test sweep_relaxations(data, prob, rand((identity, p -> -p));
                                    prepare = reach_only_info,
                                    ugroups = ug, bgroups = bg) == 0
        end
        @test tested > 0
    end
end
