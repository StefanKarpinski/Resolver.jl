# T1 preprocessing: interchangeability classes and the collapse to class
# representatives.
#
# Two oracles run throughout:
#
#   * for the *classes*, a reference partition computed straight off the raw
#     `PkgData` from the lemma's definition (same dependency set; every
#     compatibility entry in the problem, the package's own and every other
#     package's, constrains the two versions equally) — a completely different
#     computation from the bit-matrix row scan it checks;
#
#   * for the *collapse*, the ungrouped resolve. Grouping is an optimization,
#     so `resolve(…; group = true)` and `resolve(…; group = false)` must return
#     the very same answer — and the very same `nothing` — everywhere,
#     including where user constraints split classes.

using Resolver: PkgData, PkgInfo, DepsProvider, Problem, pkg_info, pkg_data,
    version_classes, class_representatives, prepare_pkg_info, is_excluded
using .TestResolver: each_potential_solution, is_valid_solution

# the class partition of package p, as a set of sets of version *values*
function class_sets(info::AbstractDict{P}, p::P) where {P}
    info_p = info[p]
    cls = info_p.classes
    vers = info_p.versions
    Set(Set(vers[i] for i in eachindex(cls) if cls[i] == c)
        for c in unique(cls))
end

# the test data's dictionaries are minimal AbstractDicts without a three-arg
# `get`, so look entries up the long way
lookup(d, k, default) = haskey(d, k) ? d[k] : default

# does version v of p forbid version w of q? (unknown packages and a package's
# entries about itself are no constraints, exactly as `pkg_info` reads them)
function ref_forbids(data, p, v, q, w)
    p == q && return false
    haskey(data, q) || return false
    c = lookup(data[p].compat, v, nothing)
    c === nothing && return false
    s = lookup(c, q, nothing)
    s === nothing && return false
    return w ∉ s
end

# the lemma's interchangeability relation, read off the raw data: same
# dependency set, and every compatibility entry in the problem constrains the
# two versions equally — p's own entries *and* every other package's entries
# that mention p, which is the half that keeps genuinely distinguishable
# versions apart
function ref_class_sets(data::AbstractDict{P}, p::P) where {P}
    others = sort!(P[q for q in keys(data) if q != p])
    depset(v) = Set(q for q in lookup(data[p].depends, v, ()) if q != p)
    conflicts(v) = [(q, w, ref_forbids(data, p, v, q, w) ||
                            ref_forbids(data, q, w, p, v))
                    for q in others for w in data[q].versions]
    key(v) = (depset(v), conflicts(v))
    groups = Dict{Any,Set}()
    for v in data[p].versions
        push!(get!(() -> Set{Any}(), groups, key(v)), v)
    end
    return Set(values(groups))
end

# the classes must be exactly the reference partition, for every package
function test_classes(data)
    info = pkg_info(data, keys(data); filter = false)
    for p in keys(data)
        @test class_sets(info, p) == ref_class_sets(data, p)
        # no class ever spans different dependency sets
        for cls in class_sets(info, p)
            @test length(Set(Set(lookup(data[p].depends, v, ())) for v in cls)) == 1
        end
    end
    return info
end

# swapping one member of a class for another preserves the validity of any
# solution — the soundness half of the lemma, checked over *every* potential
# solution rather than only the valid ones (a strictly stronger statement).
#
# The converse is deliberately not asserted: it is not a theorem. Two versions
# with different rows can still be indistinguishable in practice if every
# assignment that would tell them apart is invalid for an unrelated reason
# (say, their shared dependency conflicts with the partner whose bound differs).
# The reference partition above is what pins the class test down from the other
# side.
function test_class_swaps(data)
    info = pkg_info(data, keys(data); filter = false)
    pkgs = sort!(collect(keys(data)))
    V = eltype(info[first(pkgs)].versions)
    bad = nothing # first counterexample, if the lemma ever fails
    for (i, p) in enumerate(pkgs)
        for cls in class_sets(info, p)
            length(cls) > 1 || continue
            members = sort!(collect(cls))
            each_potential_solution(data, pkgs) do s
                s[i] === nothing && return
                t = Vector{Union{V,Nothing}}(s)
                t[i] = first(members)
                base = is_valid_solution(data, pkgs, t)
                for v in members
                    t[i] = v
                    is_valid_solution(data, pkgs, t) == base && continue
                    bad === nothing && (bad = (p, first(members), v, copy(t)))
                end
            end
        end
    end
    @test bad === nothing
end

# grouping cannot change the answer, whatever the constraints or the ordering
function test_collapse_invariance(data, prob; by::Function = identity)
    @test resolve(data, prob; by, group = true) ==
          resolve(data, prob; by, group = false)
end

# the constraints that split classes: forbid exactly one member of a
# multi-member class, or pin one member of it
function class_splitting_problems(data, reqs)
    info = pkg_info(data, keys(data); filter = false)
    probs = []
    for p in sort!(collect(keys(data)))
        vers = collect(data[p].versions)
        for cls in class_sets(info, p)
            length(cls) > 1 || continue
            for v in sort!(collect(cls))
                push!(probs, Problem(reqs;
                    compat = Dict(p => [w for w in vers if w != v])))
                push!(probs, Problem(reqs; pins = Dict(p => v)))
            end
        end
    end
    return probs
end

@testset "classes: the interchangeability test" begin
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()

    # no constraints at all: every version of every package is interchangeable
    data = Dict(
        :A => PkgData([:v3, :v2, :v1], nodeps, nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    info = test_classes(data)
    @test info[:A].classes == [1, 1, 1]
    @test info[:B].classes == [1, 1]
    @test resolve(data, [:A, :B]) == Dict(:A => :v3, :B => :v2)

    # differing dependency sets split
    data = Dict(
        :A => PkgData([:v3, :v2, :v1], Dict(:v3 => [:B], :v2 => [:B]), nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    info = test_classes(data)
    @test info[:A].classes == [1, 1, 2]

    # an *incoming* bound splits: nothing about B's own rows differs, but A@v1
    # admits only B@v2, so B's two versions are genuinely distinguishable
    data = Dict(
        :A => PkgData([:v1], Dict(:v1 => [:B]), Dict(:v1 => Dict(:B => [:v2]))),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    info = test_classes(data)
    @test info[:B].classes == [1, 2]

    # classes need not be contiguous runs: A@v3 and A@v1 share a row, A@v2
    # does not, and the merge pass has to find them
    data = Dict(
        :A => PkgData([:v3, :v2, :v1], nodeps, Dict(:v2 => Dict(:B => [:v2]))),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    info = test_classes(data)
    @test info[:A].classes == [1, 2, 1]

    # ... and the collapse keeps the best member of each of them
    keep = class_representatives(info)
    @test keep[:A] == BitVector([true, true, false])

    # a pin inside a class refines it: the pinned version becomes a class of
    # its own, so the collapse cannot drop it
    prob = Problem([:A]; pins = Dict(:A => :v1))
    @test class_representatives(info, prob)[:A] == BitVector([true, true, true])
end

@testset "classes: reference partition, complete data grids" begin
    # every dependency & compatibility pattern of the smallest grids
    for m = 1:2, n = 1:2
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for deps_bits = 0:2^d-1
            deps = make_deps(deps_bits)
            for comp_bits = 0:2^c-1
                comp = make_comp(comp_bits)
                fill_data!(m, n, deps, comp, data)
                test_classes(data)
                test_class_swaps(data)
            end
        end
    end
end

@testset "classes: reference partition, random grids" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    for (m, n) in ((2, 3), (3, 2), (2, 4), (3, 3), (4, 2), (2, 5), (5, 2))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:150
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            test_classes(data)
            test_class_swaps(data)
        end
    end
end

@testset "collapse invariance: complete data grids" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    hi = p -> p
    lo = p -> -p
    for m = 1:2, n = 1:2
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for deps_bits = 0:2^d-1
            deps = make_deps(deps_bits)
            for comp_bits = 0:2^c-1
                comp = make_comp(comp_bits)
                fill_data!(m, n, deps, comp, data)
                for reqs_bits = 0:2^m-1
                    reqs = collect(make_reqs(reqs_bits))
                    for by in (hi, lo)
                        test_collapse_invariance(data, Problem(reqs); by)
                    end
                    for _ = 1:4
                        compat, pins = random_constraints(m, n)
                        test_collapse_invariance(data, Problem(reqs; compat, pins))
                    end
                end
            end
        end
    end
end

@testset "collapse invariance: random grids" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    hi = p -> p
    lo = p -> -p
    for (m, n) in ((2, 4), (4, 2), (2, 5), (5, 2), (3, 3))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:120
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            for _ = 1:4
                compat, pins = random_constraints(m, n)
                test_collapse_invariance(data, Problem(reqs; compat, pins);
                                         by = rand((hi, lo)))
            end
        end
    end
end

@testset "collapse invariance: constraints that split classes" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    hi = p -> p
    lo = p -> -p
    splits = 0
    for (m, n) in ((2, 3), (3, 2), (3, 3), (2, 4), (4, 2), (2, 5), (5, 2))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:60
            # sparse patterns, so that classes have several members to split
            deps = make_deps(randbits(d) & randbits(d))
            comp = make_comp(randbits(c) & randbits(c) & randbits(c))
            fill_data!(m, n, deps, comp, data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            probs = class_splitting_problems(data, reqs)
            splits += length(probs)
            for prob in probs, by in (hi, lo)
                test_collapse_invariance(data, prob; by)
                # ... and the grouped answer is also the answer the same
                # universe with those versions deleted would give
                test_bake_equivalence(data, prob; by)
            end
        end
    end
    @test splits > 1000 # the sweep really did split classes
end

@testset "collapse invariance: exhaustive constraint shapes" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    hi = p -> p
    lo = p -> -p
    for (m, n) in ((2, 2), (2, 3), (3, 2))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:3
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            for shapes in Iterators.product(ntuple(_ -> CONSTRAINT_SHAPES, m)...)
                compat = Dict{Int,Vector{Int}}()
                pins   = Dict{Int,Int}()
                for (p, shape) in enumerate(shapes)
                    constrain!(compat, pins, p, shape, n)
                end
                prob = Problem(reqs; compat, pins)
                test_collapse_invariance(data, prob; by = hi)
                test_collapse_invariance(data, prob; by = lo)
            end
        end
    end
end

@testset "collapse invariance: adversarial" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    # break the solution until it is unsolvable, with constraints in force the
    # whole way, checking grouping against no grouping at every step
    for m = 2:4, n = 2:3
        (m*n)^2 ≤ 128 || continue
        make_deps, make_comp, data, d, c, bit = tiny_data_makers(m, n)
        for _ = 1:30
            deps = make_deps(0)
            comp = make_comp(0)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            compat, pins = random_constraints(m, n; mild = true)
            prob = Problem(reqs; compat, pins)
            while true
                fill_data!(m, n, deps, comp, data)
                test_collapse_invariance(data, prob)
                sol = resolve(data, prob)
                sol === nothing && break
                p = rand(collect(keys(sol)))
                v = sol[p]
                q = rand(1:m-1)
                q += q ≥ p
                w = get(sol, q, nothing)
                if isnothing(w)
                    b = bit(p, v, q)
                    @assert iszero(deps.bits & b)
                    deps = typeof(deps)(deps.bits | b)
                else
                    b = bit(p, v, q, w)
                    @assert iszero(comp.bits & b)
                    comp = typeof(comp)(comp.bits | b)
                end
            end
        end
    end
end

# reorder each package's version vector; the dependency and compatibility
# dictionaries are keyed by version *value*, so they carry over untouched
permute_versions(data::AbstractDict, perms) = Dict(
    p => PkgData(perms[p], data_p.depends, data_p.compat)
    for (p, data_p) in data)

@testset "classes: T1 is independent of the version ordering" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    for (m, n) in ((2, 3), (3, 3), (3, 2), (2, 4), (4, 2))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:20
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            base = pkg_info(data, keys(data); filter = false)
            for _ = 1:3
                perms = Dict(p => randperm(n) for p = 1:m)
                perm = permute_versions(data, perms)
                info = pkg_info(perm, keys(perm); filter = false)
                # the classes are sets of version values, unchanged
                for p = 1:m
                    @test class_sets(info, p) == class_sets(base, p)
                end
                # and the answer under the permuted ordering is the same
                # whether or not it is reached by collapsing
                test_collapse_invariance(perm, Problem(reqs))
                compat, pins = random_constraints(m, n)
                test_collapse_invariance(perm, Problem(reqs; compat, pins))
            end
        end
    end
end

@testset "T1: all requirements vs. requirement-specific" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    # `pkg_info(deps)` — every package required — is the persistable artifact.
    # Resolving any requirement set against it must give what resolving against
    # the closure of those requirements alone gives: the per-resolve filter is
    # seeded by the actual requirements and does the restriction.
    for (m, n) in ((3, 2), (3, 3), (4, 2), (2, 4))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:40
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            deps = DepsProvider(p -> data[p], collect(1:m))
            all_reqs = pkg_info(deps)          # T1 over every package
            for reqs_bits = 0:2^m-1
                reqs = collect(make_reqs(reqs_bits))
                specific = pkg_info(deps, reqs) # T1 over the closure of reqs
                compat, pins = random_constraints(m, n)
                for prob in (Problem(reqs), Problem(reqs; compat, pins))
                    @test resolve(all_reqs, prob) == resolve(specific, prob)
                    @test resolve(all_reqs, prob) == resolve(data, prob)
                    # the T1 artifacts are reusable: resolving does not
                    # consume them
                    @test resolve(all_reqs, prob) == resolve(all_reqs, prob)
                end
            end
        end
    end
end

@testset "T1: the artifact is left alone by resolving" begin
    # a caller-supplied info may be a cache: `resolve` must not shrink it
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    data = Dict(
        :A => PkgData([:v3, :v2, :v1],
                      Dict(:v3 => [:B], :v2 => [:B], :v1 => [:B]),
                      Dict(:v3 => Dict(:B => [:v2]))),
        :B => PkgData([:v2, :v1], Dict{Symbol,Vector{Symbol}}(), nocomp),
        :C => PkgData([:v1], Dict(:v1 => [:A]), nocomp),
    )
    info = pkg_info(data)
    before = deepcopy(info)
    prob = Problem([:C]; compat = Dict(:B => [:v1]))
    @test resolve(info, prob) == resolve(data, prob)
    @test info == before
    @test all(info[p].classes == before[p].classes for p in keys(info))
    # ... and preparing it explicitly yields a universe of its own
    work = prepare_pkg_info(info, prob)
    @test info == before
    @test work !== info
end
