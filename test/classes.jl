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
#   * for the *answer*, the brute force. Classes are the representation now, so
#     there is no ungrouped resolve to compare against: what the collapse has to
#     agree with is the universe itself. Every constrained resolve is checked
#     against the resolve of the *baked* data — the registry with the forbidden
#     versions deleted, which is what a constraint says — and, where the
#     enumeration is affordable, against `ref_resolve` on that data, which knows
#     nothing about matrices, classes or representatives at all.

using Resolver: PkgData, PkgInfo, DepsProvider, Problem, pkg_info, pkg_data,
    version_classes, class_ranking, prepare_pkg_info, rank_pkg_info,
    is_excluded, nclasses
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

# The answer is the universe's answer, whatever the constraints or the
# ordering. Two oracles, neither of which knows how the universe is
# represented: the same problem with the forbidden versions deleted from the
# data outright (a *different* registry, hence a different partition, resolved
# with no deactivation anywhere), and — when the potential solutions are few
# enough to enumerate — the brute force on that data.
const REF⁺ = 256 # enumeration budget for the brute-force cross-check

function test_class_space(data, prob; by::Function = identity)
    baked = bake(data, prob)
    sol = resolve(data, prob; by, diagnose = false)
    @test sol == resolve(baked, prob.reqs; by, diagnose = false)
    @test issatisfiable(data, prob) == (sol !== nothing)
    Π = prod(init = 1.0, float(length(d.versions) + 1) for d in values(baked))
    Π ≤ REF⁺ && @test sol == ref_resolve(baked, prob.reqs; by)
end

# the constraints that bear on multi-member classes: forbid exactly one member
# of one, or pin one member of it. Neither can split the class — that is the
# point — so each is a test that a class survives through the members left
function class_emptying_problems(data, reqs)
    info = pkg_info(data, keys(data); filter = false)
    probs = []
    for p in sort!(collect(keys(data)))
        vers = collect(data[p].versions)
        for cls in class_sets(info, p)
            length(cls) > 1 || continue
            for v in sort!(collect(cls))
                push!(probs, Problem(reqs;
                    compat = Dict(p => [w for w in vers if w != v])))
                push!(probs, Problem(reqs; pin = Dict(p => v)))
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

    # ... and each class stands at its best member
    reps, cperms = class_ranking(info)
    @test reps[:A] == [1, 2] # :v3 for {:v3, :v1}, :v2 for {:v2}
    @test cperms === nothing # which is the order the artifact already has

    # a pin inside a class does not *refine* it — a class is one column, and a
    # user constraint has no way to split one. All it can do is move which
    # member stands for the class, and empty the classes it admits nothing of
    prob = Problem([:A]; pin = Dict(:A => :v1))
    reps, cperms = class_ranking(info, prob)
    @test info[:A].members == [[1, 3], [2]] # the partition is untouched
    @test reps[:A] == [3, 0]  # {:v3, :v1} stands at :v1; {:v2} is empty
    @test cperms[:A] == [2, 1] # and the empty class ranks where :v2 would
    @test resolve(data, prob) == Dict(:A => :v1)
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

@testset "resolve vs. the oracles: complete data grids" begin
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
                        test_class_space(data, Problem(reqs); by)
                    end
                    for _ = 1:4
                        compat, pin = random_constraints(m, n)
                        test_class_space(data, Problem(reqs; compat, pin))
                    end
                end
            end
        end
    end
end

@testset "resolve vs. the oracles: random grids" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    hi = p -> p
    lo = p -> -p
    for (m, n) in ((2, 4), (4, 2), (2, 5), (5, 2), (3, 3))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:120
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            for _ = 1:4
                compat, pin = random_constraints(m, n)
                test_class_space(data, Problem(reqs; compat, pin);
                                 by = rand((hi, lo)))
            end
        end
    end
end

@testset "resolve vs. the oracles: constraints that empty classes" begin
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
            probs = class_emptying_problems(data, reqs)
            splits += length(probs)
            for prob in probs, by in (hi, lo)
                test_class_space(data, prob; by)
            end
        end
    end
    @test splits > 1000 # the sweep really did split classes
end

@testset "resolve vs. the oracles: exhaustive constraint shapes" begin
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
                pin    = Dict{Int,Int}()
                for (p, shape) in enumerate(shapes)
                    constrain!(compat, pin, p, shape, n)
                end
                prob = Problem(reqs; compat, pin)
                test_class_space(data, prob; by = hi)
                test_class_space(data, prob; by = lo)
            end
        end
    end
end

@testset "resolve vs. the oracles: adversarial" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    # break the solution until it is unsolvable, with constraints in force the
    # whole way, checking the answer against the oracles at every step
    for m = 2:4, n = 2:3
        (m*n)^2 ≤ 128 || continue
        make_deps, make_comp, data, d, c, bit = tiny_data_makers(m, n)
        for _ = 1:30
            deps = make_deps(0)
            comp = make_comp(0)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            compat, pin = random_constraints(m, n; mild = true)
            prob = Problem(reqs; compat, pin)
            while true
                fill_data!(m, n, deps, comp, data)
                test_class_space(data, prob)
                sol = resolve(data, prob; diagnose = false)
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
                # and the answer under the permuted ordering is the
                # universe's answer all the same
                test_class_space(perm, Problem(reqs))
                compat, pin = random_constraints(m, n)
                test_class_space(perm, Problem(reqs; compat, pin))
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
                compat, pin = random_constraints(m, n)
                for prob in (Problem(reqs), Problem(reqs; compat, pin))
                    @test resolve(all_reqs, prob; diagnose = false) ==
                          resolve(specific, prob; diagnose = false)
                    @test resolve(all_reqs, prob; diagnose = false) ==
                          resolve(data, prob; diagnose = false)
                    # the T1 artifacts are reusable: resolving does not
                    # consume them
                    @test resolve(all_reqs, prob; diagnose = false) ==
                          resolve(all_reqs, prob; diagnose = false)
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
    @test work.info !== info
    @test all(info[p].members == before[p].members for p in keys(info))
end
