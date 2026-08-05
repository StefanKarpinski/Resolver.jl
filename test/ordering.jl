# The version ordering as a resolve parameter.
#
# T1 (`pkg_info`) is built on whatever order the data comes in — the *canonical*
# order — and the preference ordering a query wants is the `order` argument to
# `resolve`. The oracle throughout is the OLD path, where the ordering was baked
# into the provider: sorting each package's version vector by a comparator and
# resolving canonically must give exactly what passing that comparator gives.
#
# That is a strictly stronger statement than the ordering-independence of T1
# (tested in classes.jl by permuting the data): here the *same* artifact serves
# every ordering, so the reordering has to happen inside `prepare_pkg_info`,
# down to the partner blocks of every conflicts matrix.

using Resolver: PkgData, PkgInfo, DepsProvider, Problem, pkg_info,
    prepare_pkg_info, version_permutations, class_ranking

# a comparator over the tiny data's integer versions, from a ranking: `ranks[v]`
# is how good v is, lower being better
ranked(ranks::AbstractVector{Int}) = (u::Int, v::Int) -> ranks[u] < ranks[v]

# `order` values: a callable package -> comparator. `nothing` is the canonical
# order, which is what the artifact carries.
const CANONICAL = nothing
reversed(m::Int, n::Int) = p -> (u::Int, v::Int) -> u > v
# one random ranking per package, so the partner blocks of each matrix are
# permuted differently from the rows
function shuffled(m::Int, n::Int)
    ranks = Dict(p => randperm(n) for p = 1:m)
    p -> ranked(ranks[p])
end
# only *some* packages reordered: exercises the mixed fast/slow paths
function partly(m::Int, n::Int)
    ranks = Dict(p => (isodd(p) ? collect(1:n) : reverse(1:n)) for p = 1:m)
    p -> ranked(collect(ranks[p]))
end

# the old baked path: the provider hands over versions already in the order the
# query wants (the version-keyed deps & compat dicts carry over untouched)
function bake_order(data::AbstractDict{P}, order) where {P}
    order === nothing && return data
    Dict(p => PkgData(sort!(collect(data_p.versions); lt = order(p)),
                      data_p.depends, data_p.compat)
         for (p, data_p) in data)
end

# passing a comparator == baking the same order into the data
function test_order_equivalence(data, prob, order; by::Function = identity)
    baked = bake_order(data, order)
    @test resolve(data, prob; by, order, diagnose = false) ==
          resolve(baked, prob; by, diagnose = false)
end

@testset "ordering: permutations and identity" begin
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    data = Dict(
        :A => PkgData([:v3, :v2, :v1], nodeps, nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    info = pkg_info(data)

    # the canonical order needs no permutations at all, and neither does a
    # comparator that agrees with it: the artifact is already laid out that way
    @test version_permutations(info, nothing) === nothing
    down = p -> (u, v) -> u > v # :v3 > :v2 > :v1, i.e. the canonical order here
    @test version_permutations(info, down) === nothing

    # a comparator that disagrees gets an entry per package it reorders
    up = p -> (u, v) -> u < v
    perms = version_permutations(info, up)
    @test perms == Dict(:A => [3, 2, 1], :B => [2, 1])

    # ... and only for the packages it reorders
    mixed = p -> p == :A ? (u, v) -> u < v : (u, v) -> u > v
    @test version_permutations(info, mixed) == Dict(:A => [3, 2, 1])

    # the representative of a class is its best member in the *active* order
    data = Dict(
        :A => PkgData([:v3, :v2, :v1], nodeps, Dict(:v2 => Dict(:B => [:v2]))),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    info = pkg_info(data, keys(data); filter = false)
    @test info[:A].classes == [1, 2, 1] # :v3 and :v1 are interchangeable
    prob = Problem([:A])
    reps, cperms = class_ranking(info, prob)
    @test reps[:A] == [1, 2] # class 1 stands at :v3, class 2 at :v2
    @test cperms === nothing # ... which is the order the artifact has
    up_perms = version_permutations(info, up)
    reps, cperms = class_ranking(info, prob, up_perms)
    @test reps[:A] == [3, 2]     # reversed: :v1 now represents the class
    @test !haskey(cperms, :A)    # ... and it still outranks :v2
    @test cperms[:B] == [2, 1]   # :B's own two classes do swap, though
    # a constraint on the best member moves the class it is in: with :v3 out,
    # class 1 competes at :v1 and falls behind the class holding :v2
    prob = Problem([:A]; compat = Dict(:A => [:v2, :v1]))
    reps, cperms = class_ranking(info, prob)
    @test reps[:A] == [3, 2]
    @test cperms[:A] == [2, 1]

    # and the answer follows the order: ascending picks the worst versions
    @test resolve(data, [:A, :B]) == Dict(:A => :v3, :B => :v2)
    @test resolve(data, [:A, :B]; order = up) == Dict(:A => :v1, :B => :v1)
end

@testset "ordering: comparator vs. baked, complete data grids" begin
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
                    for mk in (reversed, shuffled, partly), by in (hi, lo)
                        order = mk(m, n)
                        test_order_equivalence(data, Problem(reqs), order; by)
                        compat, pins = random_constraints(m, n)
                        test_order_equivalence(
                            data, Problem(reqs; compat, pins), order; by)
                    end
                end
            end
        end
    end
end

@testset "ordering: comparator vs. baked, random grids" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    hi = p -> p
    lo = p -> -p
    for (m, n) in ((2, 3), (3, 2), (3, 3), (2, 4), (4, 2), (2, 5), (5, 2))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:60
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            for mk in (reversed, shuffled, partly)
                order = mk(m, n)
                by = rand((hi, lo))
                test_order_equivalence(data, Problem(reqs), order; by)
                compat, pins = random_constraints(m, n)
                test_order_equivalence(
                    data, Problem(reqs; compat, pins), order; by)
                # ... and the ordering must not disturb the bake-equivalence of
                # the user constraints either
                prob = Problem(reqs; compat, pins)
                @test resolve(data, prob; by, order, diagnose = false) ==
                      resolve(bake(bake_order(data, order), prob), reqs;
                              by, diagnose = false)
            end
        end
    end
end

@testset "ordering: comparator vs. baked, exhaustive constraint shapes" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    hi = p -> p
    lo = p -> -p
    for (m, n) in ((2, 2), (2, 3), (3, 2))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:2
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            order = shuffled(m, n)
            for shapes in Iterators.product(ntuple(_ -> CONSTRAINT_SHAPES, m)...)
                compat = Dict{Int,Vector{Int}}()
                pins   = Dict{Int,Int}()
                for (p, shape) in enumerate(shapes)
                    constrain!(compat, pins, p, shape, n)
                end
                prob = Problem(reqs; compat, pins)
                test_order_equivalence(data, prob, order; by = hi)
                test_order_equivalence(data, prob, order; by = lo)
            end
        end
    end
end

@testset "ordering: one T1 artifact serves every ordering" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    # The universality claim, at the level the cache will be used: build the
    # all-requirements artifact once, then resolve it under several orderings
    # and several constraint sets, cross-checking each against a from-scratch
    # resolve of the data baked that way — and confirm the artifact survives.
    for (m, n) in ((3, 2), (3, 3), (4, 2), (2, 4))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:20
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            deps = DepsProvider(p -> data[p], collect(1:m))
            all_reqs = pkg_info(deps) # T1 over every package
            before = deepcopy(all_reqs)
            for reqs_bits = 1:2^m-1
                reqs = collect(make_reqs(reqs_bits))
                compat, pins = random_constraints(m, n)
                for prob in (Problem(reqs), Problem(reqs; compat, pins)),
                    order in (CANONICAL, reversed(m, n), shuffled(m, n))
                    baked = bake_order(data, order)
                    @test resolve(all_reqs, prob; order, diagnose = false) ==
                          resolve(baked, prob; diagnose = false)
                    # ... and the same artifact answers the next query too
                    @test resolve(all_reqs, prob; order, diagnose = false) ==
                          resolve(all_reqs, prob; order, diagnose = false)
                end
            end
            @test all_reqs == before
        end
    end
end
