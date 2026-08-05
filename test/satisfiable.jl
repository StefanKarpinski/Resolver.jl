# `issatisfiable`: the satisfiability verdict without the descent.
#
# The contract is agreement with `resolve`: `issatisfiable(args...)` is true
# exactly when `resolve(args...)` returns a solution. That equality is checked
# here across the entry-point family and on the edge cases, and — since it has
# to hold for every instance, not just the hand-built ones — inside
# `test_resolver` (test/setup.jl) and `test_bake_equivalence` (test/problem.jl),
# so every sweep in the suite checks it too.

using Resolver: SAT, PicoSAT, DepsProvider, PkgData, pkg_info, finalize,
    prepare_pkg_info

# The number of solves an instance has run — or `nothing` when picosat's report
# didn't reach us (see below). Picosat counts its own `picosat_sat` calls and
# reports the total in its statistics dump, which it writes to its output stream
# — stdout unless told otherwise — so the count is available without a counter
# of our own anywhere in the solving path. Two wrinkles: picosat prints the
# "N calls" line only once there have been *several* calls, so an absent line
# means at most one, which the propagation count separates from none at all (a
# solve of an instance that has clauses propagates); and reading any of it means
# capturing what the C library prints, which may not work everywhere.
function solve_count(sat::SAT)
    path, io = mktemp()
    try
        redirect_stdout(io) do
            ccall((:picosat_stats, PicoSAT.lib), Cvoid, (Ptr{Cvoid},), sat.pico)
            # picosat prints through the C library's buffered stdout, so flush
            # it before the redirect is undone
            ccall(:fflush, Cint, (Ptr{Cvoid},), C_NULL)
        end
        close(io)
        stats = read(path, String)
        props = match(r"(\d+) propagations", stats)
        props === nothing && return nothing # the dump didn't reach us
        calls = match(r"(\d+) calls", stats)
        calls === nothing && return parse(Int, props.captures[1]) > 0 ? 1 : 0
        return parse(Int, calls.captures[1])
    finally
        close(io)
        rm(path; force = true)
    end
end

@testset "issatisfiable: the entry-point family" begin
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    data = Dict(
        :A => PkgData([:v2, :v1], Dict(:v2 => [:B], :v1 => [:B]), nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    compat = Dict(:B => [:v1])
    pins = Dict(:A => :v1)

    # every shape of package data `resolve` accepts, with a `Problem` and with
    # the convenience keywords that build one
    deps = DepsProvider(p -> data[p], keys(data))
    info = pkg_info(data, Problem([:A]; compat, pins))
    for src in (data, deps, info)
        for prob in (Problem([:A]),
                     Problem([:A]; compat),
                     Problem([:A]; pins),
                     Problem([:A]; compat, pins),
                     Problem([:A]; compat = Dict(:B => Symbol[])),
                     Problem([:A]; pins = Dict(:B => :v9)),
                     Problem([:A, :B]; compat = Dict(:Z => Symbol[])),
                     Problem(Symbol[]))
            @test issatisfiable(src, prob) == (resolve(src, prob) !== nothing)
        end
        # the keyword forms are the same calls as the `Problem` ones
        @test issatisfiable(src, [:A]) == issatisfiable(src, Problem([:A]))
        @test issatisfiable(src, [:A]; compat) ==
              issatisfiable(src, Problem([:A]; compat))
        @test issatisfiable(src, [:A]; compat, pins) ==
              issatisfiable(src, Problem([:A]; compat, pins))
        # requirements may be a set as well as a vector, and default to the
        # whole universe
        @test issatisfiable(src, Set([:A])) == issatisfiable(src, [:A])
        @test issatisfiable(src) == (resolve(src) !== nothing)
    end

    # ... and a SAT instance, the shape the others build internally
    sat = SAT(prepare_pkg_info(info, Problem([:A]; compat, pins)))
    try
        @test issatisfiable(sat, [:A]) == (resolve(sat, [:A]) !== nothing)
        @test issatisfiable(sat) == (resolve(sat) !== nothing)
    finally
        finalize(sat)
    end

    # a requirement the data doesn't have is an error, exactly as for `resolve`
    @test_throws ArgumentError issatisfiable(data, [:Nope])
end

@testset "issatisfiable: edge cases" begin
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    data = Dict(
        :A => PkgData(Symbol[], nodeps, nocomp), # no versions: uninstallable
        :B => PkgData([:v1], Dict(:v1 => [:A]), nocomp),
        :C => PkgData([:v1], nodeps, nocomp),
    )
    # no requirements at all is satisfiable (`resolve` returns an empty dict)
    @test issatisfiable(data, Symbol[])
    @test resolve(data, Symbol[]) == Dict()
    # a package with no versions can't be required, directly or transitively —
    # the filter drops it, so the verdict comes from the requirement check
    @test !issatisfiable(data, [:A])
    @test !issatisfiable(data, [:B])
    @test !issatisfiable(data, [:B, :C])
    @test issatisfiable(data, [:C])
    # an unfiltered artifact reaches the same verdicts
    info = pkg_info(data, [:B]; filter = false)
    @test !issatisfiable(info, [:B])
    @test issatisfiable(info, [:C])
end

@testset "issatisfiable: one solve, no descent, no state" begin
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    # A@v2 needs B@v1 and B@v2 needs C@v1, so the best versions are not jointly
    # feasible: the descent has improvement work to do, and one solve cannot be
    # the whole story
    data = Dict(
        :A => PkgData([:v2, :v1], Dict(:v2 => [:B], :v1 => [:B]),
                      Dict(:v2 => Dict(:B => [:v1]))),
        :B => PkgData([:v2, :v1], Dict(:v2 => [:C], :v1 => [:C]),
                      Dict(:v2 => Dict(:C => [:v1]))),
        :C => PkgData([:v2, :v1], nodeps, nocomp),
    )
    reqs = [:A, :B, :C]
    prob = Problem(reqs)
    info = pkg_info(data, prob)

    # the verdict costs a single solve, and since the requirements are assumed
    # rather than asserted, the instance is left exactly as it was found
    sat = SAT(prepare_pkg_info(info, prob))
    try
        clauses = PicoSAT.clause_count(sat.pico)
        # can the solve count be read here? the clause counts below stand on
        # their own either way — asserting nothing is what leaves the instance
        # reusable, and it is also what the descent cannot do
        counts = solve_count(sat) !== nothing
        counts || @info "picosat's solve count is unreadable here: not checked"
        @test !counts || solve_count(sat) == 0
        @test issatisfiable(sat, reqs)
        @test !counts || solve_count(sat) == 1
        @test PicoSAT.clause_count(sat.pico) == clauses
        # ... so the same instance answers again, at the same price
        @test issatisfiable(sat, reqs)
        @test !counts || solve_count(sat) == 2
        @test PicoSAT.clause_count(sat.pico) == clauses
    finally
        finalize(sat)
    end

    # the descent, for comparison: it asserts the requirements and pins each
    # package it optimizes, so it both solves more than once (three times, on
    # this data) and — until the temporary clauses are rolled back — leaves
    # clauses behind
    sat = SAT(prepare_pkg_info(info, prob))
    try
        counts = solve_count(sat) !== nothing
        clauses = PicoSAT.clause_count(sat.pico)
        @test resolve(sat, reqs; restore = false) == resolve(data, reqs)
        @test !counts || solve_count(sat) > 1
        @test PicoSAT.clause_count(sat.pico) > clauses
    finally
        finalize(sat)
    end
end

@testset "issatisfiable: the verdict is order-independent" begin
    # `issatisfiable` takes no `by` and no `order` because neither can move the
    # verdict: orderings select among the valid solutions rather than deciding
    # which ones are valid. Sweep the tiny grids under both, constrained and
    # not, and check the one verdict against every `resolve` in sight.
    Random.seed!(rand(RandomDevice(), UInt64))
    hi = p -> p
    lo = p -> -p
    up = p -> (u, v) -> u > v # prefer the *lowest* version
    for (m, n) in ((2, 2), (2, 3), (3, 2), (3, 3))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:20
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            for cs in (nothing, random_constraints(m, n))
                compat, pins = cs === nothing ?
                    (Dict{Int,Vector{Int}}(), Dict{Int,Int}()) : cs
                prob = Problem(reqs; compat, pins)
                verdict = issatisfiable(data, prob)
                for by in (hi, lo), order in (nothing, up)
                    @test verdict ==
                        (resolve(data, prob; by, order) !== nothing)
                end
            end
        end
    end
end
