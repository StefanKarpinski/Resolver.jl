# `Problem`: user compat bounds & pins as resolve-time constraints.
#
# The oracle throughout is *bake-equivalence*. Forbidding a version by clause
# and deleting it from the data leave the same model set, hence — by the
# reduction-invariance theorem — the same layered answer and the same
# satisfiability verdict. So every constrained resolve is checked against the
# unconstrained resolve of the "baked" data, which is what user constraints
# used to mean.

using Resolver: Problem, PkgData, PkgInfo, pkg_info, SAT, PicoSAT,
    exclusion_masks, is_excluded, finalize, sat_assume, sat_pop, is_satisfiable

# delete the versions `prob` excludes from each package's data, old-style
function bake(
    data :: AbstractDict{P,<:PkgData{P,V}},
    prob :: Problem{P},
) where {P,V}
    Vers = Vector{V}
    Deps = Dict{V,Vector{P}}
    Comp = Dict{V,Dict{P,Any}}
    baked = Dict{P,PkgData{P,V,Any,Vers,Deps,Comp}}()
    for (p, data_p) in data
        vers = V[v for v in data_p.versions if !is_excluded(prob, p, v)]
        deps = Deps(v => P[q for q in data_p.depends[v]]
                    for v in vers if haskey(data_p.depends, v))
        comp = Comp(v => Dict{P,Any}(data_p.compat[v])
                    for v in vers if haskey(data_p.compat, v))
        baked[p] = PkgData(vers, deps, comp)
    end
    return baked
end

# random user constraints over `m` packages with `n` versions each; the m+1st
# package is never in the data, so entries for it must be no-ops. the draw
# covers every case that has to work: no entry, an allow-everything entry, an
# exclude-everything entry (an empty random subset), random subsets, pins at
# existing and non-existing versions, and compat plus a pin on one package
#
# `mild` skips the shapes that make most instances unsatisfiable outright
# (empty compat sets, pins at non-existing versions), for sweeps that need the
# instance to stay solvable for a while.
function random_constraints(m::Int, n::Int; mild::Bool = false)
    compat = Dict{Int,Vector{Int}}()
    pins   = Dict{Int,Int}()
    subset() = mild ? [v for v = 1:n if rand(Bool) || v == rand(1:n)] :
                      [v for v = 1:n if rand(Bool)]
    for p = 1:m+1
        r = rand(1:8)
        if r ≤ 3
            compat[p] = subset()
        elseif r == 4
            compat[p] = collect(1:n)  # allows everything
        elseif r == 5
            pins[p] = rand(1:n)
        elseif r == 6
            pins[p] = mild ? rand(1:n) : n + 1 # no such version
        elseif r == 7
            compat[p] = subset()
            mild || (pins[p] = rand(1:n))
        end
        # r == 8: leave p unconstrained
    end
    return compat, pins
end

# the constraint shapes worth enumerating exhaustively for one package
const CONSTRAINT_SHAPES = [
    :none,        # no entry at all
    :allow_all,   # compat admits every version (must emit nothing)
    :exclude_all, # compat admits nothing
    :best,        # compat admits only the best version
    :worst,       # compat admits only the worst version
    :middle,      # compat excludes the best version only
    :pin_best,
    :pin_worst,
    :pin_missing, # pin at a version that doesn't exist
    :pin_and_compat, # two constraint sources on one package
]

function constrain!(compat, pins, p::Int, shape::Symbol, n::Int)
    shape == :none           ? nothing :
    shape == :allow_all      ? (compat[p] = collect(1:n)) :
    shape == :exclude_all    ? (compat[p] = Int[]) :
    shape == :best           ? (compat[p] = [1]) :
    shape == :worst          ? (compat[p] = [n]) :
    shape == :middle         ? (compat[p] = collect(2:n)) :
    shape == :pin_best       ? (pins[p] = 1) :
    shape == :pin_worst      ? (pins[p] = n) :
    shape == :pin_missing    ? (pins[p] = n + 1) :
    shape == :pin_and_compat ? (compat[p] = [1]; pins[p] = n) :
    error("unknown constraint shape: $shape")
    return compat, pins
end

# `resolve(data, prob)` must agree with the unconstrained resolve of the baked
# data — including when both are `nothing`
function test_bake_equivalence(data, prob; by::Function = identity)
    sol = resolve(data, prob; by)
    @test sol == resolve(bake(data, prob), prob.reqs; by)
    # and the cheap verdict agrees with the descent's, constraints and all
    @test issatisfiable(data, prob) == (sol !== nothing)
end

@testset "Problem: construction & masks" begin
    prob = Problem([:A, :B])
    @test prob.reqs == [:A, :B]
    @test isempty(prob.compat) && isempty(prob.pins)
    prob = Problem(Set([:A]); compat = Dict(:B => [:v1]), pins = Dict(:C => :v2))
    @test prob isa Problem{Symbol,Symbol,Vector{Symbol}}
    @test !is_excluded(prob, :B, :v1)
    @test  is_excluded(prob, :B, :v2)
    @test !is_excluded(prob, :C, :v2)
    @test  is_excluded(prob, :C, :v1)
    @test !is_excluded(prob, :A, :v1) # unconstrained

    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    data = Dict(
        :A => PkgData([:v2, :v1], nodeps, nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
        :C => PkgData([:v2, :v1], nodeps, nocomp),
    )
    info = pkg_info(data, keys(data); filter = false)
    excl = exclusion_masks(info, prob)
    # only constrained packages get a mask; entries for absent packages and
    # allow-everything entries get none
    @test sort!(collect(keys(excl))) == [:B, :C]
    @test excl[:B] == BitVector([true, false]) # compat forbids :v2
    @test excl[:C] == BitVector([false, true]) # the pin forbids :v1
    @test isempty(exclusion_masks(info, Problem([:A])))
    @test isempty(exclusion_masks(info,
        Problem([:A]; compat = Dict(:A => [:v1, :v2], :Z => Symbol[]))))
end

@testset "Problem: the unconstrained case is free" begin
    # `Problem(reqs)` is on every convenience `resolve` path, so it must not
    # allocate a pair of dictionaries just to say "no constraints": the empty
    # compat and pins are one shared immutable value
    reqs = [:A, :B]
    a, b = Problem(reqs), Problem(reverse(reqs))
    @test a.compat === b.compat
    @test a.pins === b.pins
    @test isempty(a.compat) && isempty(a.pins)
    @test valtype(a.compat) == valtype(a.pins) == Any
    @test !haskey(a.compat, :A) && get(a.pins, :A, nothing) === nothing
    # the mask dictionary of an unconstrained problem is shared the same way
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    info = pkg_info(Dict(:A => PkgData([:v1], Dict{Symbol,Vector{Symbol}}(), nocomp)),
                    [:A]; filter = false)
    @test exclusion_masks(info, a) === exclusion_masks(info, b)

    # a caller's dictionary is still copied, so later mutation can't reach in
    d = Dict(:B => [:v1])
    p = Problem(reqs; compat = d)
    @test p.compat == d && p.compat !== d

    # so building one allocates the requirements vector and the struct, and
    # strictly less than the same call handed an (empty) dictionary to copy
    build(::Nothing) = Problem(reqs)
    build(c::AbstractDict) = Problem(reqs; compat = c)
    empty_dict = Dict{Symbol,Any}()
    build(nothing); build(empty_dict) # compile both
    @test @allocated(build(nothing)) < @allocated(build(empty_dict))
end

@testset "Problem: resolve's convenience keywords" begin
    # the requirements-and-keywords form of `resolve` just builds the Problem
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    data = Dict(
        :A => PkgData([:v2, :v1], Dict(:v2 => [:B], :v1 => [:B]), nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    compat = Dict(:B => [:v1])
    pins = Dict(:A => :v1)
    prob = Problem([:A]; compat, pins)
    @test resolve(data, [:A]; compat, pins) == resolve(data, prob)
    @test resolve(data, [:A]; compat) == resolve(data, Problem([:A]; compat))
    @test resolve(data, [:A]) == resolve(data, Problem([:A]))
    info = pkg_info(data, prob)
    @test resolve(info, [:A]; compat, pins) == resolve(info, prob)
end

@testset "Problem: exclusions constrain, they don't delete" begin
    # user constraints are conflict pressure in the filter, not deletions: the
    # forbidden version stays in the info (reachability keeps prefixes), and
    # the SAT clauses are what rule it out
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    data = Dict(
        :A => PkgData([:v1], Dict(:v1 => [:B]), nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    prob = Problem([:A]; compat = Dict(:B => [:v1]))
    info = pkg_info(data, prob)
    @test info[:B].versions == [:v2, :v1] # :v2 forbidden but kept
    @test resolve(data, prob) == Dict(:A => :v1, :B => :v1)
    test_bake_equivalence(data, prob)

    # every version of B forbidden: B saturates, so it keeps all its versions
    # and its dependents degrade past the versions that need it
    data = Dict(
        :A => PkgData([:v2, :v1], Dict(:v2 => [:B]), nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    prob = Problem([:A]; compat = Dict(:B => Symbol[]))
    info = pkg_info(data, prob)
    # B survives -- its uninstallability is constraint-driven, not a deletion
    # -- though redundancy collapses it, since versions that are all excluded
    # alike all dominate each other
    @test haskey(info, :B)
    @test resolve(data, prob) == Dict(:A => :v1)
    test_bake_equivalence(data, prob)
    # ... and if the dependent can't degrade, the problem is unsatisfiable
    data = Dict(
        :A => PkgData([:v1], Dict(:v1 => [:B]), nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    prob = Problem([:A]; compat = Dict(:B => Symbol[]))
    @test resolve(data, prob) === nothing
    test_bake_equivalence(data, prob)

    # an excluded version must not dominate a non-excluded one: B@v2 has no
    # constraints at all, so it would dominate B@v1 (which conflicts with
    # A@v2) if the virtual conflict column were missing -- and then pinning B
    # to v1 would have nothing to select
    data = Dict(
        :A => PkgData([:v2, :v1], nodeps,
                      Dict(:v2 => Dict(:B => [:v2]))),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    prob = Problem([:A, :B]; pins = Dict(:B => :v1))
    info = pkg_info(data, prob)
    @test info[:B].versions == [:v2, :v1]
    @test resolve(data, prob) == Dict(:A => :v1, :B => :v1)
    test_bake_equivalence(data, prob)
end

@testset "Problem: SAT selectors" begin
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    data = Dict(
        :A => PkgData([:v2, :v1], nodeps, nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    info = pkg_info(data, keys(data); filter = false)

    # no constraints at all: exactly the structural instance
    plain = SAT(info)
    same  = SAT(info, Problem([:A]))
    try
        @test isempty(plain.sels) && isempty(same.sels)
        @test PicoSAT.var_count(plain.pico) == PicoSAT.var_count(same.pico)
        @test PicoSAT.clause_count(plain.pico) == PicoSAT.clause_count(same.pico)
    finally
        finalize(plain)
        finalize(same)
    end

    # one selector per constraint source that forbids something; allow-all
    # entries and entries for absent packages get none
    prob = Problem([:A];
        compat = Dict(:A => [:v1], :B => [:v1, :v2], :Z => Symbol[]),
        pins   = Dict(:A => :v1, :Z => :v1))
    sat = SAT(info, prob)
    try
        @test Set(values(sat.sels)) == Set([(:compat, :A), (:pin, :A)])
        @test length(sat.sels) == 2
        # every mask has a selector, every selector a mask
        excl = exclusion_masks(info, prob)
        @test Set(keys(excl)) == Set(p for (_, p) in values(sat.sels))
    finally
        finalize(sat)
    end
end

@testset "Problem: popping the frame relaxes the constraints" begin
    # the user constraints are asserted as unit clauses in a push frame of
    # their own: nothing pops it in production, but one pop makes every
    # previously-forbidden assignment feasible again
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    data = Dict(
        :A => PkgData([:v2, :v1], nodeps, nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    prob = Problem([:A, :B];
        compat = Dict(:A => [:v1]), pins = Dict(:B => :v2))
    info = pkg_info(data, prob; filter = false)
    sat = SAT(info, prob)
    try
        @test length(sat.sels) == 2
        # A@v2 is forbidden by compat, B@v1 by the pin
        sat_assume(sat, :A, 1)
        @test !is_satisfiable(sat)
        sat_assume(sat, :B, 2)
        @test !is_satisfiable(sat)
        # both jointly feasible once the frame is popped
        sat_pop(sat)
        sat_assume(sat, :A, 1)
        sat_assume(sat, :B, 2)
        @test is_satisfiable(sat)
    finally
        finalize(sat)
    end
end

@testset "Problem: exclusion kinds" begin
    # `excludes` carries the *admission* knobs — "no prereleases", "no yanked
    # versions" — which are stated about versions rather than about packages, so
    # they come as `kind => predicate` pairs instead of per-package entries.
    # Semantically a kind is nothing but another constraint source: forbidding
    # exactly the versions a compat entry forbids must be the same problem.
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    data = Dict(
        :A => PkgData([:v2, :v1], Dict(:v2 => [:B], :v1 => [:B]), nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    info = pkg_info(data, keys(data); filter = false)

    odd = :test => (p, v) -> v == :v2 # forbid :v2 of every package
    prob = Problem([:A]; excludes = [odd])
    @test  is_excluded(prob, :A, :v2)
    @test !is_excluded(prob, :A, :v1)
    @test  is_excluded(prob, :Z, :v2) # a kind applies to packages too

    # the kind reaches every package in the universe, unlike compat and pins
    excl = exclusion_masks(info, prob)
    @test sort!(collect(keys(excl))) == [:A, :B]
    @test excl[:A] == excl[:B] == BitVector([true, false])

    # ... and it is the same problem as the equivalent compat bounds
    both = Problem([:A]; compat = Dict(:A => [:v1], :B => [:v1]))
    @test exclusion_masks(info, both) == excl
    @test resolve(data, prob) == resolve(data, both)
    test_bake_equivalence(data, prob)

    # one selector per source, so a kind, a compat bound and a pin on the same
    # package stay distinguishable
    prob = Problem([:A]; compat = Dict(:A => [:v1, :v2]),
                         pins = Dict(:B => :v1), excludes = [odd])
    sat = SAT(info, prob)
    try
        @test Set(values(sat.sels)) ==
              Set([(:test, :A), (:test, :B), (:pin, :B)])
    finally
        finalize(sat)
    end

    # no kinds means no cost: the shared empty vector, and an unconstrained
    # problem that still allocates nothing but its requirements
    a, b = Problem([:A]), Problem([:B])
    @test a.excludes === b.excludes
    @test isempty(a.excludes)
    @test isempty(Problem([:A]; excludes = Pair{Symbol,Any}[]).excludes)
    @test exclusion_masks(info, a) === exclusion_masks(info, b)
end

@testset "Problem: exclusion kinds are ordinary constraints" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    # A kind that forbids a random version subset per package must behave
    # exactly like the compat entries that forbid the same versions, and (by
    # bake-equivalence) exactly like deleting them from the data. Swept over the
    # tiny grids, with compat and pins in force on top.
    for (m, n) in ((2, 2), (2, 3), (3, 2), (3, 3), (2, 4), (4, 2))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:40
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            allowed = Dict(p => [v for v = 1:n if rand(Bool)] for p = 1:m+1)
            kinds = [:test => (p, v) -> !(v in get(allowed, p, 1:n))]
            for cs in (nothing, random_constraints(m, n))
                compat, pins = cs === nothing ?
                    (Dict{Int,Vector{Int}}(), Dict{Int,Int}()) : cs
                as_kind = Problem(reqs; compat, pins, excludes = kinds)
                as_compat = Problem(reqs;
                    compat = mergewith!(∩, Dict(compat), Dict(allowed)), pins)
                for by in (identity, p -> -p)
                    @test resolve(data, as_kind; by) ==
                          resolve(data, as_compat; by)
                    test_bake_equivalence(data, as_kind; by)
                end
                # ... and the collapse and the ordering stay orthogonal to it
                @test resolve(data, as_kind; group = false) ==
                      resolve(data, as_kind; group = true)
                order = p -> (u, v) -> u > v
                @test resolve(data, as_kind; order) ==
                      resolve(bake(data, as_kind), reqs; order)
            end
        end
    end
end

@testset "Problem: brute-force cross-checks" begin
    # hand-built cases: the brute-force reference resolves the baked data
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    A3B2 = Dict(
        :A => PkgData([:v3, :v2, :v1],
                      Dict(:v3 => [:B], :v2 => [:B], :v1 => [:B]), nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    conflict = Dict(
        :A => PkgData([:v3, :v2, :v1],
                      Dict(:v3 => [:B], :v2 => [:B], :v1 => [:B]),
                      Dict(:v3 => Dict(:B => [:v2]),
                           :v2 => Dict(:B => [:v1]))),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
        :C => PkgData([:v2, :v1], Dict(:v2 => [:A]), nocomp),
    )
    cases = [
        (A3B2, Problem([:A]; compat = Dict(:B => [:v1]))),
        (A3B2, Problem([:A]; pins = Dict(:B => :v2))),
        (A3B2, Problem([:A]; compat = Dict(:A => [:v2]))),
        (A3B2, Problem([:A]; compat = Dict(:B => Symbol[]))),
        (A3B2, Problem([:A]; pins = Dict(:B => :v9))),
        (A3B2, Problem([:A, :B]; compat = Dict(:Z => Symbol[]))),
        (A3B2, Problem([:A]; compat = Dict(:A => [:v1, :v2]),
                             pins = Dict(:A => :v2))),
        (conflict, Problem([:C]; pins = Dict(:B => :v2))),
        (conflict, Problem([:C]; compat = Dict(:A => [:v1, :v2]))),
        (conflict, Problem([:C, :B]; compat = Dict(:B => [:v1]))),
        (conflict, Problem([:C]; compat = Dict(:C => [:v2], :B => [:v1]))),
        (conflict, Problem([:A, :B, :C]; pins = Dict(:A => :v3, :B => :v1))),
    ]
    lo = p -> -Int(first(string(p))) # reverse alphabetical priority
    for (data, prob) in cases, by in (identity, lo)
        baked = bake(data, prob)
        @test resolve(data, prob; by) == ref_resolve(baked, prob.reqs; by)
        test_bake_equivalence(data, prob; by)
    end
    for (data, prob) in cases
        test_resolver(bake(data, prob), prob.reqs)
    end
end

@testset "Problem: bake equivalence, complete data grids" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    # every dependency & compatibility pattern of the smallest grids, with
    # random user constraints on top of each
    for m = 1:2, n = 1:2
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for deps_bits = 0:2^d-1
            deps = make_deps(deps_bits)
            for comp_bits = 0:2^c-1
                comp = make_comp(comp_bits)
                fill_data!(m, n, deps, comp, data)
                for reqs_bits = 0:2^m-1
                    reqs = collect(make_reqs(reqs_bits))
                    for _ = 1:4
                        compat, pins = random_constraints(m, n)
                        prob = Problem(reqs; compat, pins)
                        test_bake_equivalence(data, prob)
                    end
                end
            end
        end
    end
end

@testset "Problem: bake equivalence, exhaustive constraint shapes" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    hi = p -> p  # default priority: lower package id first
    lo = p -> -p # reversed priority
    # every constraint shape on every package, on random data
    for (m, n) in ((2, 2), (2, 3), (3, 2))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:3
            deps = make_deps(randbits(d))
            comp = make_comp(randbits(c))
            fill_data!(m, n, deps, comp, data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            for shapes in Iterators.product(ntuple(_ -> CONSTRAINT_SHAPES, m)...)
                compat = Dict{Int,Vector{Int}}()
                pins   = Dict{Int,Int}()
                for (p, shape) in enumerate(shapes)
                    constrain!(compat, pins, p, shape, n)
                end
                prob = Problem(reqs; compat, pins)
                test_bake_equivalence(data, prob; by = hi)
                test_bake_equivalence(data, prob; by = lo)
            end
        end
    end
end

@testset "Problem: bake equivalence, random grids" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    hi = p -> p
    lo = p -> -p
    for (m, n) in ((2, 4), (4, 2), (2, 5), (5, 2), (3, 3))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:120
            deps = make_deps(randbits(d))
            comp = make_comp(randbits(c))
            fill_data!(m, n, deps, comp, data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            for _ = 1:4
                compat, pins = random_constraints(m, n)
                prob = Problem(reqs; compat, pins)
                test_bake_equivalence(data, prob; by = rand((hi, lo)))
            end
        end
    end
end

@testset "Problem: bake equivalence, adversarial" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    # the adversarial generator from the main suite, but resolving a
    # constrained problem: break the solution until it is unsolvable, with
    # user constraints in force the whole way
    for m = 2:4, n = 2:3
        (m*n)^2 ≤ 128 || continue
        make_deps, make_comp, data, d, c, bit = tiny_data_makers(m, n)
        for _ = 1:40
            deps = make_deps(0)
            comp = make_comp(0)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            compat, pins = random_constraints(m, n; mild = true)
            prob = Problem(reqs; compat, pins)
            while true
                fill_data!(m, n, deps, comp, data)
                test_bake_equivalence(data, prob)
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
