# `Problem`: user compat bounds & pins as resolve-time constraints.
#
# The oracle throughout is *bake-equivalence*. Forbidding a version by clause
# and deleting it from the data leave the same model set, hence — by the
# reduction-invariance theorem — the same layered answer and the same
# satisfiability verdict. So every constrained resolve is checked against the
# unconstrained resolve of the "baked" data, which is what user constraints
# used to mean.

using Resolver: Problem, PkgData, PkgInfo, pkg_info, SAT, PicoSAT,
    exclusion_masks, is_excluded, finalize, sat_assume, sat_pop,
    is_satisfiable, sat_solve, rank_pkg_info, prepare_pkg_info, Constraint, relax

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
    pin    = Dict{Int,Int}()
    subset() = mild ? [v for v = 1:n if rand(Bool) || v == rand(1:n)] :
                      [v for v = 1:n if rand(Bool)]
    for p = 1:m+1
        r = rand(1:8)
        if r ≤ 3
            compat[p] = subset()
        elseif r == 4
            compat[p] = collect(1:n)  # allows everything
        elseif r == 5
            pin[p] = rand(1:n)
        elseif r == 6
            pin[p] = mild ? rand(1:n) : n + 1 # no such version
        elseif r == 7
            compat[p] = subset()
            mild || (pin[p] = rand(1:n))
        end
        # r == 8: leave p unconstrained
    end
    return compat, pin
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
    :pin_and_compat, # two constraints on one package
]

function constrain!(compat, pin, p::Int, shape::Symbol, n::Int)
    shape == :none           ? nothing :
    shape == :allow_all      ? (compat[p] = collect(1:n)) :
    shape == :exclude_all    ? (compat[p] = Int[]) :
    shape == :best           ? (compat[p] = [1]) :
    shape == :worst          ? (compat[p] = [n]) :
    shape == :middle         ? (compat[p] = collect(2:n)) :
    shape == :pin_best       ? (pin[p] = 1) :
    shape == :pin_worst      ? (pin[p] = n) :
    shape == :pin_missing    ? (pin[p] = n + 1) :
    shape == :pin_and_compat ? (compat[p] = [1]; pin[p] = n) :
    error("unknown constraint shape: $shape")
    return compat, pin
end

# `resolve(data, prob)` must agree with the unconstrained resolve of the baked
# data — including when both are `nothing`
function test_bake_equivalence(data, prob; by::Function = identity)
    # the answer, not the report an unsatisfiable resolve otherwise returns:
    # the two sides here agree on the verdict, and a report about constraints
    # only one side has is a different thing to compare (see test/diagnostics.jl)
    sol = resolve(data, prob; by, diagnose = false)
    @test sol == resolve(bake(data, prob), prob.reqs; by, diagnose = false)
    # and the cheap verdict agrees with the descent's, constraints and all
    @test issatisfiable(data, prob) == (sol !== nothing)
end

@testset "Problem: construction & masks" begin
    prob = Problem([:A, :B])
    @test prob.reqs == [:A, :B]
    @test isempty(prob.constraints)
    prob = Problem(Set([:A]); compat = Dict(:B => [:v1]), pin = Dict(:C => :v2))
    # one namespace, keyed by kind. What a constraint is made of is not part of
    # what it is: every kind is the same type, and only what it names differs
    @test prob isa Problem{Symbol}
    @test sort!(collect(keys(prob.constraints))) == [:compat, :pin]
    @test valtype(prob.constraints) == Constraint{Symbol}
    @test prob.constraints[:compat].named == Set([:B])
    @test prob.constraints[:pin].named == Set([:C])
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

@testset "Problem: every keyword is a kind, and its shape is checked" begin
    # a keyword's name *is* the kind's name, so a misspelling of a well-known
    # one is a perfectly good name for a kind of its own; what catches it is
    # that a dictionary is not a predicate. The message has to name the kind,
    # since that is the word the caller has to change — the rest of its wording
    # is not the test's business
    for bad in (:(compt = Dict(:B => [:v1])),      # a dictionary, kind unknown
                :(compat = (p, v) -> true),        # a callable where one is not
                :(pin = [:v1, :v2]),               # not a dictionary at all
                :(compat = Dict("B" => [:v1])))    # keyed by the wrong package
        kind = String(bad.args[1])
        e = @test_throws ArgumentError @eval Problem([:A]; $bad)
        @test occursin(kind, e.value.msg)
    end
    # a kind of its own takes any callable, and reaches every package
    prob = Problem([:A]; pre = (p, v) -> v == :v2)
    @test prob.constraints[:pre].named === nothing  # it speaks about all of them
    @test is_excluded(prob, :Z, :v2)
end

@testset "Problem: a constraint answers three questions" begin
    # whatever its shape, a constraint says whether it forbids a version, which
    # packages it speaks about, and what is left of it once some of them are
    # relaxed — and nothing else. Everything above is a loop over those.
    prob = Problem([:A]; compat = Dict(:A => [:v1]), pin = Dict(:B => :v2),
                         pre = (p, v) -> v == :v3)
    # relaxing a constraint for everything it names leaves nothing of it,
    # which is what makes "relaxed outright" need no convention of its own
    for (_, c) in prob.constraints
        @test relax(c, c.named) === nothing
        @test relax(c, Set{Symbol}()) !== nothing
    end
    compat = prob.constraints[:compat]
    @test compat.named == Set([:A])
    @test compat.forbids(:A, :v2) && !compat.forbids(:A, :v1)
    @test !compat.forbids(:Z, :v2) # a package it doesn't name
    pin = prob.constraints[:pin]
    @test pin.named == Set([:B])
    @test pin.forbids(:B, :v1) && !pin.forbids(:B, :v2)
    pre = prob.constraints[:pre]
    @test pre.named === nothing # a predicate cannot enumerate them
    @test pre.forbids(:Z, :v3) && pre.forbids(:Y, :v3)
    # ... but it can still be relaxed for a package, like anything else
    @test !relax(pre, [:Z]).forbids(:Z, :v3)
    @test  relax(pre, [:Z]).forbids(:Y, :v3)
end

@testset "Problem: the unconstrained case is free" begin
    # `Problem(reqs)` is on every convenience `resolve` path, so it must not
    # allocate a dictionary just to say "no constraints": the empty constraint
    # map is one shared immutable value
    reqs = [:A, :B]
    a, b = Problem(reqs), Problem(reverse(reqs))
    @test a.constraints === b.constraints
    @test isempty(a.constraints)
    @test !haskey(a.constraints, :compat)
    @test get(a.constraints, :pin, nothing) === nothing
    # the mask dictionary of an unconstrained problem is shared the same way
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    info = pkg_info(Dict(:A => PkgData([:v1], Dict{Symbol,Vector{Symbol}}(), nocomp)),
                    [:A]; filter = false)
    @test exclusion_masks(info, a) === exclusion_masks(info, b)

    # a caller's dictionary is still copied, so later mutation can't reach in
    d = Dict(:B => [:v1])
    p = Problem(reqs; compat = d)
    @test !is_excluded(p, :B, :v1) && is_excluded(p, :B, :v2)
    d[:B] = [:v2]  # mutating it afterwards cannot reach in
    @test !is_excluded(p, :B, :v1) && is_excluded(p, :B, :v2)

    # so building one allocates the requirements vector and the struct and
    # nothing else. An empty dictionary names no constraint to make, so it gets
    # the same shared map, and what it costs over the no-keyword case is the
    # keyword itself rather than a dictionary — where an entry costs one
    build(::Nothing) = Problem(reqs)
    build(c::AbstractDict) = Problem(reqs; compat = c)
    empty_dict = Dict{Symbol,Any}()
    one_entry = Dict(:B => [:v1])
    build(nothing); build(empty_dict); build(one_entry) # compile all three
    @test build(empty_dict).constraints === build(nothing).constraints
    @test @allocated(build(empty_dict)) - @allocated(build(nothing)) < 64
    @test @allocated(build(nothing)) < @allocated(build(one_entry))
end

@testset "Problem: bare requirements are the unconstrained problem" begin
    # `resolve` takes a `Problem` or bare requirements, and bare requirements
    # are `Problem(reqs)`. Constraints are the `Problem`'s alone: there is no
    # keyword here to pass them by, so `by`/`order` stay a closed set and a
    # misspelled one is a `MethodError` rather than a constraint nobody meant
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    data = Dict(
        :A => PkgData([:v2, :v1], Dict(:v2 => [:B], :v1 => [:B]), nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    info = pkg_info(data, keys(data); filter = false)
    @test resolve(data, [:A]) == resolve(data, Problem([:A]))
    @test resolve(info, [:A]) == resolve(info, Problem([:A]))
    @test issatisfiable(data, [:A]) == issatisfiable(data, Problem([:A]))
    @test_throws MethodError resolve(data, [:A]; compat = Dict(:B => [:v1]))
    @test_throws MethodError issatisfiable(data, [:A]; pin = Dict(:A => :v1))
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
    @test resolve(data, prob; diagnose = false) === nothing
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
    prob = Problem([:A, :B]; pin = Dict(:B => :v1))
    info = pkg_info(data, prob)
    @test info[:B].versions == [:v2, :v1]
    @test resolve(data, prob) == Dict(:A => :v1, :B => :v1)
    test_bake_equivalence(data, prob)
end

@testset "Problem: constraints reach SAT as deactivated classes" begin
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    data = Dict(
        :A => PkgData([:v2, :v1], nodeps, nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    # nothing distinguishes either package's versions, so each package is one
    # class of two members -- which is exactly why a constraint sparing one
    # member has nothing to say to the solver
    info = pkg_info(data, keys(data); filter = false)
    @test info[:A].classes == info[:B].classes == [1, 1]

    # no constraints at all: exactly the structural instance
    plain = SAT(info)
    same  = SAT(rank_pkg_info(info, Problem([:A])))
    nvars = nclauses = 0
    try
        @test isempty(plain.deact) && isempty(same.deact)
        nvars = PicoSAT.var_count(plain.pico)
        nclauses = PicoSAT.clause_count(plain.pico)
        @test nvars == PicoSAT.var_count(same.pico)
        @test nclauses == PicoSAT.clause_count(same.pico)
    finally
        finalize(plain)
        finalize(same)
    end

    # a constraint that leaves a member standing is likewise no clause at all:
    # :A's class still has :v1, so the instance is the structural one again
    sat = SAT(rank_pkg_info(info, Problem([:A]; compat = Dict(:A => [:v1]))))
    try
        @test isempty(sat.deact)
        @test PicoSAT.var_count(sat.pico) == nvars
        @test PicoSAT.clause_count(sat.pico) == nclauses
    finally
        finalize(sat)
    end

    # one literal per emptied class, whatever emptied it, and none for
    # allow-all entries or entries for absent packages
    prob = Problem([:A];
        compat = Dict(:A => [:v1], :B => [:v1, :v2], :Z => Symbol[]),
        pin    = Dict(:A => :v2))
    univ = rank_pkg_info(info, prob)
    sat = SAT(univ)
    try
        # compat allows only :v1, the pin only :v2: jointly nothing, so :A's
        # one class is empty. :B's is untouched
        @test univ.reps[:A] == [0]
        @test univ.reps[:B] == [1]
        @test sat.deact == [-(sat.vars[:A] + 1)]
        @test !is_satisfiable(sat, [:A])
    finally
        finalize(sat)
    end
end

@testset "Problem: popping the frame reactivates the classes" begin
    # the deactivations are asserted as unit clauses in a push frame of their
    # own: nothing pops it in production, but one pop makes every deactivated
    # class choosable again
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    # :A@v2 conflicts with :B@v2, so each package has two classes of its own
    data = Dict(
        :A => PkgData([:v2, :v1], nodeps, Dict(:v2 => Dict(:B => [:v1]))),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    prob = Problem([:A, :B];
        compat = Dict(:A => [:v1]), pin = Dict(:B => :v2))
    info = pkg_info(data, prob; filter = false)
    univ = rank_pkg_info(info, prob)
    sat = SAT(univ)
    try
        @test length(sat.deact) == 2
        # :A's :v2 class is emptied by compat, :B's :v1 class by the pin
        sat_assume(sat, :A, 1)
        @test !sat_solve(sat)
        sat_assume(sat, :B, 2)
        @test !sat_solve(sat)
        # both jointly feasible once the frame is popped
        sat_pop(sat)
        sat_assume(sat, :A, 1)
        sat_assume(sat, :B, 2)
        @test sat_solve(sat)
    finally
        finalize(sat)
    end
end

@testset "Problem: exclusion kinds" begin
    # a keyword of its own carries an *admission* knob — "no prereleases" and
    # the like — which is stated about versions rather than about packages, so
    # it comes as a predicate instead of per-package entries. Semantically a
    # kind is nothing but another constraint: forbidding exactly the versions a
    # compat entry forbids must be the same problem.
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    data = Dict(
        :A => PkgData([:v2, :v1], Dict(:v2 => [:B], :v1 => [:B]), nocomp),
        :B => PkgData([:v2, :v1], nodeps, nocomp),
    )
    info = pkg_info(data, keys(data); filter = false)

    odd = (p, v) -> v == :v2 # forbid :v2 of every package
    prob = Problem([:A]; test = odd)
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
    @test resolve(data, prob; diagnose = false) ==
          resolve(data, both; diagnose = false)
    test_bake_equivalence(data, prob)

    # kinds are indistinguishable once they have emptied a class -- and here
    # they have emptied none. Nothing tells :A's two versions apart, nor :B's,
    # so each package is one class; the kind forbids :v2 of each and the class
    # survives through :v1. Three constraints on two packages, and
    # nothing at all for any of them to say to the solver
    @test info[:A].classes == info[:B].classes == [1, 1]
    prob = Problem([:A]; compat = Dict(:A => [:v1, :v2]),
                         pin = Dict(:B => :v1), test = odd)
    univ = rank_pkg_info(info, prob)
    sat = SAT(univ)
    try
        @test univ.reps[:A] == [2] # the class of {:v2, :v1}, standing at :v1
        @test univ.reps[:B] == [2]
        @test isempty(sat.deact)
    finally
        finalize(sat)
    end
    # ... whereas a pin the kind contradicts leaves the class nothing at all
    prob = Problem([:A]; pin = Dict(:B => :v2), test = odd)
    univ = rank_pkg_info(info, prob)
    sat = SAT(univ)
    try
        @test univ.reps[:B] == [0]
        @test sat.deact == [-(sat.vars[:B] + 1)]
    finally
        finalize(sat)
    end

    # no kinds means no cost: the shared empty map, and an unconstrained
    # problem that still allocates nothing but its requirements
    a, b = Problem([:A]), Problem([:B])
    @test a.constraints === b.constraints
    @test isempty(a.constraints)
    @test Problem([:A]; compat = Dict{Symbol,Any}()).constraints ===
          a.constraints # a kind with nothing in it is a kind that isn't there
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
            test = (p, v) -> !(v in get(allowed, p, 1:n))
            for cs in (nothing, random_constraints(m, n))
                compat, pin = cs === nothing ?
                    (Dict{Int,Vector{Int}}(), Dict{Int,Int}()) : cs
                as_kind = Problem(reqs; compat, pin, test)
                as_compat = Problem(reqs;
                    compat = mergewith!(∩, Dict(compat), Dict(allowed)), pin)
                for by in (identity, p -> -p)
                    @test resolve(data, as_kind; by, diagnose = false) ==
                          resolve(data, as_compat; by, diagnose = false)
                    test_bake_equivalence(data, as_kind; by)
                end
                # ... and the ordering stays orthogonal to it
                order = p -> (u, v) -> u > v
                @test resolve(data, as_kind; order, diagnose = false) ==
                      resolve(bake(data, as_kind), reqs; order, diagnose = false)
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
        (A3B2, Problem([:A]; pin = Dict(:B => :v2))),
        (A3B2, Problem([:A]; compat = Dict(:A => [:v2]))),
        (A3B2, Problem([:A]; compat = Dict(:B => Symbol[]))),
        (A3B2, Problem([:A]; pin = Dict(:B => :v9))),
        (A3B2, Problem([:A, :B]; compat = Dict(:Z => Symbol[]))),
        (A3B2, Problem([:A]; compat = Dict(:A => [:v1, :v2]),
                             pin = Dict(:A => :v2))),
        (conflict, Problem([:C]; pin = Dict(:B => :v2))),
        (conflict, Problem([:C]; compat = Dict(:A => [:v1, :v2]))),
        (conflict, Problem([:C, :B]; compat = Dict(:B => [:v1]))),
        (conflict, Problem([:C]; compat = Dict(:C => [:v2], :B => [:v1]))),
        (conflict, Problem([:A, :B, :C]; pin = Dict(:A => :v3, :B => :v1))),
    ]
    lo = p -> -Int(first(string(p))) # reverse alphabetical priority
    for (data, prob) in cases, by in (identity, lo)
        baked = bake(data, prob)
        @test resolve(data, prob; by, diagnose = false) ==
              ref_resolve(baked, prob.reqs; by)
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
                        compat, pin = random_constraints(m, n)
                        prob = Problem(reqs; compat, pin)
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
                pin    = Dict{Int,Int}()
                for (p, shape) in enumerate(shapes)
                    constrain!(compat, pin, p, shape, n)
                end
                prob = Problem(reqs; compat, pin)
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
                compat, pin = random_constraints(m, n)
                prob = Problem(reqs; compat, pin)
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
            compat, pin = random_constraints(m, n; mild = true)
            prob = Problem(reqs; compat, pin)
            while true
                fill_data!(m, n, deps, comp, data)
                test_bake_equivalence(data, prob)
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
