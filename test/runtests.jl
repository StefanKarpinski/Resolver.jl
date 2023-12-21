include("setup.jl")

@testset "tiny tests, complete" begin
    for m = 1:2, n = 1:2
        d, c, data, make_deps, make_comp = tiny_data_makers(m, n)
        # all possible dependency patterns
        for deps_bits = 0:2^d-1
            deps = make_deps(deps_bits)
            all(deps[i].bits ≥ deps[i+1].bits for i=1:m-1) || continue
            # all possible compatibility patterns
            for comp_bits = 0:2^c-1
                comp = make_comp(comp_bits)
                fill_data!(m, n, data, deps, comp)
                # all possible requirements sets
                for reqs_bits = 0:2^m-1
                    reqs = make_reqs(reqs_bits)
                    test_resolver(data, reqs)
                end
            end
        end
    end
end

@testset "small tests, semi-full" begin
    for m = 2:3, n = 2:3
        m == n && continue # fully tested or too large
        d, c, data, make_deps, make_comp = tiny_data_makers(m, n)
        # all dependency patterns + random compatibility pattern
        for deps_bits = 0:2^d-1
            deps = make_deps(deps_bits)
            comp = make_comp(randbits(c))
            fill_data!(m, n, data, deps, comp)
            for reqs_bits = 1:2^m-1
                # @show m, n, deps_bits, reqs_bits
                # println("comp = Comp(", comp.bits, ")")
                reqs = make_reqs(reqs_bits)
                test_resolver(data, reqs)
            end
        end
        # all compatibility patterns + random dependency pattern
        for comp_bits = 0:2^c-1
            deps = make_deps(randbits(d))
            comp = make_comp(comp_bits)
            fill_data!(m, n, data, deps, comp)
            for reqs_bits = 1:2^m-1
                # @show m, n, comp_bits, reqs_bits
                # println("deps = Deps(", deps.bits, ")")
                reqs = make_reqs(reqs_bits)
                test_resolver(data, reqs)
            end
        end
    end
end

@testset "small tests, random" begin
    # chosen to hit some previously failing cases:
    Random.seed!(0x8cb0074336f2d04f)
    for m = 2:5, n = 2:5
        16 < (m*n)^2 ≤ 128 || continue
        d, c, data, make_deps, make_comp = tiny_data_makers(m, n)
        for _ = 1:1000
            # random dependencies, compatibility and requirements
            deps = make_deps(randbits(d))
            comp = make_comp(randbits(c))
            reqs = make_reqs(rand(1:2^m-1))
            fill_data!(m, n, data, deps, comp)
            test_resolver(data, reqs)
        end
    end
end

@testset "medium tests, adversarial" begin
    for m = 2:5, n = 1:5
        (m*n)^2 ≤ 128 || continue
        d, c, data, make_deps, make_comp, bit = tiny_data_makers(m, n)
        for _ = 1:10
            deps = make_deps(0) # start with no dependencies
            comp = make_comp(0) # start with no conflicts
            reqs = make_reqs(rand(1:2^m-1)) # random requirements set
            # iteratively pick a solution and break it until unsolvable
            while true
                fill_data!(m, n, data, deps, comp)
                pkgs, vers = test_resolver(data, reqs)
                @show size(vers)
                all(isnothing, vers) && break
                # pick a random solution
                k = rand(1:size(vers,2))
                # there must be some non-nothing version
                @test any(vers[i, k] !== nothing for i=1:size(vers,1))
                # pick a package with non-nothing version
                i = rand(1:size(vers,1))
                while vers[i, k] === nothing
                    i = rand(1:size(vers,1))
                end
                p = pkgs[i]
                v = vers[i, k]
                # pick a different package (version can be nothing)
                q = rand(1:m-1)
                q += q ≥ p
                j = something(findfirst(==(q), pkgs), 0)
                w = get(vers, (j, k), 0)
                # make resolved versions of p & q incompatible
                if w == 0
                    # add a dependency p@v => q
                    x = bit(p, v, q)
                    @assert iszero(deps.bits & x)
                    deps = typeof(deps)(deps.bits | x)
                else
                    # add incompatibility p@v ⊼ p@w
                    w = vers[j, k]
                    x = bit(p, v, q, w)
                    @assert iszero(deps.bits & x)
                    comp = typeof(comp)(comp.bits | x)
                end
            end
        end
    end
end
