include("setup.jl")

@testset "very small tests, full" begin
    for m = 1:2, n = 1:2
        d, c, data, make_deps, make_comp = tiny_data_makers(m, n)
        for deps_bits = 0:2^d-1
            deps = make_deps(deps_bits)
            all(deps[i].bits ≥ deps[i+1].bits for i=1:m-1) || continue
            for comp_bits = 0:2^c-1
                comp = make_comp(comp_bits)
                for i = 1:m
                    data[i] = PkgData(TinyRange(n), deps[i], comp[i])
                end
                for reqs_bits = 1:2^m-1
                    reqs = TinyVec(reqs_bits)
                    test_resolver(data, reqs)
                end
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
            deps = make_deps(randbits(d))
            comp = make_comp(randbits(c))
            reqs = TinyVec(rand(1:2^m-1))
            for i = 1:m
                data[i] = PkgData(TinyRange(n), deps[i], comp[i])
            end
            test_resolver(data, reqs)
        end
    end
end
