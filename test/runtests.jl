include("setup.jl")

@testset "small comprehensive tests" begin
    for m = 1:3, n = 1:3
        m*n ≤ 6 || continue
        @show m, n
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
                    @show deps_bits, comp_bits, reqs_bits
                    test_resolver(data, reqs)
                end
            end
        end
    end
end
