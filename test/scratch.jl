using Revise
includet("setup.jl")
includet("TinyTypes.jl")

#=
let m = 2, n = 2
    d, c, data, make_deps, make_comp = data_makers(m, n)
    for deps_bits = 0:2^d-1
        deps = make_deps(deps_bits)
        all(deps[i].bits ≥ deps[i+1].bits for i=1:m-1) || continue
        for comp_bits = 0:2^c-1
            comp = make_comp(comp_bits)
            for i = 1:m
                data[i] = PkgData(TinyRange(n), deps[i], comp[i])
            end
            for reqs_bits = 1:m^2-1
                reqs = TinyVec(reqs_bits)
                # @show deps_bits, comp_bits, reqs_bits
                test_resolver(data, reqs)
            end
        end
    end
end
=#

#=
begin
    deps = rand_deps()
    comp = rand_comp()
    data = Dict(
        i => PkgData(TinyRange(n), deps[i], comp[i]) for i = 1:m
    )
    pkgs, vers = resolve(data, [1])
    [pkgs vers]
end
=#
