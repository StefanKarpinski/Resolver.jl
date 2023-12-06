using Revise
using Resolver
includet("TinyTypes.jl")

# [(m, n, m*n*(m-1)) for m=2:5 for n=2:5 if (m*n)^2 <= 128]
const m = 3 # number of packages
const n = 2 # number of versions

@assert m*n*m*n ≤ TinyTypes.N

const deps_mask = reduce(|,
    TinyTypes.UIntN(1) << ((p-1)*n*m + (v-1)*m + (q-1))
    for v=1:n for w=1:n for p=1:m for q=1:m if p ≠ q
)
const comp_mask = reduce(|,
    TinyTypes.UIntN(1) << ((p-1)*n*m*n + (v-1)*m*n + (q-1)*n + (w-1))
    for v=1:n for w=1:n for p=1:m for q=1:m
    if isodd(p + q) ? p < q : p > q
)

const d = count_ones(deps_mask)
const c = count_ones(comp_mask)

const Deps = TinyDict{n*m, TinyDict{m, TinyVec}}
const Comp = TinyDict{n*m*n, TinyDict{m*n, TinyDict{n, TinyVec}}}

rand_deps() = deposit_bits(deps_mask, randbits(d)) |> Deps
rand_comp() = deposit_bits(comp_mask, randbits(c)) |> Comp

#=
const configs = []
for deps_bits in BinomialBits(d, d÷2)
    deps = deposit_bits(deps_mask, deps_bits) |> Deps
    all(deps[i].bits ≥ deps[i+1].bits for i=1:m-1) || continue
    for comp_bits in BinomialBits(c, c÷2)
        comp = deposit_bits(comp_mask, comp_bits) |> Comp
        for reqs_bits = 1:2^m-1
            reqs = TinyVec(reqs_bits)
            push!(configs, (deps, comp, reqs))
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
