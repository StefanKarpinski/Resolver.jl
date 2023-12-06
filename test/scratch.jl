using Revise
using Resolver
includet("TinyTypes.jl")

# [(m, n, m*n*(m-1)) for m=2:5 for n=2:5 if (m*n)^2 <= 128]
const m = 3 # number of packages
const n = 3 # number of versions

@assert m*n*m*n ≤ TinyTypes.N

const Deps = TinyDict{n*m, TinyDict{m, TinyVec}}
const Comp = TinyDict{n*m*n, TinyDict{m*n, TinyDict{n, TinyVec}}}

const deps_mask = reduce(|,
    TinyTypes.UIntN(1) << ((p-1)*n*m + (v-1)*m + (q-1))
    for v=1:n for w=1:n for p=1:m for q=1:m if p ≠ q
)
const comp_mask = reduce(|,
    TinyTypes.UIntN(1) << ((p-1)*n*m*n + (v-1)*m*n + (q-1)*n + (w-1))
    for v=1:n for w=1:n for p=1:m for q=1:m
    if isodd(p + q) ? p < q : p > q
)

const deps_ones = count_ones(deps_mask)
const comp_ones = count_ones(comp_mask)

rand_deps() = deposit_bits(deps_mask, randbits(deps_ones)) |> Deps
rand_comp() = deposit_bits(comp_mask, randbits(comp_ones)) |> Comp

#=
for _ = 1:10^6
    deps = rand_deps()
    for (p, deps_p) in deps,
        (v, deps_pv) in deps_p,
        q in deps_pv
        @assert p != q
    end
end

for _ = 1:10^6
    comp = rand_comp()
    for (p, comp_p) in comp,
        (v, comp_pv) in comp_p,
        (q, comp_pvq) in comp_pv
        @assert p != q
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
