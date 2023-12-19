include("TinyTypes.jl")

function tiny_data_makers(m::Int, n::Int)
    (m*n)^2 ≤ 128 || throw(ArgumentError("m=$m and n=$n are too big"))

    Deps = TinyDict{n*m, TinyDict{m, TinyVec}}
    Comp = TinyDict{n*m*n, TinyDict{m*n, TinyDict{n, TinyVec}}}

    deps_mask = reduce(|, init = TinyTypes.UIntN(0),
        TinyTypes.UIntN(1) << ((p-1)*n*m + (v-1)*m + (q-1))
        for v=1:n for w=1:n for p=1:m for q=1:m if p ≠ q
    )
    comp_mask = reduce(|, init = TinyTypes.UIntN(0),
        TinyTypes.UIntN(1) << ((p-1)*n*m*n + (v-1)*m*n + (q-1)*n + (w-1))
        for v=1:n for w=1:n for p=1:m for q=1:m
        if isodd(p + q) ? p < q : p > q
    )

    d = count_ones(deps_mask)
    c = count_ones(comp_mask)

    make_deps(b) = deposit_bits(deps_mask, b) |> Deps
    make_comp(b) = deposit_bits(comp_mask, b) |> Comp

    PkgD = typeof(PkgData(TinyRange(n), Deps(0)[1], Comp(0)[1]))
    data = Dict{Int,PkgD}()

    return d, c, data, make_deps, make_comp
end
