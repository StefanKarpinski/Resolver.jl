includet("TinyTypes.jl")

make_reqs(b) = TinyVec(b)

function tiny_data_makers(m::Int, n::Int)
    (m*n)^2 ≤ 128 || throw(ArgumentError("m=$m and n=$n are too big"))

    f = false
    t = true

    Deps = TinyDict{n*m, TinyDict{m, TinyVec, f}, f}
    Comp = TinyDict{n*m*n, TinyDict{m*n, TinyDict{n, TinyVec, t}, f}, f}

    𝟘 = TinyTypes.UIntN(0)
    𝟙 = TinyTypes.UIntN(1)

    bit(p, v, q) = p == q ? 𝟘 : 𝟙 << ((p-1)*n*m + (v-1)*m + (q-1))
    bit(p, v, q, w) = p == q ? 𝟘 : isodd(p + q) ⊻ (p < q) ?
        𝟙 << ((p-1)*n*m*n + (v-1)*m*n + (q-1)*n + (w-1)) :
        𝟙 << ((q-1)*n*m*n + (w-1)*m*n + (p-1)*n + (v-1))

    d_mask = reduce(|, init=𝟘, bit(p, v, q) for p=1:m, v=1:n, q=1:m)
    c_mask = reduce(|, init=𝟘, bit(p, v, q, w) for p=1:m, v=1:n, q=1:m, w=1:n)

    d = count_ones(d_mask)
    c = count_ones(c_mask)

    make_deps(b) = deposit_bits(d_mask, b) |> Deps
    make_comp(b) = deposit_bits(c_mask, b) |> Comp

    PkgD = typeof(PkgData(TinyRange(n), Deps(0)[1], Comp(0)[1]))
    data = Dict{Int,PkgD}()

    return d, c, data, make_deps, make_comp, bit
end
