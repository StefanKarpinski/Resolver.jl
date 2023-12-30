module tiny_data

export
    tiny_data_makers,
    make_reqs,
    make_data,
    fill_data!,
    Deps,
    Comp,
    Data,
    deps_mask,
    comp_mask,
    randbits,
    BinomialBits,
    deposit_bits,
    extract_bits

const UIntN = UInt128
const N = 8*sizeof(UIntN)

const 𝟘 = UIntN(0)
const 𝟙 = UIntN(1)

struct TinyDict{b,T,x} <: AbstractDict{Int,T}
    bits :: UIntN
end

struct TinyVec <: AbstractVector{Int}
    bits :: UIntN
end

struct TinyRange{n} <: AbstractUnitRange{Int}
    # no fields
end

# TinyDict

getbits(d::TinyDict{b}, k::Int) where {b} =
    (d.bits >> (b*(k-1))) & ((𝟙 << b) - 1)

val(d::TinyDict{b,T,x}, v::UIntN) where {b,T,x} =
    T(x ? ((𝟙 << b) - 1) ⊻ v : v)

Base.haskey(d::TinyDict{b,T}, k::Int) where {b,T} =
    getbits(d, k) ≠ 0

Base.getindex(d::TinyDict{b,T}, k::Int) where {b,T} =
    val(d, getbits(d, k))

function Base.iterate(d::TinyDict{b,T}, k::Int = 1) where {b,T}
    while b*k ≤ N
        v = getbits(d, k)
        v ≠ 0 && return k => val(d, v), k+1
        k += 1
    end
end

Base.length(d::TinyDict{b,T}) where {b,T} =
    count(getbits(d, k) ≠ 0 for k = 1:N÷b)

function Base.show(io::IO, d::TinyDict)
    print(io, '{')
    first = true
    for (k, v) in d
        first || print(io, ", ")
        show(io, k)
        print(io, ": ")
        let io = IOContext(io, :typeinfo => typeof(v))
            show(io, v)
        end
        first = false
    end
    print(io, '}')
end

# TinyVec

Base.size(v::TinyVec) = (count_ones(v.bits),)

function Base.getindex(v::TinyVec, i::Int)
    x = 0
    b = v.bits
    for j = 1:i
        s = trailing_zeros(b) + 1
        b >>= s
        x += s
    end
    return x
end

function Base.iterate(v::TinyVec, x::Int = 0)
    b = v.bits >> x
    x += trailing_zeros(b) + 1
    x ≤ N ? (x, x) : nothing
end

Base.in(x::Int, v::TinyVec) = isodd(v.bits >> (x - 1))

# TinyRange

TinyRange(n::Integer) = TinyRange{n}()

Base.first(r::TinyRange) = 1
Base.last(r::TinyRange{n}) where {n} = n::Int

# bit operations: deposit & extract

function deposit_bits(mask::Integer, bits::Integer)
    s = 0
    r = zero(mask)
    while !iszero(mask)
        s += t = trailing_zeros(mask) + 1
        r |= (bits & one(r)) << (s - 1)
        mask >>>= t
        bits >>>= 1
    end
    return r
end

function extract_bits(mask::Integer, bits::Integer)
    a = s = 0
    r = zero(mask)
    while !iszero(mask)
        s += t = trailing_zeros(mask) + 1
        r |= ((bits >>> (s - 1)) & one(r)) << a
        mask >>>= t
        a += 1
    end
    return r
end

# tests: make sure deposit & extract work correctly
for mask = 0x0:0xff
    bits₀ = UInt8(2^count_ones(mask)-1)
    @assert deposit_bits(mask, 0xff) == mask
    @assert extract_bits(mask, 0xff) == bits₀
    @assert extract_bits(mask, mask) == bits₀
    for bits₁ = 0x0:bits₀
        bits₂ = deposit_bits(mask, bits₁)
        @assert count_ones(bits₁) == count_ones(bits₂)
        @assert mask & bits₂ == bits₂
        @assert ~mask & bits₂ == 0
        bits₃ = extract_bits(mask, bits₂)
        @assert bits₃ == bits₁
        bits₄ = bits₂ | (rand(typeof(bits₂)) & ~mask)
        bits₅ = extract_bits(mask, bits₄)
        @assert bits₅ == bits₁
    end
end

# functions for generiting test data

using ..Resolver: PkgData

const f = false
const t = true

Deps(m, n) = TinyDict{n*m, TinyDict{m, TinyVec, f}, f}
Comp(m, n) = TinyDict{n*m*n, TinyDict{m*n, TinyDict{n, TinyVec, t}, f}, f}
Data(m, n) =
    Dict{Int,typeof(PkgData(TinyRange(n), Deps(m,n)(0)[1], Comp(m,n)(0)[1]))}

deps_bit(m, n, p, v, q) = p == q ? 𝟘 :
    𝟙 << ((p-1)*n*m + (v-1)*m + (q-1))
comp_bit(m, n, p, v, q, w) = p == q ? 𝟘 : isodd(p + q) ⊻ (p < q) ?
    𝟙 << ((p-1)*n*m*n + (v-1)*m*n + (q-1)*n + (w-1)) :
    𝟙 << ((q-1)*n*m*n + (w-1)*m*n + (p-1)*n + (v-1))

deps_mask(m, n) = reduce(|, init=𝟘,
    deps_bit(m, n, p, v, q) for p=1:m, v=1:n, q=1:m)
comp_mask(m, n) = reduce(|, init=𝟘,
    comp_bit(m, n, p, v, q, w) for p=1:m, v=1:n, q=1:m, w=1:n)

function tiny_data_makers(m::Int, n::Int)
    (m*n)^2 ≤ 128 || throw(ArgumentError("m=$m and n=$n are too big"))

    d_mask = deps_mask(m, n)
    c_mask = comp_mask(m, n)

    d = count_ones(d_mask)
    c = count_ones(c_mask)

    data = Data(m, n)()

    bit(p, v, q)    = deps_bit(m, n, p, v, q)
    bit(p, v, q, w) = comp_bit(m, n, p, v, q, w)

    make_deps(b) = deposit_bits(d_mask, b) |> Deps(m, n)
    make_comp(b) = deposit_bits(c_mask, b) |> Comp(m, n)

    return make_deps, make_comp, data, d, c, bit
end

make_reqs(b) = TinyVec(b)

function fill_data!(m, n, deps, comp, data)
    for i = 1:m
        data[i] = PkgData(TinyRange(n), deps[i], comp[i])
    end
    return data
end

make_data(m, n, deps, comp) = fill_data!(m, n, deps, comp, Data(m, n)())

# methods for generating bit patterns

randbits(b::Int) = rand(UIntN) & ((𝟙 << b) - 1)

struct BinomialBits{T<:Base.BitInteger}
    n :: T
    k :: T
    function BinomialBits{T}(n::Integer, k::Integer) where {T}
        0 ≤ n ≤ 8*sizeof(T) ||
            throw(ArgumentError("n = $n must be in [0, 8*sizeof($T)]"))
        0 ≤ k ≤ n ||
            throw(ArgumentError("k = $k must be in [0, n]"))
        n == k == 8*sizeof(T) &&
            throw(ArgumentError("n = k = 8*sizeof($T) unsupported"))
        new{T}(n, k)
    end
end
BinomialBits(n::Integer, k::Integer) = BinomialBits{typeof(n)}(n, k)

Base.eltype(b::BinomialBits{T}) where {T} = T
Base.length(b::BinomialBits) = binomial(b.n, b.k)

Base.first(b::BinomialBits) = (one(b.k) << b.k) - one(b.k)
Base.last(b::BinomialBits) = first(b) << (b.n - b.k)

function Base.iterate(
    b :: BinomialBits{T},
    s :: T = -one(T),
) where {T<:Integer}
    s == last(b) && return nothing
    z = trailing_zeros(s)
    s += one(T) << z
    d = b.k - count_ones(s)
    s += one(T) << d - one(T)
    return s, s
end

# some interesting example cases

const examples = []
let
    m, n = 3, 3
    make_deps, make_comp, data = tiny_data_makers(m, n)
    push!(examples, (
        data = make_data(m, n, make_deps(64), make_comp(8165305)),
        reqs = [1, 2, 3],
    ))
    m, n = 5, 2
    make_deps, make_comp, data = tiny_data_makers(m, n)
    push!(examples, (
        data = make_data(m, n, make_deps(524288), make_comp(85916123396)),
        reqs = [1, 2, 3, 4],
    ))
    push!(examples, (
        data = make_data(m, n, make_deps(9), make_comp(4362076160)),
        reqs = [1, 3, 4],
    ))
    push!(examples, (
        data = make_data(m, n, make_deps(17179869696), make_comp(17179869440)),
        reqs = [1, 2, 4, 5],
    ))
    push!(examples, (
        data = make_data(m, n, make_deps(8589934593), make_comp(842351773701)),
        reqs = [1, 3, 4, 5],
    ))
    push!(examples, (
        data = make_data(m, n, make_deps(2148009984), make_comp(5393925603)),
        reqs = [1, 2, 3, 4],
    ))
    push!(examples, (
        data = make_data(m, n, make_deps(2148009984), make_comp(5393991139)),
        reqs = [1, 2, 3, 4],
    ))
end

end # module

using .tiny_data
