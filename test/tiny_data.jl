module tiny_data

export tiny_data_makers, make_reqs, fill_data!
export randbits, BinomialBits

const UIntN = UInt128
const N = 8*sizeof(UIntN)

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
    (d.bits >> (b*(k-1))) & ((1 << b) - 1)

val(d::TinyDict{b,T,x}, v::UIntN) where {b,T,x} =
    T(x ? ((1 << b) - 1) ⊻ v : v)

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

using ..Resolver: PkgData

const f = false
const t = true

Deps(m, n) = TinyDict{n*m, TinyDict{m, TinyVec, f}, f}
Comp(m, n) = TinyDict{n*m*n, TinyDict{m*n, TinyDict{n, TinyVec, t}, f}, f}

function tiny_data_makers(m::Int, n::Int)
    (m*n)^2 ≤ 128 || throw(ArgumentError("m=$m and n=$n are too big"))

    𝟘 = UIntN(0)
    𝟙 = UIntN(1)

    bit(p, v, q) = p == q ? 𝟘 : 𝟙 << ((p-1)*n*m + (v-1)*m + (q-1))
    bit(p, v, q, w) = p == q ? 𝟘 : isodd(p + q) ⊻ (p < q) ?
        𝟙 << ((p-1)*n*m*n + (v-1)*m*n + (q-1)*n + (w-1)) :
        𝟙 << ((q-1)*n*m*n + (w-1)*m*n + (p-1)*n + (v-1))

    d_mask = reduce(|, init=𝟘, bit(p, v, q) for p=1:m, v=1:n, q=1:m)
    c_mask = reduce(|, init=𝟘, bit(p, v, q, w) for p=1:m, v=1:n, q=1:m, w=1:n)

    d = count_ones(d_mask)
    c = count_ones(c_mask)

    make_deps(b) = deposit_bits(d_mask, b) |> Deps(m, n)
    make_comp(b) = deposit_bits(c_mask, b) |> Comp(m, n)

    PkgD = typeof(PkgData(TinyRange(n), Deps(m,n)(0)[1], Comp(m,n)(0)[1]))
    data = Dict{Int,PkgD}()

    return d, c, data, make_deps, make_comp, bit
end

make_reqs(b) = TinyVec(b)

function fill_data!(m, n, data, deps, comp)
    for i = 1:m
        data[i] = PkgData(TinyRange(n), deps[i], comp[i])
    end
end

# methods for generating bit patterns

randbits(b::Int) = rand(UIntN) & ((1 << b) - 1)

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

end # module

using .tiny_data
