module TinyTypes

export TinyDict, TinyVec, TinyRange, BinomialBits, randbits, deposit_bits

const UIntN = UInt128
const N = 8*sizeof(UIntN)

struct TinyDict{b,T} <: AbstractDict{Int,T}
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

Base.haskey(d::TinyDict{b,T}, k::Int) where {b,T} =
    getbits(d, k) ≠ 0

Base.getindex(d::TinyDict{b,T}, k::Int) where {b,T} =
    T(getbits(d, k))

function Base.iterate(d::TinyDict{b,T}, k::Int = 1) where {b,T}
    while b*k ≤ N
        v = getbits(d, k)
        v ≠ 0 && return k => T(v), k+1
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
        show(io, v)
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

Base.in(x::Int, v::TinyVec) = (v.bits >> x) ≠ 0

# TinyRange

TinyRange(n::Integer) = TinyRange{n}()

Base.first(r::TinyRange) = 1
Base.last(r::TinyRange{n}) where {n} = n::Int

randbits(b::Int) = rand(TinyTypes.UIntN) & ((1 << b) - 1)

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

using .TinyTypes
