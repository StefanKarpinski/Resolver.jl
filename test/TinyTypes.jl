module TinyTypes

export TinyDict, TinyVec, TinyRange, randbits, deposit_bits

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

end # module

using .TinyTypes
