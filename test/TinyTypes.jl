module TinyTypes

export TinyDict, TinyVec, TinyRange

import Base: haskey, getindex, iterate, length, size, in, first, last

struct TinyDict{b,T} <: AbstractDict{Int,T}
    bits :: UInt64
end

struct TinyVec <: AbstractVector{Int}
    bits :: UInt64
end

struct TinyRange{n} <: AbstractUnitRange{Int}
    # no fields
end

# TinyDict

getbits(d::TinyDict{b}, k::Int) where {b} =
    (d.bits >> (b*(k-1))) & ((1 << b) - 1)

haskey(d::TinyDict{b,T}, k::Int) where {b,T} =
    getbits(d, k) ≠ 0

getindex(d::TinyDict{b,T}, k::Int) where {b,T} =
    T(getbits(d, k))

function iterate(d::TinyDict{b,T}, k::Int = 1) where {b,T}
    while b*k ≤ 64
        v = getbits(d, k)
        v ≠ 0 && return k => T(v), k+1
        k += 1
    end
end

length(d::TinyDict{b,T}) where {b,T} =
    count(getbits(d, k) ≠ 0 for k = 1:64÷b)

# TinyVec

size(v::TinyVec) = (count_ones(v.bits),)

function getindex(v::TinyVec, i::Int)
    x = 0
    b = v.bits
    for j = 1:i
        s = trailing_zeros(b) + 1
        b >>= s
        x += s
    end
    return x
end

function iterate(v::TinyVec, x::Int = 0)
    b = v.bits >> x
    x += trailing_zeros(b) + 1
    x ≤ 64 ? (x, x) : nothing
end

in(x::Int, v::TinyVec) = (v.bits >> x) ≠ 0

# TinyRange

TinyRange(n::Integer) = TinyRange{n}()

first(r::TinyRange) = 1
last(r::TinyRange{n}) where {n} = n::Int

end # module

using .TinyTypes
