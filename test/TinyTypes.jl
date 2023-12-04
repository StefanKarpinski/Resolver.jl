module TinyTypes

export TinyDict, TinyVec

import Base: getindex, iterate, length, size, in

struct TinyDict{b,T} <: AbstractDict{Int,T}
    bits :: UInt64
end

getbits(d::TinyDict{b}, k::Int) where {b} =
    (d.bits >> (b*(k-1))) & ((1 << b) - 1)

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

struct TinyVec <: AbstractVector{Int}
    bits :: UInt64
end

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

function iterate(v::TinyVec, x::Int = 1)
    b = v.bits >> (x - 1)
    x += trailing_zeros(b)
    x ≤ 64 ? (x, x) : nothing
end

in(x::Int, v::TinyVec) = (v.bits >> x) ≠ 0

end # module

using .TinyTypes
