# Every PkgInfo conflicts matrix is allocated with its row dimension padded
# to a whole number of 64-bit words (see padded_rows), so that column j
# occupies exactly the aligned chunk words X.chunks[(j-1)*W .+ (1:W)] where
# W = size(X, 1) >> 6. Version i of the package corresponds to bit i-1 of
# that span; the row after the last version holds the column's active flag;
# any rows beyond that are padding and must remain zero. The helpers below
# are the word-parallel column primitives the filtering passes build on.

# rows to allocate for a package with m versions:
# m version rows plus the flag row, rounded up to whole words
padded_rows(m::Int) = 64 * cld(m + 1, 64)

# number of 64-bit words per column
col_words(X::BitMatrix) = size(X, 1) >> 6

# copy column j of X into buf[1:W]
function col_copy!(buf::Vector{UInt64}, X::BitMatrix, j::Int)
    W = col_words(X)
    base = (j - 1) * W
    chunks = X.chunks
    @inbounds for w = 1:W
        buf[w] = chunks[base + w]
    end
    return buf
end

# does column j of X intersect the row set in buf?
function col_intersects(X::BitMatrix, j::Int, buf::Vector{UInt64})
    W = col_words(X)
    base = (j - 1) * W
    chunks = X.chunks
    @inbounds for w = 1:W
        iszero(chunks[base + w] & buf[w]) || return true
    end
    return false
end

# clear the bits of buf corresponding to rows above m (flag row & padding)
function clear_rows_above!(buf::Vector{UInt64}, m::Int)
    wt = (m + 63) >> 6
    @inbounds for w = wt + 1:length(buf)
        buf[w] = 0
    end
    r = m & 63
    if r != 0 && 1 ≤ wt ≤ length(buf)
        @inbounds buf[wt] &= (UInt64(1) << r) - 1
    end
    return buf
end

# copy len bits from the (word-aligned) column j of X into dest,
# starting at 0-based bit offset off
function copy_col_bits!(dest::Vector{UInt64}, off::Int, X::BitMatrix, j::Int, len::Int)
    W = col_words(X)
    base = (j - 1) * W
    src = X.chunks
    wd = (off >> 6) + 1
    sh = off & 63
    s = 0
    @inbounds while len > 0
        s += 1
        take = min(64, len)
        mask = take == 64 ? ~UInt64(0) : (UInt64(1) << take) - 1
        chunk = src[base + s] & mask
        dest[wd] = (dest[wd] & ~(mask << sh)) | (chunk << sh)
        hi = 64 - sh
        if sh != 0 && take > hi
            dest[wd + 1] = (dest[wd + 1] & ~(mask >> hi)) | (chunk >> hi)
        end
        len -= take
        wd += 1
    end
    return dest
end

# number of set bits in column j of X
function col_count(X::BitMatrix, j::Int)
    W = col_words(X)
    base = (j - 1) * W
    chunks = X.chunks
    n = 0
    @inbounds for w = 1:W
        n += count_ones(chunks[base + w])
    end
    return n
end

# is some active row of X (per the version-flag column) clear in column j?
function col_active_avoids(X::BitMatrix, j::Int)
    W = col_words(X)
    base = (j - 1) * W
    fbase = (size(X, 2) - 1) * W
    chunks = X.chunks
    @inbounds for w = 1:W
        iszero(chunks[fbase + w] & ~chunks[base + w]) || return true
    end
    return false
end

# highest set row ≤ hi in column j of X, or 0 if none
function col_max_upto(X::BitMatrix, j::Int, hi::Int)
    hi ≤ 0 && return 0
    W = col_words(X)
    base = (j - 1) * W
    chunks = X.chunks
    w = (hi + 63) >> 6
    @inbounds c = chunks[base + w]
    r = hi & 63
    r != 0 && (c &= (UInt64(1) << r) - 1)
    while true
        iszero(c) || return ((w - 1) << 6) + (64 - leading_zeros(c))
        w -= 1
        w == 0 && return 0
        @inbounds c = chunks[base + w]
    end
end

# lowest set row in lo:hi of column j of X, or 0 if none
function col_min_from(X::BitMatrix, j::Int, lo::Int, hi::Int)
    lo ≤ 0 && (lo = 1)
    hi < lo && return 0
    W = col_words(X)
    base = (j - 1) * W
    chunks = X.chunks
    wh = (hi + 63) >> 6
    w = ((lo - 1) >> 6) + 1
    @inbounds c = chunks[base + w] & (~UInt64(0) << ((lo - 1) & 63))
    while true
        if w == wh
            r = hi & 63
            r != 0 && (c &= (UInt64(1) << r) - 1)
        end
        iszero(c) || return ((w - 1) << 6) + trailing_zeros(c) + 1
        w == wh && return 0
        w += 1
        @inbounds c = chunks[base + w]
    end
end
