module CompressedSparseArrays

using SparseBase
using SparseBase: getoffset
abstract type AbstractCompressedArray{Order, Bi, Tv, Tfill, Ti, N} <: SparseBase.AbstractSparseArray{Tv, Ti, N} end
SparseBase.comptime_storageorder(::AbstractCompressedArray{O}) where O = O
SparseBase.issparse(::AbstractCompressedArray) = true
# SparseBase.levelformat!!!

Base.eltype(::AbstractCompressedArray{<:Any, <:Any, Tv, Tfill}) where {Tv, Tfill} = Union{Tv, Tfill}
mutable struct SinglyCompressedArray{Order, Bi, Tvalues, Tfill, Tindex<:Integer, 
        V<:AbstractArray{Tvalues}, I1<:AbstractVector{Tindex}, I2<:AbstractVector{Tindex}, N} <: 
        AbstractCompressedArray{Order, Bi, Tvalues, Tfill, Tindex, N}
    # Comments here reflect column major ordering. 
    # To understand as CSR simply swap references to columns with rows and vice versa.
    vlen::Int # m in CSC, n in CSR. The length of the sparse vectors.
    vdim::Int # n in CSC, m in CSR. The number of sparse vectors being stored.
    # This should be the length of ptr.

    ptr::I1{Tindex} # The pointers into i/nzval.
    idx::I2{Tindex} # the stored row(col) indices.
    v::V # the values of the stored indices.
    fill::Tfill # the fill value. A[i,j] == fill if (i,j) is compressed out.
    # COO pending insertions
    # COO pending deletions
    function SinglyCompressedArray{Order, Bi}(vlen, vdim, ptr::I, idx::I, v::V, fill::Tfill) where 
        {Order, Bi, V, I, Tfill}
        N = 2 + ndims(v) - 1
        return new{Order, Bi, eltype(V), Tfill, eltype(I), V, I, N}(vlen, vdim, ptr, idx, v, fill)
    end
end

const CSCMatrix{Bi, Tvalues, V, Tfill, Tindex, I} = SinglyCompressedArray{ColMajor(), Bi, Tvalues, V, Tfill, Tindex, I, 2} where
    {Bi, Tvalues, V<:AbstractVector, Tfill, Tindex, I}
const CSRMatrix{Bi, Tvalues, V, Tfill, Tindex, I} = SinglyCompressedArray{RowMajor(), Bi, Tvalues, V, Tfill, Tindex, I, 2} where
{Bi, Tvalues, V<:AbstractVector, Tfill, Tindex, I}


SparseBase.setfill(A::SinglyCompressedArray{Order, Bi}, f) where {Order, Bi} = 
    SinglyCompressedArray{Order, Bi}(A.vlen, A.vdim, A.ptr, A.idx, A.v, f)
function SparseBase.setfill!(A::SinglyCompressedArray, f)
    A.fill = f
    return A
end

SparseBase.nstored(A::SinglyCompressedArray) = length(A.v)

Base.size(A::SinglyCompressedArray{ColMajor(), <:Any, <:Any, <:AbstractVector}) = (A.vlen, A.vdim)
Base.size(A::SinglyCompressedArray{ColMajor(), <:Any, <:Any, <:AbstractMatrix}) = (A.vlen, A.vdim, size(A.v, 1))
Base.size(A::SinglyCompressedArray{ColMajor(), <:Any}) = (A.vlen, A.vdim, size(A.v)[1:end-1]...)

Base.size(A::SinglyCompressedArray{RowMajor(), <:Any, <:Any, <:AbstractVector}) = (A.vdim, A.vlen)
Base.size(A::SinglyCompressedArray{RowMajor(), <:Any, <:Any, <:AbstractMatrix}) = (A.vdim, A.vlen, size(A.v, 1))
Base.size(A::SinglyCompressedArray{RowMajor(), <:Any}) = (A.vdim, A.vlen, size(A.v)[1:end-1]...)

# Bi setup taken from SparseMatrixCSR.jl
SparseBase.getoffset(::SinglyCompressedArray{O, Bi}) where {O, Bi} = getoffset(Bi)

swapindices(::RowMajor, row, col) = row, col
swapindices(::ColMajor, row, col) = col, row

# indexing adapted from SparseArrays and SparseMatricesCSR:
_indexingtostored(x, ::Any, ::NTuple{N, <:Integer}) where N = x # an iso matrix.
# I wonder, will these checks slow down indexing? 
# The usual CSC/CSR single value per stored index.
_indexintostored(v::AbstractVector, k, ::NTuple{N, <:Integer}) where N = length(v) == 1 ? v[1] : v[k]
# We'll assume that A is stored colexicographic, so k will be the last index.
_indexintostored(A::AbstractArray, k, x::NTuple{N, <:Integer}) where N = length(v) == 1 ? A[1] : 
    A[x[3:end]..., k]

function Base.getindex(A::SinglyCompressedArray{Order, Bi, <:Any, <:Any, <:Any, <:Any, <:Any, N}, x::Vararg{<:Integer, N}) where 
    {Order, Bi, N}
    @boundscheck checkbounds(A, x...) # the overloaded size should take care of doing this correctly.
    row, col = x[1:2]
    i, j = swapindices(Order, row, col)
    SparseBase.nstored(A) == 0 && return A.fill # no values, return fill.
    o = getoffset(A)
    # this is all adapted directly from sparsematrix.jl ∈ SparseArrays.jl
    r1 = @inbounds A.ptr[i] + o
    r2 = A.ptr[i + 1] - Bi
    (r1 > r2) && return A.fill
    jo = j - o
    k = searchsortedfirst(A.idx, jo, r1, r2, Base.Order.Forward)
    return ((k > r2) || (A.idx[k] != jo)) ? A.fill : _indexintostored(A.v, k, x)
end

# Duplication here. Easy to split, but just two functions.
function Base.isstored(A::SinglyCompressedArray{Order, Bi, <:Any, <:Any, <:Any, <:Any, N}, x::Vararg{<:Integer, N}) where 
    {Order, Bi, N}
    @boundscheck checkbounds(A, x...) # the overloaded size should take care of doing this correctly.
    row, col = x[1:2]
    i, j = swapindices(Order, row, col)
    SparseBase.nstored(A) == 0 && return false # no values, return fill.
    o = getoffset(A)
    # this is all adapted directly from sparsematrix.jl ∈ SparseArrays.jl
    r1 = @inbounds A.ptr[i] + o
    r2 = A.ptr[i + 1] - Bi
    (r1 > r2) && return false
    jo = j - o
    k = searchsortedfirst(A.idx, jo, r1, r2, Base.Order.Forward)
    return !((k > r2) || (A.idx[k] != jo)) # no need to check dense obviously
end

function Base.similar(A::SinglyCompressedArray{Order, Bi, <:Any, <:Any, <:Any, <:Ti}, ::Type{S}, dims::Dims) where
    {Order, Bi, Ti, S}
    vdim = Order === ColMajor() ? dims[2] : dims[1]
    v = if length(dims) == 2
        similar(A.v, S, (0,))
    elseif length(dims) > 2
        similar(A.v, S, (dims[3:end]..., 0)) # this isn't terribly useful I don't think.
    else
        throw(ArgumentError("length(dims) must be >= 2"))
    end
    return SinglyCompressedArray{Order, Bi}(
        Order === ColMajor() ? dims[1] : dims[2], vdim,
        ones(Ti, vdim), # TODO: this is wrong, doesn't account for type of A vectors
        Vector{Ti}(), # TODO: this is wrong, doesn't account for type of A vectors
        v,
        A.fill
    )
end

## setindex!
end
