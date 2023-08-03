module GridTools

using Printf
using Statistics
using BenchmarkTools
using Profile
using Debugger
using Base: @propagate_inbounds

import Base.Broadcast: Extruded, Style, BroadcastStyle, ArrayStyle ,Broadcasted

export Cell, K , Edge, E2C, Field, Dimension, Connectivity, neighbor_sum, where

# Lib ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

abstract type Dimension end

# Field struct --------------------------------------------------------------------

# TODO: check for #dimension at compile time and not runtime
# TODO: <: AbstractArray{T,N} is not needed... but then we have to define our own length and iterate function for Fields
struct Field{T, N, T2 <: Tuple{Vararg{<:Dimension}}, T3 <: Tuple{Vararg{<:Dimension}}} <: AbstractArray{T,N}
    dims::T2
    data::Array{T,N}
    broadcast_dims::T3

    function Field(dims::T2, data::Array{T,N}, broadcast_dims::T3 = dims) where {T, N, T2 <: Tuple{Vararg{<:Dimension}}, T3 <: Tuple{Vararg{<:Dimension}}}
        @assert length(dims) == ndims(data)
        return new{T,N,T2,T3}(dims, data, broadcast_dims)
    end
end

Base.size(F::Field)::Tuple = size(F.data)
@propagate_inbounds Base.getindex(F::Field{T,N}, inds::Vararg{Int,N}) where {T,N} = F.data[inds...]
@propagate_inbounds Base.setindex!(F::Field{T,N}, val, inds::Vararg{Int,N}) where {T,N} = F.data[inds...] = val
Base.showarg(io::IO, F::Field, toplevel) = print(io, " Field with dims ", F.dims, " and broadcasted_dims ", F.broadcast_dims)

function (field_call::Field)(field_in::Field)::Field
    return Field(field_in.dims, map(x -> x == 0 ? 0 : getindex(field_call.data, Int.(x)), field_in.data))
end

# TODO: returns new Field. If to manipulate existing field make field mutable
function broadcast(f::Field, b_dims::D)::Field where D <: Tuple{Vararg{<:Dimension}}
    return Field(f.dims, f.data, b_dims)
end

# Connectivity struct ------------------------------------------------------------

struct Connectivity
    data::Array{Int64, 2}
    source::Tuple{Vararg{<:Dimension}}
    target::Tuple{Vararg{<:Dimension}}
    dims::Int64
end

function (conn_call::Connectivity)(neighbor::Int64 = -1)::Field
    if neighbor == -1
        return Field(conn_call.target, conn_call.data)
    else
        @assert conn_call.dims >= neighbor
        return Field((conn_call.target[neighbor],), conn_call.data[:, neighbor])
    end
end

# Built-ins ----------------------------------------------------------------------

function neighbor_sum(field_in::Field; axis::Dimension)::Field
    dim = findall(x -> x == axis, field_in.dims)[1]
    return Field((field_in.dims[1:dim-1]..., field_in.dims[dim+1:end]...), dropdims(sum(field_in.data, dims=dim), dims=dim)) 
end


@inbounds where(mask::Field, a::Field, scal::Real)::Field = Field(a.dims, ifelse(mask.data, a.data, scal))
@inbounds where(mask::Field, scal::Real, a::Field)::Field = Field(a.dims, ifelse(mask.data, a.data, scal))
@inbounds where(mask::Field, a::Field, b::Field)::Field = Field(a.dims, ifelse.(mask.data, a.data, b.data))
"""
    where(mask, tuple1, tupl2)

This function takes 3 input arguments:
    mask: a field with dtype boolean
    true branch: a tuple, a field, or a scalar
    false branch: a tuple, a field, or a scalar
The where loops over each entry of the mask and returns values corresponding to the same indexes of either the true or the false branch.
"""
where(mask::Field, t1::Tuple, t2::Tuple)::Field = map(x -> whereit(mask, x[1], x[2]), zip(t1, t2))


# Broadcast ---------------------------------------------------------------------

# Create custom broadcast style
Base.BroadcastStyle(::Type{<:Field}) = Broadcast.ArrayStyle{Field}()

# Custom similar(): Creates output object
function Base.similar(bc::Broadcasted{ArrayStyle{Field}}, ::Type{ElType}) where ElType
    output_props = axes(bc)
    
    Field(output_props[1], similar(Array{ElType}, output_props[2]), output_props[3])
end

# Custom _broadcast_getindex(): Drops additional dimensions to access array size
@propagate_inbounds _broadcast_getindex(b::Extruded{Field}, i) = b.x[f_newindex(i, b.keeps, b.defaults)]

# TODO: Could use default for index shift. Could do: return (I[1]+ Idefault[1], _newindex(tail(I), tail(keep), tail(Idefault))...)
@inline function f_newindex(I, keep, Idefault)
    if keep[1]
        return (I[1], _newindex(Base.tail(I), Base.tail(keep), Base.tail(Idefault))...)
    else
        return _newindex(Base.tail(I), Base.tail(keep), Base.tail(Idefault))
    end
end


# Custom newindexer(): Builds boolean array of dimensions to keep when iterating over all dimensions
@inline function newindexer(A::Field)
    return f_newindexer(A.dims, A.broadcast_dims, axes(A))
end

# Idefault not used atm. Use for indexshift
@inline f_newindexer(dims::Tuple, b_dims::Tuple{}, ax::Tuple) = (), ()
@inline function f_newindexer(dims::Tuple, b_dims::Tuple, ax::Tuple)
    ind1 = b_dims[1]
    keep, Idefault = f_newindexer(dims, Base.tail(b_dims), ax)
    next_keep = ind1 in dims
    (next_keep, keep...) , (next_keep ? ax[findall(x -> x == ind1, dims)][1][1] : 0, Idefault...)
end

# Custom instantiate(): Dimension check and calculation of output dimension
function Base.Broadcast.instantiate(bc::Broadcasted{ArrayStyle{Field}})
    return Broadcasted{ArrayStyle{Field}}(bc.f, bc.args, axes(bc)[2])
end

# Helper function for combine_axes
@inline get_dim_ind(dims::Tuple{}, b_dims::Tuple) =  ()
@inline get_dim_ind(dims::Tuple, b_dims::Tuple) = [findall(x -> x == dims[1], b_dims)[1] , get_dim_ind(Base.tail(dims), b_dims)...]

# Checks dimension and broadcast combatibility of all fields
@inline Base.axes(bc::Broadcasted{ArrayStyle{Field}}) =              combine_axes(bc.args...)
@inline combine_axes(A::Field, bc::Broadcasted{ArrayStyle{Field}}) = combine_axes((A.dims, axes(A), A.broadcast_dims), axes(bc))
@inline combine_axes(bc::Broadcasted{ArrayStyle{Field}}, B::Field) = combine_axes(axes(bc), (B.dims, axes(B), B.broadcast_dims))
@inline combine_axes(A::Field, B::Field) =                           combine_axes((A.dims, axes(A), A.broadcast_dims), (B.dims, axes(B), B.broadcast_dims))
@inline combine_axes(A::Field, t::Tuple{}) = (A.dims, axes(A), A.broadcast_dims)
@inline combine_axes(t::Tuple{}, A::Field) = (A.dims, axes(A), A.broadcast_dims)
@inline function combine_axes(A::Tuple, B::Tuple)
    length(A[1]) > length(B[1]) ? (A,B) = (B,A) : nothing  # swap A and B       
    
    # A,B = (dims, axes(data), broadcast_dims)
    @assert issubset(A[1], B[1]) "Dimension Mismatch between the Dimensions of the two Fields"
    @assert issubset(B[1], A[3]) "Dimension Mismatch between the broadcasted Dimensions of the two Fields"
    matching_dims = get_dim_ind(A[1], B[1])
    @assert A[2] == B[2][matching_dims] "Dimension Mismatch between the data Dimensions of the two Fields"

    return B
end

# User ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

struct Cell_ <: Dimension end
struct K_ <: Dimension end
struct Edge_ <: Dimension end
struct E2C_ <: Dimension end
Cell = Cell_()
K = K_()
Edge = Edge_()
E2C = E2C_()


end

