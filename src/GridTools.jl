module GridTools

using Printf
using Statistics
using BenchmarkTools
using Profile
using Debugger
using Base: @propagate_inbounds


import Base.Broadcast: Extruded, Style, BroadcastStyle, ArrayStyle ,Broadcasted


export Cell, K , Edge, E2C, Field, Dimension, Connectivity, neighbor_sum

# Lib --------------------------------------------------------------------------------------------------------------------------------
abstract type Dimension end

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

Base.BroadcastStyle(::Type{<:Field}) = Broadcast.ArrayStyle{Field}()

function Base.similar(bc::Broadcasted{ArrayStyle{Field}}, ::Type{ElType}) where ElType
    output_props = axes(bc)
    
    Field(output_props[1], similar(Array{ElType}, output_props[2]), output_props[3])
end

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

function (field_call::Field)(field_in::Field)::Field
    return Field(field_in.dims, map(x -> x == 0 ? 0 : getindex(field_call.data, Int.(x)), field_in.data))
end

function neighbor_sum(field_in::Field; axis::Dimension)::Field
    dim = findall(x -> x == axis, field_in.dims)[1]
    return Field((field_in.dims[1:dim-1]..., field_in.dims[dim+1:end]...), dropdims(sum(field_in.data, dims=dim), dims=dim)) 
end

# TODO: returns new Field. If to manipulate existing field make field mutable
function broadcast(f::Field, b_dims::D)::Field where D <: Tuple{Vararg{<:Dimension}}
    return Field(f.dims, f.data, b_dims)
end



@inbounds function where(mask::Field, a::Field, b::Field)::Field
    @assert size(mask.data) == size(a.data) == size(b.data)
    return Field(a.dims, ifelse.(mask.data, a.data, b.data))
end

"""
    function where(mask::Field, t1::Tuple, t2::Tuple)
        return map(x -> whereit(mask, x[1], x[2]), zip(t1, t2))
    end

Das da macht was
"""
function where(mask::Field, t1::Tuple, t2::Tuple)
    return map(x -> whereit(mask, x[1], x[2]), zip(t1, t2))
end






# User ------------------------------------------------------------------------------------------------------------------------------

struct Cell_ <: Dimension end
struct K_ <: Dimension end
struct Edge_ <: Dimension end
struct E2C_ <: Dimension end
Cell = Cell_()
K = K_()
Edge = Edge_()
E2C = E2C_()

function run_add(a::Field, b::Field)
    temp = a + b
    return b + temp
end

function nearest_cell_to_edge(cell_values::Field, e2c::Connectivity)
    return cell_values(e2c(1))
end


function sum_adjacent_cells(cells::Field, e2c::Connectivity)
    return neighbor_sum(cells(e2c()), axis=E2C)
end






# Ex --------------------------------------------------------------------------------------------------------------------------------


# a = Field((Cell, K), fill(2.0, (3,3)))
# b = Field((Cell, K), fill(3.0, (3,3)))
# c = Field((Cell, K), fill(5.0, (3,3)))


# result = Field((Cell, K), zeros(3,3))

# result = run_add(a, b)

# @printf "%f + %f = %f ± %f" 3.0 (2.0 + 3.0) mean(result) std(result)
# println()

# # ------------------------------------------------------------------------------------------------------------------------------------

# cell_values = Field((Cell,), [1.0, 1.0, 2.0, 3.0, 5.0, 8.0])


# edge_to_cell_table = [
#     [1  0];
#     [3  0];
#     [3  0];
#     [4  0];
#     [5  0];
#     [6  0];
#     [1  6];
#     [1  2];
#     [2  3];
#     [2  4];
#     [4  5];
#     [5  6]
# ]

# cell_to_edge_table = [
#     [1   7   8];
#     [8   9  10];
#     [2   3   9];
#     [4  10  11];
#     [5  11  12];
#     [6   7  12]
# ]

# E2C_offset_provider = Connectivity(edge_to_cell_table, (Cell,), (Edge, E2C), 2)

# edge_values = nearest_cell_to_edge(cell_values, E2C_offset_provider)
# println(edge_values)

# edge_values = sum_adjacent_cells(cell_values, E2C_offset_provider)

# println(edge_values.data)


# Broadcast ----------------------------------------------------------------------------------------------------------------------------------

# Drops additional dimensions to access array size
@propagate_inbounds _broadcast_getindex(b::Extruded{Field}, i) = b.x[f_newindex(i, b.keeps, b.defaults)]

# TODO: Could use default for index shift. Could do: return (I[1]+ Idefault[1], _newindex(tail(I), tail(keep), tail(Idefault))...)
@inline function f_newindex(I, keep, Idefault)
    if keep[1]
        return (I[1], _newindex(Base.tail(I), Base.tail(keep), Base.tail(Idefault))...)
    else
        return _newindex(Base.tail(I), Base.tail(keep), Base.tail(Idefault))
    end
end


# Builds boolean array of dimensions to keep when iterating over all dimensions
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



# ----------------------------------------------------------------------------------------------------------------------------------
end

