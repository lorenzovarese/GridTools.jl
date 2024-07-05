
# __precompile__(false)
module GridTools

using Printf
using Statistics
using BenchmarkTools
using Profile
using Base: @propagate_inbounds
using MacroTools
using OffsetArrays: IdOffsetRange
using Debugger

import Base.Broadcast: Extruded, Style, BroadcastStyle, ArrayStyle, Broadcasted

export Dimension,
    DimensionKind,
    HORIZONTAL,
    VERTICAL,
    LOCAL,
    Field,
    Connectivity,
    FieldOffset,
    neighbor_sum,
    max_over,
    min_over,
    where,
    concat,
    @field_operator,
    slice,
    copyfield!,
    get_dim_name,
    get_dim_kind,
    get_dim_ind

const SKIP_NEIGHBOR_INDICATOR = -1 # TODO(tehrengruber): move from atlas submodule here

# ========================================
# =============== Dimension ==============
# ========================================

abstract type DimensionKind end

struct HORIZONTAL <: DimensionKind end
struct VERTICAL <: DimensionKind end
struct LOCAL <: DimensionKind end

"""
    struct Dimension{name, kind <: DimensionKind}

Create a new `Dimension` with a specified name and kind. This struct is parameterized to allow for flexible 
representation of various dimensions, typically used in grid and spatial data structures.

## Type Parameters
- `name`: A symbol representing the name of the dimension, uniquely identifying it.
- `kind`: A subtype of `DimensionKind` which defaults to `HORIZONTAL` if not specified, indicating the 
          dimension's orientation or scope (e.g., `HORIZONTAL`, `VERTICAL`, `LOCAL`).

By default, dimensions are horizontal unless otherwise specified. This can be customized by providing 
a specific kind during creation.

## Examples
```julia
julia> Cell_ = Dimension(:Cell_)  # Defaults to horizontal
julia> Cell = Cell_()

julia> V2EDim_ = Dimension(:V2E_, LOCAL)  # Explicitly specifying local kind
julia> V2EDim = V2EDim_()
```
"""
struct Dimension{name, kind <: DimensionKind} end

# Create a new `Dimension` instance with a specified symbol and kind.
function Dimension(sym::Symbol, kind = HORIZONTAL)
    return Dimension{sym, kind}
end

# Define the length of a Dimension object as always 1. This treats Dimension as a scalar entity.
Base.length(d::Dimension)::Int64 = 1

# Enable iteration over a Dimension object. Allows a Dimension to be used once in a loop.
Base.iterate(d::Dimension, state = 1) = state == 1 ? (d, state + 1) : nothing

# Retrieve the name of the Dimension, removing the trailing underscore.
get_dim_name(d::Dimension) = string(typeof(d).parameters[1])[1:end-1]

# Obtain the kind of the Dimension, which indicates its type (e.g., HORIZONTAL, VERTICAL, LOCAL).
get_dim_kind(d::Dimension) = typeof(d).parameters[2]

# Find the first index of a specific Dimension object within a tuple of Dimensions.
get_dim_ind(source::Tuple{Vararg{Dimension}}, ind::Dimension) =
    findfirst(x -> x == ind, source)

# ========================================
# ============= FieldOffset ==============
# ========================================

"""
    FieldOffset(name::String, source::Dimension, target::Union{Dimension, Tuple{Vararg{Dimension}}})

The `FieldOffset` struct represents a transformation that can be applied to fields (or tuples of fields) over 
one domain to another domain. This transformation leverages the connectivity between the source and target domains 
to find the values of adjacent mesh elements.

## Arguments
- `name::String`: A name for the field offset.
- `source::Dimension`: The source dimension from which the transformation originates.
- `target::Union{Dimension, Tuple{Vararg{Dimension}}}`: The target dimension or a tuple of dimensions to which the 
   transformation is applied. If `target` is a tuple with more than one dimension, all but the first dimension must 
   be local dimensions.

## Returns
- A `FieldOffset` instance that represents the transformation from the source dimension to the target dimension(s).

## Examples
```julia-repl
julia> V2E = FieldOffset("V2E", source=Edge, target=(Vertex, V2EDim))
FieldOffset("V2E", Dimension{:Edge_, HORIZONTAL}(), (Dimension{:Vertex_, HORIZONTAL}(), Dimension{:V2E_, LOCAL}()))
```
"""
struct FieldOffset
    name::String
    source::Dimension
    target::Tuple{Vararg{Dimension}}

    function FieldOffset(
        name::String;
        source::Dimension,
        target::Union{Dimension, Tuple{Vararg{Dimension}}}
    )::FieldOffset
        if length(target) > 1
            @assert all(get_dim_kind.(Base.tail(target)) .== LOCAL) (
                "All but the first dimension in an offset must be local dimensions."
            )
        end
        return new(name, source, Tuple(target))
    end
end

# Returns a tuple of the FieldOffset instance and the provided index.
Base.getindex(f_off::FieldOffset, ind::Integer) = f_off, ind

# ========================================
# ================ Field =================
# ========================================

# TODO(tehrengruber): expand docstring for binary ops, e.g. +
"""
    Field{B_Dim, T, N, Dim, D} <: AbstractArray{T, N}

Represents a multidimensional data structure defined over a set of named dimensions (`dims`), capable of storing 
numeric data (`data`). This structure supports advanced indexing and data manipulation features, making it 
suitable for complex data analysis and computational science tasks.

## Fields
- `dims::Dim`: Tuple of `Dimension` objects defining the axes of the field.
- `data::D`: The multi-dimensional array storing the actual data.
- `broadcast_dims::B_Dim`: Tuple of dimensions used for broadcasting operations to align with Julia's broadcasting behavior.
- `origin::NTuple{N, Int64}`: Specifies the index origin for each dimension, allowing for offset indexing.

## Parameters
- `dims`: A tuple specifying the dimensions of the field.
- `data`: An array containing the field's data.
- `broadcast_dims` (optional): Dimensions used for broadcasting; defaults to `dims` if not specified.
- `origin`: A dictionary mapping dimensions to their respective origins for indexing, allowing shifting 
   of the data's logical view.

## Examples
Create a simple field:
```julia
julia> new_field = Field((Cell, K), fill(1.0, (3,2)))
3×2 Float64 Field with dimensions ("Cell", "K") with indices 1:3×1:2:
 1.0  1.0
 1.0  1.0
 1.0  1.0
```

Offset a dimension:
Each dimension can be offset by a certain index. This changes the index-range of the Field.
```julia
julia> new_field = Field((Cell, K), fill(1.0, (3,2)), origin = Dict(K => 1))
3×2 Float64 Field with dimensions ("Cell", "K") with indices 1:3×2:3:
 1.0  1.0
 1.0  1.0
 1.0  1.0
```

Fields also have a call operator. You can transform fields (or tuples of fields) over one domain to another 
domain by using the call operator of the source field with a field offset as argument.
As an example, you can use the field offset E2C below to transform a field over cells to a field over edges 
using edge-to-cell connectivities. The FieldOffset can take an index as argument which restricts the output 
dimension of the transform.

The call itself must happen in a field_operator since it needs the functionality of an offset_provider.

## Examples
```julia
julia> E2C = FieldOffset("E2C", source=Cell, target=(Edge, E2CDim))
julia> field = Field(Cell, ones(5))

julia> field(E2C)
julia> field(E2C[1])
```
"""
struct Field{
    B_Dim <: Tuple{Vararg{Dimension}},
    T <: Union{AbstractFloat, Integer},
    N,
    Dim <: NTuple{N, Dimension},
    D <: AbstractArray{T, N}
} <: AbstractArray{T, N}
    dims::Dim
    data::D
    broadcast_dims::B_Dim
    origin::NTuple{N, Int64}
end

# TODO: check for #dimension at compile time and not runtime
# Ensure the number of dimensions in 'dims' matches the number of dimensions in 'data'
function Field(
    dims::Dim,
    data::D,
    broadcast_dims::B_Dim = dims;
    origin::Dict{<:Dimension, Int64} = Dict{Dimension, Int64}()
) where {
    T <: Union{AbstractFloat, Integer},
    N,
    B_Dim <: Tuple{Vararg{Dimension}},
    Dim <: NTuple{N, Dimension},
    D <: AbstractArray{T, N}
}
    if ndims(data) != 0
        @assert length(dims) == ndims(data)
    end
    offsets = Tuple([get(origin, dim, 0) for dim in dims])
    return Field(dims, data, broadcast_dims, offsets)
end

# Verify that 'data' is one-dimensional if it is not empty
function Field(
    dim::Dimension,
    data::D,
    broadcast_dims::Union{Dimension, B_Dim} = dim;
    origin::Dict{<:Dimension, Int64} = Dict{Dimension, Int64}()
) where {
    T <: Union{AbstractFloat, Integer},
    N,
    B_Dim <: NTuple{N, Dimension},
    D <: AbstractArray{T, N}
}
    if ndims(data) != 0
        @assert ndims(data) == 1
    end
    return Field(Tuple(dim), data, Tuple(broadcast_dims), origin = origin)
end

# A type-stable struct representing a transformation between source and target dimensions.
struct FieldOffsetTS{
    Name,
    Source <: Dimension,
    Target <: Union{Tuple{<:Dimension, <:Dimension}, Tuple{<:Dimension}}
} end

"""
    to_type_stable_field_offset(offset::FieldOffset) -> FieldOffsetTS

Converts a `FieldOffset` instance to a type-stable `FieldOffsetTS` instance.

## Arguments
- `offset::FieldOffset`: The `FieldOffset` instance to be converted.

## Returns
- A `FieldOffsetTS` instance with type parameters derived from the `offset`.
"""
to_type_stable_field_offset(offset::FieldOffset) = FieldOffsetTS{
    Symbol(offset.name),
    typeof(offset.source),
    Tuple{map(typeof, offset.target)...}
}()

"""
    struct AllNeighbors

The `AllNeighbors` struct is used to indicate that all neighboring positions should be considered 
in a remapping operation. This can be useful in functions where a complete traversal or interaction 
with all neighboring elements is required.

## Examples
```julia
julia> nb_ind = AllNeighbors()
AllNeighbors()
```
"""
struct AllNeighbors end

"""
    remap_position(current_position::Tuple{Int64, Vararg{Int64}},
                   dims::Tuple{Dimension, Vararg{Dimension}},
                   offset::FieldOffsetTS{OffsetName, SourceDim, Tuple{TargetDim, TargetLocalDim}},
                   nb_ind::AllNeighbors,
                   conn) -> Tuple{Bool, Tuple{Int64, Vararg{Int64}}}

Remaps the given `current_position` from the `SourceDim` to the `TargetDim` using the provided `offset` and connectivity information. This function recursively processes the dimensions and their respective positions.

## Arguments
- `current_position::Tuple{Int64, Vararg{Int64}}`: The current position as a tuple of indices.
- `dims::Tuple{Dimension, Vararg{Dimension}}`: The dimensions corresponding to the current position.
- `offset::FieldOffsetTS{OffsetName, SourceDim, Tuple{TargetDim, TargetLocalDim}}`: The field offset indicating the transformation from the source dimension to the target dimensions.
- `nb_ind::AllNeighbors`: An indicator to consider all neighboring positions.
- `conn`: The connectivity information used for remapping. This should be of a type that provides connectivity data (e.g., `conn.data`).

## Returns
- A tuple `(position_exists::Bool, new_position::Tuple{Int64, Vararg{Int64}})`, where `position_exists` indicates if the remapped position is valid and `new_position` is the remapped position tuple.

## Examples
```julia
julia> current_position = (1, 2, 3)
julia> dims = (Dimension{:Edge_, HORIZONTAL}(), Dimension{:Vertex_, HORIZONTAL}(), Dimension{:V2E_, LOCAL}())
julia> offset = FieldOffsetTS{:V2E, Dimension{:Edge_, HORIZONTAL}, Tuple{Dimension{:Vertex_, HORIZONTAL}, Dimension{:V2E_, LOCAL}}}()
julia> nb_ind = AllNeighbors()
julia> conn = Connectivity([ [1,2,3], [4,5,6] ])  # Assuming `Connectivity` is a struct with appropriate `data`
julia> remap_position(current_position, dims, offset, nb_ind, conn)
(true, (new_index, new_tail_position...))
```

## Implementation Details
- The function checks if the first dimension in dims matches the TargetDim. If it matches, it performs a lookup in conn using the index and neighbor index from current_position.
- If the first dimension does not match the TargetDim, it processes the current position and dimensions recursively.
- The function ensures that if new_ind is not a SKIP_NEIGHBOR_INDICATOR and the tail position exists, the remapped position is considered valid.
""" # TODO: move to embedded.jl
function remap_position(
    current_position::Tuple{Int64, Vararg{Int64}},
    dims::Tuple{Dimension, Vararg{Dimension}},
    offset::FieldOffsetTS{OffsetName, SourceDim, Tuple{TargetDim, TargetLocalDim}},
    nb_ind::AllNeighbors,
    conn
) where {
    OffsetName,
    SourceDim <: Dimension,
    TargetDim <: Dimension,
    TargetLocalDim <: Dimension
}  # TODO: restrict conn to type ::Connectivity
    if dims[1] == TargetDim()  # since we are mapping indices not field here the target source are flipped
        ind, actual_nb_ind, tail_position... = current_position
        _, local_dim, tail_dims... = dims
        #@assert local_dim isa Dimension && string(typeof(local_dim).parameters[1]) == (string(OffsetName)*"_") && typeof(local_dim).parameters[2] == LOCAL
        new_ind = conn.data[ind, actual_nb_ind]
    else
        new_ind, tail_position... = current_position
        dim, tail_dims... = dims
    end

    tail_position_exists, new_tail_position =
        remap_position(tail_position, tail_dims, offset, nb_ind, conn)
    position_exists = (new_ind != SKIP_NEIGHBOR_INDICATOR) && tail_position_exists
    return position_exists, (new_ind, new_tail_position...)
end

# Remaps the position with a neighbor index specified as an integer.
function remap_position(
    current_position::Tuple{Int64, Vararg{Int64}},
    dims::Tuple{Dimension, Vararg{Dimension}},
    offset::FieldOffsetTS{OffsetName, SourceDim, Tuple{TargetDim, TargetLocalDim}},
    nb_ind::Int64,
    conn
) where {
    OffsetName,
    SourceDim <: Dimension,
    TargetDim <: Dimension,
    TargetLocalDim <: Dimension
}  # TODO: restrict conn to type ::Connectivity
    if dims[1] == TargetDim()  # since we are mapping indices not field here the target source are flipped
        ind, tail_position... = current_position
        _, tail_dims... = dims
        new_ind = conn.data[ind, nb_ind]
    else
        new_ind, tail_position... = current_position
        dim, tail_dims... = dims
    end

    tail_position_exists, new_tail_position =
        remap_position(tail_position, tail_dims, offset, nb_ind, conn)
    position_exists = (new_ind != SKIP_NEIGHBOR_INDICATOR) && tail_position_exists
    return position_exists, (new_ind, new_tail_position...)
end

# Override for handling the base case where `current_position` and `dims` are empty.
remap_position(
    current_position::Tuple{},
    dims::Tuple{},
    offset::FieldOffsetTS,
    nb_ind::Union{Int64, AllNeighbors},
    conn
) = (true, ())

function compute_remapped_field_info(
    field_size::Tuple{Int64, Vararg{Int64}},
    dims::Tuple{Dimension, Vararg{Dimension}},
    offset::FieldOffsetTS{OffsetName, SourceDim, Tuple{TargetDim, TargetLocalDim}},
    nb_ind::Union{Int64, AllNeighbors},
    conn
) where {
    OffsetName,
    SourceDim <: Dimension,
    TargetDim <: Dimension,
    TargetLocalDim <: Dimension
}
    length, tail_size... = field_size
    dim, tail_dims... = dims
    if dim == SourceDim()
        if nb_ind == AllNeighbors()
            new_dims_part = (TargetDim(), TargetLocalDim())
            # TODO: restrict indices to only what is needed
            new_size_part = size(conn.data)  # we just take the size of the connecitvity
        else
            new_dims_part = (TargetDim(),)
            # TODO: restrict indices to only what is needed
            new_size_part = size(conn.data)[1]  # we just take the size of the connecitvity
        end
    else
        new_dims_part = (dim,)
        new_size_part = (length,)
    end

    new_tail_dims, new_tail_size =
        compute_remapped_field_info(tail_size, tail_dims, offset, nb_ind, conn)
    new_dims = (new_dims_part..., new_tail_dims...)
    new_size = (new_size_part..., new_tail_size...)
    return new_dims, new_size
end

compute_remapped_field_info(
    field_size::Tuple{},
    dims::Tuple{},
    offset::FieldOffsetTS,
    nb_ind::Union{Int64, AllNeighbors},
    conn
) = ((), ())

function remap_broadcast_dims(
    broadcast_dims::Tuple{T, Vararg{Dimension}},
    offset::FieldOffsetTS{OffsetName, SourceDim, Tuple{TargetDim, TargetLocalDim}},
    nb_ind::Union{Int64, AllNeighbors}
) where {
    T <: Dimension,
    OffsetName,
    SourceDim <: Dimension,
    TargetDim <: Dimension,
    TargetLocalDim <: Dimension
}
    dim, dim_tail... = broadcast_dims
    if dim == SourceDim()
        if nb_ind == AllNeighbors()
            (
                TargetDim(),
                TargetLocalDim(),
                remap_broadcast_dims(dim_tail, offset, nb_ind)...
            )
        else
            (TargetDim(), remap_broadcast_dims(dim_tail, offset, nb_ind)...)
        end
    else
        (dim, remap_broadcast_dims(dim_tail, offset, nb_ind)...)
    end
end

remap_broadcast_dims(
    broadcast_dims::Tuple{},
    offset::FieldOffsetTS,
    nb_ind::Union{Int64, AllNeighbors}
) = ()

function remap_ts(
    field::Field,
    offset::FieldOffsetTS{OffsetName, SourceDim, Tuple{TargetDim, TargetLocalDim}},
    nb_ind::Union{Int64, AllNeighbors} = AllNeighbors()
)::Field where {
    OffsetName,
    SourceDim <: Dimension,
    TargetDim <: Dimension,
    TargetLocalDim <: Dimension
}
    conn = OFFSET_PROVIDER[string(OffsetName)]

    # compute new indices
    out_field_dims, out_field_size =
        compute_remapped_field_info(size(field.data), field.dims, offset, nb_ind, conn)
    #out_field = map(position -> begin
    #    neighbor_exists, new_position = remap_position(Tuple(position), out_field_dims, offset, nb_ind, conn)
    #    if neighbor_exists
    #        field.data[new_position...]
    #    else
    #        eltype(field.data)(0)
    #    end
    #end, CartesianIndices(map(len -> Base.OneTo(len), out_field_size)))
    out_field = zeros(eltype(field.data), out_field_size)
    for position in eachindex(IndexCartesian(), out_field)
        neighbor_exists, new_position =
            remap_position(Tuple(position), out_field_dims, offset, nb_ind, conn)
        if neighbor_exists
            out_field[position] = field.data[new_position...]
        end
    end
    # todo: origin
    return Field(
        out_field_dims,
        out_field,
        remap_broadcast_dims(field.broadcast_dims, offset, nb_ind)
    )
end

(field::Field)(f_off::Tuple{FieldOffset, <:Integer})::Field = field(f_off...)
function (field::Field)(
    f_off::FieldOffset,
    nb_ind::Union{Int64, AllNeighbors} = AllNeighbors()
)::Field
    result = remap_ts(field, to_type_stable_field_offset(f_off), nb_ind)
    return result
end

# Field struct interfaces
Base.axes(F::Field)::Tuple = map((i, j) -> i .+ j, axes(F.data), F.origin)
Base.size(F::Field)::Tuple = size(F.data)
Base.convert(t::Type{T}, F::Field) where {T <: Number} =
    Field(F.dims, convert.(t, F.data), F.broadcast_dims)
@propagate_inbounds function Base.getindex(
    F::Field{BD, T, N},
    inds::Vararg{Int, N}
) where {BD, T, N}
    new_inds = inds .- F.origin
    return F.data[new_inds...]
end
@propagate_inbounds function Base.setindex!(
    F::Field{BD, T, N},
    val,
    inds::Vararg{Int, N}
) where {BD, T, N}
    new_inds = inds .- F.origin
    F.data[new_inds...] = val
end
Base.showarg(io::IO, @nospecialize(F::Field), toplevel) =
    print(io, eltype(F), " Field with dimensions ", get_dim_name.(F.broadcast_dims))
function slice(F::Field, inds...)::Field
    dim_ind = findall(x -> typeof(x) <: UnitRange{Int64}, inds)
    return Field(F.dims[dim_ind], view(F.data, inds...), F.broadcast_dims)
end

# ========================================
# ============= Connectivity =============
# ========================================

"""
    Connectivity(data::Array, source::Tuple, target::Tuple, dims::Int)

# Examples
```julia
julia> new_connectivity = Connectivity(fill(1, (3,2)), Cell, Edge, 2)
Connectivity(Integer[1 1; 1 1; 1 1], Dimension{:Cell_, HORIZONTAL}(), Dimension{:Edge_, HORIZONTAL}(), 2
```
"""
struct Connectivity
    data::Array{Int64, 2}
    source::Dimension
    target::Dimension
    dims::Integer
end

# ========================================
# ============ Field Operator ============
# ========================================

struct FieldOp
    name::Symbol
    f::Function
    expr::Expr
    closure_vars::Dict
end

# Includes ------------------------------------------------------------------------------------

include("embedded/builtins.jl")
include("embedded/cust_broadcast.jl")
include("gt2py/gt2py.jl")

# Helper functions for fields -------------------------------------------------------------------

"""
    copyfield!(target::Field, source::Field)

Takes a copies the data from a source field to a target fields (equivalent to target .= source). Also works with tuples of fields.

# Examples
```julia-repl
julia> copyfields!((a, b), (c, d))
```
"""
function copyfield!(target::Tuple{Vararg{Field}}, source::Tuple{Vararg{Field}})
    for i = 1:length(target)
        target[i] .= source[i]
    end
end
copyfield!(target, source) = target .= source

# Field operator functionalities ------------------------------------------------------------

OFFSET_PROVIDER::Union{Dict{String, <:Union{Connectivity, Dimension}}, Nothing} = nothing
FIELD_OPERATORS::Dict{Symbol, PyObject} = Dict{Symbol, PyObject}()

function (fo::FieldOp)(
    args...;
    offset_provider::Dict{String, <:Union{Connectivity, Dimension}} = Dict{
        String,
        Union{Connectivity, Dimension}
    }(),
    backend::String = "embedded",
    out = nothing,
    kwargs...
)

    is_outermost_fo = isnothing(OFFSET_PROVIDER)
    if is_outermost_fo
        @assert !isnothing(out) "Must provide an out field."
        @assert typeof(out) <: Field || typeof(out) <: Tuple{Vararg{Field}} "Out argument is not a field."
        global OFFSET_PROVIDER = offset_provider
        out = backend_execution(
            Val{Symbol(backend)}(),
            fo,
            args,
            kwargs,
            out,
            is_outermost_fo
        )
        global OFFSET_PROVIDER = nothing
    else
        # TODO(tehrengruber): this breaks when fo execution fails and no cleanup is done. use try finally block
        #@assert isnothing(out)
        #@assert isempty(offset_provider)
        out = backend_execution(
            Val{Symbol(backend)}(),
            fo,
            args,
            kwargs,
            out,
            is_outermost_fo
        )
    end

    return out
end

function backend_execution(
    backend::Val{:embedded},
    fo::FieldOp,
    args,
    kwargs,
    out,
    is_outermost_fo
)
    if is_outermost_fo
        copyfield!(out, fo.f(args...; kwargs...))
        return
    else
        return fo.f(args...; kwargs...)
    end
end

function backend_execution(
    backend::Val{:py},
    fo::FieldOp,
    args,
    kwargs,
    out,
    is_outermost_fo
)
    if haskey(FIELD_OPERATORS, fo.name)
        f = FIELD_OPERATORS[fo.name]
    else
        f = py_field_operator(fo)
        FIELD_OPERATORS[fo.name] = f
    end
    p_args, p_kwargs, p_out, p_offset_provider =
        py_args.((args, kwargs, out, GridTools.OFFSET_PROVIDER))
    if is_outermost_fo
        f(p_args..., out = p_out, offset_provider = p_offset_provider; p_kwargs...)
        return
    else
        return f(p_args...; p_kwargs...)
    end
end

function get_closure_vars(expr::Expr, current_vars::Dict)::Dict
    expr_def = splitdef(expr)
    @assert all(typeof.(expr_def[:args]) .== Expr) &&
            all(typeof.(expr_def[:kwargs]) .== Expr) (
        "Field operator parameters must be type annotated."
    )

    local_vars = Set()
    closure_names = Set()
    closure_vars = Dict()

    # catch all local variables
    postwalk(expr) do x
        if @capture(x, (name_ = value_) | (name_::type_))
            if typeof(name) == Symbol
                push!(local_vars, name)
            elseif typeof(name) == Expr
                if name.head == :tuple
                    push!(local_vars, name.args...)
                elseif name.head == :call
                    push!(local_vars, name.args[1])
                else
                    throw(
                        "For the following local variable: $name we dont provide support yet. Please report."
                    ) # TODO: verify
                end
            end
        end
        return x
    end

    # catch all dimensions
    postwalk(expr.args[1]) do x
        if typeof(x) == Symbol &&
           x in keys(current_vars) &&
           typeof(current_vars[x]) == DataType &&
           current_vars[x] <: Dimension
            push!(closure_names, x)
        end
    end

    # catch all closure_variables
    prewalk(expr.args[2]) do x
        if @capture(x, name_(args__)) && !(name in local_vars) && !(name in math_ops)
            push!(closure_names, name)
            return Expr(:irgendoeppis, args...)  # small present for tehrengruber
        elseif typeof(x) == Symbol && !(x in local_vars) && !(x in math_ops)
            push!(closure_names, x)
            return x
        else
            return x
        end
    end

    # update dictionary
    for name in closure_names
        closure_vars[name] = current_vars[name]
    end

    return closure_vars
end

macro module_vars()
    return esc(
        quote
            # TODO(tehrengruber): for some reasons this was needed from some point on. cleanup
            base_vars = Dict(
                name => Core.eval(Base, name) for
                name in [:Int64, :Int32, :Float32, :Float64]
            )
            module_vars =
                Dict(name => Core.eval(@__MODULE__, name) for name in names(@__MODULE__))
            local_vars = Base.@locals
            merge(base_vars, module_vars, local_vars, GridTools.builtin_op)
        end
    )
end

"""
    @field_operator

The field_operator macro takes a function definition and creates a run environment for the function call within the GridTools package. It enables the additional argument "offset_provider", "backend" and "out.

# Examples
```julia-repl
julia> @field_operator addition(x::Int64) = x + x
```
"""
macro field_operator(expr::Expr)
    f_name = namify(expr)

    expr_dict = splitdef(expr)
    expr_dict[:name] = generate_unique_name(f_name)
    unique_expr = combinedef(expr_dict)

    return Expr(
        :(=),
        esc(f_name),
        :(FieldOp(
            namify($(Expr(:quote, expr))),
            $(esc(unique_expr)),
            $(Expr(:quote, expr)),
            get_closure_vars($(Expr(:quote, expr)), @module_vars)
        ))
    )
end

generate_unique_name(name::Symbol, value::Integer = 0) = Symbol("$(name)ᐞ$(value)")

end # End GridTools module
