



# Trying to do shit for broadcast #####################################################################


function Base.size(F::Field)::Tuple
    size_tup = [typemax(Int) for _ in 1:length(F.broadcast_dims)]
    inds_real = findall(in.(F.broadcast_dims, Ref(F.dims)))
    size_tup[inds_real] = [size(F.data)...]
    return Tuple(size_tup)
end


@propagate_inbounds function Base.getindex(F::Field{T,N}, inds::Vararg{Int,N}) where {T,N}
    inds_real = findall(in.(F.broadcast_dims, Ref(F.dims)))
    F.data[inds[inds_real]...]
    # F.data[inds...]
end
@propagate_inbounds function Base.setindex!(F::Field{T,N}, val, inds::Vararg{Int,N}) where {T,N}
    inds_real = findall(in.(F.broadcast_dims, Ref(F.dims)))
    F.data[inds[inds_real]...] = val
end


# where implementation with @generated #####################################################################

function dfs_index(exp::Expr)
    if typeof(exp.args[1]) == Symbol
        return [exp.args[2]]
    else
        return [dfs_index(exp.args[1])..., exp.args[2]]
    end 
end

## Recursive Version
# function fold_tuple(res::Expr, t::Expr, p::Expr, d::Vector, mask::Type{<:Field})::Expr
#     ind = d[1]
#     if length(d) > 1
#         @assert length(res.args) in (ind-1, ind)
#         if length(res.args) == ind-1
#             push!(res.args, Expr(:tuple))
#         end
#         res.args[ind] = fold_tuple(res.args[ind], t, p, d[2:length(d)], mask)
#     else
#         temp = :(ifelse.(mask, $t, $p))
#         push!(res.args, temp)
#     end
#     return res
# end

## Iterative version
function fold_tuple(t::Expr, p::Expr, mask::Type{<:Field})::Expr
   res = Expr(:tuple)
   for tuple_id in 1:length(t.args)             # for all tuples
        path = dfs_index(t.args[tuple_id])       # get the path of a tuple
        current_expr = res
        for path_step in 1:length(path)-1
            @assert length(current_expr.args) in (path[path_step]-1, path[path_step])
            if length(current_expr.args) == path[path_step]-1
                push!(current_expr.args, Expr(:tuple))
            end
            current_expr = current_expr.args[path[path_step]]
        end
        push!(current_expr.args, :(ifelse.(mask, $(t.args[tuple_id]), $(p.args[tuple_id]))))
   end
   return res
end

function unfold_tuple(tt::Type{<:Tuple}, tuple_name::Symbol, path::Vector = [])::Expr
    res = Expr(:tuple)
    for i in 1:length(tt.parameters)
        if tt.parameters[i] <: Tuple
            push!(res.args, unfold_tuple(tt.parameters[i], tuple_name, [path..., i]).args...)
        else
            current_expr = tuple_name
            for current_index in [path..., i]
                current_expr = Expr(:ref, current_expr, current_index)
            end
            push!(res.args, current_expr)
        end
    end
    return res
end

@inbounds @generated function where(mask::Field, tt::Tuple, pp::Tuple)   # TODO: Make @generated. Problem with naming String, vectors used like d::Vector and so on.
    @assert tt == pp
    t = unfold_tuple(tt, :tt)
    p = unfold_tuple(pp, :pp)

    ## Recursive Version
    # res = Expr(:tuple)

    # for i in 1:length(t.args)
    #     d = dfs_index(t.args[i])
    #     res = fold_tuple(res, t.args[i], p.args[i], d, mask)
    # end

    ## Iterative version
    res = fold_tuple(t, p, mask)

    return res
end


# Snippets for custom broadcast ###############################################################################

@inline function combine_axes(A::Tuple, B::Tuple)
    length(A[1]) > length(B[1]) ? (A,B) = (B,A) : nothing  # swap A and B       

    if issubset(A[1], B[1])
        # A,B = (dims, axes(data), broadcast_dims)
        @assert issubset(B[1], A[3]) "Dimension Mismatch between the broadcasted Dimensions of the two Fields"
        matching_dims = get_dim_ind(A[1], B[1])
        @assert A[2] == B[2][matching_dims] "Dimension Mismatch between the data Dimensions of the two Fields"
        return B
    elseif isempty(intersect(A[1],B[1]))
        matching_dimsA = get_dim_ind(intersect(A[1],B[1]),A[1])
        matching_dimsB = get_dim_ind(intersect(A[1],B[1]),B[1])
        @assert A[2][matching_dimsA] == B[2][matching_dimsB] "Dimension Mismatch between the data Dimensions of the two Fields"

        
        return (order(union(A[1],B[1])), )
    else
        error("Dimension Mismatch between the Dimensions of the two Fields")
    end
    
end