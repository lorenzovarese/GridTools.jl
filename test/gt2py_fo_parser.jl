using Test
using GridTools


include("mesh_definitions.jl")

struct TestFailedException <: Exception
    message::String
end

function to_py(field_operator::GridTools.FieldOp)
    try
        GridTools.py_field_operator(field_operator)
        return true
    catch e
        throw(TestFailedException("The following test: $(field_operator.name) encountered the following error: $e"))
    end
end

function test_arithmetic_operators()
    expr = @field_operator function test_b_addition(f::Field{Tuple{Cell_}, Int64}, g::Field{Tuple{Cell_}, Int64})
              return f .+ g
            end
    @test to_py(expr)
    
    expr = @field_operator function test_bit_and(f::Field{Tuple{Cell_}, Bool}, g::Field{Tuple{Cell_}, Bool})
                return f .& g
            end
    @test to_py(expr)

    expr = @field_operator function test_bit_xor(f::Field{Tuple{Cell_}, Bool}, g::Field{Tuple{Cell_}, Bool})
                return f .⊻ g
            end
    @test to_py(expr)

    expr = @field_operator function test_not(f::Field{Tuple{Cell_}, Bool}, g::Field{Tuple{Cell_}, Bool})
                return .~f
            end
    @test to_py(expr)
    
    expr = @field_operator function test_bool_and(f::Field{Tuple{Cell_}, Bool}, g::Field{Tuple{Cell_}, Bool})
                return f && g
            end
    @test to_py(expr)
    
    expr = @field_operator function test_bool_or(f::Field{Tuple{Cell_}, Bool}, g::Field{Tuple{Cell_}, Bool})
                return f || g
            end
    @test to_py(expr)
end

@testset "Testset Arithmetic Operators" test_arithmetic_operators()

@testset "Testset function header annotation" begin
    @test_throws AssertionError begin
        expr = @field_operator function test_no_annotation(inp)
            return inp
        end
        to_py(expr)
    end

    # Type annotation error
    @test_throws AssertionError begin
        expr = @field_operator function test_no_annotation_kwargs(f::Field{Tuple{Cell_}, Float64}; g)
                    return g
                end
        to_py(expr)
    end

    expr = @field_operator function test_simple_annotation(inp::Integer)
                return inp
            end
    @test to_py(expr)

    expr = @field_operator function test_kwargs(inp::Integer; inp2::AbstractFloat, inp3::AbstractFloat)
                return inp
            end
    @test to_py(expr)

    expr = @field_operator function test_return_annotation(f::Field{Tuple{Cell_}, Float64}, g::Field{Tuple{Cell_}, Float64})::Field{Tuple{Cell_}, Float64}
                return f
            end
    @test to_py(expr)

    # Type deduction error
    expr = @field_operator function test_incompatible_types(f::Field{Tuple{Cell_}, Float64}, g::Field{Tuple{Cell_}, Int64})::Field{Tuple{Cell_}, Float64}
                return f .+ g
            end
    @test_throws TestFailedException to_py(expr)
    
    # Return statement error
    expr = @field_operator function test_no_return_stmt(f::Field{Tuple{Cell_}, Float64}; g::Field{Tuple{Cell_}, Int64})
                g .+ 1
            end
    @test_throws TestFailedException to_py(expr)
end


@testset "Testset Preprocessing" begin
    expr = @field_operator function test_multi_assign(f::Field{Tuple{Cell_}, Float64})
                a, b = 5., 6.
                return f .+ a .- b
            end
    @test to_py(expr)

    expr = @field_operator function test_renaming_variables(f::Field{Tuple{Cell_}, Float64})
                tmp = 1.
                tmp = 2. + tmp
                tmp = 3. + tmp
                return tmp
            end
    @test to_py(expr)

    expr = @field_operator function test_ternary_expr(f::Field{Tuple{Cell_}, Float64}, g::Field{Tuple{Cell_}, Float64})::Field{Tuple{Cell_}, Float64}
                return 1 < 2 ? f :  g
            end
    @test to_py(expr)

    expr = @field_operator function test_renaming_conditional(f::Field{Tuple{Cell_}, Int32})
                if 1. .< 10.
                    tmp = f .+ 1
                elseif 1. < 10.
                    tmp = f .- 1
                else
                    tmp = f .+ 1
                end
                return tmp
            end
    @test to_py(expr)

    expr = @field_operator function test_renaming_nested_conditional(f::Field{Tuple{Cell_}, Int32}, g::Field{Tuple{Cell_}, Int32})::Field{Tuple{Cell_}, Int32}
                tmp = f
                if 1. .< 10.
                    tmp = f .+ 1
                    if 1. .< 10.
                        tmp = tmp .- 1
                        tmp = tmp .* 1
                    elseif 1. .< 10.
                        tmp = 1. .< 10. ? tmp : tmp ./ 1
                    else 
                        tmp = tmp .* 1
                    end
                    tmp = tmp .+ 1
                elseif 1. .< 10.
                    tmp = f .- 1
                else
                    tmp = tmp .* 1
                    tmp = tmp .+ 1
                    tmp = tmp ./ 1
                end
                return tmp
            end
    @test to_py(expr)

    expr = @field_operator function test_compair_chain(f::Field{Tuple{Cell_}, Int32}, g::Field{Tuple{Cell_}, Int32})::Field{Tuple{Cell_}, Bool}
                return 1 < 10 .< f .< g .< 50 > -1
            end
    @test to_py(expr)

    # Only constants in condition, Should this actually fail?  TODO
    expr = @field_operator function test_if_with_field(f::Field{Tuple{Cell_}, Float64}, g::Field{Tuple{Cell_}, Float64})::Field{Tuple{Cell_}, Float64}
                return f .< 2 ? f :  g
            end
    @test_throws TestFailedException to_py(expr)
end


@testset "Test built-ins and closure variables" begin

    expr = @field_operator function test_where(mask::Field{Tuple{Cell_}, Bool}, f::Field{Tuple{Cell_}, Float64}, g::Field{Tuple{Cell_}, Float64})::Field{Tuple{Cell_}, Float64}
                return where(mask, f, g)
            end
    @test to_py(expr)

    # Error due to non boolean mask argument
    expr = @field_operator function test_where(mask::Field{Tuple{Cell_}, Bool}, f::Field{Tuple{Cell_}, Float64}, g::Field{Tuple{Cell_}, Float64})::Field{Tuple{Cell_}, Float64}
                return where(f, f, g)
            end
    @test_throws TestFailedException to_py(expr)


    expr = @field_operator function test_maxover(f::Field{Tuple{Cell_, K_}, Float64})::Field{Tuple{Edge_}, Float64}
                return max_over(f(E2C[1]), axis=K)
            end
    @test to_py(expr)

    expr = @field_operator function test_neighborsum(f::Field{Tuple{Cell_, K_}, Float64})::Field{Tuple{Edge_, K_}, Float64}
                return neighbor_sum(f(E2C), axis=E2CDim)
            end
    @test to_py(expr)

    expr = @field_operator function test_broadcast(f::Field{ Tuple{Cell_}, Float64})
                return broadcast(f, (Cell, K))
            end
    @test to_py(expr)

    expr = @field_operator function test_julia_builtin(f::Field{Tuple{Cell_, K_}, Float64})::Field{Tuple{Cell_, K_}, Float64}
                return sin.(f)
            end
    @test to_py(expr)

end

