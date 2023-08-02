using Statistics

# Setup ------------------------------------------------------------------------------------------------------------------------------
struct Cell_ <: Dimension end
struct K_ <: Dimension end
struct Edge_ <: Dimension end
struct E2C_ <: Dimension end

Cell = Cell_()
K = K_()
Edge = Edge_()
E2C = E2C_()

num_cells = 5
num_layers = 6
grid_shape = (num_cells, num_layers)
a_value = 2.0
b_value = 3.0

a = Field((Cell, K), fill(a_value, grid_shape))
b = Field((Cell, K), fill(b_value, grid_shape))
err = Field((Cell,), fill(b_value, 4))


cell_values = Field((Cell,), [1.0, 1.0, 2.0, 3.0, 5.0, 8.0])
edge_to_cell_table = [
    [1  0];
    [3  0];
    [3  0];
    [4  0];
    [5  0];
    [6  0];
    [1  6];
    [1  2];
    [2  3];
    [2  4];
    [4  5];
    [5  6]
]

E2C_offset_provider = Connectivity(edge_to_cell_table, (Cell,), (Edge, E2C), 2)
C2E_offset_provider = Connectivity(edge_to_cell_table, (Edge, E2C), (Cell,), 2)

# Tests ------------------------------------------------------------------------------------------------------------------------------

@testset "Testset arithmetic broadcast" begin
	
	addition = a .+ b
	@test a_value+b_value == mean(addition.data)
	@test 0 == std(addition.data)

    subtraction = a .- b
    @test a_value-b_value == mean(subtraction.data)
	@test 0 == std(subtraction.data)

    mult = a .* b
    @test a_value*b_value == mean(mult.data)
	@test 0 == std(mult.data)

    div = a ./ b
    @test a_value/b_value == mean(div.data)
	@test 0 == std(div.data)

	@test isa(addition, Field)
    @test isa(subtraction, Field)
    @test isa(mult, Field)
    @test isa(div, Field)
	@test_throws AssertionError (a .+ err)
end

@testset "Testset field_call" begin
	
end

@testset "Testset conn_call" begin
	
end


@testset "Testset neighbor_sum" begin
	
	edge_values = neighbor_sum(cell_values(E2C_offset_provider()), axis=E2C)

	@test edge_values.data == [1.,  2.,  2.,  3.,  5.,  8.,  9.,  2.,  3.,  4.,  8., 13.]

	@test isa(edge_values, Field)
end