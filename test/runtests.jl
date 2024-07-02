# push!(LOAD_PATH, "/Users/jeffreyzweidler/Documents/CSCS/GridTools/src")dlopen(PyCall.libpython, RTLD_GLOBAL)

using GridTools
using Test
using OffsetArrays

#include("../src/gt2py/gt2py.jl")

# Setup --------------------------------------

@testset "GridTools tests" begin

    @testset "Embedded Testsuit" begin
        include("embedded_test.jl")
    end

    @testset "GT2Py Testsuit" begin
        include("gt2py_fo_parser.jl")
        include("gt2py_fo_exec.jl")
    end
end