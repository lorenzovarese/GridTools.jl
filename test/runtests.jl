using .GridTools
using Test

@testset "GridTools tests" begin

    @testset "Dummy Testsuit" begin
        include("dummy_test.jl")
    end

end



# using .Example: greet, simple_add, type_multiply
# oder
# import Example: greet, simple_add, type_multiply  # obwohl das eigentlich ned funktioniere sött... ohni punkt isch nur für offizielli packages



# unterschied zwüsched import und using isch dass wenn im modul wo mit "using" gadded wird, d funktione exportiert worde sind, denn chasch sie normal ufrüefe mit greet().
# Mit import muesch sie immer mit Example.greet() au wenn sie im modul exportiert worde sind. Also
# using Test
# import Test
