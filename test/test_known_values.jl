using Test
using CSV

include("../src/diagrams.jl")

function on_known_values()
    for row in CSV.Rows("diagrams_known_values.csv";comment="#",delim=";",types=[String,Bool,Bool,Float64,String],ignoreemptylines=true)
        println("$(row.graph_path):")
        @test is_compact_respectively_finvol("../graphs/"*row.graph_path) == (row.compact,row.finvol)
    end
end


@testset "Diagrams" begin

    @testset "Known compactness/finite volume values" begin
        on_known_values() 
    end

    @testset "somethin else" begin
        @test true == true 
    end
end

