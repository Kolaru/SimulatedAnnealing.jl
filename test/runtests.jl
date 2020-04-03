using Test

@testset "Travelling salesman" begin
    include("../example/travelling_salesman.jl")

    @test tour_length == 12.0
end 