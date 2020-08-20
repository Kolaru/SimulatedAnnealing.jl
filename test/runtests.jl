using SimulatedAnnealing
using Statistics
using Test

@testset "Travelling salesman" begin
    include("../example/travelling_salesman.jl")

    energies = energy.(samples)
    decrement_rules = Dict(
        "constant" => ConstantDecrementRule,
        "AVL" => AVLDecrementRule)
    stop_criteria = Dict(
        "OVG" => () -> OVGCriterion(energy_reference=mean(energies)),
        "SSV" => SSVCriterion)

    for (rulename, rule) in decrement_rules
        for (critname, criterion) in stop_criteria
            @testset "Decrement: $rulename, Stop: $critname" begin
                sa_prob = AnnealingOptimization(
                    energy,
                    propose_candidate,
                    criterion(),
                    rule(),
                    std(energies),
                    samples[1],
                    n*(n - 1))
                
                best_tour, tour_length = simulated_annealing(sa_prob)
                @test tour_length == 12
            end
        end
    end
end 