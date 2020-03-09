module SimulatedAnnealing

using Statistics

export SimulatedAnnealing, AnnealingState
export SA_optimize


include("interface.jl")

"""
    SimulatedAnnealing{C}

Object representing an optimization done using a Simulated Annealing
algorithm.

Implement single stage homogeneous Simulated Annealing as (very well) described
in
Varanelli, 1996, *On the acceleration of simulated annealing*, Chapter 2.

Notation and naming convention is consistent with that work.

A configuration of the system being optimized is represented by the parametric
type `C`. The configurations to optimize can be of any type, but the functions
`energy` and `propose_candidate` must be defined for it (see their
documentation for details).

The object is iterable to allow access to the internal state during the
optimization process. The iteration stops when the stop measure given by
Otten-van Ginneken criterion (eq. 2.47 in Varanelli) falls below the
`stop_parameter`.


Fields
------
- `initial_temperature`
- `neighborhood_size::Int`: Number of states that can be accessed from a given one.
- `distance_parameter`: Parameter controlling the speed of the temperature decrement.
- `stop_parameter`: Parameter determining the stop criterion.
"""
struct SimulatedAnnealing{C}
    initial_temperature::Float64
    initial_configuration::C
    energy_reference::Float64
    neighborhood_size::Int
    distance_parameter::Float64
    stop_parameter::Float64
end


"""
    SimulatedAnnealing(samples, neighborhood_size, distance_parameter, stop_parameter)

Given samples of the configuration space, determine the initial temperature
as well as the energy reference automatically. The first sample is used as
initial configuration.

Initial temperature is the standard deviation of the energy in the sample,
according to White criterion (see eq. 2.36 in Varanelli).
"""
function SimulatedAnnealing(samples::Vector{C}, ns, dp, sp) where C
    energies = energy.(samples)
    return SimulatedAnnealing{C}(
        std(energies),
        samples[1],
        mean(energies),
        ns, dp, sp)
end

"""
    AnnealingState{C}

Struct representing the sate of the Simulated Annealing algorithm at the end
of a Markov chain step.

Fields
------
- `temperature`
- `energy_reference`: Energy reference (μ0) used in Otten-van Ginneken stop criterion.
- `stop_measure`: Value of the Otten-van Ginneken stop criterion.
- `current_configuration`
- `current_energy`
- `bsf_configuration`: Best configuration encountered So Far
- `bsf_energy`: Energy of the Best So Far configuration.
- `energies`: All energies encountered during the previous Markov chain.
"""
struct AnnealingState{C}
    temperature::Float64
    current_configuration::C
    current_energy::Float64
    bsf_configuration::C
    bsf_energy::Float64
    energies::Vector{Float64}
    stop_measure::Float64
end

function AnnealingState(temperature, configuration)
    E = energy(configuration)
    return AnnealingState(temperature,
                          configuration, E,
                          configuration, E,
                          zeros(0), Inf)
end


Base.eltype(::Type{SA}) where {C, SA <: SimulatedAnnealing{C}} = AnnealingState{C}
Base.IteratorSize(::Type{SA}) where {SA <: SimulatedAnnealing} = Base.SizeUnknown()

function Base.iterate(
    search::SimulatedAnnealing{C},
    state=AnnealingState(C, search.initial_temperature, search.initial_configuration)) where C

    state.stop_measure < search.stop_parameter && return nothing

    # Initialize all internal loop variables
    configuration = state.current_configuration
    E = state.current_energy

    bsf = state.bsf_configuration
    bsf_energy = state.bsf_energy
    
    energies = zeros(search.neighborhood_size)

    # Markov chain at constant temperature
    for k in 1:search.neighborhood_size
        candidate, dE = propose_candidate(configuration)

        if dE < 0 || rand() < exp(dE/state.temperature)
            configuration = candidate
            E += dE

            if E < state.bsf_energy
                bsf = configuration
                bsf_energy = E
            end
        end

        energies[k] = E
    end
    
    μ = mean(energies)
    σ = std(energies)

    μ0 = state.energy_reference
    dμ = μ0 - μ

    # Otten-van Ginneken stop criterion (eq. 2.47 in Varanelli)
    if dμ >= 0
        stop_measure = σ^2/(state.temperature*dμ)
    else
        stop_measure = Inf
    end
    
    T = decrement_rule(state.temperature, σ, search.distance_parameter)

    new_state = AnnealingState(T, μ0, stop_measure, configuration, E,
                               bsf, bsf_energy, energies)

    return new_state, new_state
end


"Aarts and van Laarhoven temperature decrement rule (eq. 2.41 in Varanelli)."
decrement_rule(T, σ, δ) = T/(1 + T*log(1 + δ)/3σ)


function SA_optimize(samples::Vector{C},
                     neighborhood_size ;
                     n_thermal_sample=1000,
                     distance_parameter=0.085,
                     stop_parameter=0.0001) where C

    search = SimulatedAnnealing(samples, neighborhood_size,
                                distance_parameter, stop_parameter)
    
    local state = nothing  # Avoid being shadowed inside the loop

    for s in search  # Go to the end of the search
        state = s
    end

    return state.bsf_configuration, state.bsf_energy
end

end