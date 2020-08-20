module SimulatedAnnealing

using Parameters
using Statistics

export 
    AnnealingOptimization, AnnealingState,
    OVGCriterion, SSVCriterion,
    ConstantDecrementRule, AVLDecrementRule,
    simulated_annealing


"""
    AnnealingOptimization

Object representing an optimization done using a Simulated Annealing
algorithm.

Implement single stage homogeneous Simulated Annealing as (very well) described
in
Varanelli, 1996, *On the acceleration of simulated annealing*, Chapter 2.

Notation and naming convention is consistent with that work.

A configuration of the system being optimized is represented by the parametric
type `C`. The configurations to optimize can be of any type.

The object is iterable to allow access to the internal state during the
optimization process.

Parameters
==========
energy: callable(::C)::T
    Callable returning the energy of a configuration.
propose_candiate: callable(::C)::Tuple{::C, ::T}
    Callable returning a new configuration together based on an existing one
    togetehr with the energy difference between the two.
stop_criterion: callable(::AnnealingState)::Bool
    Callable returning whether the annealing must stop.
    Default is Otten-van Ginneken adaptive stop criterion.
decrement_rule: callable(::AneealingState)::T
    Callable returning the new temperature to use at the end of a Markov step.
    Default is Aarts and van Laarhoven decrement rule.
initial_temperature: Real
initial_configuration: C
neighborhood_size: integer
    Number of states that can be accessed from a given one.
"""
struct AnnealingOptimization{T, C, EFUNC, PFUNC, SFUNC, TFUNC}
    energy::EFUNC
    propose_candidate::PFUNC
    stop_criterion::SFUNC
    decrement_rule::TFUNC
    initial_temperature::T
    initial_configuration::C
    neighborhood_size::Int
end

function AnnealingOptimization(
        energy,
        propose_candidate,
        initial_configuration,
        initial_temperature::Real,
        neighborhood_size::Integer)

    return AnnealingOptimization(
        energy,
        propose_candidate,
        SSVCriterion(),
        AVLDecrementRule(),
        initial_temperature,
        initial_configuration,
        neighborhood_size)
end


"""
    AnnealingOptimization(energy, propose_candidate, samples, neighborhood_size)

Given samples of the configuration space, determine the initial temperature.
The first sample is used as initial configuration.

Initial temperature is the standard deviation of the energy in the sample,
according to White criterion (see eq. 2.36 in Varanelli).
"""
function AnnealingOptimization(energy, propose_candidate, samples::Vector,
                               neighborhood_size::Integer)

    energies = energy.(samples)
    return AnnealingOptimization(
        energy,
        propose_candidate,
        samples[1],
        std(energies),
        neighborhood_size)
end

"""
    AnnealingState

Struct representing the sate of the Simulated Annealing algorithm at the end
of a Markov chain step.

Fields
======
temperature
current_configuration
current_energy
bsf_configuration
    Best configuration found so far.
bsf_energy
    Energy of the best configuratin found so far.
energies
    Array of all energies encountered during the last Markov chain step.
"""
struct AnnealingState{T, C}
    temperature::T
    current_configuration::C
    current_energy::T
    bsf_configuration::C
    bsf_energy::T
    energies::Vector{T}
end

function initial_state(search::AnnealingOptimization{T}) where T
    config = search.initial_configuration
    E = search.energy(config)
    return AnnealingState(search.initial_temperature,
                          config, E,
                          config, E,
                          zeros(T, 0))
end

include("decrement_rule.jl")
include("stop_criterion.jl")

Base.eltype(::Type{A}) where {T, C, A <: AnnealingOptimization{T, C}} = AnnealingState{T, C}
Base.IteratorSize(::Type{A}) where {A <: AnnealingOptimization} = Base.SizeUnknown()

function Base.iterate(
        search::AnnealingOptimization{T},
        state=initial_state(search)) where T
    length(state.energies) > 0 && search.stop_criterion(state) && return nothing

    # Initialize all internal loop variables
    configuration = state.current_configuration
    E = state.current_energy
    bsf = state.bsf_configuration
    bsf_energy = state.bsf_energy
    energies = zeros(T, search.neighborhood_size)

    # Markov chain at constant temperature
    for k in 1:search.neighborhood_size
        candidate, dE = search.propose_candidate(configuration)

        if dE < 0 || rand() < exp(-dE/state.temperature)
            configuration = candidate
            E += dE

            if E < state.bsf_energy
                bsf = configuration
                bsf_energy = E
            end
        end

        energies[k] = E
    end

    new_state = AnnealingState(search.decrement_rule(state),
                               configuration, E,
                               bsf, bsf_energy,
                               energies)

    return new_state, new_state
end

function simulated_annealing(search::AnnealingOptimization)
    local state = nothing  # Avoid being shadowed inside the loop

    for s in search  # Go to the end of the search
        state = s
    end

    return state.bsf_configuration, state.bsf_energy
end

end