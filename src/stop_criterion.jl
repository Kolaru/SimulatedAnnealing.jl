"""
    OVGCriterion

Otten-van Ginneken adaptive stop criterion.

Parameters
==========
threshold

Reference
=========
Varanelli, eq. 2.47.
"""
@with_kw struct OVGCriterion{T}
    energy_reference::T
    threshold::T = 0.001
end

function (criterion::OVGCriterion)(state::AnnealingState)
    σ = std(state.energies)
    dμ = criterion.energy_reference - mean(state.energies)

    σ == 0 && return true
    
    if dμ >= 0
        stop_measure = σ^2/(state.temperature*dμ)
    else
        stop_measure = Inf
    end

    return stop_measure < criterion.threshold
end


"""
    SSVCriterion

Sechen and Sangiovanni-Vincentelli stop criterion. Stops if the same energy
appears a given amount of time at the end of Markov chains.

Parameters
==========
maximum_repeat: maximum number of time the same BSF energy can be repeated.

Reference
=========
Varanelli, p. 26.
"""
@with_kw mutable struct SSVCriterion{T}
    maximum_repeat::Int = 3
    last_energy::T = Inf
    repeat_count::Int = 1
end

function (criterion::SSVCriterion)(state::AnnealingState)
    if state.bsf_energy == criterion.last_energy
        criterion.repeat_count += 1
    else
        criterion.last_energy = state.bsf_energy
        criterion.repeat_count = 1
    end
    return criterion.repeat_count == criterion.maximum_repeat
end