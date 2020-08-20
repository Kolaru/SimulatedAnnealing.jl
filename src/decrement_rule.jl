"""
    ConstantDecrementRule

Fixed temperature decrement rule. The temperature is decrease by a constant
factor after each Markov chain.

Parameters
==========
factor: the constant factor multiplying the temperature.
"""
@with_kw struct ConstantDecrementRule{T}
    factor::T = 0.9
end

function (decrement_rule::ConstantDecrementRule)(state::AnnealingState)
    return decrement_rule.factor * state.temperature
end


"""
    AVLDecrementRule

Aarts and van Laarhoven temperature decrement rule.

Parameters
==========
distance_parameter

Reference
=========
Varanelli, eq. 2.42.
"""
@with_kw struct AVLDecrementRule{T}
    distance_parameter::T = 0.085
end

function (decrement_rule::AVLDecrementRule)(state::AnnealingState)
    T = state.temperature
    σ = std(state.energies)
    δ = decrement_rule.distance_parameter

    return T/(1 + T*log(1 + δ)/3σ)
end