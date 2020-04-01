export energy, propose_candidate, sample

"""
    energy(configuration)

Compute the energy of a configuration in a SA computation.
"""
function energy end

"""
    propose_candidate(configuration)

Return tuple containing a neighboring configuration (a new candidate) and
the energy difference between the candidate and the original configuration.
The energy difference should be equivalent to `energy(candidate) -
energy(configuration)`
"""
function propose_candidate end
