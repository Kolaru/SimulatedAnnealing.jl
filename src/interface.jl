export energy, propose_candidate, sample

"""
    energy(configuration)

Compute the energy of a configuration in a SA computation.
"""
function energy end

"""
    propose_candidate(configuration)

Return a candidate configuration based on the current configuration given.
Return both the candidate configuration and the energy difference between it
and the current configuration.
"""
function propose_candidate end
