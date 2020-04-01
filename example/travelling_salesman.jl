using SimulatedAnnealing
using Random

"""
    Tour

Struct representing a tour visiting cities.

Fields
------
`order::Vector{Int}` order in which each city is visited.
`D::Matrix` Matrix of distances between each pair of cities. This is stored
to be able to pass it down to new candidates without the need of
recomputing it. 
"""
struct Tour
    order::Vector{Int}
    D::Matrix{Float64}
end

"""
    Tour(D::Matrix)

Create a tour with random ordering of the cities.
"""
Tour(D::Matrix) = Tour(randcycle(size(D, 1)), D)

"""
    energy(tour::Tour)

Compute the energy of a `Tour` by summing the distances of successive cities
in the tour.
"""
function SimulatedAnnealing.energy(tour::Tour)
    s = 0.0

    for (k, c1) in enumerate(tour.order)
        c2 = tour.order[mod1(k + 1, length(tour.order))]
        s += tour.D[c1, c2]
    end

    return s
end

"""
    propose_candidate(tour::Tour)

Propose a new candidate `Tour` by reversing the order of cities between two
random indices.
"""
function SimulatedAnnealing.propose_candidate(tour::Tour)
    n = length(tour.order)
    i1, i2 = rand(1:n, 2)  # Indices between which the order is reversed
    order = copy(tour.order)
    
    # Make sure that i1 < i2
    if i1 > i2
        i1, i2 = i2, i1
    end

    if i1 == 1 && i2 == 12
        return Tour(order, tour.D), 0.0
    end

    city_a1 = order[i1]
    city_a2 = order[i2]

    city_b1 = order[mod1(i1 - 1, n)]  # Index following i1
    city_b2 = order[mod1(i2 + 1, n)]  # Index preceding i2


    dE = (tour.D[city_a1, city_b2] + tour.D[city_a2, city_b1]
          - tour.D[city_a1, city_b1] - tour.D[city_a2, city_b2])

    order[i1:i2] = reverse(order[i1:i2])

    return Tour(order, tour.D), dE
end


cities = [
    [0.0 1.0 2.0 3.0 3.0 3.0 3.0 2.0 1.0 0.0 0.0 0.0] ;
    [0.0 0.0 0.0 0.0 1.0 2.0 3.0 3.0 3.0 3.0 2.0 1.0]
]

# Distance matrix, common for all Tour objects
D = sqrt.((cities[1, :] .- cities[1, :]').^2 +
          (cities[2, :] .- cities[2, :]').^2)
n = size(D, 1)

# Initial sampling of the configuration space
samples = [Tour(D) for _ in 1:1000]

# Run the algorithm with default parameters
best_tour, tour_length = simulated_annealing(samples, 10n*(n - 1))

println("Shortest tour length $tour_length, for the order:")
println("$(best_tour.order)")