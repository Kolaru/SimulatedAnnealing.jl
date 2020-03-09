using SimulatedAnnealing
using Random


struct Tour
    order::Vector{Int}
    D::Matrix{Float64}
end

function following(array, k)
    if k == length(array)
        return array[1]
    end

    return array[k + 1]
end

function preceding(array, k)
    if k == 1
        return array[end]
    end

    return array[k - 1]
end

function energy(tour::Tour)
    s = 0.0

    for (k, c1) in enumerate(tour.order)
        c2 = following(tour.order, k)

        s += tour.D[c1, c2]
    end

    return s
end

function propose_candidate(tour::Tour)
    n = length(tour.order)
    i1, i2 = rand(1:n, 2)

    order = copy(tour.order)
    
    if i1 > i2
        reverse!(order)
        i1, i2 = i2, i1
    end

    j1 = preceding(order, i1)
    j2 = following(order, i2)

    dE = tour.D[j1, i2] + tour.D[i1, j2] - tour.D[j1, i1] - tour.D[i2, j2]

    order[i1:i2] = reverse(order[i1:i2])

    return Tour(order, tour.D), dE
end

function sample(::Type{Tour}, N, D)
    n, n = size(D)
    return [Tour(randperm(n), D) for _ in 1:n]
end


cities = [
    [0.0 1.0 2.0 3.0 3.0 3.0 3.0 2.0 1.0 0.0 0.0 0.0] ;
    [0.0 0.0 0.0 0.0 1.0 2.0 3.0 3.0 3.0 3.0 2.0 1.0]
]

D = sqrt.((cities[1, :] .- cities[1, :]').^2 +
          (cities[2, :] .- cities[2, :]').^2)
