# SimulatedAnnealing.jl

This package aims is to provide a simple and flexible implementation of the simulated annealing algorithm in Julia.

The implementation is fully based on the following PhD thesis:
> Varanelli, 1996, *On the acceleration of simulated annealing*

Currently, nothing fancy is done and only the basic homogenous algorithm as describe in Chapter 2 of this work is implemented (in particular the whole part about accelerating simulated annealing is ignored).

## Basics

Simulated annealing try to minimize the energy of a configuration by successively probing neighboring configurations. In this package, no restriction is imposed on the type of configuration being optimize.

Before starting the algorithm we need several things:

- An `energy` function, that takes a configuration as an argument and return its energy, the quantity that the algorithm minimizes.
- A `propose_candidate` function that takes a configuration and return a candidate configuration and the energy difference between the two. This is the main function used in the algorithm and the performance mostly depends on it. Returning the energy difference is required as it is often possible to make this computation very efficient compared to explicitly computing the energy of the two configurations and taking their difference. In fact if there is no efficient way of computing this difference, it may indicate simulated annealing is not the correct algorithm to use for the problem.
- An unifrom sampling of the space of configuration. This is used to determine the initial temperature as well as an energy reference for the stop criterion. Around 1000 samples should be enough.
- The size of the neighborhood of a configuration. This is needed to have sufficient sampling and thus ensure (in probability) convergence.

When all that is set, an optimization problem can be created and the algorithm can be run on it, yielding the best configuration and its energy.

```julia
sa_problem = AnnealingOptimization(energy, propose_candidate, samples, neighborhood_size)
best_configuration, best_energy = simulated_annealing(sa_problem)
```

Please refer to the example folder for a complete example.

## Customization

Various temperature decrement rules and stop criterion can be chosen. New ones can also be created. Please refer to the docstrings for more information, or open an issue if something is not clear (I am too lazy to write a full documentation before seeing if there is interest for the package ^^).