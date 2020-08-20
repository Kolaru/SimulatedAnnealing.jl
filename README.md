# SimulatedAnnealing.jl

This package aims is to provide a simple and flexible implementation of the simulated annealing algorithm in Julia.

The implementation is fully based on the following PhD thesis:
> Varanelli, 1996, *On the acceleration of simulated annealing*

Currently, nothing fancy is done and only the basic homogenous algorithm as describe in Chapter 2 of this work is implemented (in particular the whole part about accelerating simulated annealing is ignored).

## Basics

Simulated annealing try to minimize the energy of a configuration by successively probing neighboring configurations. In this package, no restriction is imposed on the type of configuration being optimize.

TODO: More explanation and maybe even an example

Before starting the algorithm we need two more things:

- A unifrom sampling of the space of configuration. This is used to determine the initial temperature as well as an energy reference for the stop criterion. Around 1000 samples should be enough.
- The size of the neighborhood of a configuration. This is needed to have sufficient sampling and thus ensure (in probability) convergence.


## `SimulatedAnnealing` iterator

TODO: Nice exemple of using ProgressMeter to show progress