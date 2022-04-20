#!/usr/bin/env julia1
push!(LOAD_PATH, pwd())
using Distributed, SpinDynamics, BenchmarkTools, UserInputs
using InitialCondition, EffectiveField, Energy, EffectiveSize, Dipolar

allParams = UserInputs.paramObjs

if userParams.cp.parallel == 1

  addprocs(userParams.cp.numCores)
  println("Initializing parallel cores. Number of processors = ", nprocs(),
        ", number of workers = ", nworkers())

    @everywhere module Test
        push!(LOAD_PATH, pwd())

        using SpinDynamics, Distributed
        include("Params.jl")

        function parallelize(allParams, rp)
          # All params is a list of all combinations of the requested inputs
          pmap(launchcomputation, allParams)
        end

    end

  using .Test
  Test.parallelize(userParams,rp)

else
    map(launchcomputation, allParams)
end
