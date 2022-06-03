#!/usr/bin/env julia1
push!(LOAD_PATH, pwd())
using Distributed, UserInputs, SpinDynamics

allParams = UserInputs.paramObjs

if allParams[1].cp.parallel == 1

  addprocs(allParams[1].cp.numCores)
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
