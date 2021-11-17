#!/usr/bin/env julia1
push!(LOAD_PATH, pwd())
include("Parameters.jl")
using Distributed, SpinDynamics, BenchmarkTools, UserInputs
using InitialCondition, EffectiveField, Energy, EffectiveSize, Dipolar

userParams = Parameters.buildUserInputParam()
rp = Parameters.paramRanges()#UserInputs.Ranges.r)#, UserInputs.Ranges.pbc)

userParams.ic.type = "prec"

if userParams.cp.parallel == 1

  addprocs(userParams.cp.numCores)
  println("Initializing parallel cores. Number of processors = ", nprocs(),
        ", number of workers = ", nworkers())

    @everywhere module Test
        push!(LOAD_PATH, pwd())

        using SpinDynamics, Distributed
        include("Parameters.jl")

        function parallelize(userParams, rp)
          # All params is a list of all combinations of the requested inputs
          allParams = []
          Parameters.getparamlist!(userParams, allParams, rp)
          pmap(launchcomputation, allParams)
        end

    end

  using .Test
  Test.parallelize(userParams,rp)

else

    allParams = []
    Parameters.getparamlist!(userParams, allParams, rp)

    map(launchcomputation, allParams)

end
