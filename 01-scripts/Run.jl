#!/usr/bin/env julia1
push!(LOAD_PATH, pwd())
include("Parameters.jl")
using Distributed, SpinDynamics, BenchmarkTools, UserInputs
using InitialCondition, EffectiveField, Energy, EffectiveSize, Dipolar

# Only build the following matrices if the dipolar interaction is nonzero
if UserInputs.Material.ed!=0.0
    v = vdmatrices(UserInputs.Material.nx, UserInputs.Material.ny,
        UserInputs.Material.nz, 1.0==UserInputs.Material.pbc)
else
    v = [[0.0 0.0; 0.0 0.0],[0.0 0.0; 0.0 0.0],[0.0 0.0; 0.0 0.0],
        [0.0 0.0; 0.0 0.0]]
end

# Create all the structs to be used by the rest of the program.
mp = Parameters.materialParams(UserInputs.Material.j,
    UserInputs.Material.h, UserInputs.Material.a, UserInputs.Material.dz,
    UserInputs.Material.ed, UserInputs.Material.nx, UserInputs.Material.ny,
    UserInputs.Material.nz, UserInputs.Material.pbc, v)

cp = Parameters.compParams(UserInputs.Computation.tMax, UserInputs.Computation.dt,
    UserInputs.Computation.nSteps, UserInputs.Computation.tol,
    UserInputs.Computation.damping, UserInputs.Computation.T,
    UserInputs.Computation.parallel, UserInputs.Computation.numCores,
    UserInputs.Computation.runRelaxation, UserInputs.Computation.sMax,
    UserInputs.Computation.faConst,UserInputs.Computation.nRot,
    UserInputs.Computation.faTol, UserInputs.Computation.rChirality)

ic = Parameters.icParams(UserInputs.InitialCondition.type,
    UserInputs.InitialCondition.r, UserInputs.InitialCondition.chirality,
    UserInputs.InitialCondition.icx, UserInputs.InitialCondition.icy)

pp = Parameters.pinningParams(UserInputs.Pinning.pinField)

dp = Parameters.defectParams(UserInputs.Defect.defType,
    UserInputs.Defect.defStrength, UserInputs.Defect.defWidth,
    UserInputs.Defect.defX, UserInputs.Defect.defY, [])

cc = Parameters.currentParams(UserInputs.Current.xCurr,
    UserInputs.Current.yCurr, UserInputs.Current.tOff)

sp = Parameters.saveChoices(UserInputs.SaveChoices.totalE,
    UserInputs.SaveChoices.exchE, UserInputs.SaveChoices.zeemE,
    UserInputs.SaveChoices.dmiE, UserInputs.SaveChoices.pmaE,
    UserInputs.SaveChoices.ddiE, UserInputs.SaveChoices.magn,
    UserInputs.SaveChoices.size, UserInputs.SaveChoices.charge,
    UserInputs.SaveChoices.loc, UserInputs.SaveChoices.fieldDuring,
    UserInputs.SaveChoices.chirality)

rp = Parameters.paramRanges(UserInputs.Ranges.r,UserInputs.Ranges.pbc)

testParams = Parameters.params(mp, cp, ic, pp, dp, cc, sp)

if testParams.cp.parallel == 1

  addprocs(testParams.cp.numCores)
  println("Initializing parallel cores. Number of processors = ", nprocs(),
        ", number of workers = ", nworkers())

    @everywhere module Test
        push!(LOAD_PATH, pwd())

        using Parameters, SpinDynamics, Distributed

        function parallelize(testParams, rp)
          # All params is a list of all combinations of the requested inputs
          allParams = []
          Parameters.getparamlist!(testParams, allParams, rp)
          pmap(launchcomputation, allParams)
        end

    end

  using .Test
  Test.parallelize(testParams,rp)

else

allParams = []
Parameters.getparamlist!(testParams, allParams, rp)

for i in 1:length(allParams)
    println()
    println(allParams[i])
end

# map(launchcomputation, allParams)

end
