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
    UserInputs.Computation.faTol)

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
    UserInputs.SaveChoices.loc, UserInputs.SaveChoices.fieldDuring)

rp = Parameters.paramRanges(UserInputs.Ranges.rRange)

testParams = Parameters.params(mp, cp, ic, pp, dp, cc, sp, rp)#, hp)

println("Build skyrmion")
s0 = buildinitial(testParams.ic, testParams.mp)

print("Get material parameter fields")
@btime [getfield(testParams.mp, x)
    for x in fieldnames(typeof(testParams.mp))]
j,h,a,dz,ed,nx,ny,nz,pbc,vd = [getfield(testParams.mp, x)
    for x in fieldnames(typeof(testParams.mp))]

println("\nNx=Ny=", nx, " lattice")

println("\nCompute energy")
print("Exchange")
@btime exchange_energy(s0, j, pbc)

print("Zeeman")
@btime zeeman_energy(s0, h)

print("DMI")
@btime dmi_energy(s0, a, pbc)

print("PMA")
@btime pma_energy(s0, dz)

if ed!=0.0
    print("DDI")
    @btime ddi_energy(s0, ed, pbc, vd)
end

Heff = zeros(3,nx,ny)
println("\nCompute effective field")
print("Exchange")
@btime exchangefield!(Heff, s0, j, pbc, testParams.defect)

print("Zeeman")
@btime zeemanfield!(Heff, s0, h)

print("DMI")
@btime dmifield!(Heff, s0, a, pbc)

print("PMA")
@btime pmafield!(Heff, s0, dz)

if ed!=0.0
    print("DDI")
    @btime ddifield(s0, ed, pbc, vd)
end

print("Full effective field")
@btime effectivefield(s0, testParams)
