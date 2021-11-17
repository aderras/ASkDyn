#!/usr/bin/env julia1
push!(LOAD_PATH, pwd())
include("Params.jl")
using SpinDynamics, BenchmarkTools, UserInputs,InitialCondition, EffectiveField,
Energy, EffectiveSize, Dipolar, Chirality, LLequation, RungeKutta

testParams = Params.buildUserInputParam()

println("Build skyrmion")
s0 = buildinitial(testParams.ic, testParams.mp)

print("Compute chirality")
stmp = zeros(3)
@btime Chirality.computeGamma(stmp, s0,testParams.cp.rChirality)

print("Get material parameter fields")
@time [getfield(testParams.mp, x)
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

# Initialize arrays for storage
p, m, n = size(s0)
Heff = zeros(p,m,n)
SDotH = zeros(1,m,n)

# Initialize arrays for the RK steps
K1 = zeros(p,m,n)
K2 = zeros(p,m,n)
K3 = zeros(p,m,n)
K4 = zeros(p,m,n)
tmp = zeros(p,m,n);

println("\nCompute effective field")
print("Exchange")
@btime exchangefield!(Heff, s0, j, pbc)

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
@btime effectivefield!(Heff, s0, testParams)

Runge Kutta
effectivefield!(Heff, s0, testParams)
@btime rk4!(0.0, s0, RHS!, testParams, [Heff, SDotH],
        [K1, K2, K3, K4, tmp])
