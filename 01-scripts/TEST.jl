#!/usr/bin/env julia1
push!(LOAD_PATH, pwd())
include("Parameters.jl")
using Distributed, SpinDynamics, BenchmarkTools, UserInputs
using Energy, TopologicalCharge, LLequation, Magnetization, RungeKutta,
Normalize, NoiseRotation, EffectiveSize, InitialCondition, UserInputs,
EffectiveField, Dipolar, BenchmarkTools, HDF5

function h5overwrite(name, array)
    fid = h5open(name, "w")
    write(fid, "Dataset1", array)
    close(fid)
end

j = 1.0
h = parse(Float64,ARGS[1])
a = parse(Float64,ARGS[2])
dz = parse(Float64,ARGS[3])
ed = parse(Float64,ARGS[4])
nx = ny = parse(Int,ARGS[8])
nz = 10
pbc = parse(Float64,ARGS[7])

tMax = 1
dt = 0.1
nSteps = 10
tol = 10^-6

type = "skyrmion"
r = 10
chirality = pi/2
icx = parse(Int,ARGS[8])/2
icy = parse(Int,ARGS[8])/2

# Only build the following matrices if the dipolar interaction is nonzero
if parse(Float64,ARGS[4])!=0.0
    v = vdmatrices(nx,ny,nz,1.0==pbc)
else
    v = [[0.0 0.0; 0.0 0.0],[0.0 0.0; 0.0 0.0],[0.0 0.0; 0.0 0.0],
        [0.0 0.0; 0.0 0.0]]
end

# Create all the structs to be used by the rest of the program.
mp = Parameters.materialParams(j,h,a,dz,ed,nx,ny,nz,pbc,v)

cp = Parameters.compParams(tMax, dt, nSteps, tol, parse(Float64,ARGS[6]),
    parse(Float64,ARGS[5]), 0, 1, 0, 1, 0.3, 10, tol)

ic = Parameters.icParams(type, r, chirality, icx, icy)

pp = Parameters.pinningParams(0.0)
dp = Parameters.defectParams(0.0,0.0,0.0,0,0,0)
cc = Parameters.currentParams(0.0, 0.0, 0)
sp = Parameters.saveChoices(1,0,0,0,0,0,0,1,1,1,0)
rp = Parameters.paramRanges([])

p = Parameters.params(mp, cp, ic, pp, dp, cc, sp, rp)

dataDir = string(dirname(pwd()),"/02-data/test/")
println("Saving results in ", dataDir)

# Build the initial spin field and save.
s0 = buildinitial(p.ic, p.mp)
h5overwrite(string(dataDir,"/initial-spin-field.h5"),s0)

if ed != 0.0
    println("\nSaving VDD matrices")
    h5overwrite(string(dataDir,"/vddxx.h5"), p.mp.v[1])
    h5overwrite(string(dataDir,"/vddyy.h5"), p.mp.v[2])
    h5overwrite(string(dataDir,"/vddzz.h5"), p.mp.v[3])
    h5overwrite(string(dataDir,"/vddxy.h5"), p.mp.v[4])

    println("\nSaving FHD result")
    h5overwrite(string(dataDir,"/fhdResult.h5"), fhd(s0, p.mp.v, pbc))

    println("\nTiming calcDdiField")
    @btime ddifield(s0, p.mp.ed, pbc, p.mp.v)

    println("\nSaving convolution")
    convAns = convfft(s0[1,:,:], p.mp.v[2], pbc)
    h5overwrite(string(dataDir,"/convolution-result.h5"), convAns)
end

println("\nSaving energy")
energy(s0, p)
h5overwrite(string(dataDir, "/energy.h5"), energy(s0, p))

println("\nSaving effective field")
effField = effectivefield(s0, p)
@btime effectivefield(s0, p)
h5overwrite(string(dataDir, "/effectiveField.h5"), effField)

println("\nSaving spin field after 5 steps")
rk4!(s0, RHS!, p)
# rk4!(s0, RHS!, p)
h5overwrite(string(dataDir, "/afterRk.h5"), s0)
