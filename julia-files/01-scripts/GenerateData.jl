#!/usr/bin/env julia1
push!(LOAD_PATH, "/home/amel/Documents/skyrmion-micromagnetics/julia-files/")

using Parameters, Distributed, SpinDynamics, BenchmarkTools

using InitialCondition, EffectiveField, Energy, EffectiveSize

j = 1.0
h = 0.0
a = 0.0
beta = 1.5
dz = 0.02
ed = beta/(4*pi*beta)
nx = ny = 128
nz = 1
pbc = 2.0 # 0.0=nothing, 1.0=pbc, 2.0=linear extrapolation

mpList = [j, h, a, dz, ed, nx, ny, nz, pbc] # material parameters list

tMax = 5
dt = 0.2
nSteps = 1
tol = 10^-5
damping = 0.0
T = 0.0
numRuns = 1
parallel = 0
numCores = 1
runRelaxation = 0 # 0=skip, 1=LLG with damping set to 1, 2=FA relaxation

# If using FA relaxation, the following few parameters also are required
constrained = true
targetR = 7.0


llgList = [tMax, dt, nSteps, tol, damping, T, numRuns, parallel, numCores,
    runRelaxation]

faList = [tMax, 0.3, 10, tol, numRuns, parallel, numCores]

type = "skyrmion"
r = 7
chirality = pi/2
icx = nx/2
icy = ny/2

icList = [type, r, chirality, icx, icy]

pinField = 0.0

pinList = [pinField]

defType = 0.0 # 0.0=None, 1.0=point, 2.0=Gaussian
defStrength = 0.0
defWidth = 0.0
defX = 0
defY = 0

defList = [defType, defStrength, defWidth, defX, defY]

saveE = 0
saveExc = 0
saveZee = 0
saveDmi = 0
savePma = 0
saveDdi = 0
saveMag = 0
saveSize = 0
saveQ = 0
savePos = 0
saveField = 0

saveList = [saveE, saveExc, saveZee, saveDmi, savePma, saveDdi, saveMag,
    saveSize, saveQ, savePos, saveField]

xCurr = 0.0
yCurr = 0.0
tOff = 0

currList = [xCurr, yCurr, tOff]

jRange = []
hRange = []
aRange = []
dzRange = []
edRange = []
nxRange = []
nyRange = []
nzRange = []
pbcRange = []
xCurrRange = []
yCurrRange = []
pxRange = []

rangeList = [jRange, hRange, aRange, dzRange, edRange, nxRange, nyRange, nzRange,
    pbcRange, xCurrRange, yCurrRange, pxRange]

# Generating data for different skyrmion locations and initial radii
# Storing results in 50 timestep increments
# Exporting each look in these increments in case something crashes

using RungeKutta, LLequation, TopologicalCharge, HDF5

nsteps=50

for r in 10:10
    for icx in nx/2:2:nx/2
        for icy in ny/2:2:ny/2
            icx = round(Int,icx)
            icy = round(Int,icy)

            icList = [type, r, chirality, icx, icy]
            params = buildparam(mpList, llgList, faList, icList, pinList,
                defList, saveList, rangeList, currList)
            s = buildinitial(params.ic, params.mp)

            datachunk = zeros(3,nx,ny,nsteps)

            for i in 1:nsteps
                rk4!(s, RHS!, params, false)
                datachunk[:,:,:,i] = s
                if i==nsteps
                    println("r = ",r,", Final Q = ", calcQ(s))
                end
            end

            datadir = "/media/amel/easystore1/data/"
            filename = string("rinit=",r,"_xinit=",icx,"_yinit=",icy,"_j=",j,"_h=",h,"_a=",a,"_dz=",dz,"_beta=",beta,".h5")

            h5write(string(datadir,filename),"in",datachunk)
        end
    end
end
