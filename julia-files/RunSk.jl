#!/usr/bin/env julia1
push!(LOAD_PATH, pwd())

using UserInputs, Distributed, SpinDynamics, BenchmarkTools

using InitialCondition, EffectiveField

j = 1.0
h = -0.0065
a = 0.04
ed = 0.001
dz = 4*pi*ed
nx = ny = 64
nz = 10
pbc = 2.0 # 0.0=nothing, 1.0=pbc, 2.0=linear extrapolation

mpList = [j, h, a, dz, ed, nx, ny, nz, pbc]

tMax = 5
dt = 0.2
nSteps = 5
tol = 10^-5
damping = 0.1
T = 0.0
numRuns = 1
parallel = 0
numCores = 1
runRelaxation = 1 # 0=skip, 1=LLG with damping set to 1, 2=FA relaxation

llgList = [tMax, dt, nSteps, tol, damping, T, numRuns, parallel, numCores,
    runRelaxation]

faList = [tMax, 0.3, 10, tol, numRuns, parallel, numCores]

type = "skyrmion"
r = 10
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

saveE = 1
saveExc = 0
saveZee = 0
saveDmi = 0
savePma = 0
saveDdi = 0
saveMag = 1
saveSize = 1
saveQ = 1
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

testParams = buildparam(mpList, llgList, faList, icList, pinList,
    defList, saveList, rangeList, currList)

confirm = getconfirmation(testParams)
if confirm == false
    exit()
end

if testParams.llg.parallel == 1

  addprocs(testParams.llg.numCores)
  println("Initializing parallel cores. Number of processors = ", nprocs(),
        ", number of workers = ", nworkers())

    @everywhere module Test
    push!(LOAD_PATH, pwd())

    using UserInputs, SpinDynamics, Distributed

    function dosomething(testParams)

      # All params is a list of all combinations of the requested inputs
      allParams = []
      getparamlist!(testParams, allParams)

      pmap(launchcomputation, allParams)

    end

  end

  using .Test
  Test.dosomething(testParams)

else

  allParams = []
  getparamlist!(testParams, allParams)

  map(launchcomputation, allParams)

end
