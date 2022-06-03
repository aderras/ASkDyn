module UserInputs

    import Params
    import Dipolar, DefectFunctions

    function itemizedList(params...)
      #=
        Given an arbitrary list of items and lists of items, return a list of tuples
        containing combinations of every argument and each item in the input list.
      =#
      return collect(Iterators.product(params...))
    end

    function buildStructFromTuple(structName, inputTuple::Tuple)
      #=
        Given a struct and a tuple containing all elements of the struct, return
        the filled-in struct
      =#
      return structName(inputTuple...)
    end

    function buildItemizedStruct(structName, itemizedArgs)
      #=
        Given a struct name and a list of tuples each containing the arguments of
        the struct, return a list of structs of type structName.
      =#
      if length(itemizedArgs)==1 && !(typeof(itemizedArgs) <: Vector)
          return [buildStructFromTuple.(structName, itemizedArgs)]
      else
          return buildStructFromTuple.(structName, itemizedArgs)
      end
    end

    j = 1.0
    h = 0.0 # Negative field is into the plane. +z points out of page
    a = 0.0
    ed = 0.0
    dz = 0.0 #4*pi*ed # Negative dz is into the plane
    nx = 256
    ny = nx
    nz = 1
    bc = 0.0  # OPTIONS: Float between 1.0 and 6.0,
              # 0.0=nothing, 1.0=pbc,
              # 2.0/3.0 = 1st order backward euler, linear/quadratic interp
              # 4.0/5.0 = 4th order backward euler, linear/quadratic interp
              # 6.0/7.0 = 4th central difference, linear/quadratic interp
    if ed!=0.0
      v = [vdmatrices(UserInputs.Material.nx, UserInputs.Material.ny,
          UserInputs.Material.nz, 1.0==UserInputs.Material.pbc)]
    else
      v = [[[0.0 0.0; 0.0 0.0],[0.0 0.0; 0.0 0.0],[0.0 0.0; 0.0 0.0],
          [0.0 0.0; 0.0 0.0]]]
    end

    # Computation parameters
    solver = 0      # 0=RK4, 1=RK5
    maxSteps = 20    # This is the maximum number of computation steps
    dt = 0.2
    nSteps = 10
    tol = 10^-5
    damping = 0.0
    T = 0.0
    parallel = 0
    numCores = 2
    runRelaxation = 0  # 0=import data if it exists, otherwise build initial
                       # condition and run LLG, 1=LLG with damping set to
                       # 1.0, 2=FA relaxation
    sMax = 1
    faConst = 0.3
    nRot = 10
    faTol = 10^-4

    rChirality = 20    # MAKE SURE THIS MAKES SENSE

    type = ["skyrmion"] # OPTIONS: "skyrmion", "domainWall"
    r = 5.0
    chirality = pi/2
    icx = nx/2-10
    icy = ny/2

    pinField = 0.0

    defType = 2.0 # 0.0=None, 1.0=point, 2.0=Gaussian
    defStrength = -0.7
    defWidth = 1.0
    defX = nx/2
    defY = ny/2
    if defType!=0.0
        jmat = [DefectFunctions.buildJmats(nx,ny,defType,defStrength,defWidth,
                defX, defY)]
    else
        jmat=[[[],[]]]
    end

    xCurr = 0.0
    yCurr = 0.0
    tOff = 0

    # Parameter to save
    totalE = 1
    exchE = 0
    zeemE = 0
    dmiE = 0
    pmaE = 0
    ddiE = 0
    magn = 0
    size = 1
    charge = 1
    loc = 1
    fieldDuring = 0
    chirality = 0

    # All the paths for importing and exporting files
    module Paths
        relaxation = string(dirname(pwd()),"/02-data/relaxation/")
        output = string(dirname(pwd()),"/02-data/out/")
        initialConds = string(dirname(pwd()),"/02-data/store/")
    end

    module Filenames
        # timestampString = string(Dates.format(Dates.now(),"HH_MM_SS"),
        #     trunc(Int,rand()*10000))
        using Params, Helpers

        # Choose labels and parameters here. Make sure order and elements match
        # elements of values(p).
        labels = ["damp","t","bc","nx","ny","nz","j","h","a","dz","ed","r","dt"]
        labels = [string("_",l) for l in labels]

        function values(p)
            return [p.cp.damp, p.cp.temp, p.mp.bc, string(p.mp.nx), p.mp.ny,
                    p.mp.nz, p.mp.j, p.mp.h, p.mp.a, p.mp.dz, p.mp.ed, p.ic.r,
                    p.cp.dt*p.cp.nn]
        end

        function outputSuffix(p)
            return zipFlatten(labels,values(p))
        end

        # If you would like to import pre-computed relaxation data and run
        # dynamics, change the filename here.
        function inputName(p)
            return string("s_final",zipFlatten(labels,values(p)),".h5")
        end

        exch = "exchange"
        zee = "zeeman"
        dmi = "dmi"
        pma = "pma"
        ddi = "ddi"
        magn = "magnetization"
        size = "size"
        loc = "loc"
        energy = "en"
        charge = "charge"
        finalField = "s_final"
        chirality = "gamma"

        # The names in this list have to match saveChoices options
        allFilenames = [energy,exch,zee,dmi,pma,ddi,magn,size,charge,loc,finalField,chirality]
    end


    ## END USER INPUTS ########################################################

    # Create structs to be used by the program
    mpItems = itemizedList(j, h, a, ed, dz, nx, ny, nz, bc, v)
    cpItems = itemizedList(solver, maxSteps, dt, nSteps, tol, damping,
        T, parallel, numCores, runRelaxation, sMax, faConst, nRot, faTol,
        rChirality)
    icItems = itemizedList(type, r, chirality, icx, icy)
    pinItems = itemizedList(pinField)
    defectItems = itemizedList(defType, defStrength, defWidth, defX, defY, jmat)
    currItems = itemizedList(xCurr, yCurr, tOff)
    saveItems = itemizedList(totalE, exchE, zeemE, dmiE, pmaE, ddiE, magn, size, charge,
                loc, fieldDuring, chirality)

    mp = buildItemizedStruct(Params.materialParams, mpItems)
    cp = buildItemizedStruct(Params.compParams, cpItems)
    ic = buildItemizedStruct(Params.icParams, icItems)
    pin = buildItemizedStruct(Params.pinningParams, pinItems)
    dp = buildItemizedStruct(Params.defectParams, defectItems)
    curr = buildItemizedStruct(Params.currentParams, currItems)
    save = buildItemizedStruct(Params.saveChoices, saveItems)

    paramItems = itemizedList(mp, cp, ic, pin, dp, curr, save)
    paramObjs = buildItemizedStruct(Params.params, paramItems)
    if length(paramObjs)==1 && (typeof(paramObjs[1]) <: Array)
        paramObjs = paramObjs[1]
    end
end
