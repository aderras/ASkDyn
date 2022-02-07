module UserInputs

    module Material
        j = 1.0
        h = 0.0 # Negative field is into the plane. +z points out of page
        a = 0.0
        ed = 0.0
        dz = 0.0 #4*pi*ed # Negative dz is into the plane
        nx = 1024
        ny = nx
        nz = 1
        pbc = 0.0 # OPTIONS: Float between 1.0 and 6.0,
                  # 0.0=nothing, 1.0=pbc,
                  # 2.0/3.0 = 1st order backward euler, linear/quadratic interp
                  # 4.0/5.0 = 4th order backward euler, linear/quadratic interp
                  # 6.0/7.0 = 4th central difference, linear/quadratic interp
    end

    module Computation
        # Computation parameters
        solver = 0      # 0=RK4, 1=RK5
        maxSteps = 500    # This is the maximum number of computation steps
        dt = 0.2
        nSteps = 10
        tol = 10^-5
        damping = 0.0
        T = 0.0
        parallel = 0
        numCores = 1
        runRelaxation = 0  # 0=import data if it exists, otherwise build initial
                           # condition and run LLG, 1=LLG with damping set to
                           # 1.0, 2=FA relaxation
        sMax = 1
        faConst = 0.3
        nRot = 10
        faTol = 10^-4

        rChirality = 20    # MAKE SURE THIS MAKES SENSE
    end

    module InitialCondition
        import UserInputs.Material
        type = "skyrmion" # OPTIONS: "skyrmion", "domainWall"
        r = 12.0
        chirality = pi/2
        icx = Material.nx/2
        icy = Material.ny/2
    end

    module Pinning
        pinField = 0.0
    end

    module Defect
        defType = 0.0 # 0.0=None, 1.0=point, 2.0=Gaussian
        defStrength = 0.0
        defWidth = 0.0
        defX = 0
        defY = 0
    end

    module Current
        xCurr = 0.0
        yCurr = 0.0
        tOff = 0
    end

    module SaveChoices
        totalE = 0
        exchE = 0
        zeemE = 0
        dmiE = 0
        pmaE = 0
        ddiE = 0
        magn = 0
        size = 0
        charge = 0
        loc = 0
        fieldDuring = 0
        chirality = 0
    end

    #=
        If user would like to run the script for a range of values, store
        those values in the following struct. Elements of this struct are
        1d arrays.

        The name of each field in paramRanges must exactly match the name of
        the field in the original struct. All structs are located in Params.jl
    =#
    Base.@kwdef mutable struct paramRanges
        r = [5.0:20.0;]
        # damp = [0.01,1.0]
        # pbc = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
        # dz = [0.00001,0.0001,0.001,0.01,0.1]
        # h = [-0.0001,-0.001,-0.01,-0.1,0.01,0.001,0.0001]
    end

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
        labels = ["bc","nx","ny","nz","j","h","a","dz","ed","r","dt"]
        labels = [string("_",l) for l in labels]

        function values(p)
            return [p.mp.pbc, string(p.mp.nx), p.mp.ny, p.mp.nz, p.mp.j, p.mp.h,
                    p.mp.a, p.mp.dz, p.mp.ed, p.ic.r, p.cp.dt*p.cp.nn]
        end

        function outputSuffix(p)
            return string(zipFlatten(labels,values(p)),".h5")
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
        loc = "location"
        energy = "en"
        charge = "charge"
        finalField = "s_final"
        chirality = "gamma"


        # The names in this list have to match saveChoices options
        allFilenames = [energy,exch,zee,dmi,pma,ddi,magn,size,charge,loc,finalField,chirality]

    end

end
