module UserInputs

    module Material
        j = 1.0
        h = 0.0
        a = 0.0
        ed = 0.0
        dz = 4*pi*ed
        nx = 512
        ny = nx
        nz = 1
        pbc = 5.0 # OPTIONS: Float between 1.0 and 6.0,
                  # 0.0=nothing, 1.0=pbc,
                  # 2.0/3.0 = 1st order backward euler, linear/quadratic interp
                  # 4.0/5.0 = 4th order backward euler, linear/quadratic interp
                  # 6.0/7.0 = 4th central difference, linear/quadratic interp
    end

    module Computation
        # Computation paramters
        tMax = 10
        dt = 0.2
        nSteps = 10
        tol = 10^-5
        damping = 0.0
        T = 0.0
        parallel = 0
        numCores = 1
        runRelaxation = 0   # 0=import data if it exists, otherwise build initial
                            # condition and run LLG, 1=LLG with damping set to
                            # 1.0, 2=FA relaxation

        sMax = 1
        faConst = 0.3
        nRot = 10
        faTol = 10^-4

        rChirality = 20    #MAKE SURE THIS MAKES SENSE
    end

    module InitialCondition
        import UserInputs.Material
        type = "skyrmion"
        r = 7.0
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

    module Ranges
        r = []
        pbc = []
    end

    module SaveChoices
        totalE = 1
        exchE = 0
        zeemE = 0
        dmiE = 0
        pmaE = 0
        ddiE = 0
        magn = 1
        size = 1
        charge = 1
        loc = 0
        fieldDuring = 0
        chirality = 1
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

        # If you would like to import pre-computed relaxation data and run
        # dynamics, change the filename here
        function outputSuffix(p)
            return string("_BC=",p.mp.pbc,"_NX=",p.mp.nx,"_NY=",p.mp.ny,"_NZ=",
            p.mp.nz,"_J=",p.mp.j,"_H=",p.mp.h,"_A=",p.mp.a,"_DZ=",
            round(p.mp.dz,digits=5),"_ED=",p.mp.ed,"_R=",p.ic.r,"_DT=",
            p.cp.nn*p.cp.dt,".h5")
        end
        function inputName(p)
            return string("S_CONSTRAINED_NX=",p.mp.nx,"_NY=",p.mp.ny,"_NZ=",
            p.mp.nz,"_J=",p.mp.j,"_H=",p.mp.h,"_A=",p.mp.a,"_DZ=",
            round(p.mp.dz,digits=5),"_ED=",p.mp.ed,"_R=",p.ic.r,".h5")
        end

        exch = "EXCHANGE"
        zee = "ZEEMAN"
        dmi = "DMI"
        pma = "PMA"
        ddi = "DDI"
        magn = "MAGNETIZATION"
        size = "SIZE"
        loc = "LOCATION"
        energy = "EN"
        charge = "CHARGE"
        finalField = "S_FINAL"
        chirality = "GAMMA"

        # The names in this list have to match saveChoices options
        allFilenames = [energy,exch,zee,dmi,pma,ddi,magn,size,charge,loc,finalField,chirality]

    end

end
