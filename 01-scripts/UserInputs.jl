module UserInputs

    module Material
        j = 1.0
        h = 0.0
        a = 0.0
        ed = 0.0
        dz = 4*pi*ed
        nx = ny = 512
        nz = 1
        pbc = 0.0 # 0.0=nothing, 1.0=pbc, 2.0=linear extrapolation
    end

    module Computation
        # Computation paramters
        tMax = 1
        dt = 0.2
        nSteps = 10
        tol = 10^-5
        damping = 0.1
        T = 0.0
        parallel = 0
        numCores = 1
        runRelaxation = 0   # 0=import data if it exists, otherwise build initial
                            # condition and run LLG, 1=LLG with damping set to
                            # 1.0, 2=FA relaxation

        sMax = 1000
        faConst = 0.3
        nRot = 10
        faTol = 10^-4

        rChirality = 100
    end

    module InitialCondition
        import UserInputs.Material
        type = "skyrmion"
        r = 12
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
        rRange=[8.0:12.0;]
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

        starting = string(dirname(pwd()),"/02-data/store/")
    end

    module Filenames
        # timestampString = string(Dates.format(Dates.now(),"HH_MM_SS"),
        #     trunc(Int,rand()*10000))

        # If you would like to import pre-computed relaxation data and run
        # dynamics, change the filename here
        function outputSuffix(p)
            return string("_NX=",p.mp.nx,"_NY=",p.mp.ny,"_NZ=",p.mp.nz,"_J=",p.mp.j,
            "_H=",p.mp.h,"_A=",p.mp.a,"_DZ=",round(p.mp.dz,digits=5),
            "_BETA=",p.mp.ed,"_R=",p.ic.r,".h5")
        end
        function inputName(p)
            return string("S_1.0",outputSuffix(p))
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
