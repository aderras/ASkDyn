#=

    This module contains structs and functions related to parameters that are
    passed to the spin dynamics computation functions.

    Params are grouped by catagory and stored in structs. E.g. the
    materialParams struct contains the values for exchange interaction,
    external field, shape of the array, etc.

    All of the parameters are then group into one parent struct of the type
    "params." This parent struct is passed to different functions, and can
    be extended to more specific problems by including additional structs.
    E.g. the parameters for a defect in the lattice are stored in the
    defectParams struct, which is an element of params.

    Below the structs are the functions

    IMPROVEMENTS:
        + add warning if user requests RK discretization that doesn't make sense
        + also warn if the range of parameters requested doesn't make sense

=#
module Params

    using Dipolar, EffectiveField, UserInputs
    export getParams, buildparam, params, getparamlist!,
        setfield_nestedstruct!, getconfirmation
    export materialParams, llgParams, faParams, icParams, pinningParams,
        defectParams, currentParams, saveChoices, params

    # Material parameters
    mutable struct materialParams

        j::Float64      # Exchange interaction
        h::Float64      # Zeeman field
        a::Float64      # Dzayaloshinskii-Moriya interaction
        dz::Float64     # Perpendicular magnetic anisotropy
        ed::Float64     # Dipole-dipole constant
        nx::Int         # Lattice size in x
        ny::Int         # Lattice size in y
        nz::Int         # Lattice size in z
        pbc::Float64    # Periodic boundary conditions

        v               # Matrices to be used to compute DDI. It
                        # is advantageous to compute these once and
                        # store them because this is slow.
    end

    mutable struct compParams
        maxSteps::Int         # Maximum number of steps to make
        dt::Float64       # Step size
        nn::Float64       # Skip this many steps to save
        tol::Float64      # Tolerance for convergence
        damp::Float64     # Damping in llg
        temp::Float64     # Temperature
        parallel::Int     # Parallel (1 or 0)
        numCores::Int     # Number of parallel cores
        relax             # Run relaxation if data doesn't exist

        # These are for FA relaxation, if using
        sMax::Int         # Maximum number of steps to make
        faconst::Float64  # Field alignment param
        nRot::Float64     # Number of rotations to make
        faTol::Float64    # Tolerance for convergence

        rChirality::Int
    end

    # Initial condition parameters
    mutable struct icParams

        type::String        # Type of spin structure
        r::Float64          # Radius (only makes sense for skyrmion)
        gamma::Float64      # Chirality (only skyrmion)
        px         # Position in x direction (0 < px < nx)
        py         # Position in y direction (0 < py < ny)

    end

    # Pinning parameters
    mutable struct pinningParams
        hPin::Float64       # Pinning field strength
    end

    # Defect parameters
    mutable struct defectParams

        t                   # Type of defect (1 = Point, 2 = Gaussian)
        strength::Float64   # Strength of the defect
        width::Float64      # Width of impact
        dx::Int           # x position (0 < dx < nx)
        dy::Int           # y position (0 < dy < ny)
        jMat
    end

    mutable struct currentParams
        jx::Float64          # Current in the x direction
        jy::Float64          # y direction
        tf::Float64          # Cutoff time for current application
    end

    # Choices for data saving
    mutable struct saveChoices

        totE::Int       # Total energy
        excE::Int       # Exchange energy
        zeeE::Int       # Zeeman energy
        dmiE::Int       # DMI energy
        pmaE::Int       # PMA energy
        ddiE::Int       # DDI energy
        magn::Int       # Magnetization
        size::Int       # Effective size (skyrmion)
        qCharge::Int    # Topological charge
        location::Int   # Position (skyrmion)
        spinField::Int  # Spin field (The whole matrix. Data intensive.)
        chir::Int       # Chirality
    end

    # All parameters are stored in this struct.
    mutable struct params

        mp          # materialParams
        cp         # computationParams
        ic          # icParams
        pin         # pinningParams
        defect      # defectParams
        current     # current
        save        # saveChoices

    end

    # getTestParams creates a params struct to use for testing, so that the user
    # does not need to navigate through user input every time. This function is
    # also used to create a params struct of the same structure that is used in
    # the program. It is useful in some areas when copying between structs has
    # to be done.
    function buildparam()

        myMatParams = materialParams(0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, [])
        myCompParams = compParams(0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0,0.0,0,0.0,0)
        myICParams = icParams("", 0, 0.0, 0, 0)
        myPinningParams = pinningParams(0.0)
        myDefectParams = defectParams(0, 0.0, 0.0, 0, 0, [])
        mySaveChoices = saveChoices(0,0,0,0,0,0,0,0,0,0,0,0)
        myCurrentChoices = currentParams(0.0, 0.0, 100)

        # Now put all of the user choices into one struct
        allParams = params(myMatParams, myCompParams, myICParams,
            myPinningParams,myDefectParams, myCurrentChoices,
            mySaveChoices)

        return allParams

    end

    # Given a field and value, this searches the nested struct for approapriate
    # location of the value. Works recursively. And slowly.
    function setfield_nestedstruct!(nestedStruct, fieldname, value)

        names = fieldnames(typeof(nestedStruct))

        # If there are no fields, you reached the end of the nested structure
        if length(names) == 0
            return
        # If the field is in this list, then set the value to what you want
        elseif fieldname in names
            setfield!(nestedStruct, fieldname, value)
            return
        # If the field is not in this level, it means there is another level to
        # recurse through.
        else
            for f in names
                setfield_nestedstruct!(getfield(nestedStruct, f), fieldname,
                    value)
            end
        end
    end

    # Deep copies struct of structs. Works recursively. Source and dest are
    # structs with the same structure.
    #
    # There is probably a native Julia way to do this. Couldn't find it.
    function copystruct!(source, dest)
        fields = fieldnames(typeof(source))

        # If there are no fields, you reached the end of the nested structure
        if length(fields) == 0
            return
        # If the parent module of this field is not Core, there is another
        # level of nested structs to dig through. Call copystruct! again.
        elseif parentmodule(typeof(getfield(source, fields[1]))) != Core
            # For all the fields, recursively call copyStruct
            for f in fields
                copystruct!(getfield(source, f), getfield(dest, f))
            end

        # If the parent module of this level is Core, you've reached the values
        # to copy, so do it.
        else
            # For all the fields, copy values from source to destination
            for f in fields
                setfield!(dest, f, getfield(source, f))
            end
        end
    end

    # Set field of a nested struct given multiple assignments
    function setfield_multi!(paramStruct, fields, values)
        for fv in collect(zip(fields,values))
            setfield_nestedstruct!(paramStruct, fv[1], fv[2])
        end
    end

    # Builds a list of all the compuatation parameters requested by the user. It
    # does this by reursively sesarching through the user requests to make
    # combinations of all possible  Stores all user requests in the
    # input array, paramList
    function getparamlist!(allParams, paramList, rangeParams, iters=[],
        iCurr=[], count=1)

        # Get field names of the range struct
        fields = fieldnames(typeof(rangeParams))
        vals=[getfield(rangeParams,f) for f in fields]
        if all(length.(vals).==0) push!(paramList,allParams); return end
        iters=[1:length(vec) for vec in vals]


        for i in iters[count]
            if count==1 iCurr=[] else iCurr=iCurr[1:(count-1)] end
            push!(iCurr, i)

            if count==length(iters)
                values=[vals[k][iCurr[k]] for k in 1:length(iCurr)]
                newParamElem = buildparam()
                copystruct!(allParams, newParamElem)
                setfield_multi!(newParamElem, fields, values)
                push!(paramList, newParamElem)
            else
                getparamlist!(allParams, paramList, rangeParams, iters, iCurr, count+1)
            end
        end
    end

    function buildUserInputParam()

        # Only build the following matrices if the dipolar interaction is nonzero
        if UserInputs.Material.ed!=0.0
            v = vdmatrices(UserInputs.Material.nx, UserInputs.Material.ny,
                UserInputs.Material.nz, 1.0==UserInputs.Material.pbc)
        else
            v = [[0.0 0.0; 0.0 0.0],[0.0 0.0; 0.0 0.0],[0.0 0.0; 0.0 0.0],
                [0.0 0.0; 0.0 0.0]]
        end

        # Create all the structs to be used by the rest of the program.
        mp = materialParams(UserInputs.Material.j,
            UserInputs.Material.h, UserInputs.Material.a, UserInputs.Material.dz,
            UserInputs.Material.ed, UserInputs.Material.nx, UserInputs.Material.ny,
            UserInputs.Material.nz, UserInputs.Material.pbc, v)

        cp = compParams(UserInputs.Computation.maxSteps, UserInputs.Computation.dt,
            UserInputs.Computation.nSteps, UserInputs.Computation.tol,
            UserInputs.Computation.damping, UserInputs.Computation.T,
            UserInputs.Computation.parallel, UserInputs.Computation.numCores,
            UserInputs.Computation.runRelaxation, UserInputs.Computation.sMax,
            UserInputs.Computation.faConst,UserInputs.Computation.nRot,
            UserInputs.Computation.faTol, UserInputs.Computation.rChirality)

        ic = icParams(UserInputs.InitialCondition.type,
            UserInputs.InitialCondition.r, UserInputs.InitialCondition.chirality,
            UserInputs.InitialCondition.icx, UserInputs.InitialCondition.icy)

        pp = pinningParams(UserInputs.Pinning.pinField)

        dp = defectParams(UserInputs.Defect.defType,
            UserInputs.Defect.defStrength, UserInputs.Defect.defWidth,
            UserInputs.Defect.defX, UserInputs.Defect.defY, [])

        cc = currentParams(UserInputs.Current.xCurr,
            UserInputs.Current.yCurr, UserInputs.Current.tOff)

        sp = saveChoices(UserInputs.SaveChoices.totalE,
            UserInputs.SaveChoices.exchE, UserInputs.SaveChoices.zeemE,
            UserInputs.SaveChoices.dmiE, UserInputs.SaveChoices.pmaE,
            UserInputs.SaveChoices.ddiE, UserInputs.SaveChoices.magn,
            UserInputs.SaveChoices.size, UserInputs.SaveChoices.charge,
            UserInputs.SaveChoices.loc, UserInputs.SaveChoices.fieldDuring,
            UserInputs.SaveChoices.chirality)


        userParams = params(mp, cp, ic, pp, dp, cc, sp)

        return userParams

    end

end
