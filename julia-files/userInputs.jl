#=

    This module contains structs and functions related to parameters that are 
    passed to the spin dynamics computation functions. 

    Parameters are grouped by catagory and stored in structs. E.g. the 
    materialParams struct contains the values for exchange interaction, 
    external field, shape of the array, etc. 

    All of the parameters are then group into one parent struct of the type 
    "params." This parent struct is passed to different functions, and can 
    be extended to more specific problems by including additional structs. 
    E.g. the parameters for a defect in the lattice are stored in the 
    defectParams struct, which is an element of params.

    Below the structs are the functions 

    Next steps:
        + add warning if user requests RK discretization that doesn't make sense
        + also warn if the range of parameters requested doesn't make sense 

=#
module userInputs

    using dipoleDipole
    export getUserParams, getTestParams, params, copyStruct!,
    getParamList!

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

        vdd             # Matrices to be used to compute DDI. It 
                        # is advantageous to compute these once and 
                        # store them because this is slow.

        # Create the struct given an array of appropriate length 
        function materialParams( arr )

            if arr[5] != 0.0
                v = vddMatrices( Int(arr[6]) , Int(arr[7]), Int(arr[8]), 1.0==arr[9])
                return new( arr[1], arr[2], arr[3], arr[4], arr[5], 
                    Int(arr[6]), Int(arr[7]), Int(arr[8]), arr[9], v )
            else
                return new( arr[1], arr[2], arr[3], arr[4], arr[5], 
                    Int(arr[6]), Int(arr[7]), Int(arr[8]), arr[9], [])       
            end 
        end

    end

    # If user would like to run the script for a range of values, store 
    # those values in the following struct. Elements of this struct are 
    # lists of the same physical values as labeled in "materialParams"
    mutable struct paramRanges

        j
        h
        a
        dz
        ed
        nx
        ny
        nz
        pbc
        jx  # current in x
        jy  # current in y

        function paramRanges( arr )
            return new( arr[1], arr[2], arr[3], arr[4], arr[5], arr[6],
                arr[7], arr[8], arr[9], arr[10], arr[11] )
        end


    end

    # Parameters for the Landau-Lifshitz-Gilbert calculation
    mutable struct llgParams

        tMax::Int64         # Maximum number of steps to make
        dt::Float64         # Step size 
        saveEvery::Float64  # Skip this many steps to save 
        tol::Float64        # Tolerance for convergence
        damp::Float64       # Damping in llg
        temp::Float64       # Temperature 
        nRuns::Int64        # Number of times to run 
        parallel::Int64     # Parallel (1 or 0)
        numCores::Int64     # Number of parallel cores 

        function llgParams( arr )
            return new( Int(arr[1]), arr[2], arr[3], arr[4], arr[5], arr[6],
                Int(arr[7]), Int(arr[8]), Int(arr[9]) )
        end

    end

    # Parameters for field-alignment algorithm
    mutable struct faParams

        sMax::Int64         # Maximum number of steps to make 
        param::Float64      # Field alignment param 
        tol::Float64        # Tolerance for convergence 
        temp::Float64       # Temperature
        nRuns::Float64      # Number of times to run 
        parallel::Int64     # Parallel (1 or 0)
        numCores::Int64     # Number of parallel cores 

        function faParams( arr )
            return new( Int(arr[1]), arr[2], arr[3], arr[4], arr[5], Int(arr[6]), 
                Int(arr[7]) )
        end

    end

    # Initial condition parameters 
    mutable struct icParams

        type::String        # Type of spin structure 
        r::Float64          # Radius (only makes sense for skyrmion)
        gamma::Float64      # Chirality (only skyrmion)
        px::Int64           # Position in x direction (0 < px < nx)
        py::Int64           # Position in y direction (0 < py < ny)
        
        function icParams( arr )
            return new( arr[1], arr[2], arr[3], arr[4], arr[5] )
        end
    end

    # Pinning parameters
    mutable struct pinningParams

        hPin::Float64       # Pinning field strength
        px::Int64           # Position of pinning field in x (0 < px < nx)
        py::Int64           # Position of pinning field in y (0 < py < ny)
        
        function pinningParams( arr )
            return new( arr[1], arr[2], arr[3] )
        end

    end

    # Defect parameters
    mutable struct defectParams

        t::Float64          # Type of defect (1 = Point, 2 = Gaussian)
        strength::Float64   # Strength of the defect
        width::Float64      # Width of impact
        dx::Int64           # x position (0 < dx < nx)
        dy::Int64           # y position (0 < dy < ny)

        function defectParams( arr )
            return new( arr[1], arr[2], arr[3], arr[4], arr[5] )
        end

    end

    mutable struct currentParams

        jx::Float64          # Current in the x direction
        jy::Float64          # y direction
        tf::Float64          # Cutoff time for current application

        function currentParams( arr )
            return new( arr[1], arr[2], arr[3] )
        end

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

        function saveChoices( arr )
            return new( arr[1], arr[2], arr[3], arr[4], arr[5],
                arr[6], arr[7], arr[8], arr[9], arr[10], arr[11] )
        end

    end

    # All parameters are stored in this struct.
    mutable struct params

        mp          # materialParams
        llg         # llgParams
        fa          # faParams
        ic          # icParams
        pin         # pinningParams
        defect      # defectParams
        save        # saveChoices
        range       # paramRanges
        current     # current

        function params( arr )
            return new( arr[1], arr[2], arr[3], arr[4], arr[5], arr[6], arr[7], arr[8], arr[9] )
        end
    end

    #######################################################################################
    # Begin functions used in building the parameter struct 

    # getUserInput asks for user input and tests whether it's of the right type
    #
    # in: T = type of input to request, msg = message to send to user 
    # out: Value input of type T
    function getUserInput( T=String, msg="" )
        
        print( "$msg" )

        # If requesting a string, just return it. 
        if T == String
            return readline()

        # If requesting something else, try to parse the answer. If there's a problem, 
        # try again. 
        else
            try
                return parse( T, readline() )
            catch
                println("I could not interpret your answer. Please try again. ")
                getUserInput(T, msg)
            end
        end
    end


    # Request a series of parameters with the following function.
    #
    # in: paramType = type of parameter struct you will be constructing,
    # paramList = array of strings containing parameters to request, 
    # defaultList = default values ot paramType. Types of defaultList must 
    # reflect the type required in the parameter struct 
    # out: struct paramType with either user inputs or default values
    function getComputationParams( paramType, paramList::Array{String,1}, defaultList)

        ans = getUserInput( String, string("\nEnter paramters ", [i for i in paramList], "? y/n ")  )
        
        # If user wants to enter parameters, get them
        if ans == "y"
            paramVals = zeros( length(paramList) )

            for i in 1:length(paramList)

                paramVals[i] = getUserInput( typeof(defaultList[i]), string( paramList[i], " = ") )

            end
            return paramType( paramVals )

        # Otherwise set default values
        else
            
            println( string("Setting default values ", [j for j in defaultList] ) )
            return paramType( defaultList )

        end      

    end

    # In order to set ranges for parameters, have to modify the rangeParams struct. Do that 
    # with the following function. 
    function modifyMatParam!(mat, val, num)

        if num == 1
            mat.j = val
        elseif num == 2
            mat.h = val
        elseif num == 3
            mat.a = val
        elseif num == 4
            mat.dz = val
        elseif num == 5
            mat.ed = val
        elseif num == 6
            mat.nx = val
        elseif num == 7
            mat.ny = val
        elseif num == 8
            mat.nz = val
        elseif num == 9
            mat.pbc = val
        elseif num == 10
            mat.jx = val
        elseif num == 11
            mat.jy = val
        else
            println(string("Error modifying material parameter struct. No ", num, "th value found. Quitting."))
            quit()
        end

    end


    # getTestParams creates a params struct to use for testing, so that the user does not 
    # need to navigate through user input every time. This function is also used to create a 
    # params struct of the same structure that is used in the program. It is useful in some 
    # areas when copying between structs has to be done. 
    function getTestParams()

        myMatParams = materialParams([1.0, -0.03, 0.05, 0.0, 0.0, 64, 64, 1, true] )
        myLlgParams = llgParams([1000, 0.2, 5.0, 10^-6, 0.1, 0.0, 1, 0, 0])
        myFaParams = []
        myICParams = icParams( ["skyrmion", 5, pi/2 , myMatParams.nx/2, myMatParams.ny/2] )
        myPinningParams = pinningParams( [0.0, 32, 32] )
        myDefectParams = defectParams( [0.0, 0.0, 0.0, 0, 0] )
        mySaveChoices = saveChoices( [1,0,0,0,0,0,0,1,1,1,0] )
        myParamRanges = paramRanges( [ [], [], [], [], [], [], [], [], [], [], [], [] ] )
        myCurrentChoices = currentParams( [ 0.2, -0.1, 100 ] )

        # Now put all of the user choices into one struct 
        allParams = params( [myMatParams, myLlgParams, myFaParams, myICParams, myPinningParams, 
            myDefectParams, mySaveChoices, myParamRanges, myCurrentChoices] )

        return allParams

    end

    # This function prompts the user for different computational values, filling the structs at the top of 
    # this script.
    function getUserParams()
        
        # Ask for material parameters 
        myMatParams = getComputationParams( materialParams, 
            ["J", "H", "A", "Dz", "Ed", "Nx", "Ny", "Nz","pbc" ], 
            [1.0, -0.2, 0.1, 0.01, 0.0, 64, 64, 1, true] )

        ans = getUserInput( String, 
        string("\nWould you like to run the computation over a range of parameters? 
        Enter n to choose 'no' and run for only selected choices. y/n ")  )

        myParamRanges = paramRanges( [ [], [], [], [], [], [], [], [], [] ] )

        if ans == "y"
            ans = getUserInput( String, 
            string("\nSelect which values you would like to vary:
                1 = Exchange interaction
                2 = Zeeman Field
                3 = DMI constant
                4 = PMA constant
                5 = DDI constant
                6 = Nx (Lattice size in x)
                7 = Ny (Lattice size in y)
                8 = Nz
                9 = pbc
                10 = Current in x direction
                11 - Current in y direction
            
                Enter your choice(s) separated by commas: ")  )

            varyChoices = parse.(Int64, split(ans,",") )

            for var in varyChoices

                ans = getUserInput( String, 
                    string( "\nEnter the start, step size, and end value of choice ", 
                        var, " separated by commas " ) )
                
                paramRanges = parse.( Float64, split(ans,",") )
                paramVals = [ paramRanges[1]:paramRanges[2]:paramRanges[3]; ]

                modifyMatParam!( myParamRanges, paramVals, var )

            end
        elseif ans == "n" || ans == ""
            # Keep going
        else 
            println( "I did not understand your request. Let's try again." )
            return getUserInput()
        end


        ans = getUserInput( String, string( "\nComputation type? 
            1 = Landau-Lifshitz-Gilbert 
            2 = Field Alignment relaxation ") )

        if ans == "1"

            myLlgParams = getComputationParams( llgParams, [ "Max iterations", "Step size", 
                "Skip size", "Tolerance", "Damping", "Temperature", "Number of runs", 
                "Parallelize (1 or 0)", "If parallel, number of cores"], 
                [10000, 0.2, 5.0, 10^-6, 1.0, 0.0, 1, 0, 1] )
            myFaParams = []

        elseif ans == "2"

            myFaParams = getComputationParams( faParams, ["Max iterations", "FA Param","Num rotations", 
                "Temperature", "Number of runs", "Parallelize (1 or 0)"], [10000, 0.03, 10, 0.0, 1, 0] )
            myLlgParams = []

        else
            println("Invalid input. Choosing default, LLG.")
            myLlgParams = getComputationParams( llgParams, [ "Max iterations", "Step size", 
                "Skip size", "Tolerance", "Damping", "Temperature", "Number of runs", 
                "Parallelize (1 or 0)", "If parallel, number of cores"], 
                [10000, 0.2, 5.0, 10^-6, 1.0, 0.0, 1, 0, 1] )
            myFaParams = []

        end

        # Initial condition
        myICParams = getComputationParams( icParams, 
            ["Initial condition", "lambda", "gamma", "x position", "y position"], 
            ["skyrmion", 5, pi/2 , myMatParams.nx/2, myMatParams.ny/2] )


        ################################################
        # Problem-specific options start here
        # Ask if the user wants to pin the skyrmion
        ans = getUserInput( String, string("\nPin skyrmion? y/n ")  )
            
        # Get params. Default set pin field at center of the lattice
        if ans == "y"
            myPinningParams = getComputationParams( pinningParams, ["H_pin", "X", "Y"], 
                [0.01, myMatParams.nx/2 , myMatParams.ny/2] )
        else
            myPinningParams = []
        end

        # Is there a defect in the lattice
        ans = getUserInput( String, string("\nPlace defect in the lattice? y/n ")  )

        # If user wants to enter defect parameters, get them
        if ans == "y"
            myDefectParams = getComputationParams( defectParams, 
                ["type", "strength (relative to J)", "width", "x position", "y position"], 
                ["gaussian", -0.5, 1.0, 0, 0] )
        else
            myDefectParams = []
        end


        # Is there a defect in the lattice
        ans = getUserInput( String, string("\nIntroduce current in lattice? y/n ")  )

        # If user wants to introduce a current, get values 
        if ans == "y"
            myCurrentParams = getComputationParams( defectParams, 
                ["Ix", "Iy"], 
                [0.2, -0.1] )
        else
            myCurrentParams = []
        end


        # Choose what to save
        ans = getUserInput( String, string("\nSelect which values you would like to save",
            " during the computation. (Note that choosing too many may increase computation",
            " time.) The options are: 
            1 = Total energy
            2 = Exchange Energy
            3 = Zeeman energy
            4 = DMI energy
            5 = PMA energy
            6 = DDI energy
            7 = Magnetization
            8 = Effective Size
            9 = Topological charge
            10 = Postion (only for skyrmion initial condition)
            11 = Spin field
            
            Enter your choices separated by commas: ")  )

        # If user does not select anything, move on
        if ans == ""
            println("\nComputation will not save data.")
            mySaveChoices = []

        else 

            # If user makes multiple choices, parse them into an Int array
            if occursin( ",",ans )
                userChoices = parse.(Int64, split(ans,",") )

            # If user makes 1 choice, store it in a 1 element array
            else
                userChoices = [ parse(Int64,ans) ]
            end

            userChoicesBool = zeros(Int64, 11)
            [userChoicesBool[i] = 1 for i in userChoices]
            mySaveChoices = saveChoices( userChoicesBool )

        end

        # Now put all of the user choices into one struct 
        allParams = params( [myMatParams, myLlgParams, myFaParams, myICParams, 
            myPinningParams, myDefectParams, mySaveChoices, myParamRanges] )

        return allParams

    end


    # Builds a list of all the compuatation parameters requested by the user. It 
    # does this by reursively sesarching through the user requests to make 
    # combinations of all possible parameters. Stores all user requests in the 
    # input array, allParams
    function getParamList!( par, allParams )

        # Get all the ranges requested. These are a stored in a struct. 
        rangeParams = par.range 

        # Get field names of the range struct
        fields = fieldnames( typeof(rangeParams) )

        # Get the length of all the range params to see if user selected any ranges
        lenRangeParams = [ length( getfield( rangeParams, x ) ) for x in fields ]

        if sum( lenRangeParams ) > 0    

        # find the first nonzero range 
        rangePos = findall(x->x != 0, lenRangeParams)[1]

        currentRange = getfield( rangeParams, fields[ rangePos ] )

        # Make a copy of the params struct. Call copyStruct to deep copy
        copyParams = getTestParams()
        copyStruct!( par, copyParams )

        # Set the range in the copied struct to empty
        setfield!( copyParams.range, fields[rangePos], [] )

        # For all the values in the range, set the field and call run comp 
        # again to recurse through the rest of the ranges. 
        for i in range( 1, stop=length( currentRange ) )
            
            setfield!( copyParams.mp, fields[rangePos], currentRange[i] )
        
            # Then recurseRun again, until all of the ranges have 
            # been cycled through.
            getParamList!( copyParams, allParams )

        end

        # If there are no ranges, then store the params 
        else

            # Need to copy the struct first so that it won't pass by ref 
            copyParams = getTestParams()
            copyStruct!( par, copyParams )

            push!( allParams, copyParams )

        end
    end

    # Deep copies struct of structs. Works recursively. Source and dest are
    # structs with the same structure. 
    #
    # There is probably a native Julia way to do this. Couldn't find it. 
    function copyStruct!( source, dest )

        fields = fieldnames( typeof( source ) )

        # If there are no fields, you reached the end of the nested structure
        if length( fields ) == 0
            
            return

        # If the parent module of this field is not Core, there is another
        # level of nested structs to dig through. Call copyStruct! again.
        elseif parentmodule( typeof( getfield( source, fields[1] ) ) ) != Core 

            # For all the fields, recursively call copyStruct
            for f in fields 

                copyStruct!( getfield( source, f), getfield( dest, f) )

            end

        # If the parent module of this level is Core, you've reached the values 
        # to copy, so do it. 
        else

            # For all the fields, copy values from source to destination
            for f in fields 

                setfield!( dest, f, getfield( source, f)  )

            end

        end
    end
 
end
