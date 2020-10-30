
module userInputs

    using dipoleDipole
    export getUserParams, getTestParams, params, copyStruct!,
    getParamList!

    mutable struct materialParams

        j::Float64
        h::Float64
        a::Float64
        dz::Float64
        ed::Float64
        nx::Int
        ny::Int
        nz::Int
        pbc::Float64
        vdd

        function materialParams( arr )

            if arr[5] != 0.0
                v = vddMatrices( Int(arr[6]) , Int(arr[7]), 1.0==arr[9])
                return new( arr[1], arr[2], arr[3], arr[4], arr[5], 
                    Int(arr[6]), Int(arr[7]), Int(arr[8]), arr[9], v )
            else
                return new( arr[1], arr[2], arr[3], arr[4], arr[5], 
                    Int(arr[6]), Int(arr[7]), Int(arr[8]), arr[9], [])       
            end 
        end

    end

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

        function paramRanges( arr )
            return new( arr[1], arr[2], arr[3], arr[4], arr[5], arr[6],
                arr[7], arr[8], arr[9] )
        end


    end


    mutable struct llgParams

        tMax::Int64
        dt::Float64
        saveEvery::Float64
        tol::Float64
        damp::Float64
        temp::Float64
        nRuns::Int64
        parallel::Int64
        numCores::Int64

        function llgParams( arr )
            return new( Int(arr[1]), arr[2], arr[3], arr[4], arr[5], arr[6],
                arr[7], arr[8], arr[9] )
        end

    end

    mutable struct faParams

        sMax::Int64
        param::Float64
        tol::Float64
        temp::Float64
        nRuns::Float64
        parallel::Float64

        function faParams( arr )
            return new( arr[1], arr[2], arr[3], arr[4], arr[5], arr[6] )
        end

    end

    mutable struct icParams

        type::String
        r::Float64
        gamma::Float64
        px::Int64
        py::Int64
        
        function icParams( arr )
            return new( arr[1], arr[2], arr[3], arr[4], arr[5] )
        end
    end

    mutable struct pinningParams

        hPin::Float64
        px::Int64
        py::Int64
        
        function pinningParams( arr )
            return new( arr[1], arr[2], arr[3] )
        end

    end

    mutable struct defectParams

        t::Float64
        strength::Float64
        width::Float64
        dx::Int64
        dy::Int64

        function defectParams( arr )
            return new( arr[1], arr[2], arr[3], arr[4], arr[5] )
        end

    end

    mutable struct saveChoices

        totE::Int
        excE::Int
        zeeE::Int
        dmiE::Int
        pmaE::Int
        ddiE::Int
        magn::Int
        size::Int
        qCharge::Int
        location::Int
        spinField::Int

        function saveChoices( arr )
            return new( arr[1], arr[2], arr[3], arr[4], arr[5],
                arr[6], arr[7], arr[8], arr[9], arr[10], arr[11] )
        end

    end

    # All parameters stored in this struct. Pass this struct 
    # to all the functions. 
    mutable struct params

        mp 
        llg
        fa 
        ic
        pin 
        defect 
        save 
        range

        function params( arr )
            return new( arr[1], arr[2], arr[3], arr[4], arr[5], arr[6], arr[7], arr[8] )
        end
    end

    function getCompuationParams( inputType, paramList::Array{String,1}, defaultList)

        ans = getUserInput( String, string("\nEnter paramters ", [i for i in paramList], "? y/n ")  )
        
        # If user wants to enter parameters, get them
        if ans == "y"
            paramVals = zeros( length(paramList) )

            for i in 1:length(paramList)

                paramVals[i] = getUserInput( typeof(defaultList[i]), string( paramList[i], " = ") )

            end

            return inputType( paramVals )

        # Otherwise set default values
        else
            
            println( string("Setting default values ", [j for j in defaultList] ) )
            return inputType( defaultList )

        end      

    end

    function getUserInput( T=String, msg="" )
        
        print( "$msg" )

        if T == String
            return readline()
        else
            try
                return parse(T,readline())
            catch
                println("Sorry, I could not interpret your answer. Please try again")
                getUserInput(T,msg)
            end
        end

    end

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
        else
            println(string(" Error modifying material parameter struct. No ", num, "th value found. Quitting."))
            quit()
        end

    end



    function getTestParams()

        myMatParams = materialParams([1.0, -0.2, 0.1, 0.01, 0.0, 64, 64, 1, true] )
        myLlgParams = llgParams([10, 0.2, 5.0, 10^-6, 1.0, 0.0, 1, 0, 0])
        myFaParams = []
        myICParams = icParams( ["skyrmion", 5, pi/2 , myMatParams.nx/2, myMatParams.ny/2] )
        myPinningParams = pinningParams( [0.1, 32, 32] )
        myDefectParams = defectParams( [1.0, -0.7, 1.0, 0, 0] )
        mySaveChoices = saveChoices( [1,0,0,0,0,0,0,1,1,0,1] )
        myParamRanges = paramRanges( [ [],[], [], [], [], [], [], [], [], [] ] )

        # Now put all of the user choices into one struct 
        allParams = params( [myMatParams, myLlgParams, myFaParams, myICParams, myPinningParams, 
            myDefectParams, mySaveChoices, myParamRanges] )

        return allParams

    end

    function getUserParams()
        
        ##################################################################################################
        # Begin user prompts
        #
        # Ask user to input data
        # The user may want to run the computation for various different parameters 
        myMatParams = getCompuationParams(materialParams, ["J", "H", "A", "Dz", "Ed", "Nx", "Ny", "Nz","pbc" ], 
            [1.0, -0.2, 0.1, 0.01, 0.0, 64, 64, 1, true] )

        ans = getUserInput( String, string("\nWould you like to run over a range of parameters? Enter n to choose 'no' and run for only selected choices. y/n ")  )
        myParamRanges = paramRanges( [ [], [], [], [], [], [], [], [], [] ] )

        if ans == "y"


            ans = getUserInput( String, string("\nSelect which values you would like to vary:
                1 = Exchange interaction
                2 = Zeeman Field
                3 = DMI constant
                4 = PMA constant
                5 = DDI constant
                6 = Nx (Lattice size in x)
                7 = Ny (Lattice size in y)
                8 = Nz
                9 = pbc
            
                Enter your choice(s) separated by commas: ")  )

            varyChoices = parse.(Int64, split(ans,",") )

            for var in varyChoices
                ans = getUserInput( String, string( "\nEnter the start, step size, and end value of choice ", var, " separated by commas " ) )
                
                paramRanges = parse.( Float64, split(ans,",") )
                paramVals = [ paramRanges[1]:paramRanges[2]:paramRanges[3]; ]

                modifyMatParam!( myParamRanges, paramVals, var )

            end

        end

        ans = getUserInput( String, string( "\nComputation type? 
            1 = Landau Lifshitz, 
            2 = Field Alignment relaxation ") )

        if ans == "1"

            myLlgParams = getCompuationParams( llgParams, [ "Max iterations", "Step size", "Skip size", 
                "Tolerance", "Damping", "Temperature", "Number of runs", "Parallelize (1 or 0)", "If parallel, number of cores"], 
                [10000, 0.2, 5.0, 10^-6, 1.0, 0.0, 1, 0, 1] )
            myFaParams = []

        elseif ans == "2"

            myFaParams = getCompuationParams( faParams, ["Max iterations", "FA Param","Num rotations", 
                "Temperature", "Number of runs", "Parallelize (1 or 0)"], [10000, 0.03, 10, 0.0, 1, 0] )
            myLlgParams = []

        else
            println("Invalid input. Quitting.")
            exit()
        end

        # Initial condition
        myICParams = getCompuationParams( icParams, 
            ["Initial condition", "lambda", "gamma", "x position", "y position"], 
            ["skyrmion", 5, pi/2 , myMatParams.nx/2, myMatParams.ny/2] )


        ################################################
        # Problem-specific options start here

        # Ask if the user wants to pin the skyrmion
        ans = getUserInput( String, string("\nPin skyrmion? y/n ")  )
            
        # Get params. Default set pin field at center of the lattice
        if ans == "y"
            myPinningParams = getCompuationParams( pinningParams, ["H_pin", "X", "Y"], 
                [0.01, myMatParams.nx/2 , myMatParams.ny/2] )
        else
            myPinningParams = []
        end

        # Is there a defect in the lattice
        ans = getUserInput( String, string("\nPlace defect in the lattice? y/n ")  )

        # If user wants to enter defect parameters, get them
        if ans == "y"
            myDefectParams = getCompuationParams( defectParams, 
                ["type", "strength (relative to J)", "width", "x position", "y position"], 
                ["gaussian", -0.5, 1.0, 0, 0] )
        else
            myDefectParams = []
        end

        # Choose what to save
        ans = getUserInput( String, string("\nSelect which values you would like to save during the computation. (Note that choosing too many may increase computation time.) The options are: 
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
            
            Enter your choices in a comma-separated list: ")  )

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
