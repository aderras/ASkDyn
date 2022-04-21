# This module contains functions that call the numerical solver.
module SpinDynamics

    using HDF5, Energy, TopologicalCharge, LLequation, Magnetization,
    RungeKutta, LLequation, Normalize, NoiseRotation, Dates, EffectiveSize,
    InitialCondition, Distributed, FieldAlignment, UserInputs, Chirality,
    Helpers

    export launchcomputation, rundynamics!, runRelaxation!, computeLL!

    # Given a struct of parameters, compute the spin dynamics
    function launchcomputation(p)

        # First run relaxation if the user requested. If they didn't, check
        # the relaxation folder for data there. If there's no data,
        # print a warning and continue
        if p.cp.relax == 1
            println("Computing relaxation with LLG and high damping.")
            s0 = buildinitial(p.ic, p.mp)
            rundynamics!(s0, p, true)
        elseif p.cp.relax == 2
            println("Computing relaxation with field-alignment algorithm.")
            s0 = buildinitial(p.ic, p.mp)
            runfieldalignment!(s0, p)
        else
            # Check for starting data
            inputDir = UserInputs.Paths.initialConds
            filename = string(inputDir, UserInputs.Filenames.inputName(p))

            if isfile(filename)
                println("Importing file ", filename)
                s0 = h5read(filename,"Dataset1")
                # If importing initial conditions from MM, have to reverse dims
                if size(s0)[1]!=3 s0=permutedims(s0,[3,2,1]) end
                println("Dimensions of imported data: ", size(s0))
            else
                println("No relaxation data with name ", filename, " found. ",
                    "Continuing evaluation without relaxation.")
                s0 = buildinitial(p.ic, p.mp)
            end
        end
        rundynamics!(s0, p)
        println("Completed eval on worker ", myid())
    end

    # in: mat = (N,N,3) initial condition, params struct, optional relaxation
    # boolean determines whether to set damping to 1.0 or keep user input.
    # out: returns nothing. Modifies mat and exports data to disk
    function rundynamics!(mat::Array{Float64,3}, params, relaxation=false)

        maxLoop = params.cp.maxSteps

        j,h,a,dz,ed,nx,ny,nz,pbc,vdd =
            [getfield(params.mp, x) for x in fieldnames(typeof(params.mp))]
        α = params.cp.damp

        # Directory where results are saved.
        reldir = UserInputs.Paths.output

        # Suffix used to distinguish files. Append "relaxation" if applicable
        filesuffix = UserInputs.Filenames.outputSuffix(params)

        timestamp = Dates.format(Dates.DateTime(Dates.now()), "_yyyy-m-d_HMSsss")

        if relaxation filesuffix = string("relaxation", filesuffix) end

        # If temperature is nonzero, divide tFree by 2 so that we apply noise
        # rotation halfway through (see computLL!() below)
        if params.cp.temp !=0.0 && relaxation==false
            params.cp.nn%2==0 || error("compParams.nn must be even.")
            println("Halving the number of steps for nonzero temperature.")
            params.cp.nn = round(Int64,params.cp.nn/2)
        end

        ########################################################################
        # Initialize storage

        # allArrays is an arrays of arrays. Each element is a 1D array storing
        # parameter values throughout the computation.
        allArrays = [[0.0] for i in 1:length(fieldnames(typeof(params.save)))]

        # Initialize arrays for all parameters user asked to save. All others
        # remain [0.0].
        count = 1
        for x in fieldnames(typeof(params.save))
            if getfield(params.save, x)==1 allArrays[count] = zeros(maxLoop) end
            count = count+1
        end

        # Always measure the energy. Also always measure Q if this is a skyrmion
        # Initialize arrays for both.
        if params.save.totE == 0 allArrays[1] = zeros(maxLoop) end
        if params.ic.type == "skyrmion" allArrays[9] = zeros(maxLoop) end

        p,m,n = size(mat)

        # Initialize matrices for the LLG equation calculation
        Heff = zeros(p,m,n)
        SDotH = zeros(m,n)

        # Initialize arrays for the Runge Kutta solver
        tmp = zeros(p,m,n)
        K1 = zeros(p,m,n)
        K2 = zeros(p,m,n)
        K3 = zeros(p,m,n)
        K4 = zeros(p,m,n)
        if params.cp.solver==1 # If user requested 6th order solver
            K5 = zeros(p,m,n)
            K6 = zeros(p,m,n)
            Kvec = [K1, K2, K3, K4, K5, K6, tmp]
        else
            Kvec = [K1, K2, K3, K4, tmp]
        end

        stmp1 = zeros(p)
        stmp2 = zeros(p)

        # Create arrays of parameter values to be passed into the LL compuataion
        mpValues = [j, h, a, ed, dz, vdd, pbc]
        cpValues = [params.cp.dt, params.cp.nn]

        ########################################################################
        # Begin computation

        # Stopping criteria: energy converges to within some tolerance
        # or topological charge becomes negative
        en = 10.0
        enPrev = 0.0

        for i in 1:maxLoop

            computeLL!(mat, mpValues, cpValues, params, [Heff, SDotH],
                Kvec, [stmp1, stmp2], relaxation, α)

            # According to some julia discussion boards, it may be useful to
            # call the garbage collector when you modify elements of an array
            # several times. This slows down the solver though.
            # GC.gc()
            # GC.gc()
            # GC.gc()
            # GC.gc()

            en = energy(mat, mpValues)
            allArrays[1][i] = en
            if params.save.excE == 1.0
                allArrays[2][i] = exchange_energy(mat, j, pbc, params.defect)
            end
            if params.save.zeeE == 1.0
                allArrays[3][i] = zeeman_energy(mat, h)
            end
            if params.save.dmiE == 1.0
                allArrays[4][i] = dmi_energy(mat, a, pbc)
            end
            if params.save.pmaE == 1.0
                allArrays[5][i] = pma_energy(mat, dz)
            end
            if params.save.ddiE == 1.0
                allArrays[6][i] = ddi_energy(mat, ed, pbc, vdd)
            end
            if params.save.magn == 1.0
                allArrays[7][i] = magnetization(mat)
            end
            if params.save.size == 1.0
                allArrays[8][i] = effectivesize(mat)
            end
            if params.ic.type == "skyrmion" allArrays[9][i] = calcQ(mat) end

            # THE FOLLOWING LOCATION SAVER IS INCORRECT
            if params.save.location == 1.0
                # @views maxPos = argmax(mat[3,:,:])
                # locArray[(i*2-1)] = maxPos[1]
                # locArray[i*2] = maxPos[2]
                allArrays[10][i] = 0.0
            end
            if params.save.spinField == 1.0
                if i%1 ==0
                    h5overwrite(string(reldir,"S_",
                        round(i*params.cp.dt*params.cp.nn, digits=2),
                        filesuffix),mat)
                end
            end
            if params.save.chir == 1.0
                allArrays[12][i] = Chirality.computeGamma(stmp1, mat,
                    params.cp.rChirality)
            end

            # # Print if debugging
            # println("i = ", i, ", E = ", allArrays[1][i], ", |Delta(E)| = ",
            #     abs(en-enPrev), ", Q = ", allArrays[9][i])

            # Check for collapse if this is a skyrmion
            if params.ic.type == "skyrmion"
                if abs(allArrays[9][i]) < 0.1 break end
            end
            # Check for convergence with energy if this is relaxation
            if relaxation && i > 10 && abs(en-enPrev) < params.cp.tol break end
            # If the energy is NaN, something went wrong
            if isnan(en) break end

            enPrev = en
        end

        # Save with a timestamp and random number generator if you are running
        # the same computation multiple times.
        count = 1
        for x in fieldnames(typeof(params.save))
            if getfield(params.save, x)==1
                h5overwrite(string(reldir,
                    UserInputs.Filenames.allFilenames[count], filesuffix),
                    filter(x->x!=0.0,allArrays[count]))
            end
            count = count+1
        end

        # Always save the final spin field
        h5overwrite(string(reldir,UserInputs.Filenames.finalField,filesuffix)
            ,mat)

    end

    function runRK!(s::Array{Float64,3}, mpValues::Array{Any,1},
        cpValues::Array{Float64,1}, params, fArgs, rkMatrices, relaxation=false,
        α=0.0)
        if params.cp.solver == 0
            rk4!(s, 0.0, RHS!, fArgs, mpValues, cpValues, params, rkMatrices,
                relaxation, α)
        elseif params.cp.solver == 1
            rk5!(s, 0.0, RHS!, fArgs, mpValues, cpValues, params, rkMatrices,
                relaxation, α)
        end
    end

    # Calculates Landau Lifshitz equation
    function computeLL!(s::Array{Float64,3}, mpValues::Array{Any,1},
        cpValues::Array{Float64,1}, params, fArgs,
        rkMatrices::Vector{Array{Float64,3}},
        tmps::Vector{Vector{Float64}}, relaxation=false,
        α=0.0)

        # This is the pulse-noise algorithm
        if params.cp.temp == 0.0
            runRK!(s, mpValues, cpValues, params, fArgs, rkMatrices, relaxation,
                    α)
            normalizelattice!(s)
        else
            runRK!(s, mpValues, cpValues, params, fArgs, rkMatrices, relaxation,
                    α)
            noisyrotate!(s,tmps, cpValues, α, params.cp.temp)
            runRK!(s, mpValues, cpValues, params, fArgs, rkMatrices, relaxation,
                    α)
            normalizelattice!(s)
        end
    end
end
