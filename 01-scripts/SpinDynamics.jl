# This module contains functions that call the numerical solver.
module SpinDynamics

    using HDF5, Energy, TopologicalCharge, LLequation, Magnetization, RungeKutta,
    LLequation, Normalize, NoiseRotation, Dates, EffectiveSize,
    InitialCondition, Distributed, FieldAlignment, UserInputs, Chirality

    using BenchmarkTools

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

        ########################################################################

        j,h,a,dz,ed,nx,ny,nz,pbc,vdd =
            [getfield(params.mp, x) for x in fieldnames(typeof(params.mp))]

        # Directory where results are saved.
        reldir = UserInputs.Paths.output

        # Suffix used to distinguish files. Append "relaxation" if applicable
        filesuffix = UserInputs.Filenames.outputSuffix(params)

        if relaxation filesuffix = string("_RELAXATION", filesuffix) end

        # If temperature is nonzero, divide tFree by 2 so that we apply noise
        # rotation halfway through (see computLL!() below)
        if params.cp.temp !=0.0 && relaxation==false
            params.cp.nSteps = round(Int64,p.cp.nSteps/2)
        end

        allArrays = [[0.0] for i in 1:length(fieldnames(typeof(params.save)))]

        count = 1
        for x in fieldnames(typeof(params.save))
            if getfield(params.save, x)==1 allArrays[count] = zeros(maxLoop) end
            count = count+1
        end
        # Always measure the energy. Also always measure Q if this is a skyrmion
        allArrays[1] = zeros(maxLoop)
        if params.ic.type == "skyrmion" allArrays[9] = zeros(maxLoop) end


        ########################################################################
        # Begin computation

        # Stopping criteria: energy converges to within some tolerance
        # or topological charge becomes negative
        en = 10.0
        enPrev = 0.0

        p,m,n = size(mat)

        Heff = zeros(p,m,n)
        SDotH = zeros(m,n)

        # Initialize arrays for the RK steps
        K1 = zeros(p,m,n)
        K2 = zeros(p,m,n)
        K3 = zeros(p,m,n)
        K4 = zeros(p,m,n)
        tmp = zeros(p,m,n);

        stmp = zeros(p)

        for i in 1:maxLoop

            computeLL!(mat, params, [Heff, SDotH], [K1, K2, K3, K4, tmp],
                relaxation)

            en = energy(mat, params)
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
                        "_",filesuffix),mat)
                end
            end
            if params.save.chir == 1.0
                allArrays[12][i] = Chirality.computeGamma(stmp, mat,
                    params.cp.rChirality)
            end

            # Print if debugging
            println("i = ", i, ", E = ", allArrays[1][i], ", |Delta(E)| = ",
                abs(en-enPrev))
            # println("s[10,10] = ", mat[1,10,10],mat[3,10,10],mat[2,10,10])

            # Check for collapse if this is a skyrmion
            if params.ic.type == "skyrmion"
                if abs(allArrays[9][i]) < 0.1 break end
            end
            # Check for convergence with energy if this is relaxation
            if relaxation && i > 10 && abs(en-enPrev) < params.cp.tol break end
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

    function h5overwrite(name, array)
        fid = h5open(name, "w")
        write(fid, "Dataset1", array)
        close(fid)
    end

    # Calculates Landau Lifshitz equation
    function computeLL!(s::Array{Float64,3}, params, fArgs, rkMatrices,
        relaxation=false)

        # This is the pulse-noise algorithm
        if params.cp.temp == 0.0
            rk4!(0.0, s, RHS!, params, fArgs, rkMatrices, relaxation)
            normalizelattice!(s)
        else
            rk4!(0.0, s, RHS!, params, fArgs, rkMatrices, relaxation)
            noisyrotate!(s, params.cp)
            rk4!(0.0, s, RHS!, params,  fArgs, rkMatrices, relaxation)
            normalizelattice!(s)
        end
    end
end
