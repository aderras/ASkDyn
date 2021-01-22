# This module contains functions that call the numerical solver. Useful
# function is runComputation(), which runs the runge kutta solver and
# saves the resulting data.
module SpinDynamics

    using HDF5, Energy, TopologicalCharge, LLequation, Magnetization, RungeKutta,
    LLequation, Normalize, NoiseRotation, Dates, EffectiveSize,
    InitialCondition, Distributed, FieldAlignment

    export launchcomputation, rundynamics!, runRelaxation!, compute_ll!

    # Given a struct of parameters, compute the spin dynamics
    function launchcomputation(p)

        # First run relaxation if the user requested. If they didn't, check
        # the folder "relaxation-data/" for data there. If there's no data,
        # print a warning and continue
        if p.llg.relax == 1

            println("Computing relaxation with LLG and high damping. ")
            s0 = buildinitial(p.ic, p.mp)
            rundynamics!(s0, p, true)

        elseif p.llg.relax == 2

            println("Computing relaxation with field-alignment algorithm.")
            s0 = buildinitial(p.ic, p.mp)
            runfieldalignment!(s0, p)

        else

            # Check for relaxation data
            relaxDir = string(dirname(pwd()),"/relaxation-data/")

            j,h,a,dz,ed = [getfield(p.mp, x)
                for x in fieldnames(typeof(p.mp))[1:5]]

            filesuffix = string("relaxation_T=",p.llg.temp,"_H=",h,"_A=",
                a,"_DZ=",round(dz,digits=5),"_ED=",round(ed,digits=5),"_ICX=",
                p.ic.px,"_.h5")

            filename = string(relaxDir,"final-spin-field_", filesuffix)

            if isfile(filename)

                println("Importing relaxation file ", filename)
                s0 = h5read(filename,"Dataset1")

            else

                println("No relaxation data with name ", filename, " found. ",
                    "Continuing evaluation without relaxation.")
                s0 = buildinitial(p.ic, p.mp)

            end

        end

        # Then dynamics
        rundynamics!(s0, p)

        println("Completed eval on worker ", myid())

    end

    # in: mat = (N,N,3) initial condition, params struct
    # out: returns nothing. Modifies mat and exports data to disk
    function rundynamics!(mat::Array{Float64,3}, params, relaxation = false)

        maxLoop = params.llg.tMax

        ###########################################################################
        # Make arrays to store data based on what user requested
        # Can't think of a smarter way to do this....
        j,h,a,dz,ed,nx,ny,nz,pbc,vdd =
            [getfield(params.mp, x) for x in fieldnames(typeof(params.mp))]

        # Directory where results are saved. By default, saves to
        # parentDirectory() + "/data/"
        reldir = string(pwd(),"/data/")

        filesuffix = string("T=",params.llg.temp,"_H=",h,"_A=",
            a,"_DZ=",round(dz,digits=5),"_ED=",round(ed,digits=5),"_ICX=",
            params.ic.px,"_.h5")

        if relaxation
            filesuffix = string("relaxation_", filesuffix)
        end

        if params.save.excE == 1.0
            excArray = zeros(maxLoop)
        end
        if params.save.zeeE == 1.0
            zeeArray = zeros(maxLoop)
        end
        if params.save.dmiE == 1.0
            dmiArray = zeros(maxLoop)
        end
        if params.save.pmaE == 1.0
            pmaArray = zeros(maxLoop)
        end
        if params.save.ddiE == 1.0
            ddiArray = zeros(maxLoop)
        end
        if params.save.magn == 1.0
            magnArray = zeros(maxLoop)
        end
        if params.save.size == 1.0
            sizeArray = zeros(maxLoop)
        end
        if params.save.location == 1.0
            locArray = zeros(maxLoop*2)
        end

        # Always store energy and topological charge, so make arrays for those
        enArray = zeros(maxLoop)
        qArray = zeros(maxLoop)

        ##########################################################################
        # Begin computation
        # First, run relaxation on the spin lattice to find the stable skyrmion
        # configuration
        # println("Running relaxation on the spin lattice. ")
        # @time runRelaxation!(mat, params)
        # h5write(string(reldir,"initial-spin-field_", filesuffix), "Dataset1", mat)

        # Stopping criteria: energy converges to within some tolerance
        # or topological charge becomes negative
        for i in 1:maxLoop

            compute_ll!(mat, params, relaxation)

            if params.save.excE == 1.0
                excArray[i] = exchange_energy(mat, j, pbc==1.0, params.defect)
            end
            if params.save.zeeE == 1.0
                zeeArray[i] = zeeman_energy(mat, h)
            end
            if params.save.dmiE == 1.0
                dmiArray[i] = dmi_energy(mat, a, pbc==1.0)
            end
            if params.save.pmaE == 1.0
                pmaArray[i] = pma_energy(mat, dz)
            end
            if params.save.ddiE == 1.0 
                ddiArray[i] = ddi_energy(mat, ed, pbc, vdd)
            end
            if params.save.magn == 1.0
                magnArray[i] = magnetization(mat)
            end
            if params.save.size == 1.0
                sizeArray[i] = effectivesize(mat)
            end
            if params.save.location == 1.0
                @views maxPos = argmax(mat[3,:,:])
                locArray[(i*2-1)] = maxPos[1]
                locArray[i*2] = maxPos[2]
            end
            if params.save.spinField == 1.0 && i%10==0
                h5write(string(reldir,"spin-field-i=",i,"_", filesuffix), "Dataset1",mat)
            end

            enArray[i] = energy(mat, params)
            qArray[i] = calcQ(mat)

            # Print if debugging
            #println("i = ", i, ", Q = ", qArray[i], ", E = ", enArray[i])

            # Check for collapse
            if abs(qArray[i]) < 0.1
                break
            end

            if relaxation && i > 40 &&
                    abs(enArray[i]-enArray[i-1]) < params.llg.tol
                break
            end
        end

        # Save with a timestamp and random number generator if you are running the same computation
        # multiple times. (Not using right now because still testing and want to overwrite saved
        # files instead of filling up the folder with distinct data.)
        timestampString = string(Dates.format(Dates.now(),"HH_MM_SS"),
                trunc(Int,rand()*10000))

        # Define file suffix to distinguish the files The following contains defect info
        # filesuffix = string("T=",params.llg.temp,"_H=",h,"_A=",a,"_DZ=",dz,"_ED=",round(ed,digits=5),"_DEFECT=",
        #         params.defect.t,"_D-X=",params.defect.dx,"_D-Y=",params.defect.dy,"_D-STRENGTH=",
        #         params.defect.strength,"_D-WIDTH=",params,defect,width,"_",timestampString,"_.h5")

        if params.save.excE == 1.0
            h5write(string(reldir,"exchange-energy_", filesuffix), "Dataset1",
            filter(x->x!=0.,excArray))
        end
        if params.save.zeeE == 1.0
            h5write(string(reldir,"zeeman-energy_", filesuffix), "Dataset1",
            filter(x->x!=0.,zeeArray))
        end
        if params.save.dmiE == 1.0
            h5write(string(reldir,"dmi-energy_", filesuffix), "Dataset1",
            filter(x->x!=0.,dmiArray))
        end
        if params.save.pmaE == 1.0
            h5write(string(reldir,"pma-energy_", filesuffix), "Dataset1",
            filter(x->x!=0.,pmaArray))
        end
        if params.save.ddiE == 1.0
            h5write(string(reldir,"ddi-energy_", filesuffix), "Dataset1",
            filter(x->x!=0.,ddiArray))
        end
        if params.save.magn == 1.0
            h5write(string(reldir,"magnetization_", filesuffix), "Dataset1",
            filter(x->x!=0.,magnArray))
        end
        if params.save.size == 1.0
            h5write(string(reldir,"effective-size_", filesuffix), "Dataset1",
            filter(x->x!=0.,sizeArray))
        end
        if params.save.location == 1.0
            h5write(string(reldir,"max-spin-location_", filesuffix), "Dataset1",
            filter(x->x!=0.,locArray))
        end
        if params.save.totE == 1.0
            h5write(string(reldir,"total-energy_", filesuffix), "Dataset1",
            filter(x->x!=0.,enArray))
        end
        if params.save.qCharge == 1.0
            h5write(string(reldir,"topological-charge_", filesuffix), "Dataset1",
            filter(x->x!=0.,qArray))
        end

        # Always save the final spin field
        h5write(string(reldir,"final-spin-field_", filesuffix), "Dataset1", mat)

    end

    # Calculates Landau Lifshitz equation
    function compute_ll!(s::Array{Float64,3}, params, relaxation = false)

        # This is the pulse-noise algorithm
        if params.llg.temp == 0.0
            rk4!(s, RHS!, params, relaxation)
            rk4!(s, RHS!, params, relaxation)
            normalizelattice!(s)
        else
            rk4!(s, RHS!, params, relaxation)
            noisyrotate!(s, params.llg)
            rk4!(s, RHS!, params, relaxation)
            normalizelattice!(s)
        end
    end
end
