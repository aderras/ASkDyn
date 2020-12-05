# This module contains functions that call the numerical solver. Useful
# function is runComputation(), which runs the runge kutta solver and
# saves the resulting data.
module spinDynamics

    using modifyFiles, energy, topologicalCharge, LLequation,
    magnetization, rungeKutta, LLequation, normalize, noiseRotation, Dates,
    effectiveSize, initialCondition, Distributed

    export evaluateSpinDynamics, evaluateLL!, runRelaxation!, pulseNoiseStep!

    # Given a struct of parameters, compute the spin dynamics
    function evaluateSpinDynamics( p )

        s0 = buildInitial( p.ic, p.mp )

        # First run relaxation
        evaluateLL!( s0, p, true )

        # Then dynamics
        evaluateLL!( s0, p )

        println("Completed eval on worker ", myid())

    end

    # in: mat = (N,N,3) initial condition, params struct
    # out: returns nothing. Modifies mat and exports data to disk
    function evaluateLL!( mat::Array{Float64,3}, params, relaxation = false )

        maxLoop = params.llg.tMax

        ###########################################################################
        # Make arrays to store data based on what user requested
        # Can't think of a smarter way to do this....
        j,h,a,dz,ed,nx,ny,nz,pbc,vdd =
            [ getfield( params.mp, x ) for x in fieldnames( typeof(params.mp) ) ]

        # Directory where results are saved. By default, saves to
        # parentDirectory() + "/data/"
        reldir = string(dirname(pwd()),"/data/")

        if relaxation
            filesuffix = string("relaxation_T=",params.llg.temp,"_H=",h,"_A=",
                a,"_DZ=",round(dz,digits=5),"_ED=",round(ed,digits=5),"_JX=",
                params.current.jx,"_JY=",params.current.jy,"_.h5")
        else
            filesuffix = string("T=",params.llg.temp,"_H=",h,"_A=",a,"_DZ=",
                round(dz,digits=5),"_ED=",round(ed,digits=5),"_JX=",
                params.current.jx,"_JY=",params.current.jy,"_.h5")
        end

        if params.save.excE == 1.0
            excArray = zeros( maxLoop )
        end
        if params.save.zeeE == 1.0
            zeeArray = zeros( maxLoop )
        end
        if params.save.dmiE == 1.0
            dmiArray = zeros( maxLoop )
        end
        if params.save.pmaE == 1.0
            pmaArray = zeros( maxLoop )
        end
        if params.save.ddiE == 1.0
            ddiArray = zeros( maxLoop )
        end
        if params.save.magn == 1.0
            magnArray = zeros( maxLoop )
        end
        if params.save.size == 1.0
            sizeArray = zeros( maxLoop )
        end
        if params.save.location == 1.0
            locArray = zeros( maxLoop*2 )
        end

        # Always store energy and topological charge, so make arrays for those
        enArray = zeros( maxLoop )
        qArray = zeros( maxLoop )

        ##########################################################################
        # Begin computation
        # First, run relaxation on the spin lattice to find the stable skyrmion
        # configuration
        # println( "Running relaxation on the spin lattice. " )
        # @time runRelaxation!( mat, params )
        # writeDataH5(string(reldir,"initial-spin-field_",filesuffix),mat)

        # Stopping criteria: energy converges to within some tolerance
        # or topological charge becomes negative
        for i in 1:maxLoop

            enArray[i]  = calcEnergy( mat, params )
	        qArray[i]	= calcQ(mat)

            pulseNoiseStep!(mat, params, relaxation)

            if params.save.excE == 1.0
                excArray = exchangeEnergy( mat, j, pbc==1.0, params.defect)
            end
            if params.save.zeeE == 1.0
                zeeArray = zeemanEnergy( mat, h )
            end
            if params.save.dmiE == 1.0
                dmiArray = dmiEnergy( mat, a, pbc==1.0 )
            end
            if params.save.pmaE == 1.0
                pmaArray = pmaEnergy( mat, dz )
            end
            if params.save.ddiE == 1.0
                ddiArray = ddiEnergy( mat, ed, vdd )
            end
            if params.save.magn == 1.0
                magnArray[i] = calcM( mat )
            end
            if params.save.size == 1.0
                sizeArray[i] = calcEffectiveSize( mat )
            end
            if params.save.location == 1.0
                @views maxPos = argmax( mat[3,:,:] )
                locArray[(i*2-1)] = maxPos[1]
                locArray[i*2] = maxPos[2]
            end
            if params.save.spinField == 1.0
                writeDataH5(string(reldir,"spin-field-i=",i,"_",filesuffix),mat)
            end

            enArray[i] = calcEnergy( mat, params )

            qArray[i] = calcQ( mat )

            # println("i = ", i, ", Q = ", qArray[i], ", E = ", enArray[i] )

            # Check for collapse. Currently only checking that topological charge decreases
            if abs(qArray[i]) < 0.1
                break
            end

            if (relaxation && i > 40 &&
                    abs(enArray[i]-enArray[i-1]) < params.llg.tol)
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
            writeDataH5( string(reldir,"exchange-energy_",filesuffix),filter(x->x!=0.,excArray) )
        end
        if params.save.zeeE == 1.0
            writeDataH5(string(reldir,"zeeman-energy_",filesuffix),filter(x->x!=0.,zeeArray))
        end
        if params.save.dmiE == 1.0
            writeDataH5(string(reldir,"dmi-energy_",filesuffix),filter(x->x!=0.,dmiArray))
        end
        if params.save.pmaE == 1.0
            writeDataH5(string(reldir,"pma-energy_",filesuffix),filter(x->x!=0.,pmaArray))
        end
        if params.save.ddiE == 1.0
            writeDataH5(string(reldir,"ddi-energy_",filesuffix),filter(x->x!=0.,ddiArray))
        end
        if params.save.magn == 1.0
            writeDataH5(string(reldir,"magnetization_",filesuffix),filter(x->x!=0.,magnArray))
        end
        if params.save.size == 1.0
            writeDataH5(string(reldir,"effective-size_",filesuffix),filter(x->x!=0.,sizeArray))
        end
        if params.save.location == 1.0
            writeDataH5(string(reldir,"max-spin-location_",filesuffix),filter(x->x!=0.,locArray))
        end
        if params.save.totE == 1.0
            writeDataH5(string(reldir,"total-energy_",filesuffix),filter(x->x!=0.,enArray))
        end
        if params.save.qCharge == 1.0
            writeDataH5(string(reldir,"topological-charge_",filesuffix),filter(x->x!=0.,qArray))
        end

        # Alwasys save the final spin field
        writeDataH5(string(reldir,"final-spin-field_",filesuffix),mat)

    end

    # runs a single pulse noise step
    function pulseNoiseStep!( s::Array{Float64,3}, params, relaxation = false )

        # This is the pulse-noise algorithm
        if params.llg.temp == 0.0
            rk4!(s, RHS!, params, relaxation)
            rk4!(s, RHS!, params, relaxation)
            normalizeSpins!(s)
        else
            rk4!(s, RHS!, params, relaxation)
            rotateSpins!(s, params.llg)
            rk4!(s, RHS!, params, relaxation)
            normalizeSpins!(s)
        end

    end

end
