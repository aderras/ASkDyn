# This module contains functions that call the numerical solver. Useful
# function is runComputation(), which runs the runge kutta solver and 
# saves the resulting data. 
module spinDynamics

    using modifyFiles, energy, topologicalCharge, LLequation,
    magnetization, rungeKutta, LLequation, normalize, noiseRotation, Dates,
    effectiveSize
    
    export evaluateLL!, runRelaxation!, pulseNoiseStep!

    # The following function evaluates the LLG equation. In the comments you'll find
    # two options: pulseNoise or rk4. For zero temperature relaxation, it doesn't matter
    # which you use. For nonzero temperature, must use pulse noise. 
    #
    # inputs: mat = (N,N,3) initial condition, params 
    function evaluateLL!( mat::Array{Float64,3}, params )


        maxLoop = params.llg.tMax

        #########################################################
        # Make arrays to store data based on what user requested
        # Can't think of a smarter way to do this....
        j,h,a,dz,ed,nx,ny,nz,pbc,vdd = 
            [ getfield( params.mp, x ) for x in fieldnames( typeof(params.mp) ) ]

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

        # Have to store energy and topological charge. 
        enArray = zeros( maxLoop )
        qArray = zeros( maxLoop )

        ###############################################################
        # Begin computation
        # Stopping criteria: energy converges to within some tolerance
        # or topological charge becomes negative
        for i in 1:maxLoop

            enArray[i]  = calcEnergy( mat, params.mp )
	        qArray[i]	= calcQ(mat)
	  
            pulseNoiseStep!(mat, params)

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

            # If the user asked to save all the energies 
            # independently, then just sum those numbers.
            # Otherwise calculate the total energy
            if params.save.excE == 1.0 &&
                params.save.zeeE == 1.0 &&
                params.save.dmiE == 1.0 &&
                params.save.pmaE == 1.0 &&
                params.save.ddiE == 1.0

                enArray[i] = excArray[i] + zeeArray[i] +
                    dmiArray[i] + pmaArray[i] + ddiArray[i]

            else 
                enArray[i] = calcEnergy( mat, params.mp, params.defect )
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

            qArray[i] = calcQ( mat )

            # Check for collapse
            if abs(qArray[i]) < 0.1
                break 
            end

        end

        # For loop ended, so save data. 
        reldir = "/data/"
        timestampString = string(Dates.format(Dates.now(),"HH_MM_SS"),
                trunc(Int,rand()*10000))

        # Define file suffix to distinguish the files 
        # filesuffix = string("T=",params.llg.temp,"_H=",h,"_A=",a,"_DZ=",dz,"_ED=",round(ed,digits=5),"_DEFECT=",
        #         params.defect.t,"_D-X=",params.defect.dx,"_D-Y=",params.defect.dy,"_D-STRENGTH=",
        #         params.defect.strength,"_D-WIDTH=",params,defect,width,"_",timestampString,"_.h5")

        filesuffix = string("T=",params.llg.temp,"_H=",h,"_A=",a,"_DZ=",dz,"_ED=",round(ed,digits=5),"_.h5")

        if params.save.excE == 1.0
            writeDataH5(string(reldir,"exchange-energy_",filesuffix),excArray)
        end
        if params.save.zeeE == 1.0
            writeDataH5(string(reldir,"zeeman-energy_",filesuffix),zeeArray)
        end
        if params.save.dmiE == 1.0
            writeDataH5(string(reldir,"dmi-energy_",filesuffix),dmiArray)
        end
        if params.save.pmaE == 1.0
            writeDataH5(string(reldir,"pma-energy_",filesuffix),pmaArray)
        end
        if params.save.ddiE == 1.0
            writeDataH5(string(reldir,"ddi-energy_",filesuffix),ddiArray)
        end
        if params.save.magn == 1.0
            writeDataH5(string(reldir,"magnetization_",filesuffix),magnArray)
        end
        if params.save.size == 1.0
            writeDataH5(string(reldir,"effective-size_",filesuffix),sizeArray)
        end
        if params.save.location == 1.0
            writeDataH5(string(reldir,"max-spin-location_",filesuffix),locArray)
        end
        if params.save.totE == 1.0
            writeDataH5(string(reldir,"total-energy_",filesuffix),enArray)
        end
        if params.save.spinField == 1.0
            writeDataH5(string(reldir,"final-spin-field_",filesuffix),mat)
        end
        if params.save.qCharge == 1.0
            writeDataH5(string(reldir,"topological-charge_",filesuffix),qArray)
        end


    end

    # runs a single pulse noise step
    function pulseNoiseStep!( s::Array{Float64,3}, params )

        # This is the pulse-noise algorithm
        rk4!(s, RHS!, params)
        rotateSpins!(s, params.llg)
        rk4!(s, RHS!, params)
        normalizeSpins!(s)

    end

    function runRelaxation!( mat, params )

        maxLoop = params.llg.tMax
        
        prevEnergy = 0.0
        currEnergy = 0.0
        diff = 10.0
        
        @time for i in 1:maxLoop

            rk4!( mat, RHS!, params )
            
            currEnergy = calcEnergy( mat, params.mp );
            diff = currEnergy - prevEnergy

            if abs(diff) < 0.00004
                break
            else
                prevEnergy = currEnergy
            end

        end

    end

end