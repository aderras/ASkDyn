#=

    This module contains functions to be used to calculate the minimum
    energy configuration using the field alignment relaxation. Details
    on the numerical method can be found at
      https://arxiv.org/abs/1311.7211

    IMPROVEMENTS
      + randomdirection!() is clunky because it modifies 1 element arrays. I
      wasn't able to modify regular variables in place. You probably don't need
      to though.

=#
module FieldAlignment

    using EffectiveField, Energy, TopologicalCharge, EffectiveSize,
    LinearAlgebra, Normalize, HDF5
    export runfieldalignment!, modifyspin!

    # This function runs modifyspin! until the magnetization settles
    # to some value within a specified tolerance. The input matrix,
    # mat is modified during the computation.
    #
    # in: mat = the spin matrix, params = struct of all the compuataion
    # parameters
    #
    # out: nothing
    function runfieldalignment!(mat::Array{Float64,3}, params)

        maxLoop, alpha, nrot, tol =
            [getfield(params.fa, x) for x in fieldnames(typeof(params.fa))]

        j,h,a,dz,ed,nx,ny,nz,pbc,v =
            [getfield(params.mp, x) for x in fieldnames(typeof(params.mp))]

        reldir = string(pwd(),"/data/")

        filesuffix = string("fa_H=",h,"_A=",a,"_DZ=",round(dz,digits=5),"_ED=",
            round(ed,digits=5),"_ICX=",params.ic.px,"_.h5")

        p, m, n = size(mat)

        # Initialize parameter storage
        if params.save.totE == 1.0
            enArray = zeros(maxLoop)
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
        end
        if params.save.location == 1.0
            locArray = zeros(maxLoop*2)
        end

        sizeArray = zeros(maxLoop)
        qArray = zeros(maxLoop)

        ddiField = zeros(p, m, n)

        # Run the random spin relaxation
        for i in 1:maxLoop

            if ed != 0.0
                ddiField = ddifield(mat, ed, pbc==1.0, v)
            end

            for i in 1:nrot modifyspin!(mat, ddiField, params) end
            normalizelattice!(mat)

            if params.save.totE == 1.0
                enArray[i] = energy(mat, params)
            end
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
                ddiArray[i] = ddi_energy(mat, ed, pbc, v)
            end
            if params.save.magn == 1.0
                magnArray[i] = magnetization(mat)
            end
            if params.save.location == 1.0
                @views maxPos = argmax(mat[3,:,:])
                locArray[(i*2-1)] = maxPos[1]
                locArray[i*2] = maxPos[2]
            end

            sizeArray[i] = effectivesize(mat)
            qArray[i] = calcQ(mat)

            # println("i = ", i, ", Q = ", qArray[i], ", R = ", sizeArray[i])

            # If the skyrmion collapses or the effective field converges,
            # break the loop
            if abs(qArray[i]) < 0.1 || (i > 10 &&
                abs(sizeArray[i-1] - size[i]) < tol)
                break
            end

        end


        if params.save.totE == 1.0
            h5write(string(reldir,"total-energy_", filesuffix), "Dataset1",
            filter(xx->xx!=0.0,enArray))
        end
        if params.save.excE == 1.0
            h5write(string(reldir,"exchange-energy_", filesuffix), "Dataset1",
            filter(xx->xx!=0.,excArray))
        end
        if params.save.zeeE == 1.0
            h5write(string(reldir,"zeeman-energy_", filesuffix), "Dataset1",
            filter(xx->xx!=0.,zeeArray))
        end
        if params.save.dmiE == 1.0
            h5write(string(reldir,"dmi-energy_", filesuffix), "Dataset1",
            filter(xx->xx!=0.,dmiArray))
        end
        if params.save.pmaE == 1.0
            h5write(string(reldir,"pma-energy_", filesuffix), "Dataset1",
            filter(xx->xx!=0.,pmaArray))
        end
        if params.save.ddiE == 1.0
            h5write(string(reldir,"ddi-energy_", filesuffix), "Dataset1",
            filter(xx->xx!=0.,ddiArray))
        end
        if params.save.magn == 1.0
            h5write(string(reldir,"magnetization_", filesuffix), "Dataset1",
            filter(xx->xx!=0.,magnArray))
        end
        if params.save.location == 1.0
            h5write(string(reldir,"max-spin-location_", filesuffix), "Dataset1",
            filter(xx->xx!=0.,locArray))
        end
        if params.save.qCharge == 1.0
            h5write(string(reldir,"topological-charge_", filesuffix), "Dataset1",
            filter(xx->xx!=0.,qArray))
        end

        # Always save the final spin field
        h5write(string(reldir,"effective-size_", filesuffix), "Dataset1",
            filter(xx->xx!=0.,sizeArray))
        h5write(string(reldir,"final-spin-field_", filesuffix), "Dataset1",  mat)

    end


    # The function changes the input matrix itself. (This means allocating less
    # space for the computation.) Running this function once means running through
    # the entire spin field once and changing the spins according to some randomness
    # parameter, α.

    # inputs: mat = spin matrix, α = randomness parameter, N = number of spins
    # outputs: nothing
    function modifyspin!(mat::Array{Float64,3}, ddiField::Array{Float64,3},
            params)

        p,m,n = size(mat)

        r = 0.0
        nhcurrent = 0.0
        pref =  0.0

        startx = [1]
        starty = [1]
        endx = [m]
        endy = [n]
        dx = [1]
        dy = [1]

        # initialize these as (3,1) arrays
        # scurrent = a copy of the spin element we're looking at
        # hcurrent = the effective field at scurrent
        # sreplace = the spin to change into
        scurrent = Array{Float64}(undef,3)
        hcurrent = Array{Float64}(undef,3)
        sreplace = Array{Float64}(undef,3)

        randomdirection!(startx,starty,endx,endy,dx,dy)

        # If including DDI, calculate the effective field first. This
        # is because this particular computation takes a long time and is
        # prohibitively slow if computed N x N times.
        for i in startx[1]:dx[1]:endx[1]

            for j in starty[1]:dy[1]:endy[1]

                # scurrent is a spin x,y,x array at a certain i,j
                for k in 1:3 scurrent[k] = mat[k,i,j] end
                # println("nx = ", i, ", ny = ", j,", s0 = ", scurrent)

                # hcurrent has to start as zero
                hcurrent[1] = hcurrent[2] = hcurrent[3] = 0.0

                # this function changes hcurrent so that its new values
                # are the effective field at i,j
                effectivefieldelem!(hcurrent,mat, i, j, params)
                for k in 1:3 hcurrent[k] = hcurrent[k] +  ddiField[k,i,j] end

                if params.mp.ed != 0.0
                    for k in 1:3 hcurrent[k] += ddiField[k, i, j] end
                end
                # println("Heff = ", hcurrent)

                # for some reason norm allocates a lot of space, so
                # manually computed the norm, nhcurrent
                nhcurrent = sqrt(dot(hcurrent, hcurrent))

                r = rand()

                # If the random number is less than alpha, then replace
                # the spin at that point with the Normalized h field
                if r < params.fa.param
                    for k in 1:3 sreplace[k] = hcurrent[k]/nhcurrent end
                    # Otherwise, replace the spin with the following
                else
                    # calculate prefactor
                    pref = 2(dot(scurrent,hcurrent))/nhcurrent^2
                    for k in 1:3 sreplace[k] = pref*hcurrent[k]-scurrent[k] end
                end

                # store results
                for k in 1:3 mat[k,i,j] = sreplace[k] end

            end
        end

        nothing

    end

    # Randomly decide whether to traverse the array forward or backward, toward
    # the top or toward the bottom.
    #
    # in: start and end values of the randomness
    #
    # out: begining and end range in both x and y in a randomnly chosen
    # direction. Returns (4,1) array
    function randomdirection!(startx::Array{Int,1},starty::Array{Int,1},
        endx::Array{Int,1},endy::Array{Int,1},dx::Array{Int,1},
        dy::Array{Int,1})

        if rand(0:1) == 0
            # do nothing to x range
        else
            startx[1] = endx[1]
            endx[1] = 1
            dx[1] = -1
        end

        if rand(0:1) == 0
            # do nothing to y range
        else
            starty[1] = endy[1]
            endy[1] = 1
            dy[1] = -1
        end

    end

end
