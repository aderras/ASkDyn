# This module contains functions to be used to calculate the minimum
# energy configuration using the field alignment relaxation. Details 
# on the numerical method can be found at
#   https://arxiv.org/abs/1311.7211
#
module fieldAlignment

    using effectiveField, energy, topologicalCharge,
    effectiveSize, modifyFiles, LinearAlgebra, normalize
    export runFA!, modifyspin!

    # This function runs modifyspin! until the magnetization settles
    # to some value within a specified tolerance. The input matrix,
    # mat is modified along the trajectory.
    #
    # inputs: mat = the spin matrix, matParams = array of material parameters,
    # evalParams = array evaluation parameters, phiMatrices = array of 
    # matrices used to callculate DDI effective field and energy, 
    # defectParams = array of params defining defect values
    #
    # outputs: nothing
    function runFA!(mat::Array{Float64,3}, matParams::Array{Float64,1},
        evalParams::Array{Float64,1}, 
        phiMatrices::Array{Array{Float64,2},1} = [zeros(2,2)],
        defectParams::Array{Float64,1} = zeros(10))

        J, H, DMI, Anis, ed, pbc = matParams
        alpha, nrot, tol = evalParams

        # Initialize all relaxation parameters
        count = 0
        maxLoop = 10000
        p, m, n = size(mat)

        # Initialize parameter storage
        varCalc         = Array{Float64}(undef,maxLoop)
        reffcalc        = Array{Float64}(undef,maxLoop)
        topchargecalc   = Array{Float64}(undef,maxLoop)

        ddiField = zeros(p,m,n)

        var = 0.0
        varNew = 10.0

        # run the random spin relaxation
        while abs(var-varNew) > tol && count < maxLoop
            
            count+=1

            var = varNew

            if ed != 0.0
                ddiField = calculateDdiField(mat, ed, phiMatrices, pbc==1.0)
            end

            for i in 1:nrot modifyspin!(mat,matParams,evalParams,defectParams,ddiField) end
            normalizeSpins!(mat)

            varNew = calcEffectiveSize(mat)
            varCalc[count] = varNew

        end

        return varCalc

    end

    # inputs: mat = spin matrix, α = randomness parameter, N = number of spins
    # outputs: nothing
    #
    # The function changes the input matrix itself. (This means allocating less
    # space for the computation.) Running this function once means running through 
    # the entire spin field once and changing the spins according to some randomness
    # parameter, α. 
    function modifyspin!(mat::Array{Float64,3}, params::Array{Float64,1},
        evalParams::Array{Float64,1},defectParams::Array{Float64,1},
        ddiField::Array{Float64,3})
        
        J, H, DMI, Anis, ed, pbc = params
        alpha, nrot, tol = evalParams

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

        randomizeRunDirection(startx,starty,endx,endy,dx,dy)

        # If including DDI, calculate the effective field first. This
        # is because this particular computation takes a long time and is
        # prohibitively slow if computed N x N times. 

        for i in startx[1]:dx[1]:endx[1]
            
            for j in starty[1]:dy[1]:endy[1] 

                # scurrent is a spin x,y,x array at a certain i,j
                for k in 1:3 scurrent[k] = mat[k,i,j] end
                # println("nx = ", i, ", ny = ", j,", s0 = ", scurrent)

                # hcurrent has to start as zero
                hcurrent[1] = hcurrent[2] = hcurrent[3] = 0.0;

                # this function changes hcurrent so that its new values
                # are the effective field at i,j
                effectivefield!(mat,hcurrent,i,j,params,defectParams);
                if ed != 0.0
                    for k in 1:3 hcurrent[k] += ddiField[k, i, j] end
                end
                # println("Heff = ", hcurrent)

                # for some reason norm allocates a lot of space, so 
                # manually computed the norm, nhcurrent
                nhcurrent = sqrt(dot(hcurrent,hcurrent))

                r = rand()

                # If the random number is less than alpha, then replace 
                # the spin at that point with the normalized h field
                if r < alpha
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

    # inputs: start and end values of the randomness
    # outputs: begining and end range in both x and y in a randomnly
    # chosen direction. Returns (4,1) array
    function randomizeRunDirection(startx::Array{Int64,1},starty::Array{Int64,1},
        endx::Array{Int64,1},endy::Array{Int64,1},dx::Array{Int64,1},
        dy::Array{Int64,1})
        
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