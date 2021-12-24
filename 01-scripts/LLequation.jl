#=

    This module contains functions that define the Landau Lifshitz Gilbert
    equations with current.

           dS/dt = S x H_eff - lambda * S x (S x H_eff) - [Current Term]

           [Current Term] = Sum_{nu} j_{nu} S x (dS/dnu x S)

    All the 'x' symbols are cross products. The sum over nu is the sum over
    dimensions [x,y]. The effective field is calculated in the eponymous module.

=#
module LLequation

    using EffectiveField
    using LoopVectorization

    # The following are imported to test out interpolating ds/dt at the edges
    import BoundaryConditions, Helpers

    export RHS!

    # RHS computes the right side of ds/dt = ... from LLG by updating 'mat'
    #
    # in: t = current timestep, mat = current spin array, params = all the
    # parameters, relax = relaxation == true/false (default false)
    #
    # out: nothing


    # With @avx macro, this loop speeds up by 10x, but it
    # introduces an error on the order of 1e-16. Interesting
    function dotSum!(dest, mat1, mat2)
        p, m, n = size(mat1)
        # Used to have @avx, but now returns errors
        for i in 1:m, j in 1:n
            dest[i,j] = 0
            for q in 1:p dest[i,j] += mat1[q,i,j]*mat2[q,i,j] end
        end
    end

    # Improved version
    function RHS!(matRHS::Array{Float64,3},
            t::Float64,
            mat::Array{Float64,3},
            rhsArgs,
            mpValues::Array{Any,1},
            params,
            relax=false,
            damping=0.0)

        # If running relaxation, set damping to 1. Otherwise use user input.
        if relax lambda = 1.0 else lambda = damping end

        # Calculate effective field.
        effectivefield!(rhsArgs[1], mat, mpValues, params)

        # Calculate S dot H if damping is nonzero
        if lambda!=0.0 dotSum!(rhsArgs[2], mat, rhsArgs[1]) end

        fillRHS!(matRHS, mat, rhsArgs[1], rhsArgs[2], lambda)

        # Only add current if nonzero and this is dynamics (not relaxation).
        if (params.current.jx != 0.0 && relax==false) ||
            (params.current.jy != 0.0 && relax==false)
            addCurrent!(mat, matRHS, params.current)
        end

        ## DELETE WHEN DONE ###################################################
        # Instead of interpolating the effective field for fictitous spins,
        # we can interpolate the time derivative of the edge spins based on the
        # inner ones, taking care to interpolate only once at the corners.
        # pbc = mpValues[end]
        # if pbc>0
        #     bcInt = round(Int64,pbc)
        #     Helpers.interpEdges!(matRHS, BoundaryConditions.extrap[bcInt])
        # end

    end

    # LLG implemented here.
    #
    # in: mat = (3,m,n) spin array, Heff = (3,m,n) effective field
    # matrix, SDotH = (m,n,1) array of SDotH, matRHS = (3,m,n) right
    # side of LLG equation (modifies this matrix)
    #
    # out: nothing
    function fillRHS!(matRHS::Array{Float64,3}, mat::Array{Float64,3},
        Heff::Array{Float64,3}, SDotH::Array{Float64,2}, lambda::Float64)

        p,m,n = size(mat)

        for i in 1:m, j in 1:n
            matRHS[1,i,j] = mat[2,i,j]*Heff[3,i,j] - mat[3,i,j]*Heff[2,i,j]
            matRHS[2,i,j] = mat[3,i,j]*Heff[1,i,j] - mat[1,i,j]*Heff[3,i,j]
            matRHS[3,i,j] = mat[1,i,j]*Heff[2,i,j] - mat[2,i,j]*Heff[1,i,j]
        end
        if lambda != 0.0
            for i in 1:m, j in 1:n, k in 1:p
                matRHS[k,i,k] += lambda*(Heff[k,i,j]-mat[k,i,j]*SDotH[i,j])
            end
        end
    end

    # Current implemented here.
    #
    # in: s = (3,m,n) spin array, matRHS = (3,m,n) right side of LLG equation
    # (matRHS is updated in the function), current = struct from Params
    # containing current parameters
    #
    # out: nothing
    function addCurrent!(matRHS::Array{Float64,3}, s::Array{Float64,3}, current)

        jx = current.jx
        jy = current.jy

        p, m, n = size(s)

        for i in 1:m, j in 1:n

            # Periodic BC always used. Define neighbors here.
            iNext =  i%m + 1
            if i==1
                iPrev = m
            else
                iPrev = i-1
            end
            jNext = j%n + 1

            if j==1
                jPrev = n
            else
                jPrev = j-1
            end

            # First order finite difference used for deltaX and Delta Y
            # (It may be faster to compute this matrix first, rather
            # than in every iteration.)
            deltaX = (1/2) * (s[:,iNext,j] - s[:,iPrev,j])
            deltaY = (1/2) * (s[:,i,jNext] - s[:,i,jPrev])

            # Cross product of deltaX and s_{i,j}
            deltaXCrossS = [ deltaX[2]*s[3,i,j] - deltaX[3]*s[2,i,j],
                                -deltaX[1]*s[3,i,j] + deltaX[3]*s[1,i,j],
                                deltaX[1]*s[2,i,j] - deltaX[2]*s[1,i,j] ]

            deltaYCrossS = [ deltaY[2]*s[3,i,j] - deltaY[3]*s[2,i,j],
                                -deltaY[1]*s[3,i,j] + deltaY[3]*s[1,i,j],
                                deltaY[1]*s[2,i,j] - deltaY[2]*s[1,i,j] ]

            matRHS[1,i,j] += -jx* (s[2,i,j]*deltaXCrossS[3] -
                    s[3,i,j]*deltaXCrossS[2]) - jy* (s[2,i,j]*deltaYCrossS[3] -
                    s[3,i,j]*deltaYCrossS[2])
            matRHS[2,i,j] += -jx* (-s[1,i,j]*deltaXCrossS[3] +
                    s[3,i,j]*deltaXCrossS[1]) - jy* (-s[1,i,j]*deltaYCrossS[3] +
                    s[3,i,j]*deltaYCrossS[1])
            matRHS[3,i,j] += -jx* (s[1,i,j]*deltaXCrossS[2] -
                    s[2,i,j]*deltaXCrossS[1]) - jy* (s[1,i,j]*deltaYCrossS[2] -
                    s[2,i,j]*deltaYCrossS[1])
        end

    end

end
