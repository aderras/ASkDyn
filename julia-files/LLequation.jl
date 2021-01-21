#=

    This module contains functions that define the Landau Lifshitz Gilbert
    equations with current.

           dS/dt = S x H_eff - lambda * S x (S x H_eff) - [Current Term]

           [Current Term] = Sum_{nu} j_{nu} S x (dS/dnu x S)

    All the 'x' symbols are cross products. The sum over nu is the sum over
    dimensions [x,y]. The effective field is calculated in the eponymous module.

=#
module LLequation

    import EffectiveField
    export RHS!

    # RHS computes the right side of ds/dt = ... from LLG by updating 'mat'
    #
    # in: t = current timestep, mat = current spin array, params = all the
    # parameters, relax = relaxation == true/false (default false)
    #
    # out: nothing
    function RHS!(t::AbstractFloat, mat::Array{AbstractFloat,3},
        matRHS::Array{AbstractFloat,3}, params, relax=false)

        p, m, n = size(mat)

        # If running relaxation, set damping to 1. Otherwise use user input.
        if relax
            lambda = 1.0
        else
            lambda = params.llg.damp
        end

        Heff = Array{AbstractFloat}(undef,p,m,n)

        # Calculate effective field.
        Heff = EffectiveField.effectivefield(mat, params)

        # Calculate S dot H.
        SDotH = zeros(1,m,n)
        SDotH .= sum(mat.*Heff, dims=1)

        fillRHS!(mat, Heff, SDotH, matRHS, lambda)

        # Only add current if nonzero and this is dynamics (not relaxation).
        if (params.current.jx != 0.0 && relax==false) ||
            (params.current.jy != 0.0 && relax==false)

            addCurrent!(mat, matRHS, params.current)

        end

    end

    # LLG implemented here.
    #
    # in: mat = (3,m,n) spin array, Heff = (3,m,n) effective field
    # matrix, SDotH = (m,n,1) array of SDotH, matRHS = (3,m,n) right
    # side of LLG equation (modifies this matrix)
    #
    # out: nothing
    function fillRHS!(mat::Array{AbstractFloat,3}, Heff::Array{AbstractFloat,3},
        SDotH::Array{AbstractFloat,3}, matRHS::Array{AbstractFloat,3},
        lambda::AbstractFloat)

        p, m, n = size(mat)

        for i in 1:m, j in 1:n
            matRHS[1,i,j] = mat[2,i,j]*Heff[3,i,j] - mat[3,i,j]*Heff[2,i,j] +
                lambda*(Heff[1,i,j]-mat[1,i,j]*SDotH[1,i,j])
            matRHS[2,i,j] = mat[3,i,j]*Heff[1,i,j] - mat[1,i,j]*Heff[3,i,j] +
                lambda*(Heff[2,i,j]-mat[2,i,j]*SDotH[1,i,j])
            matRHS[3,i,j] = mat[1,i,j]*Heff[2,i,j] - mat[2,i,j]*Heff[1,i,j] +
                lambda*(Heff[3,i,j]-mat[3,i,j]*SDotH[1,i,j])
        end

    end

    # Current implemented here.
    #
    # in: s = (3,m,n) spin array, matRHS = (3,m,n) right side of LLG equation
    # (matRHS is updated in the function), current = struct from UserInputs
    # containing current parameters
    #
    # out: nothing
    function addCurrent!(s::Array{AbstractFloat,3},
        matRHS::Array{AbstractFloat,3}, current)

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
