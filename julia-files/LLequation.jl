# This module contains functions that define the Landau Lifshitz Gilbert equations
#
module LLequation

    import effectiveField
    export RHS!

    # RHS computes the right side of ds/dt = ... 
    #
    # inputs: t = time of current step, mat = spin array, matRHS = right hand
    # side of the LL equation. This array is modified when calling the function
    # lambda = damping constant, matParams = material parameters, phiMatrices = 
    # matrices used to compute DDI field
    #
    # outputs: nothing
    function RHS!( t::Float64, mat::Array{Float64,3}, matRHS::Array{Float64,3}, params, flag = true )
       
        p, m, n = size(mat)

        llgParams = params.llg 

        tMax, hStep, nn, tol, lambda, T, nRuns, par = 
            [ getfield( llgParams, x ) for x in fieldnames( typeof(llgParams) ) ]


        Heff = Array{Float64}(undef,p,m,n)

        # calculate effective field 
        Heff = effectiveField.getFullEffField!( mat, params )

        # # calculate S dot H. 
        SDotH = zeros(1,m,n)       
        SDotH .= sum( mat.*Heff, dims=1 ) 
        
        fillRHS!( mat, Heff, SDotH, matRHS, lambda )


        if (params.current.jx!=0 && flag==true ) || (params.current.jy!=0 && flag==true)

            addCurrent!( mat, matRHS, params.current )

        end

    end
    
    # Right side of LL equation
    #
    # inputs: mat = (3, m, n) spin array, Heff = (3, m, n) effective field matrix, 
    # SDotH = (m, n, 1) array of SDotH, matRHS = (3, m, n) right hand side of 
    # LL equation. Updates this value
    #
    # outputs: nothing
    function fillRHS!( mat::Array{Float64,3}, Heff::Array{Float64,3}, 
        SDotH::Array{Float64,3}, matRHS::Array{Float64,3}, lambda::Float64)
    
        p, m, n = size(mat)

        for i in 1:m, j in 1:n
            matRHS[1,i,j] = mat[2,i,j]*Heff[3,i,j] - mat[3,i,j]*Heff[2,i,j] + lambda*(Heff[1,i,j]-mat[1,i,j]*SDotH[1,i,j])
            matRHS[2,i,j] = mat[3,i,j]*Heff[1,i,j] - mat[1,i,j]*Heff[3,i,j] + lambda*(Heff[2,i,j]-mat[2,i,j]*SDotH[1,i,j])
            matRHS[3,i,j] = mat[1,i,j]*Heff[2,i,j] - mat[2,i,j]*Heff[1,i,j] + lambda*(Heff[3,i,j]-mat[3,i,j]*SDotH[1,i,j])
        end
	
    end

    function addCurrent!( s::Array{Float64,3}, matRHS::Array{Float64,3}, current )

        jx = current.jx
        jy = current.jy

        p, m, n = size(s)

        for i in 1:m, j in 1:n

            # Periodic BC always used
            iNext = (i % m) + 1
            if i==1
                iPrev = m
            else
                iPrev = i-1 
            end
            jNext = (j+1) % n

            if j==1
                jPrev = n 
            else 
                jPrev = j-1
            end
            # First order finite difference used for DeltaX and Delta Y
            # It would probably be faster to compute this matrix first, rather
            # than in every iteration
            DeltaX = (1/2) * ( s[:,iNext,j] - s[:,iPrev,j] )
            DeltaY = (1/2) * ( s[:,i,jPrev] - s[:,i,jPrev] )

            # Cross product of DeltaX and s_{i,j}
            DeltaXCrossS = [ DeltaX[2]*s[3,i,j] - DeltaX[3]*s[2,i,j], 
                                -DeltaX[1]*s[3,i,j] + DeltaX[3]*s[1,i,j], 
                                DeltaX[1]*s[2,i,j] - DeltaX[2]*s[1,i,j] ]

            DeltaYCrossS = [ DeltaY[2]*s[3,i,j] - DeltaY[3]*s[2,i,j], 
                                -DeltaY[1]*s[3,i,j] + DeltaY[3]*s[1,i,j], 
                                DeltaY[1]*s[2,i,j] - DeltaY[2]*s[1,i,j] ]

            matRHS[1,i,j] += -jx* ( s[2,i,j]*DeltaXCrossS[3] - s[3,i,j]*DeltaXCrossS[2] ) - 
                              jy* ( s[2,i,j]*DeltaYCrossS[3] - s[3,i,j]*DeltaYCrossS[2] )
            matRHS[2,i,j] += -jx* ( -s[1,i,j]*DeltaXCrossS[3] + s[3,i,j]*DeltaXCrossS[1] ) - 
                              jy* ( -s[1,i,j]*DeltaYCrossS[3] + s[3,i,j]*DeltaYCrossS[1] )
            matRHS[3,i,j] += -jx* ( s[1,i,j]*DeltaXCrossS[2] - s[2,i,j]*DeltaXCrossS[1] ) - 
                              jy* ( s[1,i,j]*DeltaYCrossS[2] - s[2,i,j]*DeltaYCrossS[1] )
        end

    end

end
