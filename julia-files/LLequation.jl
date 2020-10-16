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
    function RHS!( t::Float64, mat::Array{Float64,3}, matRHS::Array{Float64,3}, params )
       
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

end
