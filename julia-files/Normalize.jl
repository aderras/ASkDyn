#=
    This module contains the function to normalize a 2D spin lattice
=#
module Normalize

    export normalizelattice!

    # in: 3D array
    #
    # out: nothing
    function normalizelattice!(mat::Array{AbstractFloat,3})

        normVal = 0.0

        p,m,n = size(mat)

        for i in 1:m, j in 1:n

            normVal = 0.0
            for l in 1:p normVal+=(mat[l,i,j]^2) end
            normVal = sqrt(normVal)

            for k in 1:p
                mat[k,i,j] = mat[k,i,j]/normVal
            end
        end

    end
end
