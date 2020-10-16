# This module contains the normalizeSpins function
#
#

module normalize
    
    export normalizeSpins!
    
    function normalizeSpins!(mat::Array{Float64,3})
        
        normVal = 0.0

        p,m,n = size(mat)

        for i in 1:m
            for j in 1:n
                
                normVal = 0.0
                for l in 1:p normVal+=(mat[l,i,j]^2) end
                normVal = sqrt(normVal)

                for k in 1:p
                    mat[k,i,j] = mat[k,i,j]/normVal
                end
            end
        end
    
    end

end
