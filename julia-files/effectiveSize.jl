# This module contains a function that computes the effective size
# of a skyrmion. 
#
module effectiveSize

    export calcEffectiveSize

    # The following function calculates the effective size of
    # a single skyrmion in the spin array.
    #
    # inputs: spin array (3, m, n)
    # outputs: float representing effecive size of skyrmion
    function calcEffectiveSize(mat::Array{Float64,3})
        
        p,m,n = size(mat)

        lsum = 0.0
        
        for i in 1:m
            for j in 1:n
                lsum += (mat[3,i,j])
            end
        end
        
        lsum = lsum/(m*n)
        
        return sqrt( m*n*(lsum+1)/(2*pi) )

    end

end
