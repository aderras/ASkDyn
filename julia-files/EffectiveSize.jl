module EffectiveSize

    export effectivesize

    # The following function calculates the effective size of
    # a single skyrmion in the spin array.
    #
    # in: spin array (3, m, n)
    # out: float
    function effectivesize(mat::Array{AbstractFloat,3})

        p, m, n = size(mat)

        lsum = 0.0

        for i in 1:m, j in 1:n
            lsum += (mat[3,i,j])
        end

        lsum = lsum/(m*n)

        return sqrt( m*n*(lsum+1)/(2*pi) )

    end

end
