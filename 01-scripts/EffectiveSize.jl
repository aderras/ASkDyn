module EffectiveSize

    export effectivesize

    # The following function calculates the effective size of
    # a single skyrmion in the spin array.
    #
    # in: spin array (3, m, n)
    # out: float
    function effectivesize(mat::Array{Float64,3})

        p,m,n = size(mat)
        lsum = 0.0
        mLambda=4

        for i in 1:m, j in 1:n lsum += (mat[3,i,j]+1)^mLambda end

        return sqrt((mLambda-1)*lsum/(2^mLambda*pi))
    end

end
