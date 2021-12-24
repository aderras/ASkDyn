module Helpers

    using Printf    # Use this package to control output numbers
    using HDF5
    export h5overwrite, formatFloat2Str

    function h5overwrite(name::String, array)
        fid = h5open(name, "w")
        write(fid, "Dataset1", array)
        close(fid)
    end

    @inline convFloatToDecStr(n) = @sprintf("%10.10f",n)

    function deleteTrailingZeros!(s)
        if s[end]=='0' s=s[1:end-1]; deleteTrailingZeros!(s) else return s end
    end

    @inline formatFloat2Str(n) = deleteTrailingZeros!(convFloatToDecStr(n))

    ################## DELETE WHEN DONE
    # Probably won't need the following two functions because we're
    # approximating the effective field with additional spins at edges. 
    function interpCorners!(dest)
        p,m,n = size(dest)

        # Interpolate top right corner
        for k in 1:p dest[k,m,n] = -dest[k,m-1,n-1]+dest[k,m,n-1]+dest[k,m-1,n] end
        for k in 1:p dest[k,1,1] = -dest[k,2,2]+dest[k,2,1]+dest[k,1,2] end

        # You should check the following two interpolations again.
        for k in 1:p dest[k,1,n] = dest[k,2,n-1]-dest[k,2,n]+dest[k,1,n-2] end
        for k in 1:p dest[k,m,1] = dest[k,m-1,2]+dest[k,m,2]-dest[k,m-1,1] end
    end

    # This function uses the inner N values to interpolate the edges of an
    # array. N is the length of BC.
    #
    # in: dest is the array you would like to interpolate the edges of. bc is
    # the array of coefficients to be used in the interpolation.
    #
    # The corners should probably be dealt with more carefully.
    function interpEdges!(dest, bc)

        p,m,n = size(dest)
        ll = length(bc)

        bcVal = zero(eltype(dest))

        for j in 1:n, k in 1:p
            bcVal = 0.0
            for q in 1:ll bcVal += bc[q]*dest[k,q,j] end
            dest[k,1,j] = bcVal

            bcVal = 0.0
            for q in 1:ll bcVal += bc[ll-q+1]*dest[k,m-ll+q,j] end
            dest[k,m,j] = bcVal
        end

        for i in 1:m, k in 1:p
            bcVal = 0.0
            for q in 1:ll bcVal += bc[q]*dest[k,i,q] end
            dest[k,i,1] = bcVal

            bcVal = 0.0
            for q in 1:ll bcVal += bc[ll-q+1]*dest[k,i,n-ll+q] end
            dest[k,i,n] = bcVal
        end

        interpCorners!(dest)
    end

end
