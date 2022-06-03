module Helpers

    using Printf    # Use this package to control output numbers
    using HDF5
    export h5overwrite, formatFloat2Str, zipFlatten, h5overwriteMulti

    # Export a file to hdf5 format
    #
    # in: name=full path of the export object, obj=array to export
    # out: nothing
    function h5overwrite(name::String, obj, fieldname="Dataset1")
        fid = h5open(name, "w")
        write(fid, fieldname, obj)
        close(fid)
    end

    function h5overwriteMulti(filename::String, fieldnames, objs)
        if length(objs)!=length(fieldnames)
            print("Error saving data: length(objs)!=length(fieldnames)")
        end

        fid = h5open(filename, "w")
        for i in 1:length(objs)
            write(fid, fieldnames[i], objs[i])
        end
        close(fid)

    end

    # Convert a float in scientific notation to a string in standard form
    #
    # in: n=number to convert formatted XeY
    # out: string representing standard form of the input
    @inline convFloatToDecStr(n) = @sprintf("%10.10f",n)

    # Remove ending zeros from a string
    #
    # in: s=string respresenting a number
    # out: string representing the same number with trailing zeros removed
    function deleteTrailingZeros!(s)
        if s[end]=='0' s=s[1:end-1]; deleteTrailingZeros!(s) else return s end
    end

    # Format a float which normally would be written in scientific notation
    # to one in standard form
    #
    # in: n=float
    # out: n as a string formatted in standard form.
    @inline formatFloat2Str(n) = deleteTrailingZeros!(convFloatToDecStr(n))

    # Make one string from values of vectors a and b. Used to generate a
    # filename.
    #
    # in: a, b are vectors of the same length
    # out: c is a string made up of elements of a and b such that
    # c = "_a1b1_a2b2_..._aNbN"
    function zipFlatten(a, b)
        length(a)==length(b) || error("Length mismatch")
        c = ""
        for (x, y) in zip(a, b) c = string(c, x, y) end
        return c
    end


    #=
        The following functions are used to compute convolution with a
        nonuniform Fourier transform. May delete if that idea doesn't work.
    =#
    function skipGrid1dIndices(N::Int, skip::Int=2)
        indices = Array{Int64,1}(undef,N)
        nMax = round(Int64,N/skip)
        indices = [skip*n for n in 1:nMax]
        return indices
    end

    function skipGrid2dIndices(Nx::Int, Ny::Int, skip::Int=2)
        colValues = skipGrid1dIndices(Nx, skip)
        rowValues = skipGrid1dIndices(Ny, skip)

        rowcols = Array{Int,3}(undef,2, length(rowValues), length(colValues))
        for i in 1:length(rowValues)
            for j in 1:length(colValues)
                rowcols[1,i,j] = rowValues[i]
                rowcols[2,i,j] = colValues[j]
            end
        end
        return rowcols
    end
end
