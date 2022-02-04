module Helpers

    using Printf    # Use this package to control output numbers
    using HDF5
    export h5overwrite, formatFloat2Str, zipFlatten

    # Export a file to hdf5 format
    #
    # in: name=full path of the export object, obj=array to export
    # out: nothing
    function h5overwrite(name::String, obj)
        fid = h5open(name, "w")
        write(fid, "Dataset1", obj)
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
end
