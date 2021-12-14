module Helpers

    using HDF5
    export h5overwrite

    function h5overwrite(name::String, array)
        fid = h5open(name, "w")
        write(fid, "Dataset1", array)
        close(fid)
    end

end
