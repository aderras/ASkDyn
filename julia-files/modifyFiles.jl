# The following module contains functions related to file saving and reading.

module modifyFiles

    using HDF5
    export getH5data,writeDataH5

    # Imports spin data
    #
    # inputs: filename = string of file path relative to code directory
    #
    # outputs: mat = (N,N,3) array
    # 
    function getH5data( filename::String )
        
        # import initial data. The following are HDF5 objects
        rawfile = h5open(pwd()*filename,"r")
        latticeData = rawfile["Dataset1"]

        # extract matrix size
        p, m, n = size(latticeData)

        # initialize storage so I can use regular matrix
        mat = Array{Float64}(undef,p,m,n)
        mat[:,:,:] = latticeData[:,:,:] #for loop wasn't working for some reason

        close(rawfile)

        return mat

    end

    # Export data in H5 format
    #
    # inputs: filename = string of relative file path. data = object you
    # want to save. 
    #
    # outputs: nothing
    #
    function writeDataH5( filename::String, data )

        file = h5open(pwd()*filename,"w")
        write(file,"Dataset1",data)
        close(file)

    end

end
