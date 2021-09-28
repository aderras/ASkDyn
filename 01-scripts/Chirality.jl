
module Chirality

    # in: (3,Nx,Ny) matrix containing spin values everywhere
    # float: chirality value
    function computeGamma(mat, params)

        p,m,n = size(mat)

        xPos = round(Int64, params.cp.rChirality + m/2)
        yPos = round(Int64, n/2)

        sChoice = mat[:,xPos,yPos]

        return acos(sChoice[1]/sqrt(1-sChoice[3]^2))

    end

end
