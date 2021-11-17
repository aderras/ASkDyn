
module Chirality

    # in: (3,Nx,Ny) matrix containing spin values everywhere
    # float: chirality value
    function computeGamma(sChoice, mat, rMeasure)

        p,m,n = size(mat)

        xPos = round(Int64, rMeasure + m/2)
        yPos = round(Int64, n/2)

        for i in 1:p sChoice[i] = mat[i,xPos,yPos] end

        return acos(sChoice[1]/sqrt(1-sChoice[3]^2))

    end

end
