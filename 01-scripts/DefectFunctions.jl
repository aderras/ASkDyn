module DefectFunctions

    # Computes the exchange energy of mat
    #
    # in: mat = spin matrix, J = exchange constant, pbc = periodic boundary
    # conditions, defParams = struct containing information on the exchange
    # defect
    #
    # out: float
    function exchange_energy(mat::Array{Float64,3}, jMat::Vector{Array{Float64,2}}, pbc)

        p, m, n = size(mat)
        en = 0.0

        # Exchange energy computation differs for the type of defect. If there
        # is not a defect or if defExists == 1.0 (point defect) then compute
        # # the following
        # if defType == 1.0
        #     for j in 1:n, i in 1:m-1, k in 1:p
        #         en += mat[k,i,j] * mat[k,i+1,j]
        #     end
        #     for j in 1:n-1, i in 1:m, k in 1:p
        #         en += mat[k,i,j] * mat[k,i,j+1]
        #     end
        #     if pbc # periodic boundary conditions
        #         for j in 1:n, k in 1:p
        #             en += mat[k,1,j] * mat[k,m,j]
        #         end
        #         for i in 1:m, k in 1:p
        #             en += mat[k,i,1] * mat[k,i,n]
        #         end
        #         en -= (2 * m * n)
        #     else
        #         en -= (2 * m * n - m - n)
        #     end
        #     return -J*en
        #     # If there's a point defect, modify the total energy by changing
        #     # the contribution of the nearest dJ bonds
        #     # Find the number of sites affected by the defect
        #     nsites = 1 + 2*dJ
        #
        #     # The exchange energy at these sites is modified by aJ amount
        #     en = (en + aJ*nsites)
        #
        #     return -J * en
        #
        # # If there is a gaussian-type defect, compute the total energy using an
        # # exponential
        # elseif defType == 2.0
            for j in 1:n, i in 1:m-1, k in 1:p
                en += -jMat[4][i,j] * mat[k,i,j] * mat[k,i+1,j]
            end
            for j in 1:n-1, i in 1:m, k in 1:p
                en += -jMat[2][i,j] * mat[k,i,j] * mat[k,i,j+1]
            end
            en += (2 * m * n - m - n)
            return en
        # end
    end

    # Compute the exchange field of the entire spin array, mat, that contains
    # an exchange-modifying defect. Nearest neighbors considered. We precompute
    # the exchange-modification array and store it in params.defect.jMat[k],
    # where k distinguishes between left neighbor, right neighbor, etc. We
    # distinguis all four because the exchange reduction occurs at the bond
    # between lattice points, so to compute the effective field at lattice site
    # (nx,ny), we need the locations of all the bonds connected to it. This is
    # time-intensive, so precomputing improves speed.
    # NOTE: EXTRAPOLATED BC NOT IMPLEMENTED
    #
    # in: mat = (3,m,n) spin array, params = struct of all material params,
    # Heff = (3,m,n) array of the effective field which is modified to store
    # the result
    #
    # out: nothing
    function exchangefield!(Heff::Array{Float64,3},
        mat::Array{Float64,3},jmat, bc)

        p,m,n = size(mat)

        for j in 1:n, i in 1:m-1, k in 1:p
            Heff[k,i,j] += jmat[2][i,j]*mat[k,i+1,j]
        end
        for j in 1:n, i in 2:m, k in 1:p
            Heff[k,i,j] += jmat[1][i,j]*mat[k,i-1,j]
        end
        for j in 1:n-1, i in 1:m, k in 1:p
            Heff[k,i,j] += jmat[4][i,j]*mat[k,i,j+1]
        end
        for j in 2:n, i in 1:m, k in 1:p
            Heff[k,i,j] += jmat[3][i,j]*mat[k,i,j-1]
        end

        if bc==1.0
            for j in 1:n, k in 1:p
                Heff[k,m,j] += jmat[2][1,j]*mat[k,1,j]
                Heff[k,1,j] += jmat[1][m,j]*mat[k,m,j]
            end
            for i in 1:m, k in 1:p
                Heff[k,i,1] += jmat[3][i,n]*mat[k,i,n]
                Heff[k,i,n] += jmat[4][i,1]*mat[k,i,1]
            end
        elseif bc>1.0
            # bcInt = round(Int64,bc)
            # addbc!(Heff, mat, BoundaryConditions.extrap[bcInt])
        end
    end

    # Calculates the gaussian-like modification for each type of bond (left,
    # right, top, bottom) at every point in the lattice.
    #
    # in: aJ = strength of modification, dJ = widfh of modification, jx =
    # position of the defect in x, jy = position of the defect in y.
    #
    # out: [A1,A2,A3,A4] where each A is an (nx,ny) matrix
    function buildJmats(nx, ny, defType, aJ, dJ, jx, jy)

        left = zeros(nx,ny)
        right = zeros(nx,ny)
        top = zeros(nx,ny)
        bott = zeros(nx,ny)

        for i in 1:nx, j in 1:ny
            left[i,j] = (1+aJ*exp(-((i-jx-1/2)^2 + (j-jy)^2)/dJ^2))
            right[i,j] = (1+aJ*exp(-((i-jx+1/2)^2 + (j-jy)^2)/dJ^2))
            top[i,j] = (1+aJ*exp(-((i-jx)^2 + (j-jy-1/2)^2)/dJ^2))
            bott[i,j] = (1+aJ*exp(-((i-jx)^2 + (j-jy+1/2)^2)/dJ^2))
        end

        return [left, right, top, bott]
    end

end
