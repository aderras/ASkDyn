module DefectFunctions

    export exchangefield_defect!, exchangefieldelem_defect!,
    exchange_nergy_defect

    # Computes the exchange energy of mat
    #
    # in: mat = spin matrix, J = exchange constant, pbc = periodic boundary
    # conditions, defParams = struct containing information on the exchange
    # defect
    #
    # out: float
    function exchange_energy_defect(mat::Array{Float64,3}, J::Float64,
        pbc::Float64, dParams)

        p, m, n = size(mat)
        en = 0.0


        # If there is a defect and the point of interest is within the range
        # of the defect, change the total exchange energy
        defType,aJ,dJ,jx,jy = [getfield(dParams, x)
            for x in fieldnames(typeof(dParams))]

        # Exchange energy computation differs for the type of defect. If there
        # is not a defect or if defExists == 1.0 (point defect) then compute
        # the following
        if defType == 0.0 || defType == 1.0

            for j in 1:n, i in 1:m-1, k in 1:p
                en += mat[k,i,j] * mat[k,i+1,j]
            end

            for j in 1:n-1, i in 1:m, k in 1:p
                en += mat[k,i,j] * mat[k,i,j+1]
            end

            if pbc==1.0
                for j in 1:n, k in 1:p
                    en += mat[k,1,j] * mat[k,m,j]
                end
                for i in 1:m, k in 1:p
                    en += mat[k,i,1] * mat[k,i,n]
                end

                en -= (2 * m * n)

            elseif pbc==2.0
                for j in 1:n, k in 1:p
                    en += mat[k,m,j] * rightN[k,j]
                end
                for i in 1:m, k in 1:p
                    en += mat[k,i,n] * bottomN[k,i]
                end
                en -= (2 * m * n)
            else
                en -= (2 * m * n - m - n)
            end

            # If there's a point defect, modify the total energy by changing
            # the contribution of the nearest dJ bonds
            if 1.0 == defType

                # Find the number of sites affected by the defect
                nsites = 1 + 2*dJ

                # The exchange energy at these sites is modified by aJ amount
                en = (en + aJ*nsites)

            end

            return -J * en

        # If there is a gaussian-type defect, compute the total energy using an
        # exponential
        elseif defType == 2.0

            for j in 1:n, i in 1:m-1, k in 1:p
                en += -J * dParams.jMat[4][i,j] * mat[k,i,j] * mat[k,i+1,j]
            end

            for j in 1:n-1, i in 1:m, k in 1:p
                en += -J * dParams.jMat[2][i,j] * mat[k,i,j] * mat[k,i,j+1]
            end

            en += (2 * m * n - m - n)

            return en
        end
    end

    # This function modifies effField to store the exchange field of the array,
    # mat, which contains a Gaussian-type defect. Nearest neighbor spins are
    # considered.
    #
    # in: mat = (3, m, n), effField = (3, 1) array used to store answer,
    # nx, ny = positions in mat to calculate exchange field, params = struct
    # of all computation parameters
    #
    # out: nothing
    function exchangefieldelem_defect!(effField::Array{Float64,1},
        mat::Array{Float64,3}, nx::Int, ny::Int, params)

        defType,aJ,dJ,jx,jy = [getfield(params.defect, x)
            for x in fieldnames(typeof(params.defect))]
        J = params.mp.j
        pbc = params.mp.pbc

        # This is the Gaussian-type of defect we consider. The exchange
        # constant is modified by a maximum of aJ at the location of the
        # defect, (jx,jy).
        @inline Jmod(i,j) = J*(1 + aJ*exp(-((i - jx)^2 + (j -jy)^2)/dJ^2))

        p, m, n = size(mat)

        if nx > 1 && nx < m
            for k in 1:3
                effField[k] = effField[k] +
                    Jmod(nx-1/2, ny) * (mat[k,nx-1,ny] + mat[k,nx+1,ny])
            end
        elseif nx == 1
            if pbc==1.0
                for k in 1:3
                    effField[k] = effField[k] + Jmod(nx-1/2, ny) * mat[k,m,ny]
                end
            end

            for k in 1:3
                effField[k] = effField[k] + Jmod(nx+1/2, ny) * mat[k,nx+1,ny]
            end
        elseif nx == m
            if pbc==1.0
                for k in 1:3
                    effField[k] = effField[k] + Jmod(nx+1/2, ny) * mat[k,1,ny]
                end
            end

            for k in 1:3
                effField[k] = effField[k] + Jmod(nx+1/2, ny) * mat[k,nx-1,ny]
            end
        end


        if ny > 1 && ny < n
            for k in 1:3
                effField[k] = effField[k] +
                    Jmod(nx, ny-1/2) * (mat[k,nx,ny-1] + mat[k,nx,ny+1])
            end
        elseif ny == 1
            if pbc==1.0
                for k in 1:3
                    effField[k] = effField[k] + Jmod(nx, ny-1/2) * mat[k,nx,n]
                end
            end
            for k in 1:3
                effField[k] = effField[k] + Jmod(nx, ny+1/2) * mat[k,nx,ny+1]
            end

        elseif ny == n
            if pbc==1.0
                for k in 1:3
                    effField[k] = effField[k] + Jmod(nx, ny+1/2) * mat[k,nx,1]
                end
            end
            for k in 1:3
                effField[k] = effField[k] + Jmod(nx, ny+1/2) * mat[k,nx,ny-1]
            end
        end
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
    function exchangefield_defect!(Heff::Array{Float64,3},
        mat::Array{Float64,3}, params, edges=[])

        p, m, n = size(mat)
        pbc = params.mp.pbc

        # There's a tiny improvement in benchmark speed by writing the for
        # loop this way. Not sure if it's real improvement, but I'll take what
        # I can get.
        for k in 1:3
            for ny in 1:m
                for nx in 1:n

                    if pbc==1.0 || 2.0
                        nxNext = nx%m + 1
                        nyNext = ny%n + 1

                        if nx == 1
                            nxPrev = m
                        else
                            nxPrev = nx-1
                        end
                        if ny == 1
                            nyPrev = n
                        else
                            nyPrev = ny-1
                        end

                        Heff[k,nx,ny] = Heff[k,nx,ny] +
                            params.defect.jMat[1][nx,ny] * mat[k,nxPrev,ny] +
                            params.defect.jMat[2][nx,ny] * mat[k,nxNext,ny] +
                            params.defect.jMat[3][nx,ny] * mat[k,nx,nyPrev] +
                            params.defect.jMat[4][nx,ny] * mat[k,nx,nyNext]
                    else
                        if nx > 1
                            Heff[k,nx,ny] = Heff[k,nx,ny] +
                                params.defect.jMat[1][nx,ny]*mat[k,nx-1,ny]
                        end
                        if ny > 1
                            Heff[k,nx,ny] = Heff[k,nx,ny] +
                                params.defect.jMat[3][nx,ny]*mat[k,nx,ny-1]
                        end
                        if nx < m
                            Heff[k,nx,ny] = Heff[k,nx,ny] +
                                params.defect.jMat[2][nx,ny]*mat[k,nx+1,ny]
                        end
                        if ny < n
                            Heff[k,nx,ny] = Heff[k,nx,ny] +
                                params.defect.jMat[4][nx,ny]*mat[k,nx,ny+1]
                        end
                    end
                end
            end
        end
    end

    # Calculates the gaussian-like modification for each type of bond (left,
    # right, top, bottom) at every point in the lattice.
    #
    # in: aJ = strength of modification, dJ = widfh of modification, jx =
    # position of the defect in x, jy = position of the defect in y.
    #   NOTE: This only works with the rest of the program if the defect is at
    #   the center of the lattice. Need to generalize this.
    #
    # out: [A1,A2,A3,A4] where each A is an (m,n) matrix
    function defectarray(aJ, dJ, jx, jy)

        nx = round(Int, jx*2)
        ny = round(Int, jy*2)

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
