#=

    This module contains functions that calculate the energy of a
    spin lattice.

    IMPROVEMENTS
      + It would be more logical if defect params were an optional input for
      exchange_energy

=#
module Energy

    import DipoleDipole, InitialCondition
    export energy, exchange_energy, zeeman_energy, dmi_energy,
        pma_energy, ddi_energy

    # Calculates the total energy of matrix mat for material parameters
    # in param
    #
    # in: spin array (3, m, n), matParams = material parameters,
    # defectParams = optional array argument containing information about
    # a defect in the lattice, phiMatrices = array of matrices used to
    # compute DDI (optional argument)
    #
    # out: Float
    function energy(mat::Array{Float64,3}, params)

        j,h,a,dz,ed,nx,ny,nz,pbc,vd = [getfield(params.mp, x)
            for x in fieldnames(typeof(params.mp))]

        extrap = []
        if pbc==2.0 extrap = InitialCondition.extrapolateEdges(mat) end

        en = exchange_energy(mat, j, pbc, params.defect, extrap)

        if h != 0.0
            en += zeeman_energy(mat, h)
        end
        if a != 0.0
            en += dmi_energy(mat, a, pbc, extrap)
        end
        if dz != 0.0
            en += pma_energy(mat, dz)
        end
        if ed != 0.0
            en += ddi_energy(mat, ed, pbc, vd)
        end

        return en
    end

    # Computes the zeeman energy of mat for external field Hext
    #
    # in: mat = spin matrix, h = external field
    #
    # out: float
    function zeeman_energy(mat::Array{Float64,3}, h::Float64)

        p, m, n = size(mat)

        en = 0.0

        for i in 1:n
            for j in 1:m
                en += -h * mat[3,i,j]
            end
        end

        return en

    end

    # Computes the dmi energy of mat for DMI constant and optional defect
    #
    # in: mat = spin matrix, DMI = DMI constant, pbc = periodic boundary
    # conditions
    #
    # out: float
    function dmi_energy(mat::Array{Float64,3}, a::Float64, pbc::Float64,
        edges=[])

        p, m, n = size(mat)
        en = 0.0

        if pbc==2.0 leftN, rightN, topN, bottomN = edges end

        # The following is for nx > 1
        for i in 2:m
            for j in 1:n
                en += mat[2,i,j]*mat[3,i-1,j] - mat[3,i,j]*mat[2,i-1,j]
            end
        end

        # The following is for ny > 1
        for i in 1:m
            for j in 2:n
                en+=mat[3,i,j]*mat[1,i,j-1] - mat[1,i,j]*mat[3,i,j-1]
            end
        end

        # Deal with the edge of the matrix. If periodic boundary
        if pbc==1.0
            for i in 1:n
                en += (mat[2,1,i]*mat[3,m,i] - mat[3,1,i]*mat[2,m,i])
            end
            for i in 1:m
                en += (mat[3,i,1]*mat[1,i,n] - mat[1,i,1]*mat[3,i,n])
            end

        elseif pbc==2.0
            for i in 1:n
                en += bottomN[2,i]*mat[3,m,i] - bottomN[3,i]*mat[2,m,i]
            end
            for i in 1:m
                en += rightN[3,i]*mat[1,i,n] - rightN[1,i]*mat[3,i,n]
            end
        end

        en = -a * en

        return en

    end


    # Computes the exchange energy of mat
    #
    # in: mat = spin matrix, J = exchange constant, pbc = periodic boundary
    # conditions, defParams = struct containing information on the exchange
    # defect
    #
    # out: float
    function exchange_energy(mat::Array{Float64,3}, J::Float64,pbc::Float64,
        dParams, edges=[])

        p, m, n = size(mat)
        en = 0.0

        if pbc==2.0 leftN, rightN, topN, bottomN = edges end

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


    # Computes the DDI energy of mat
    #
    # in: mat = spin matrix, ed = DDI constant, pbc = periodic boundary
    # conditions phiMatrices = array of matrices used to compute DDI
    #
    # out: float
    function ddi_energy(mat::Array{Float64,3}, ed::Float64,
        pbc::Float64, phiMatrices::Array{Array{Float64,2},1})

        if ed == 0
            return 0.0
        end

        p, m, n     = size(mat)
        ddiEnArray  = Array{Float64}(undef,p,m,n)

        ddiEnArray  = DipoleDipole.fhd(mat,phiMatrices,pbc)

        tot = 0.0

        for i in eachindex(view(mat,:,:,:))
            tot += mat[i] * 0.5 * ddiEnArray[i]
        end

        return -ed * tot

    end

    # Computes the anisotropy energy of mat
    #
    # in: mat = spin matrix, dz = PMA constant
    #
    # out: float
    function pma_energy(mat::Array{Float64,3}, dz::Float64)

        en = 0.0
        p, m, n = size(mat)

        for i in 1:m
            for j in 1:n
                en -= 0.5 * dz * (mat[3,i,j] * mat[3,i,j] - 1)
            end
        end

        return en
    end

end
