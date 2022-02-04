#=

    This module contains functions that calculate the energy of a
    spin lattice.

=#
module Energy

    import Dipolar, InitialCondition
    import DefectFunctions
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
    function energy(mat::Array{Float64,3}, mpValues)


        j,h,a,ed,dz,vd,pbc = mpValues
        en = exchange_energy(mat, j, pbc)

        if h != 0.0 en += zeeman_energy(mat, h) end
        if a != 0.0 en += dmi_energy(mat, a, pbc) end
        if dz != 0.0 en += pma_energy(mat, dz) end
        if ed != 0.0 en += Energy.ddi_energy(mat, ed, pbc, vd) end

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
    function dmi_energy(mat::Array{Float64,3}, a::Float64, pbc::Float64)

        p, m, n = size(mat)
        en = 0.0

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
        end

        en = -a * en

        return en

    end

    # Computes the exchange energy of mat
    #
    # in: mat = spin matrix, J = exchange constant, pbc = periodic boundary
    # conditions
    #
    # out: float
    function exchange_energy(mat::Array{Float64,3}, J::Float64, pbc::Float64)

        p, m, n = size(mat)
        en = 0.0

        for j in 1:n, i in 1:m-1, k in 1:p
            en += mat[k,i,j] * mat[k,i+1,j]
        end

        for j in 1:n-1, i in 1:m, k in 1:p
            en += mat[k,i,j] * mat[k,i,j+1]
        end

        if pbc==1.0 # periodic boundary conditions
            for j in 1:n, k in 1:p
                en += mat[k,1,j] * mat[k,m,j]
            end
            for i in 1:m, k in 1:p
                en += mat[k,i,1] * mat[k,i,n]
            end
            en -= (2 * m * n)
        else
            en -= (2 * m * n - m - n)
        end
        return -J*en
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

        ddiEnArray  = Dipolar.fhd(mat,phiMatrices,pbc)

        tot = 0.0

        for i in eachindex(view(mat,:,:,:))
            tot += mat[i] * 0.5 * ddiEnArray[i]
        end

        return -ed * tot

    end

    # Computes the anisotropy energy of mat. dz>0 is uniaxial anisotropy in the
    # z direction. (dz<0 would mean easy plane anisotropy in the xy plane)
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
