# This module contains functions that calculate the energy of a 
# spin lattice. 
#
# Useful functions: 
#   calcEnergy computes the total energy of the lattice 
module energy

    import dipoleDipole
    export calcEnergy, exchangeEnergy, zeemanEnergy, dmiEnergy,
    pmaEnergy, ddiEnergy


    # Calculates the total energy of matrix mat for material parameters
    # in param
    #
    # inputs: spin array (3, m, n), matParams = material parameters,
    # defectParams = optional array argument containing information about
    # a defect in the lattice, phiMatrices = array of matrices used to 
    # compute DDI (optional argument)
    # 
    # outputs: Float
    #
    function calcEnergy( mat::Array{Float64,3}, params, dParams=[] )

        j,h,a,dz,ed,nx,ny,nz,pbc,vdd = [ getfield( params, x ) for x in fieldnames( typeof(params) ) ]

        totalenergy = exchangeEnergy(mat, j, pbc==1.0, dParams)

        if h != 0.0
            totalenergy += zeemanEnergy(mat, h) 
        end
        if a != 0.0
            totalenergy += dmiEnergy(mat, a, pbc==1.0)
        end
        if dz != 0.0
            totalenergy += pmaEnergy(mat, dz)
        end
        if ed != 0.0
            totalenergy += ddiEnergy(mat, ed, pbc, vdd)
        end
        
        return totalenergy
    end

    # Computes the zeeman energy of mat for external field Hext
    #
    # inputs: mat = spin matrix, Hext = external field
    #
    # outputs: float
    #
    function zeemanEnergy( mat::Array{Float64,3}, Hext::Float64 )

        p, m, n = size(mat)

        en = 0.0
        
        for i in 1:n 
            for j in 1:m
                en += -Hext * mat[3,i,j]
            end
        end
        
        return en
    
    end 

    # Computes the dmi energy of mat for DMI constant and optional defect
    #
    # inputs: mat = spin matrix, DMI = DMI constant, pbc = periodic boundary conditions
    # defParams = array containing information on the dmi defect
    #
    # outputs: float
    #
    function dmiEnergy( mat::Array{Float64,3}, DMI::Float64, pbc::Bool,
        defParams::Array{Float64}=zeros(5) )

        p, m, n = size(mat)
        en = 0.0

        # The following is for nx > 1
        for i in 2:m
            for j in 1:n
                en += mat[2,i,j] * mat[3,i-1,j] - mat[3,i,j]*mat[2,i-1,j]
            end
        end
        
        # The following is for ny > 1
        for i in 1:m
            for j in 2:n
                en+=mat[3,i,j]*mat[1,i,j-1] - mat[1,i,j]*mat[3,i,j-1]
            end
        end
        
        # Deal with the edge of the matrix. If periodic boundary 
        # conditions, run the following. Otherwise do nothing
        if pbc
            for i in 1:n 
                en += (mat[2,1,i]*mat[3,m,i] - mat[3,1,i]*mat[2,m,i])
            end
            
            for i in 1:m 
                en += (mat[3,i,1]*mat[1,i,n] - mat[1,i,1]*mat[3,i,n])
            end
        end

        # If there is a defect and the point of interest is within the range
        # of the defect, change the total dmi energy
        defExists,jx,jy,aDmi,dDmi = defParams
        if 1.0 == defExists 
            
            # Find the number of sites affected by the defect
            nsites = 1 + 2*dDmi

            # The exchange energy at these sites is modified by aDmi amount
            en = en + (aDmi*nsites)
        end

        en = -DMI * en
        
        return en

    end


    # Computes the exchange energy of mat 
    #
    # inputs: mat = spin matrix, J = exchange constant, pbc = periodic boundary conditions
    # defParams = array containing information on the exchange defect (optional)
    #
    # outputs: float
    #
    function exchangeEnergy(mat::Array{Float64,3}, J::Float64, pbc::Bool,
        dParams) 

        p, m, n = size(mat)
        en = 0.0
        
        # If there is a defect and the point of interest is within the range
        # of the defect, change the total exchange energy
        if dParams != []
            defType,aJ,dJ,jx,jy = [ getfield( dParams, x ) for x in fieldnames( typeof(dParams) ) ]
        else
            defType,aJ,dJ,jx,jy = zeros(5)
        end

        # Exchange energy computation differes for the type of defect. If there is not
        # a defect or if defExists == 1.0 (point defect) then compute the following
        if defType == 0.0 || defType == 1.0

            for j in 1:n       
                for i in 1:m-1
                    for k in 1:p
                        en += mat[k,i,j] * mat[k,i+1,j]
                    end
                end
            end

            for j in 1:n-1    
                for i in 1:m   
                    for k in 1:p
                        en += mat[k,i,j] * mat[k,i,j+1]
                    end
                end
            end

            if pbc
                for j in 1:n
                    for k in 1:p
                        en += mat[k,1,j] * mat[k,m,j]
                    end
                end
                for i in 1:m
                    for k in 1:p
                        en += mat[k,i,1] * mat[k,i,n]
                    end
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

            for j in 1:n       
                for i in 1:m-1
                    for k in 1:p
                        en += -J * (1 + aJ * exp( -((i-jx + 1/2)^2 + (j-jy)^2)/dJ^2 ) ) * mat[k,i,j] * mat[k,i+1,j]
                    end
                end
            end

            for j in 1:n-1    
                for i in 1:m   
                    for k in 1:p
                        en += -J * (1 + aJ * exp( -((i-jx)^2 + (j-jy - 1/2)^2)/dJ^2 ) ) * mat[k,i,j] * mat[k,i,j+1]
                    end
                end
            end

            en += (2 * m * n - m - n)

            return en
        end

        
    end


    # Computes the DDI energy of mat
    #
    # inputs: mat = spin matrix, ed = DDI constant, pbc = periodic boundary conditions
    # phiMatrices = array of matrices used to compute DDI
    #
    # outputs: float
    #
    function ddiEnergy(mat::Array{Float64,3}, ed::Float64, pbc::Float64, 
        phiMatrices::Array{Array{Float64,2},1})

        p, m, n     = size(mat)
        ddiEnArray  = Array{Float64}(undef,p,m,n)

        ddiEnArray  = dipoleDipole.FHD(mat,phiMatrices,pbc==1.0)
        
        tot = 0.0

        for i in eachindex(view(mat,:,:,:))
            tot += mat[i] * 0.5 * ddiEnArray[i]
        end

        return -ed * tot

    end


    # Computes the anisotropy energy of mat
    #
    # inputs: mat = spin matrix, AnisZ = PMA constant
    #
    # outputs: float
    #
    function pmaEnergy(mat::Array{Float64,3}, AnisZ::Float64)
        
        en = 0.0
        p, m, n = size(mat)

        for i in 1:m 
            for j in 1:n
                en -= 0.5 * AnisZ * (mat[3,i,j] * mat[3,i,j] - 1)
            end
        end

        return en
    end

end
