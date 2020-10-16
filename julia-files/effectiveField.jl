# This module contains functions that calculate the effective field at a certain
# point of the spin lattice, and the function that builds the entire effective 
# field matrix.
#
# Useful functions: 
#   effectivefield! computes effective field at a point
#   getFullEffField! computes effective field of the whole lattice
# Note: effectivefieldd!() does not compute DDI field. This is because the DDI
# field must be computed for the entire lattice at once. (It's possible to 
# compute element-wise DDI by executing FFT, but it would be prohibitively slow.)

module effectiveField

    import dipoleDipole
    export effectivefield!, getFullEffField!, calculateDdiField

    # effectivefield! modifies the (3, 1) array to equal the effective field 
    # at some point nx,ny. Prealocating this way improves speed. 
    #
    # inputs: mat = (3, m, n) spin matrix, effField = (3, 1) arbitrary array 
    # used to store answer, nx,ny = the position of the spin component of interest,
    # params =  [J, H, DMI, PMA, ED, PBC, PIN, POS] material parameters, 
    # defectParams = defect parameters (optional argument)
    #
    # outputs: nothing
    #
    function effectivefield!( mat::Array{Float64,3}, effField::Array{Float64,1},
        nx::Int64, ny::Int64, params )

        matParams = params.mp
        defParams = params.defect
        pinParams = params.pin

        j,h,a,dz,ed,n,m,nz,pbc,vdd = 
            [ getfield( matParams, x ) for x in fieldnames( typeof( matParams) ) ]

        if defParams != []
            dParamsList = [ getfield( defParams, x ) for x in fieldnames( typeof(defParams) ) ]
        else 
            dParamsList = zeros(5)
        end

        if pinParams != []
            hPin, px, py = [ getfield( pinParams, x ) for x in fieldnames( typeof(pinParams) ) ]
        end

        # Compute zeeman field
        effField[1] = 0.
        effField[2] = 0.
        effField[3] = h

        # Compute exchange field
        exchangefield!( mat, effField, nx, ny, j, pbc==1.0, dParamsList )

        # Compute DMI and PMA if they're not set to zero
        if a != 0.0
            dmiblochfield!( mat, effField, nx, ny, a, pbc==1.0 )
        end
        if dz != 0.0
            anisotropyfield!( mat, effField, nx, ny, dz )
        end

	    # To pin the skyrmion to the center add field to the pinning point
        if pinParams != [] && nx == px && ny == py
            effField[3] = effField[3] + hPin
        end

        nothing

    end

    # This function modifies effField to include the exchange field of the spin array,
    # mat. Nearest neighbor spins are considered. 
    #
    # inputs: mat = (3, m, n), effField = (3, 1) array used to store answer,
    # nx, ny = positions in mat to calculate exchange field, J = exchange 
    # constant, pbc = periodic boundary conditions, defectParams = array 
    # containing details of the exchange defect if one is present. Last 
    # argument is optional and default computation is for an array without
    # a defect.
    # 
    # outputs: nothing
    #
    function exchangefield!( mat::Array{Float64,3}, effField::Array{Float64,1},
        nx::Int, ny::Int, J::Float64, pbc::Bool, 
        defectParams=zeros(5) )
        
        p, m, n = size(mat)

        if nx > 1 
            for k in 1:3 effField[k] = effField[k] + mat[k,nx-1,ny] end
        else
            if pbc
                for k in 1:3 effField[k] = effField[k] + mat[k,m,ny] end
            end
        end

        if nx < m
            for k in 1:3 effField[k] = effField[k] + mat[k,nx+1,ny] end
        else
            if pbc
                for k in 1:3 effField[k] = effField[k] + mat[k,1,ny] end
            end
        end

        if ny > 1
            for k in 1:3 effField[k] = effField[k] + mat[k,nx,ny-1] end
        else
            if pbc
                for k in 1:3 effField[k] = effField[k] + mat[k,nx,n] end
            end 
        end

        if ny < n
            for k in 1:3 effField[k] = effField[k] + mat[k,nx,ny+1] end
        else
            if pbc
                for k in 1:3 effField[k] = effField[k] + mat[k,nx,1] end
            end 
        end
        
        # If there is a defect and the point of interest is within the range
        # of the defect, change the effective exchange field.
        # 
        # defExists = 1.0 means a point defect which modifies the exchange 
        # interaction within some circle of radius dJ as a step function. 
        #
        # defExists = 2.0 means create a gaussian-type defect which affects
        # neighboring atoms according to an exponential function.
        defType,aJ,dJ,jx,jy = defectParams

        if (1.0 == defType) && sqrt((nx-jx)^2+(ny-jy)^2) <= dJ
 
            for i in 1:3 effField[i] = J*effField[i]*(1+aJ) end
 
        elseif (2.0 == defType) 
 
            for i in 1:3 
                effField[i] = J*effField[i]*(1+aJ * 
                        exp( -((nx-jx)^2 + (ny-jy)^2)/dJ^2 ) ) 
            end

        else

            for i in 1:3 effField[i] = J*effField[i]  end

        end
        
        nothing

    end

    # This function modifies effField to add the dmi contribution to the 
    # effective field at some point nx, ny.
    # 
    # inputs: mat = (3, m, n), effField = (3, 1) array used to store answer,
    # nx, ny = positions in mat to calculate exchange field, J = exchange 
    # constant, pbc = periodic boundary conditions, defectParams = array 
    # containing details of the dmi defect if one is present. Last 
    # argument is optional and default computation is for an array without
    # a defect.
    #
    # outputs: nothing  
    #
    function dmiblochfield!( mat::Array{Float64,3}, effField::Array{Float64,1}, 
        nx::Int, ny::Int, dmi::Float64, pbc::Bool, 
        defParams::Array{Float64,1} = zeros(5) )

        p,m,n = size(mat)

        if ny==1
            
            if pbc
                effField[1]-=dmi*mat[3,nx,n];
                effField[3]+=dmi*mat[1,nx,n];
            end

            effField[1]+=dmi*mat[3,nx,ny+1];
            effField[3]-=dmi*mat[1,nx,ny+1];

        elseif ny==n

            if pbc
                effField[1]+=dmi*mat[3,nx,1];
                effField[3]-=dmi*mat[1,nx,1];
            end

            effField[1]-=dmi*mat[3,nx,ny-1];
            effField[3]+=dmi*mat[1,nx,ny-1];

        else

            effField[1]+=dmi*(mat[3,nx,ny+1]-mat[3,nx,ny-1]);
            effField[3]+=dmi*(mat[1,nx,ny-1]-mat[1,nx,ny+1]);

        end

        if nx==1

            if pbc

                effField[2]+=dmi*mat[3,m,ny];
                effField[3]-=dmi*mat[2,m,ny];

            end

            effField[2]-=dmi*mat[3,nx+1,ny];
            effField[3]+=dmi*mat[2,nx+1,ny];

        elseif nx==m

            if pbc
            
                effField[2]-=dmi*mat[3,1,ny];
                effField[3]+=dmi*mat[2,1,ny];
            
            end
            
            effField[2]+=dmi*mat[3,nx-1,ny];
            effField[3]-=dmi*mat[2,nx-1,ny];
        
        else
        
            effField[2]+=dmi*(mat[3,nx-1,ny]-mat[3,nx+1,ny]);
            effField[3]+=dmi*(mat[2,nx+1,ny]-mat[2,nx-1,ny]);
        
        end

        # If there is a defect and the point of interest is within the range
        # of the defect, change the effective exchange field
        defExists,kx,ky,aDmi,dDmi = defParams
        if (1.0 == defExists) && sqrt((nx-kx)^2+(ny-ky)^2) <= dDmi  
            for i in 1:3 effField[i] = (1 + aDmi)*effField[i] end
        end

    end

    # anisotropyfield! modifies effField to include the PMA contribution at 
    # point (nx, ny)
    #
    # inputs: mat = spin matrix, effField = effective field at some nx, ny,
    # nx & ny are coordinates of in the spin matrix, pma = anisotropy constant
    #
    # outputs: nothing. This function changes effField to the effective field.
    #
    function anisotropyfield!( mat::Array{Float64,3}, effField::Array{Float64,1},
        nx::Int, ny::Int, pma::Float64 )

        effField[3] += pma*mat[3,nx,ny]

    end

  
    # Computes the effective field for the entire spin array. Returns (3, m, n)
    # array of the matrix 
    #
    # inputs: mat = (3, m, n) spin array, effField = (3,1) effective field at a
    # particular point, params = material parameters, defectParams = optional 
    # argument containing information about a defect, phiMatrices = optional 
    # argument of matrices used to compute DDI. 
    #
    # outputs: (3, m, n) array defining effective field at every point in the 
    # input array, mat
    #
    function getFullEffField!( mat::Array{Float64,3},  params) 
        
        matParams = params.mp
        j,h,a,dz,ed,n,m,nz,pbc,vdd = 
            [ getfield( matParams, x ) for x in fieldnames( typeof( matParams) ) ]


        Heff = Array{Float64}(undef, 3, m, n)
        effField = zeros(3)

        for j in 1:n
            for i in 1:m
                effectivefield!(mat, effField, i, j, params)

                for k in 1:3 Heff[k, i, j] = effField[k] end
            end 
        end

        # If DDI exists compute the effective field
        if ed != 0.0
            dipField    = Array{Float64}(undef, 3, m, n)
            dipField    = calculateDdiField(mat, ed, pbc==1.0, vdd)

            return Heff + dipField
        end

        return Heff
    end
  
    # Compute the DDI field of a spin array, mat.
    #
    # inputs: mat = (3, m, p) array of spins, ed = DDI constant, pbc =
    # boolean defining periodic boundary conditions, phi matrices = 
    # array of matrices used to compute DDI. 
    #
    # outputs: field = (3, m, n) array containing values of DDI field
    # at every point of mat. 
    #
    function calculateDdiField( mat::Array{Float64,3}, ed::Float64, 
        pbc::Bool, phiMatrices::Array{Array{Float64,2},1} )

        p, m, n = size(mat)
        field = Array{Float64}(undef, p, m, n)

        # Rewrite this code to only compute vddMatrices at the beginning
        # of computation. This is not necessary.
        field = ed * dipoleDipole.FHD(mat, phiMatrices, pbc)

        return field

    end
    
end
