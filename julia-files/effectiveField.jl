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

        j,h,a,dz,ed,n,m,nz,pbc,vdd =
            [ getfield( matParams, x ) for x in fieldnames( typeof( matParams) ) ]

        # Compute zeeman field
        effField[1] = 0.
        effField[2] = 0.
        effField[3] = h

        # # # Compute exchange field
        # if params.defect.t==0
        #     exchangefield!( mat, effField, nx, ny, j, pbc==1.0 )
        # elseif params.defect.t==2
        #     # run exchange interaction with defect
        #     exFieldGaussianDefect!( mat, effField, nx, ny, params  )
        # end

        # Compute DMI and PMA if they're not set to zero
        if a != 0.0
            dmiblochfield!( mat, effField, nx, ny, a, pbc==1.0 )
        end
        # if dz != 0.0
        #     anisotropyfield!( mat, effField, nx, ny, dz )
        # end

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
        nx::Int, ny::Int, J::Float64, pbc::Bool )

        p, m, n = size(mat)

        if nx > 1 && nx < m
            for k in 1:3
                effField[k] = effField[k] + J*(mat[k,nx-1,ny] + mat[k,nx+1,ny])
            end
        elseif nx==1
            if pbc
                for k in 1:3 effField[k] = effField[k] + J * mat[k,m,ny] end
            end
            for k in 1:3 effField[k] = effField[k] + J * mat[k,nx+1,ny] end
        elseif nx == m
            if pbc
                for k in 1:3 effField[k] = effField[k] + J * mat[k,1,ny] end
            end
            for k in 1:3 effField[k] = effField[k] + J * mat[k,nx-1,ny] end
        end

        if ny > 1 && ny < n
            for k in 1:3
                effField[k] = effField[k] + J*(mat[k,nx,ny-1] + mat[k,nx,ny+1])
            end
        elseif ny==1
            if pbc
                for k in 1:3 effField[k] = effField[k] + J * mat[k,nx,n] end
            end
            for k in 1:3 effField[k] = effField[k] + J * mat[k,nx,ny+1] end
        elseif ny==n
            if pbc
                for k in 1:3 effField[k] = effField[k] + J * mat[k,nx,1] end
            end
            for k in 1:3 effField[k] = effField[k] + J * mat[k,nx,ny-1] end
        end
    end

    function exFieldGaussianDefect!(mat, effField, nx, ny, params)
        # If there is a defect and the point of interest is within the range
        # of the defect, change the effective exchange field.
        #
        # defExists = 1.0 means a point defect which modifies the exchange
        # interaction within some circle of radius dJ as a step function.
        #
        # defExists = 2.0 means create a gaussian-type defect which affects
        # neighboring atoms according to an exponential function.
        defType,aJ,dJ,jx,jy = [ getfield( params.defect, x )
            for x in fieldnames( typeof( params.defect ) ) ]
        J = params.mp.j
        pbc = params.mp.pbc==1.0

        @inline Jmod(i,j) = J*(1 - aJ*exp( -(( i - jx)^2 + (j -jy)^2)/dJ^2 ) )

        p, m, n = size(mat)

        if nx > 1 && nx < m
            for k in 1:3
                effField[k] = effField[k] +
                    Jmod( nx-1/2, ny ) * (mat[k,nx-1,ny] + mat[k,nx+1,ny])
            end
        elseif nx == 1
            if pbc
                for k in 1:3
                    effField[k] = effField[k] + Jmod( nx-1/2, ny ) * mat[k,m,ny]
                end
            end

            for k in 1:3
                effField[k] = effField[k] + Jmod( nx+1/2, ny ) * mat[k,nx+1,ny]
            end
        elseif nx == m
            if pbc
                for k in 1:3
                    effField[k] = effField[k] + Jmod( nx+1/2, ny ) * mat[k,1,ny]
                end
            end

            for k in 1:3
                effField[k] = effField[k] + Jmod( nx+1/2, ny ) * mat[k,nx-1,ny]
            end
        end


        if ny > 1 && ny < n
            for k in 1:3
                effField[k] = effField[k] +
                    Jmod( nx, ny-1/2 ) * ( mat[k,nx,ny-1] + mat[k,nx,ny+1])
            end
        elseif ny == 1
            if pbc
                for k in 1:3
                    effField[k] = effField[k] + Jmod( nx, ny-1/2 ) * mat[k,nx,n]
                end
            end
            for k in 1:3
                effField[k] = effField[k] + Jmod( nx, ny+1/2 ) * mat[k,nx,ny+1]
            end

        elseif ny == n
            if pbc
                for k in 1:3
                    effField[k] = effField[k] + Jmod( nx, ny+1/2 ) * mat[k,nx,1]
                end
            end
            for k in 1:3
                effField[k] = effField[k] + Jmod( nx, ny+1/2 ) * mat[k,nx,ny-1]
            end
        end
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
        if (2.0 == defExists) && sqrt((nx-kx)^2+(ny-ky)^2) <= dDmi
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
            [ getfield( matParams, x ) for x in fieldnames( typeof(matParams) ) ]

        Heff = zeros(3, m, n)
        effField = zeros(3)

        #Exchange effective field
        calcExchangeField!( mat, params, Heff )

        if params.defect.t == 2
            Jmat = zeros(3,m,n)
            Jmat = gaussianDefectExchange( params )
            Heff = Heff .* Jmat
        end

        calcZeemanField!( mat, params, Heff )

        # DMI effective field
        if a != 0.0
            calcDmiField!( mat, params, Heff )
        end
        if dz != 0.0
            calcPmaField!( mat, params, Heff )
        end

        # If DDI exists compute the effective field
        if ed != 0.0
            dipField    = Array{Float64}(undef, 3, m, n)
            dipField    = calcDdiField(mat, ed, pbc==1.0, vdd)

            return Heff + dipField
        end

        # If there is a pinning field, add the field at the point
        # To pin the skyrmion to the center add field to the pinning point
        hPin, px, py = [ getfield( params.pin, x )
            for x in fieldnames( typeof( params.pin ) ) ]
        if hPin != 0.0
            Heff[3,px,py] = Heff[3,px,py] + hPin
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
    function calcDdiField( mat::Array{Float64,3}, ed::Float64,
        pbc::Bool, phiMatrices::Array{Array{Float64,2},1} )

        p, m, n = size(mat)
        field = Array{Float64}(undef, p, m, n)

        # Rewrite this code to only compute vddMatrices at the beginning
        # of computation. This is not necessary.
        field = ed * dipoleDipole.FHD(mat, phiMatrices, pbc)

        return field
    end

    function calcExchangeField!( mat, params, Heff )

        p, m, n = size(mat)

        J = params.mp.j
        pbc = params.mp.pbc == 1.0

        for nx in 1:m, ny in 1:n

            if pbc
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

                for k in 1:3
                    Heff[k,nx,ny] = Heff[k,nx,ny] + J *
                        ( mat[k,nxPrev,ny] + mat[k,nxNext,ny] +
                        mat[k,nx,nyPrev] + mat[k,nx,nyNext] )
                end

            else

                if nx > 1
                    for k in 1:3
                        Heff[k,nx,ny] = Heff[k,nx,ny] + J * mat[k,nx-1,ny]
                    end
                end
                if ny > 1
                    for k in 1:3
                        Heff[k,nx,ny] = Heff[k,nx,ny] + J * mat[k,nx,ny-1]
                    end
                end
                if nx < m
                    for k in 1:3
                        Heff[k,nx,ny] = Heff[k,nx,ny] + J * mat[k,nx+1,ny]
                    end
                end
                if ny < n
                    for k in 1:3
                        Heff[k,nx,ny] = Heff[k,nx,ny] + J * mat[k,nx,ny+1]
                    end
                end
            end
        end

    end

    function gaussianDefectExchange( params )

        nx = params.mp.nx
        ny = params.mp.ny

        defType,aJ,dJ,jx,jy = [ getfield( params.defect, x )
            for x in fieldnames( typeof( params.defect ) ) ]
        J = params.mp.j

        exMat = zeros(3,nx,ny)

        for i in 1:nx, j in 1:ny

            for k in 1:3
                exMat[k, i,j] = (1 - aJ*exp( -(( i - jx)^2 + (j -jy)^2)/dJ^2 ) )
            end
        end

    end

    function calcDmiField!(mat, params, Heff)

        p,m,n = size(mat)

        dmi = params.mp.a
        pbc = params.mp.pbc == 1.0

        for nx in 1:m, ny in 1:n
            if ny==1

                if pbc
                    Heff[1,nx,ny] = Heff[1,nx,ny] - dmi*mat[3,nx,n]
                    Heff[3,nx,ny] = Heff[3,nx,ny] + dmi*mat[1,nx,n]
                end

                Heff[1,nx,ny] = Heff[1,nx,ny] + dmi*mat[3,nx,ny+1]
                Heff[3,nx,ny] = Heff[3,nx,ny] - dmi*mat[1,nx,ny+1]

            elseif ny==n

                if pbc
                    Heff[1,nx,ny] = Heff[1,nx,ny] + dmi*mat[3,nx,1]
                    Heff[3,nx,ny] = Heff[3,nx,ny] - dmi*mat[1,nx,1]
                end

                Heff[1,nx,ny] = Heff[1,nx,ny] - dmi*mat[3,nx,ny-1]
                Heff[3,nx,ny] = Heff[3,nx,ny] + dmi*mat[1,nx,ny-1]

            else

                Heff[1,nx,ny] = Heff[1,nx,ny]+dmi*(mat[3,nx,ny+1]-mat[3,nx,ny-1])
                Heff[3,nx,ny] = Heff[3,nx,ny]+dmi*(mat[1,nx,ny-1]-mat[1,nx,ny+1])
            end


            if nx==1

                if pbc

                    Heff[2,nx,ny] = Heff[2,nx,ny] + dmi*mat[3,m,ny]
                    Heff[3,nx,ny] = Heff[3,nx,ny] - dmi*mat[2,m,ny]

                end

                Heff[2,nx,ny] = Heff[2,nx,ny] - dmi*mat[3,nx+1,ny]
                Heff[3,nx,ny] = Heff[3,nx,ny] + dmi*mat[2,nx+1,ny]

            elseif nx==m

                if pbc

                    Heff[2,nx,ny] = Heff[2,nx,ny] - dmi*mat[3,1,ny]
                    Heff[3,nx,ny] = Heff[3,nx,ny] + dmi*mat[2,1,ny]

                end

                Heff[2,nx,ny] = Heff[2,nx,ny] + dmi*mat[3,nx-1,ny]
                Heff[3,nx,ny] = Heff[3,nx,ny] - dmi*mat[2,nx-1,ny]

            else

                Heff[2,nx,ny]= Heff[2,nx,ny]+dmi*(mat[3,nx-1,ny]-mat[3,nx+1,ny])
                Heff[3,nx,ny]= Heff[3,nx,ny]+dmi*(mat[2,nx+1,ny]-mat[2,nx-1,ny])

            end

        end

    end

    function calcZeemanField!(mat, params, Heff)


        p, m, n = size(mat)

        H = params.mp.h

        for i in 1:m, j in 1:n

            Heff[3,i,j] = Heff[3,i,j] + H
        end

    end

    function calcPmaField!(mat, params, Heff)

        p, m, n = size(mat)

        Dz = params.mp.dz

        for i in 1:m, j in 1:n

            Heff[3,i,j] = Heff[3,i,j] + Dz * mat[3, i, j]
        end

    end


end
