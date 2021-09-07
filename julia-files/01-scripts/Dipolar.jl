#=

This module contains functions used to compute the dipolar energy of thin film.
We compute the interaction with convolution, which is why the module requires
FFTW.

The dipolar energy of the thin film computed by assuming that the spins are
approximately constant along the thickness of the film. This allows us to
calculate the interaction between columns of spins. For details of the
method, see doi:10.1103/PhysRevB.100.014432.


IMPROVEMENTS

    + slicematrix(M) deepcopies a (3,M,N) into 3 separate arrays, one for each
    of the 3 dimensions. There is probably a native way to do this, but the ones
    I tried were slower.

    + If there are periodic boundary conditions, the four vdmatrices are
    symmetric about the x and y axes. (You should quantify this.) This should
    mean that we can compute DDI with convolution of (M,N) matrices instead of
    (2M, 2N) matrices.

    + There is a parameter in the implementation of the periodic boundary
    condition which I set to 2000. (See phipbc() function.) You should
    generalize this to any system size according to the paper:
    https://doi.org/10.1016/j.commatsci.2010.04.024

=#
module Dipolar

    using FFTW, PaddedViews
    export fhd, convfft, buildvd, vdmatrices

    # We precompute matrices used to calculate the dipole-dipole interaction
    # with this function. It essentially contains information about the
    # distance of every spin to every other spin.
    #
    # in: Nx, Ny, Nz = dimensions of the system, pbc = periodic boundary
    #
    # out: Array of 4 (Nx,Ny) matrices
    function vdmatrices(Nx::Int, Ny::Int, Nz::Int, pbc)

        phixx = Array{Float64}(undef,2*Nx,2*Ny)
        phiyy = Array{Float64}(undef,2*Nx,2*Ny)
        phizz = Array{Float64}(undef,2*Nx,2*Ny)
        phixy = Array{Float64}(undef,2*Nx,2*Ny)

        phixx = buildvd("xx",Nx,Ny,Nz,pbc)
        phiyy = buildvd("yy",Nx,Ny,Nz,pbc)
        phizz = buildvd("zz",Nx,Ny,Nz,pbc)
        phixy = buildvd("xy",Nx,Ny,Nz,pbc)

        return [phiyy,phixx,phizz,phixy]

    end

    # Builds one of the 4 types of vd matrices based on input "stype"
    #
    # in: stype = the dimension of the matrix we're working on (Options: "xx",
    # "yy","zz","xy"), Nx,Ny,Nz are dimensions of the problem, PBC = periodic
    # boundary conditions.
    #
    # out: (Nx,Ny) matrix
    function buildvd(stype::String, Nx::Int, Ny::Int, Nz::Int, pbc)

        if pbc == 1.0

            if Nx%2 == Ny%2 == 0
                VD = zeros(2*Nx,2*Ny)
            else
                VD = zeros(2*Nx-1,2*Ny-1)
            end
            for i in -Nx+1:Nx-1, j in -Ny+1:Ny-1
                VD[j+Ny,i+Nx] = phipbc(i,j,stype,Nx,Ny,Nz)
            end

        else

            if Nx%2 == Ny%2 == 0
                VD = zeros(2*Nx,2*Ny)
            else
                VD = zeros(2*Nx-1,2*Ny-1)
            end

            for i in -Nx+1:Nx-1, j in -Ny+1:Ny-1
                VD[j+Ny,i+Nx] = phi(i,j,stype,Nz)
            end

        end

        return VD
    end

    # The following functions come from the equation for dipole-dipole
    # interaction
    @inline phixx(a,b) = if a==b==0 0. else
        (2.0* b^2-a^2)/denominator(a,b,0.) end
    @inline phiyy(a,b) = if a==b==0 0. else
        (2.0* a^2-b^2)/denominator(a,b,0.) end
    @inline phizz(a,b) = if a==b==0 0. else
        (-1)/(sqrt((a^2 + b^2)^2 * (a^2 + b^2))) end
    @inline phixy(a,b) = if a==b==0 0. else
        (3.0* a * b)/denominator(a,b,0.) end
    @inline denominator(a::Float64,b::Float64,c::Float64) =
        sqrt((a^2 + b^2 + c^2)^2*(a^2 + b^2 + c^2)^2*(a^2 + b^2 + c^2))

    @inline Fxx(nx,ny,Tx,Ty) = nx*ny/(4*pi*Tx*Ty*(nx^2)*sqrt(nx^2 + ny^2 ))
    @inline Fzz(nx,ny,Tx,Ty) = -nx*ny*(nx^2+ny^2)/
        (4*pi*Tx*Ty*(nx^2)*(ny^2)*sqrt(nx^2 + ny^2))
    @inline Fxy(nx,ny,Tx,Ty) = -1.0/(4*pi*Tx*Ty*sqrt(nx^2 + ny^2 ))

    # Provides an element of the VD matrix with periodic boundary conditions
    #
    # in: nx,ny = position of VD matrix, stype = type of matrix, Nx, Ny, and Nz
    # are dimensions of the system
    #
    # out: variable 'sum' which is element nx, ny of VD matrix of stype
    function phipbc(nx::Int,ny::Int,stype::String, Nx, Ny, Nz)

        # In implementing periodic boundary conditions, there is a parameter
        # that determines how far in the sum you go. Here, I choose the
        # value of 2000 which worked well for the system size of 128.
        paramPbcX = round(Int,2000/float(Nx))
        paramPbcY = round(Int,2000/float(Ny))

        Xc = paramPbcX + 0.5
        Yc = paramPbcY + 0.5

        Tx = Nx
        Ty = Ny

        dim = Nz
        dimFloat = float(dim)
        sum = 0.

        nxfloat = nyfloat = AbstractFloat
        nxfloat = float(nx)
        nyfloat = float(ny)

        for ix in -paramPbcX:paramPbcX
            for iy in -paramPbcY:paramPbcY
               sum = sum + phi(nx + ix*Nx, ny + iy*Ny, stype, Nz)
            end
        end

        if stype == "xy"
            sum +=  Fxy(nxfloat + Xc, nyfloat - Yc,Tx,Ty) +
                    Fxy(nxfloat - Xc, nyfloat + Yc,Tx,Ty) -
                    Fxy(nxfloat + Xc, nyfloat + Yc,Tx,Ty) -
                    Fxy(nxfloat - Xc, nyfloat - Yc,Tx,Ty)
        elseif stype =="xx" || stype == "yy"
            sum +=  Fxx(nxfloat + Xc, nyfloat - Yc,Tx,Ty) +
                    Fxx(nxfloat - Xc, nyfloat + Yc,Tx,Ty) -
                    Fxx(nxfloat + Xc, nyfloat + Yc,Tx,Ty) -
                    Fxx(nxfloat - Xc, nyfloat - Yc,Tx,Ty)
        elseif stype == "zz"
            sum +=  Fzz(nxfloat + Xc, nyfloat - Yc,Tx,Ty) +
                    Fzz(nxfloat - Xc, nyfloat + Yc,Tx,Ty) -
                    Fzz(nxfloat + Xc, nyfloat + Yc,Tx,Ty) -
                    Fzz(nxfloat - Xc, nyfloat - Yc,Tx,Ty)
        end

	    return sum

    end

    # Returns the nx-th and ny-th matrix element of VD of type stype
    #
    # in: nx, ny are coordinates of interest, stype is the kind of VD matrix
    # being built, dim is the dimension of the system in the z direction.
    #
    # out: variable 'sum,' which is the nx,ny element of VD
    function phi(nx::Int, ny::Int, stype::String, dim)

        dimFloat = float(dim)
        sum = 0.

        nxfloat = nyfloat = AbstractFloat
        nxfloat = float(nx)
        nyfloat = float(ny)

        #implementing the so called smart formulas
        if stype == "xy"
            sum += phixy(nxfloat,nyfloat) + 2/dim*sumphixy(nxfloat,nyfloat,dim)
        elseif stype == "xx"
            sum += phixx(nxfloat,nyfloat) + 2/dim*sumphixx(nxfloat,nyfloat,dim)
        elseif stype == "yy"
            sum += phiyy(nxfloat,nyfloat) + 2/dim*sumphiyy(nxfloat,nyfloat,dim)
        elseif stype == "zz"
            sum += phizz(nxfloat,nyfloat) + 2/dim*sumphizz(nxfloat,nyfloat,dim)
        end

        return sum
    end


    function sumphiyy(a,b,N)
        tot = 0.
        for i in 1:N-1
            tot += (N-i)*(2. *a^2-b^2-i^2)/denominator(a,b,float(i))
        end
        return tot
    end

    function sumphixx(a,b,N)
        tot = 0.
        for i in 1:N-1
            tot += (N-i)*(2. *b^2-a^2-i^2)/denominator(a,b,float(i))
        end
        return tot
    end

    function sumphizz(a,b,N)
        tot = 0.
        for i in 1:N-1
            tot += (N-i)*(2. *i^2-b^2-a^2)/denominator(a,b,float(i))
        end
        return tot
    end

    function sumphixy(a,b,N)
        tot = 0.
        for i in 1:N-1
            tot += (N-i)*(3. *a*b)/denominator(a,b,float(i))
        end
        return tot
    end


    # Computes the dipolar field (without multiplying by dipolar constant)
    # by convolving and summing approproate dimensions of the spin matrix,
    # mat, and the VD matrices, called phi here.
    #
    # in: mat = spin field, phi = VD matrices, pbc = periodic boundary cond.
    #
    # out: (3,Nx,Ny) array representing dipolar interaction of spin field mat
    # at every point.
    function fhd(mat, phi, pbc=false)

        p,m,n = size(mat)

        matx = Array{Float64}(undef,m,n)
        maty = Array{Float64}(undef,m,n)
        matz = Array{Float64}(undef,m,n)
        matx,maty,matz = slicematrix(mat)

        return permutedims(cat(dims=3,
            convfft(matx,phi[2],pbc)+convfft(maty,phi[4],pbc),
            convfft(maty,phi[1],pbc)+convfft(matx,phi[4],pbc),
            convfft(matz,phi[3],pbc)),[3,1,2])


    end

    # Convolution of two arrays using FFT
    #
    # in: matrices a and b, where a has dimensions (Nx,Ny), pbc = periodic
    # boundary conditions. Need to know pbc because when true, we do not need
    # to pad the input spin array. When pbc = false, input spin array is padded
    # with zeros so that its dimensions match that of the vd matrix
    #
    # out = (Nx,Ny) matrix
    function convfft(a, b, pbc=false)

        ax,ay = size(a)
        bx,by = size(b)

        res     = Array{Float64}(undef,bx,by)

        apadded = PaddedView(0.0, a, (Base.OneTo(bx), Base.OneTo(by)))
        res = real(ifft(fft(apadded).*fft(b)))

        return res[bx-ax:bx-1,by-ay:by-1]

    end

    # This function splits a (3,Nx,Ny) array into a vector of 3 (Nx,Ny) arrays.
    # This was implemented because there was too much memory allocation when
    # using [1,:,:] operators in the fhd function. (Even with the @views macro.)
    # There is probably a better way to do this wih native julia functions, but
    # this works for now.
    #
    # in: (3, Nx, Ny) array
    #
    # out: [Bx,By,Bz] where each B is an (Nx,Ny) array
    function slicematrix(A::Array{Float64,3})

        p,m,n = size(A)

        B1 = Array{Float64}(undef,m,n)
        B2 = Array{Float64}(undef,m,n)
        B3 = Array{Float64}(undef,m,n)

        for i in 1:m
            for j in 1:n
                B1[i,j] = A[1,i,j]
                B2[i,j] = A[2,i,j]
                B3[i,j] = A[3,i,j]
            end
        end

        return [B1,B2,B3]
    end


end
