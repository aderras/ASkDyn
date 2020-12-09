#=
This module contains functions used to compute the dipolar energy of the
spin lattice. We compute the interaction with convolution, which is why the
module requires FFTW.
=#
module dipoleDipole

    using FFTW, PaddedViews
    export FHD, convfft, buildVD, vdMatrices

    # We precompute matrices used to calculate the dipole-dipole interaction
    # with this function. It's essentially the distance part of the equation.
    #
    # in: Nx, Ny, Nz = dimensions of the system, pbc = periodic boundary
    # conditions == true
    #
    # out: Array of 4 (Nx,Ny) matrices
    function vdMatrices(Nx::Int,Ny::Int,Nz::Int,pbc::Bool=false)

        phixx = Array{Float32}(undef,2*Nx, 2*Ny)
        phiyy = Array{Float32}(undef,2*Nx, 2*Ny)
        phizz = Array{Float32}(undef,2*Nx, 2*Ny)
        phixy = Array{Float32}(undef,2*Nx, 2*Ny)

        phixx = buildVD("xx",Nx,Ny,Nz,pbc)
        phiyy = buildVD("yy",Nx,Ny,Nz,pbc)
        phizz = buildVD("zz",Nx,Ny,Nz,pbc)
        phixy = buildVD("xy",Nx,Ny,Nz,pbc)

        return [phixx,phiyy,phizz,phixy]

    end

    # Builds one of the 4 types of VD matrices based on input "stype"
    #
    # in: stype = the dimension of the matrix we're working on (Options: "xx",
    # "yy","zz","xy"), Nx,Ny,Nz are dimensions of the problem, PBC = periodic
    # boundary conditions.
    #
    # out: (Nx,Ny) matrix
    function buildVD( stype::String, Nx, Ny, Nz, pbc )

        if pbc == true

            VD = zeros(Nx,Ny)

            for i in -Nx+1:0, j in -Ny+1:0
                VD[j+Ny,i+Nx] = Phipbc(i,j,stype,Nx,Ny,Nz)
            end

        else

            if (Nx % 2 == Ny % 2 == 0)
                VD = zeros(2*Nx,2*Ny)
            else
                VD = zeros(2*Nx-1,2*Ny-1)
            end

            for i in -Nx+1:Nx-1, j in -Ny+1:Ny-1
                VD[j+Ny,i+Nx] = Phi(i,j,stype,Nz)
            end

        end

        return VD
    end

    # The following functions come from the equation for dipole-dipole interaction
    @inline phixx(a,b) = if a==b==0 0. else (2.0* b^2-a^2)/denominator(a,b,0.) end
    @inline phiyy(a,b) = if a==b==0 0. else (2.0* a^2-b^2)/denominator(a,b,0.) end
    @inline phizz(a,b) = if a==b==0 0. else (-1)/(sqrt((a^2 + b^2)^2 * (a^2 + b^2))) end
    @inline phixy(a,b) = if a==b==0 0. else (3.0* a * b)/denominator(a,b,0.) end
    @inline denominator(a::Float64,b::Float64,c::Float64) =  sqrt((a^2 + b^2 + c^2)^2*(a^2 + b^2 + c^2)^2*(a^2 + b^2 + c^2))

    @inline Fxx(nx,ny,Tx,Ty) = nx*ny/(4*pi*Tx*Ty*(nx^2)*sqrt(nx^2 + ny^2 ))
    @inline Fzz(nx,ny,Tx,Ty) = -nx*ny*(nx^2+ny^2)/(4*pi*Tx*Ty*(nx^2)*(ny^2)*sqrt(nx^2 + ny^2))
    @inline Fxy(nx,ny,Tx,Ty) = -1.0/(4*pi*Tx*Ty*sqrt(nx^2 + ny^2 ))

    # Provides an element of the VD matrix with periodic boundary conditions
    #
    # in: nx,ny = position of VD matrix, stype = type of matrix, NX,Ny, and Nz
    # are dimensions of the system
    #
    # out: variable 'sum,' which is element nx,ny of VD matrix of stype
    function Phipbc(nx::Int,ny::Int,stype::String, Nx, Ny, Nz)

        paramPbcX = round(Int,2000/float(Nx))
        paramPbcY = round(Int,2000/float(Ny))

        Xc = paramPbcX + 0.5
        Yc = paramPbcY + 0.5

        Tx = Nx
        Ty = Ny

        dim = Nz
        dimFloat = float(dim)
        sum = 0.

        nxfloat = nyfloat = Float64
        nxfloat = float(nx)
        nyfloat = float(ny)

        for ix in -paramPbcX:paramPbcX
            for iy in -paramPbcY:paramPbcY
               sum = sum + Phi(nx + ix*Nx, ny + iy*Ny, stype, Nz)
            end
        end

        if stype == "xy"
            sum += Fxy(nxfloat + Xc, nyfloat - Yc,Tx,Ty) + Fxy(nxfloat - Xc, nyfloat + Yc,Tx,Ty) -
                Fxy(nxfloat + Xc, nyfloat + Yc,Tx,Ty) - Fxy(nxfloat - Xc, nyfloat - Yc,Tx,Ty)
        elseif stype =="xx" || stype == "yy"
            sum += Fxx(nxfloat + Xc, nyfloat - Yc,Tx,Ty) + Fxx(nxfloat - Xc, nyfloat + Yc,Tx,Ty) -
                Fxx(nxfloat + Xc, nyfloat + Yc,Tx,Ty) - Fxx(nxfloat - Xc, nyfloat - Yc,Tx,Ty)
        elseif stype == "zz"
            sum += Fzz(nxfloat + Xc, nyfloat - Yc,Tx,Ty) + Fzz(nxfloat - Xc, nyfloat + Yc,Tx,Ty) -
                Fzz(nxfloat + Xc, nyfloat + Yc,Tx,Ty) - Fzz(nxfloat - Xc, nyfloat - Yc,Tx,Ty)
        end

	    return sum

    end

    # Returns the nx-th and ny-th matrix element of VD of type stype
    #
    # in: nx, ny are coordinates of interest, stype is the kind of VD matrix
    # being built, Nz is the dimension of the system in the z direction.
    #
    # out: variable 'sum,' which is the nx,ny element of VD
    function Phi(nx::Int, ny::Int, stype::String, Nz)

        dim = Nz
        dimFloat = float(dim)
        sum = 0.

        nxfloat = nyfloat = Float64
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
    function FHD(mat,phi,pbc::Bool=false)

        p,m,n = size(mat)

        matx = Array{Float64}(undef,m,n)
        maty = Array{Float64}(undef,m,n)
        matz = Array{Float64}(undef,m,n)
        matx,maty,matz = slicematrix(mat)

	return permutedims(cat(dims=3, convfft(matx,phi[1],pbc)+convfft(maty,phi[4],pbc),
        convfft(maty,phi[2],pbc)+convfft(matx,phi[4],pbc),convfft(matz,phi[3],pbc)),[3,1,2])

    end

    # Convolution of two arrays using FFT
    #
    # in: matrices a and b, where a has dimensions (Nx,Ny), pbc = periodic
    # boundary conditions. Need to know pbc because when true, we do not need
    # to pad the input spin array. When pbc = false, input spin array is padded
    # with zeros so that its dimensions match that of the vd matrix
    #
    # out = (Nx,Ny) matrix
    function convfft(a,b, pbc::Bool=false)

        ax,ay = size(a)
        bx,by = size(b)

        res     = Array{Float64}(undef,bx,by)

        if pbc == false
            apadded = PaddedView(0.0, a, (Base.OneTo(bx), Base.OneTo(by)))
            res = real(ifft(fft(apadded).*fft(b)))

            return res[bx-ax:bx-1,by-ay:by-1]

        else
            res = real(ifft(fft(a).*fft(b)))
            return res
        end

    end

    # This function splits a (3,Nx,Ny) array into a vector of 3 (Nx,Ny) arrays.
    # This was implemented because there was too much memory allocation when
    # using [1,:,:] operators in the FHD function. (Even with the @views macro.)
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
