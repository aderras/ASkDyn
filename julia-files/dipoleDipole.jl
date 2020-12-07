module dipoleDipole

    using FFTW, PaddedViews
    export FHD, convfft, buildVD, vddMatrices

    function vddMatrices(Nx::Int,Ny::Int,Nz::Int,pbc::Bool=false)

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

    @inline phixx(a,b) = if a==b==0 0. else (2.0* b^2-a^2)/denominator(a,b,0.) end
    @inline phiyy(a,b) = if a==b==0 0. else (2.0* a^2-b^2)/denominator(a,b,0.) end
    @inline phizz(a,b) = if a==b==0 0. else (-1)/(sqrt((a^2 + b^2)^2 * (a^2 + b^2))) end
    @inline phixy(a,b) = if a==b==0 0. else (3.0* a * b)/denominator(a,b,0.) end
    @inline denominator(a::Float64,b::Float64,c::Float64) =  sqrt((a^2 + b^2 + c^2)^2*(a^2 + b^2 + c^2)^2*(a^2 + b^2 + c^2))

    @inline Fxx(nx,ny,Tx,Ty) = nx*ny/(4*pi*Tx*Ty*(nx^2)*sqrt(nx^2 + ny^2 ))
    @inline Fzz(nx,ny,Tx,Ty) = -nx*ny*(nx^2+ny^2)/(4*pi*Tx*Ty*(nx^2)*(ny^2)*sqrt(nx^2 + ny^2))
    @inline Fxy(nx,ny,Tx,Ty) = -1.0/(4*pi*Tx*Ty*sqrt(nx^2 + ny^2 ))

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

    # Returns the xth and yth matrix element of phi
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
    function convfft(a,b, pbc::Bool=false)#,pfor,planifft)

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
