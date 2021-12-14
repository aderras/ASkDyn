#=
    This module contains functions used to build different spin configurations

    IMPROVEMENTS:
        + Throw warnings if a user selects something unrealistic. E.g. skyrmion
        of radius 100 on a (64,64) lattice
=#
module InitialCondition

    export buildinitial

    # Use the following to build the initial condition using the complex
    # to real plane transformation.
    @inline xcomp(z::Complex) = if isnan(real(z)/(1+(abs(z)^2)/4)) 0.
        else real(z)/(1+(abs(z)^2)/4) end
    @inline ycomp(z::Complex) = if isnan(imag(z)/(1+(abs(z)^2)/4)) 0.
        else imag(z)/(1+(abs(z)^2)/4) end
    @inline zcomp(z::Complex) = if isnan((-1+(abs(z)^2)/4)/(1+(abs(z)^2)/4)) 1.
        else -(-1+(abs(z)^2)/4)/(1+(abs(z)^2)/4) end

    @inline skyr(x,y,z0,l) = (x+y*1im-z0)/l
    @inline antiskyr(x,y,z0,l) = l/(x+y*1im-z0)
    @inline skyrAnti(x,y,z1,z2,l) = (x+y*1im-z1)*(x-y*1im-z2)/l^2

    # in: t = string defining the type of initial condition, l = size
    # of the topological object, z0 = center of the topological object,
    # Nx, Ny = dimensions of the initial condition array.
    #
    # out: mat = (3, Nx, Ny) array containing the initial condition
    function buildinitial(ic, mp)

        t = ic.type
        l = ic.r
        chi = ic.gamma
        xstart = ic.px
        ystart = ic.py

        Nx = mp.nx
        Ny = mp.ny

        mat = Array{Float64}(undef,3,Nx,Ny)
        xx = 0.
        yy = 0.

        dx = 1
        dy = 1

        vec = zeros(3)
        if t == "skyrmion"
            for i in eachindex(view(mat,1,1:Nx,1:Ny))

               skyrmioncomponent!(vec,i[1]-(xstart+1/2),i[2]-(ystart+1/2),l)

                mat[1,i] = -vec[1]
                mat[2,i] = -vec[2]
                mat[3,i] = vec[3]

            end
        elseif t == "domainWall"
            for i in eachindex(view(mat,1,1:Nx,1:Ny))

                # If easy-axis anisotropy constant is not zero, the domain wall
                # width is sqrt(J/Dz). If it is, set δ equal to user input
                if mp.dz!=0.0 δ = sqrt(mp.j/mp.dz) else δ=l end

                domainwallcomponent!(vec,δ,round(Int,i[2]-(ystart+1/2)))

                mat[1,i] = vec[1]
                mat[2,i] = vec[2]
                mat[3,i] = vec[3]

            end
        elseif t == "skyrmionAntiskyrmion"
            for i in eachindex(view(mat,1,1:Nx,1:Ny))

                xx = -xstart + dx*(i[1])
                yy = -ystart + dy*(i[2])

                z = skyrAnti(xx,yy,ic.px,-ic.px,l)
                mat[1,i] = xcomp(z)
                mat[2,i] = ycomp(z)
                mat[3,i] = zcomp(z)

            end
        elseif t == "antiskyrmion"
            for i in eachindex(view(mat,1,1:Nx,1:Ny))

                xx = -xstart + dx*(i[1] + 1)
                yy = -ystart + dy*(i[2] + 1)

                z = antiskyr(xx,yy,ic.px,l)
                mat[1,i] = xcomp(z)
                mat[2,i] = ycomp(z)
                mat[3,i] = zcomp(z)

            end
        elseif t == "AFMskyrmion"
            for i in 1:Nx, j in 1:Ny
                xx = -xstart + dx*(i + 1)
                yy = -ystart + dy*(j + 1)

                z = skyr(xx,yy,z0,l)

                if i%2==0 && j%2==1 ||
                    j%2==0 && i%2==1
                    mat[1,i,j] = -xcomp(z)
                    mat[2,i,j] = -ycomp(z)
                    mat[3,i,j] = -zcomp(z)
                else
                    mat[1,i,j] = xcomp(z)
                    mat[2,i,j] = ycomp(z)
                    mat[3,i,j] = zcomp(z)
                end
            end
        elseif t == "AFMskyrmionAntiskyrmion"
            for i in 1:Nx, j in 1:Ny
                xx = -xstart + dx*(i + 1)
                yy = -ystart + dy*(j + 1)

                z = skyrAnti(xx,yy,z0,-z0,l)

                if i%2==0 && j%2==1 ||
                    j%2==0 && i%2==1
                    mat[1,i,j] = -xcomp(z)
                    mat[2,i,j] = -ycomp(z)
                    mat[3,i,j] = -zcomp(z)
                else
                    mat[1,i,j] = xcomp(z)
                    mat[2,i,j] = ycomp(z)
                    mat[3,i,j] = zcomp(z)
                end
            end
        elseif t == "ferromagnet"
            for i in 1:Nx, j in 1:Ny
                mat[3,i,j] = 1.0
            end
        elseif t == "prec"
            for i in 1:Nx, j in 1:Ny
                mat[3,i,j] = 1.0
            end
            ang = pi/4
            mat[2,round(Int,Nx/2),round(Int,Ny/2)] = sin(ang)
            mat[3,round(Int,Nx/2),round(Int,Ny/2)] = cos(ang)
        else
            println("Initial condition selected is not an option in
            InitialCondition.jl.")
            exit()
        end

        return mat

    end

    # in: x,y = position on the plane of the centrosymmetic skyrmion,
    # λ = the radius, Q = topological charge, ϕ = chirality
    #
    # out: {x,y,z} vector representing skyrmion components at x,y
    function skyrmioncomponent!(val::Array{Float64,1}, x::Float64, y::Float64,
        λ::Float64)

        Q = 1.
        ϕ = pi/2

        rsquared = sz = sperp = ang = 0.

        rsquared = x^2+y^2
        sz = -(rsquared-λ^2)/(rsquared+λ^2)
        sperp = sqrt(1-sz^2)

        if y>0
            ang = acos(x/sqrt(rsquared))
        else
            ang = 2*pi-acos(x/sqrt(rsquared))
        end

        # for some reason allocating this way leads to smaller allocation time
        val[1] = sperp*cos(Q*ang+ϕ)
        val[2] = sperp*sin(Q*ang+ϕ)
        val[3] = sz

    end

    # The following function can be either plus or minus. I'm choosing plus for
    # now. This means θ(y=-∞) = 0 and θ(y = ∞)
    @inline θ(y, δ) = 2.0*atan(exp(y/δ))

    # The following defines a Bloch domain wall
    #   Mx(y) = M_0 sin(θ(y))
    #   My(y) = 0
    #   Mz(y) = M_0 cos(θ(y))
    # For normalization, M_0 = 1
    # The result is a 3x1 vector stored in the `val` input variable
    # δ is the domain wall width, y is the position in the plane
    function domainwallcomponent!(val::Array{Float64,1}, δ::Float64, y::Int64)

        val[1] = sin(θ(y, δ))
        val[2] = 0.0
        val[3] = cos(θ(y, δ))

    end

end
