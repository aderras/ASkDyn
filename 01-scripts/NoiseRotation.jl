#=
    This module contains the function that randomly rotates spins in a lattice
    to model temperature.
=#
module NoiseRotation

    using Normalize, LinearAlgebra
    export noisyrotate!

    function noisyrotate!(mat::Array{Float64,3}, tmps::Vector{Vector{Float64}},
        cpValues, damping, T)

        p, m, n = size(mat)

        ϕ = tmps[1]
        s = tmps[2]

        sdotn = 0.0

        hStep, nStep = cpValues
        tFree = hStep*nStep*2.0

        param = sqrt(2.0*damping*T*tFree)   # noise rotation parameter

        for i in 1:m, j in 1:n

            ϕnorm = 0.0
            for k in 1:3
                s[k] = mat[k,i,j]
                ϕ[k] = param*rand()

                ϕnorm += ϕ[k]^2
            end
            ϕnorm = sqrt(ϕnorm) + 10.0^(-13)

            for k in 1:3 ϕ[k] = ϕ[k]/ϕnorm end

            sinϕ = sin(ϕnorm)
            cosϕ = cos(ϕnorm)

            sdotn = 0.0
            for i in 1:3 sdotn+=s[i]*ϕ[i]*(1-cosϕ) end

            # replace values in spin array with randomly rotated ones
            mat[1,i,j] = s[1]*cosϕ+sinϕ*(ϕ[2]*s[3]-ϕ[3]*s[2])+ϕ[1]*sdotn
            mat[2,i,j] = s[2]*cosϕ+sinϕ*(ϕ[3]*s[1]-ϕ[1]*s[3])+ϕ[2]*sdotn
            mat[3,i,j] = s[3]*cosϕ+sinϕ*(ϕ[1]*s[2]-ϕ[2]*s[1])+ϕ[3]*sdotn

        end
    end

end
