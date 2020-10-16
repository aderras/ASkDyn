# This module contains the function that randomly rotates spins in a lattice to 
# model temperature. 
#
# 
module noiseRotation

    using Distributions, normalize, LinearAlgebra
    export rotateSpins!

    function rotateSpins!( mat::Array{Float64,3}, llgParams )

        p, m, n = size(mat)

        Φ = Array{Float64}(undef,p,m,n)
        ϕ = Array{Float64}(undef,3)
        s = Array{Float64}(undef,3)
        nn = Array{Float64}(undef,3)

        sdotn = 0.0

        tMax, hStep, nStep, tol, damping, T, nRuns, par = 
            [ getfield( llgParams, x ) for x in fieldnames( typeof(llgParams) ) ]

        tFree = hStep*nStep*2

        param = sqrt(2*damping*T*tFree)   # noise rotation parameter

        Φ = param*rand(Normal(),p,m,n)
        
        for i in 1:m
            for j in 1:n

                ϕnorm = 0.0
                for k in 1:3  
                    s[k] = mat[k,i,j] 
                    ϕ[k] = Φ[k,i,j]
                
                    ϕnorm += ϕ[k]^2
                end
                ϕnorm = sqrt(ϕnorm) + 1.0*10.0^-13

                for k in 1:3 nn[k] = ϕ[k]/ϕnorm end
                
                sinϕ = sin(ϕnorm)
                cosϕ = cos(ϕnorm)

                sdotn = dot(s,nn) * (1.0-cosϕ)

                # replace values in spin array with the new randomly rotated ones
                mat[1,i,j] = s[1]*cosϕ+sinϕ*(nn[2]*s[3]-nn[3]*s[2])+nn[1]*sdotn
                mat[2,i,j] = s[2]*cosϕ+sinϕ*(nn[3]*s[1]-nn[1]*s[3])+nn[2]*sdotn
                mat[3,i,j] = s[3]*cosϕ+sinϕ*(nn[1]*s[2]-nn[2]*s[1])+nn[3]*sdotn

            end
        end
    end

end
