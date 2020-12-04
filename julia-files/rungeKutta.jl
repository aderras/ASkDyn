# This module contains the Runge-Kutta stepping function.
module rungeKutta

    using LLequation
    export rk4!

    # Calling this function once leads to 'n' iterations of the
    # fourth order runge-kutta solver. (See vars in main.jl for what
    # the values mean.)
    # Rather than return an array of the past n solutions, this modifies
    # the X0 array to constantly contain most recent solution.
    function rk4!( X0::Array{Float64,3}, f::Function, params, flag = true )

        p, m, n = size(X0)

        llgParams = params.llg

        tMax, hStep, nn, tol, lambda, T, nRuns, par =
            [ getfield( llgParams, x ) for x in fieldnames( typeof(llgParams) ) ]

        # If running relaxation, use high damping
        if flag
            lambda = 1.0
        end

        # Initialize arrays for the RK steps
        K1 = Array{Float64}(undef,p,m,n)
        K2 = Array{Float64}(undef,p,m,n)
        K3 = Array{Float64}(undef,p,m,n)
        K4 = Array{Float64}(undef,p,m,n)
        X = Array{Float64}(undef,p,m,n)

        temp = zeros(p,m,n)

        # Initialize vars
        t0 = 0.0

        X .= X0 # starting with specified

        for k in 1:nn

            # In this particular instance the equation were solving
            # has no t in RHS, so the t value is arbitrary.
            t = t0 + k*hStep

            f(t, X, K1, params, flag) # function modifies last argument to equal RK Coeff
            K1 .= hStep*K1

            f(t + 0.5*hStep, ( X .+ 0.5 .* K1 ) , K2, params, flag)
            K2 .= hStep*K2

            f(t + 0.5*hStep, ( X .+ 0.5 .* K2 ), K3, params, flag)
            K3 .= hStep*K3

            f(t + hStep, ( X .+ K3 ), K4, params, flag)
            K4 .= hStep*K4

            X .= X .+ (1/6).*K1 .+ (1/3).*K2 .+ (1/3).*K3 .+ (1/6).*K4
        end

        X0 .= X # rewrite input variable with the result

    end

end
