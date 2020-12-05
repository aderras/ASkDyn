# This module contains the Runge-Kutta stepping function.
module rungeKutta

    using LLequation
    export rk4!

    # Calling this function once leads to 'n' iterations of the
    # fourth order runge-kutta solver. (See vars in main.jl for what
    # the values mean.)
    # Rather than return an array of the past n solutions, this modifies
    # the X0 array to constantly contain most recent solution.
    function rk4!( X0::Array{Float64,3}, f::Function, params, flag = false )

        p, m, n = size(X0)

        llgParams = params.llg

        tMax, hStep, nn, tol, lambda, T, nRuns, par =
            [ getfield( llgParams, x ) for x in fieldnames(typeof(llgParams)) ]

        # If running relaxation, use high damping
        if flag
            lambda = 1.0
        end

        # Initialize vars
        t0 = 0.0

        # Initialize arrays for the RK steps
        K1 = Array{Float64}(undef,p,m,n)
        K2 = Array{Float64}(undef,p,m,n)
        K3 = Array{Float64}(undef,p,m,n)
        K4 = Array{Float64}(undef,p,m,n)
        X = Array{Float64}(undef,p,m,n)

        X .= X0 # starting with specified

        for k in 1:nn

            # In this particular instance the equation were solving
            # has no t in RHS, so the t value is arbitrary.
            t = t0 + k*hStep

            f(t, X, K1, params, flag) # function modifies last argument to equal RK Coeff

            f(t + 0.5*hStep, ( X .+ (0.5*hStep).* K1 ) , K2, params, flag)

            f(t + 0.5*hStep, ( X .+ (0.5*hStep) .* K2 ), K3, params, flag)

            f(t + hStep, ( X .+ hStep*K3 ), K4, params, flag)

            X .= X .+ (1/6*hStep).*K1 .+ (1/3*hStep).*K2 .+
                (1/3*hStep).*K3 .+ (1/6*hStep).*K4
        end

        X0 .= X # rewrite input variable with the result

    end

end
