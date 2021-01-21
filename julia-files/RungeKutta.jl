#=
    This module contains the Runge-Kutta discretization scheme. It is used to
    propagate the spins in time.
=#
module RungeKutta

    using LLequation
    export rk4!

    # Calling this function once leads to 'n' iterations of the
    # fourth order runge-kutta solver. (See vars in main.jl for what
    # the values mean.)
    # Rather than return an array of the past n solutions, this modifies
    # the X0 array to constantly contain most recent solution.
    #
    # in: X0 = spin array, f = time evolution function obeyed by the vectors in
    # X0, params = struct of evaluation parameters, relax = relaxation (T/F)
    #
    # out = nothing
    function rk4!(X0::Array{AbstractFloat,3}, f::Function, params, relax=false)

        p, m, n = size(X0)

        llgParams = params.llg

        tMax, hStep, nn, tol, lambda, T, nRuns, par =
            [getfield(llgParams, x) for x in fieldnames(typeof(llgParams))]

        # Initialize vars
        t0 = 0.0

        # Initialize arrays for the RK steps
        K1 = Array{AbstractFloat}(undef,p,m,n)
        K2 = Array{AbstractFloat}(undef,p,m,n)
        K3 = Array{AbstractFloat}(undef,p,m,n)
        K4 = Array{AbstractFloat}(undef,p,m,n)
        X = Array{AbstractFloat}(undef,p,m,n)

        X .= X0 # starting with specified

        for k in 1:nn

            # In this implentation the equation were solving
            # has no t in RHS, so the t value is arbitrary.
            t = t0 + k*hStep

            f(t, X, K1, params, relax)

            f(t + 0.5*hStep, (X .+ (0.5*hStep).* K1) , K2, params, relax)

            f(t + 0.5*hStep, (X .+ (0.5*hStep) .* K2), K3, params, relax)

            f(t + hStep, (X .+ hStep*K3), K4, params, relax)

            X .= X .+ (1/6*hStep).*K1 .+ (1/3*hStep).*K2 .+
                (1/3*hStep).*K3 .+ (1/6*hStep).*K4
        end

        X0 .= X # rewrite input variable with the result

    end

end
