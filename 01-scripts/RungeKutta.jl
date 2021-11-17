#=
    This module contains the Runge-Kutta discretization scheme. It is used to
    propagate the spins in time.
=#
module RungeKutta

    using LLequation
    using LoopVectorization
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
    function elemSum!(dest::Array{Float64,3}, A::Array{Float64,3},
        B::Array{Float64,3}, a::Float64, b::Float64)
        p,m,n = size(A)
        @avx for i in 1:m, j in 1:m, k in 1:p
           dest[k,i,j] = a*A[k,i,j] + b*B[k,i,j]
        end
    end

    function elemSum!(dest::Array{Float64,3}, A::Array{Float64,3},
            B::Array{Float64,3}, C::Array{Float64,3},D::Array{Float64,3},
            E::Array{Float64,3}, a::Float64, b::Float64, c::Float64, d::Float64,
            e::Float64)
        p,m,n = size(A)
        @avx for i in 1:m, j in 1:m, k in 1:p
           dest[k,i,j] = a*A[k,i,j] + b*B[k,i,j] + c*C[k,i,j] +
               d*D[k,i,j] + e*E[k,i,j]
        end
    end

    function rk4!(t0, X::Array{Float64,3}, f::Function,
            params,
            fargs,
            K,
            relax=false)

        p, m, n = size(X)
        tmp = K[end] # last element is a dummy array

        t = t0
        hStep = params.cp.dt

        for k in 1:params.cp.nn
            t += k*params.cp.dt
            f(t, X, K[1], fargs, params, relax)
            elemSum!(tmp, X, K[1], 1.0, 0.5*hStep)
            f(t + 0.5*hStep, tmp , K[2], fargs, params, relax)
            elemSum!(tmp, X, K[2], 1.0, 0.5*hStep)
            f(t + 0.5*hStep, tmp, K[3], fargs, params, relax)
            elemSum!(tmp, X, K[3], 1.0, hStep)
            f(t + hStep, tmp, K[4], fargs, params, relax)
            elemSum!(tmp, X, K[1], K[2], K[3], K[4],
                1.0, (1/6*hStep), (1/3*hStep), (1/3*hStep), (1/6*hStep))
            X .= tmp
        end
    end
end
