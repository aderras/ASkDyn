#=
    This module contains the Runge-Kutta discretization scheme. It is used to
    propagate the spins in time.
=#
module RungeKutta

    using LLequation
    using LoopVectorization
    export rk4!, rk5!

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
        # Used to have @avx

        for i in 1:m, j in 1:m, k in 1:p
           dest[k,i,j] = a*A[k,i,j] + b*B[k,i,j]
        end
    end

    function elemSum!(dest::Array{Float64,3}, A::Array{Float64,3},
            B::Array{Float64,3}, C::Array{Float64,3}, a::Float64, b::Float64,
            c::Float64)
        p,m,n = size(A)
        # Used to have @avx
        for i in 1:m, j in 1:m, k in 1:p
           dest[k,i,j] = a*A[k,i,j] + b*B[k,i,j] + c*C[k,i,j]
        end
    end

    function elemSum!(dest::Array{Float64,3}, A::Array{Float64,3},
            B::Array{Float64,3}, C::Array{Float64,3},D::Array{Float64,3},
            E::Array{Float64,3}, a::Float64, b::Float64, c::Float64, d::Float64,
            e::Float64)
        p,m,n = size(A)
        # Used to have @avx
        for i in 1:m, j in 1:m, k in 1:p
           dest[k,i,j] = a*A[k,i,j] + b*B[k,i,j] + c*C[k,i,j] +
               d*D[k,i,j] + e*E[k,i,j]
        end
    end

    function elemSum!(dest::Array{Float64,3}, A::Array{Float64,3},
            B::Array{Float64,3}, C::Array{Float64,3},D::Array{Float64,3},
            E::Array{Float64,3}, F::Array{Float64,3}, a::Float64, b::Float64, c::Float64, d::Float64,
            e::Float64, f::Float64)
        p,m,n = size(A)
        # Used to have @avx
        for i in 1:m, j in 1:m, k in 1:p
           dest[k,i,j] = a*A[k,i,j] + b*B[k,i,j] + c*C[k,i,j] +
               d*D[k,i,j] + e*E[k,i,j] + f*F[k,i,j]
        end
    end

    function rk4!(X::Array{Float64,3}, t0::Float64, f::Function,
        fargs,
        matParams::Array{Any,1},
        cpParams::Array{Float64,1},
        params,
        K,
        relax=false,
        damping=0.0)

        p,m,n = size(X)
        tmp = K[end] # last element is a dummy array

        t = t0
        hStep, nn = cpParams

        for k in 1:nn

            t += k*hStep
            f(K[1], t, X, fargs, matParams, params, relax, damping)

            elemSum!(tmp, X, K[1], 1.0, 0.5*hStep)
            f(K[2], t+0.5*hStep, tmp, fargs, matParams, params, relax, damping)

            elemSum!(tmp, X, K[2], 1.0, 0.5*hStep)
            f(K[3], t+0.5*hStep, tmp, fargs, matParams, params, relax, damping)

            elemSum!(tmp, X, K[3], 1.0, hStep)
            f(K[4], t+hStep, tmp, fargs, matParams, params, relax, damping)

            elemSum!(tmp, X, K[1], K[2], K[3], K[4],
                1.0, (1/6*hStep), (1/3*hStep), (1/3*hStep), (1/6*hStep))
            X .= tmp

        end
    end

    function rk5!(X::Array{Float64,3}, t0::Float64, f::Function,
        fargs,
        matParams::Array{Any,1},
        cpParams::Array{Float64,1},
        params,
        K,
        relax=false,
        damping=0.0)

        p,m,n = size(X)
        tmp = K[end] # last element is a dummy array

        t = t0
        hStep, nn = cpParams

        for k in 1:nn

            t += k*hStep
            f(K[1], t, X, fargs, matParams, params, relax, damping)

            elemSum!(tmp, X, K[1], 1.0, (1/4)*hStep)
            f(K[2], t+0.25*hStep, tmp, fargs, matParams, params, relax, damping)

            elemSum!(tmp, X, K[1], K[2], 1.0, (1/8)*hStep, (1/8)*hStep)
            f(K[3], t+0.25*hStep, tmp, fargs, matParams, params, relax, damping)

            elemSum!(tmp, X, K[2], K[3], 1.0, -(1/2)*hStep, hStep)
            f(K[4], t+0.5*hStep, tmp, fargs, matParams, params, relax, damping)

            elemSum!(tmp, X, K[1], K[4], 1.0, (3/16)*hStep, (9/16)*hStep)
            f(K[5], t+0.75*hStep, tmp, fargs, matParams, params, relax, damping)

            elemSum!(tmp, X, K[1], K[2], K[3], K[4], K[5], 1.0, (-3/7)*hStep,
                (2/7)*hStep, (12/7)*hStep, (-12/7)*hStep, (8/7)*hStep)
            f(K[6], t+0.75*hStep, tmp, fargs, matParams, params, relax, damping)

            elemSum!(tmp, X, K[1], K[3], K[4], K[5], K[6],
                1.0, (7/90)*hStep, (32/90)*hStep, (12/90)*hStep, (32/90)*hStep,
                (7/90)*hStep)
            X .= tmp

        end
    end
end
