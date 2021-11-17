module BoundaryConditions
    using StaticArrays

    extrap0 = SVector{1}([0.0]) # Need this as a placeholder

    # First order backward euler, linear and quadratic interp
    extrap1 = SVector{2}([2.0, -1.0])
    extrap2 = SVector{3}([3.0, -3.0, 1.0])

    # 4th order backward difference, linear and quadratic interp
    extrap3 = SVector{5}([37.0, -48.0, 36.0, -16.0, 3.0]/12.0)
    extrap4 = SVector{6}([41.0, -101.0, 125.0, -86.0, 32.0, -5.0]/6.0)

    # 4th order central difference, linear and quadratic interp
    extrap5 = SVector{5}([-0.25, 2.0, 1.0, -2.0, 0.25])
    extrap6 = SVector{5}([-2.75, 14.0, -21.5, 10.0, -0.5])

    extrap = [extrap0, extrap1, extrap2, extrap3, extrap4, extrap5, extrap6]

end
