#=

    This module contains structs and functions related to parameters that are
    passed to the spin dynamics computation functions.

    Params are grouped by catagory and stored in structs. E.g. the
    materialParams struct contains the values for exchange interaction,
    external field, shape of the array, etc.

    All of the parameters are then group into one parent struct of the type
    "params." This parent struct is passed to different functions, and can
    be extended to more specific problems by including additional structs.
    E.g. the parameters for a defect in the lattice are stored in the
    defectParams struct, which is an element of params.

    Below the structs are the functions

    IMPROVEMENTS:
        + add warning if user requests RK discretization that doesn't make sense
        + also warn if the range of parameters requested doesn't make sense

=#
module Params

    # Material parameters
    mutable struct materialParams

        j::Float64      # Exchange interaction
        h::Float64      # Zeeman field
        a::Float64      # Dzayaloshinskii-Moriya interaction
        dz::Float64     # Perpendicular magnetic anisotropy
        ed::Float64     # Dipole-dipole constant
        nx::Int         # Lattice size in x
        ny::Int         # Lattice size in y
        nz::Int         # Lattice size in z
        bc::Float64     # Periodic boundary conditions

        v               # Matrices to be used to compute DDI. It
                        # is advantageous to compute these once and
                        # store them because this is slow.
    end

    mutable struct compParams
        solver::Int       # Float
        maxSteps::Int     # Maximum number of steps to make
        dt::Float64       # Step size
        nn::Float64       # Skip this many steps to save
        tol::Float64      # Tolerance for convergence
        damp::Float64     # Damping in llg
        temp::Float64     # Temperature
        parallel::Int     # Parallel (1 or 0)
        numCores::Int     # Number of parallel cores
        relax             # Run relaxation if data doesn't exist

        # These are for FA relaxation, if using
        sMax::Int         # Maximum number of steps to make
        faconst::Float64  # Field alignment param
        nRot::Float64     # Number of rotations to make
        faTol::Float64    # Tolerance for convergence

        rChirality::Int
    end

    # Initial condition parameters
    mutable struct icParams

        type::String        # Type of spin structure
        r::Float64          # Radius (only makes sense for skyrmion)
        gamma::Float64      # Chirality (only skyrmion)
        px         # Position in x direction (0 < px < nx)
        py         # Position in y direction (0 < py < ny)

    end

    # Pinning parameters
    mutable struct pinningParams
        hPin::Float64       # Pinning field strength
    end

    # Defect parameters
    mutable struct defectParams

        t                   # Type of defect (1 = Point, 2 = Gaussian)
        strength::Float64   # Strength of the defect
        width::Float64      # Width of impact
        dx::Int           # x position (0 < dx < nx)
        dy::Int           # y position (0 < dy < ny)
        jmat
    end

    mutable struct currentParams
        jx::Float64          # Current in the x direction
        jy::Float64          # y direction
        tf::Float64          # Cutoff time for current application
    end

    # Choices for data saving
    mutable struct saveChoices

        totE::Int       # Total energy
        excE::Int       # Exchange energy
        zeeE::Int       # Zeeman energy
        dmiE::Int       # DMI energy
        pmaE::Int       # PMA energy
        ddiE::Int       # DDI energy
        magn::Int       # Magnetization
        size::Int       # Effective size (skyrmion)
        qCharge::Int    # Topological charge
        location::Int   # Position (skyrmion)
        spinField::Int  # Spin field (The whole matrix. Data intensive.)
        chir::Int       # Chirality
    end

    # All parameters are stored in this struct.
    mutable struct params

        mp          # materialParams
        cp          # computationParams
        ic          # icParams
        pin         # pinningParams
        defect      # defectParams
        current     # current
        save        # saveChoices

    end

end
