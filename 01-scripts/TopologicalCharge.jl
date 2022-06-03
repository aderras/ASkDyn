# This module contains functions used to calculate the topological charge
# of a spin field
#
#
module TopologicalCharge

    using ShiftedArrays
    export calcQ

    # The following function calculates the topological charge of
    # the entire lattice.
    #
    # in: mat = NxNx3 spin array
    # out: q = float
    function calcQ(mat::Array{Float64,3})

        p,m,n = size(mat)
        q = 0.0

        dxsQ = Array{Float64}(undef,p,m,n)
        dysQ = Array{Float64}(undef,p,m,n)

		dxsQ .= (-1 .*ShiftedArrays.circshift(mat,(0,2,0)).+
	           8 .* ShiftedArrays.circshift(mat,(0,1,0)).-
	           8 .* ShiftedArrays.circshift(mat,(0,-1,0)).+
	            ShiftedArrays.circshift(mat,(0,-2,0)))./12
		dysQ .= (-1 .* ShiftedArrays.circshift(mat,(0,0,2)).+
		    8 .* ShiftedArrays.circshift(mat,(0,0,1)).-
		    8 .* ShiftedArrays.circshift(mat,(0,0,-1)).+
		    ShiftedArrays.circshift(mat,(0,0,-2)))./12

		for i in 1:m, j in 1:n
		    q += mat[1,i,j]*dxsQ[2,i,j]*dysQ[3,i,j] +
		    mat[2,i,j]*dxsQ[3,i,j]*dysQ[1,i,j] +
		    mat[3,i,j]*dxsQ[1,i,j]*dysQ[2,i,j] -
		    mat[3,i,j]*dxsQ[2,i,j]*dysQ[1,i,j] -
		    mat[2,i,j]*dxsQ[1,i,j]*dysQ[3,i,j] -
		    mat[1,i,j]*dxsQ[3,i,j]*dysQ[2,i,j]
		end
		q = q/(4*pi)
		return q

    end

end
