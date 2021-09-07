#=
	Functions used to calculate total magnetization of the
	spin lattice
=#
module Magnetization

export magnetization

	function magnetization(mat)

		p, m, n = size(mat)
		mz = 0.

		for i in 1:m, j in 1:n
			mz = mz + mat[3,i,j]
		end

		return mz/(m*n)

	end

end
