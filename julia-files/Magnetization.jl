#=
	Functions used to calculate total magnetization of the
	spin lattice
=#
module Magnetization

export magnetization

	function magnetization(mat)

		p,m,n = size(mat)
		mz = 0.

		for i in eachindex(view(mat,1,1:m,1:n))
			mz = mz + mat[i[1],i[2],3]
		end

		return mz/(m*n)

	end

end
