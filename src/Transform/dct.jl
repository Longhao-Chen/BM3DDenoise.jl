"""
type-II discrete cosine transform (DCT)
"""
function dct!(data::Array{T, 3}) where {T <:AbstractFloat}
	dctmat = genmat(size(data, 1))
	dctmat1 = genmat(size(data, 2))'
	@views @inbounds @fastmath for i in 1:size(data, 3)
		data[:, :, i] .= dctmat * data[:, :, i] * dctmat1
	end

	dctmat = genmat(size(data, 3))
	@views @inbounds @fastmath for j in 1:size(data, 2), i in 1:size(data, 1)
		data[i, j, :] .= dctmat * data[i, j, :]
	end
end

"""
inverse discrete cosine transform (DCT)
"""
function idct!(data::Array{T, 3}) where {T <:AbstractFloat}
	dctmat = genmat(size(data, 3))'
	@views @inbounds @fastmath for j in 1:size(data, 2), i in 1:size(data, 1)
		data[i, j, :] .= dctmat * data[i, j, :]
	end

	dctmat = genmat(size(data, 1))'
	dctmat1 = genmat(size(data, 2))
	@views @inbounds @fastmath for i in 1:size(data, 3)
		data[:, :, i] .= dctmat * data[:, :, i] * dctmat1
	end
end

function genmat(h::Int)::Array{Float64, 2}
	mat = Array{Float64, 2}(undef, h, h)
	mat[1, :] .= sqrt(1 / h)
	@inbounds @fastmath for j in 1:h
		@simd for i in 2:h
			mat[i, j] = cos(π * (i - 1) * (j - 0.5) / h) * sqrt(2 / h)
		end
	end
	mat
end
