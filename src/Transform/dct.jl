"""
type-II discrete cosine transform (DCT)
"""
function dct!(data::AbstractArray{<:AbstractFloat,3})
	dctmat = genmat(size(data, 1))
	dctmat1 = genmat(size(data, 2))'
	@views @inbounds for i in 1:size(data, 3)
		data[:, :, i] .= dctmat * data[:, :, i] * dctmat1
	end
	dctmat = genmat(size(data, 3))
	@views @inbounds for j in 1:size(data, 2), i in 1:size(data, 1)
		data[i, j, :] .= dctmat * data[i, j, :]
	end
end

"""
inverse discrete cosine transform (DCT)
"""
function idct!(data::AbstractArray{<:AbstractFloat,3})
	dctmat = genmat(size(data, 3))'
	@views @inbounds for j in 1:size(data, 2), i in 1:size(data, 1)
		data[i, j, :] .= dctmat * data[i, j, :]
	end
	dctmat = genmat(size(data, 1))'
	dctmat1 = genmat(size(data, 2))
	@views @inbounds for i in 1:size(data, 3)
		data[:, :, i] .= dctmat * data[:, :, i] * dctmat1
	end
end

function genmat(h)
	mat = zeros(h, h)
	@inbounds for j in 1:h, i in 1:h
		mat[i, j] = cos(π * (i - 1) * (j - 0.5) / h)
	end
	mat[1, :] ./= sqrt(h)
	mat[2:end, :] .*= sqrt(2 / h)
	mat
end
