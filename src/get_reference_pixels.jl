"""
	get_reference_pixels(imgSize::Vector{Int64},
					patchSize::CartesianIndex{2},
					searchStride::CartesianIndex{2})

Get locations of reference pixels, the upper-left corner of each patch.

imgSize - size of input image, including border
patchSize - size of patch in pixels
searchStride - step in x and y between reference pixels
"""
function get_reference_pixels(imgSize::Vector{Int64},
				patchSize::CartesianIndex{2},
				searchStride::CartesianIndex{2})

	imgH, imgW = imgSize

	if (imgH - patchSize[1]) % searchStride[1] != 0
		tmp_H = 1:searchStride[1]:(imgH - patchSize[1] + 1)
		H = Vector{UInt}(undef, length(tmp_H) + 1)
		H[1:end - 1] .= tmp_H
		H[end] = imgH - patchSize[1] + 1
	else
		H = 1:searchStride[1]:(imgH - patchSize[1] + 1)
	end

	if (imgW - patchSize[2] ) % searchStride[2] != 0
		tmp_W = 1:searchStride[2]:(imgW - patchSize[2] + 1)
		W = Vector{UInt}(undef, length(tmp_W) + 1)
		W[1:end - 1] .= tmp_W
		W[end] = imgW - patchSize[2] + 1
	else
		W = 1:searchStride[2]:(imgW - patchSize[2] + 1)
	end

	refIndex = Array{CartesianIndex{2}}(undef, length(H), length(W))

	for j in 1:length(W), i in 1:length(H)
		refIndex[i, j] = CartesianIndex(H[i], W[j])
	end
	return refIndex
end