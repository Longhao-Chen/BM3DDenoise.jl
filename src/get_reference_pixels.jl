"""
	get_reference_pixels(imgSize::Vector{Int64},
					patchSize::Vector{Int64},
					stepSize::Vector{Int64},
					nBorder::Vector{Int64})

Get locations of reference pixels, the upper-left corner of each patch.

imgSize - size of input image, including border
patchSize - size of patch in pixels
stepSize - step in x and y between reference pixels
nBorder - number of border pixels on each side in x and y
"""
function get_reference_pixels(imgSize::Vector{Int64},
				patchSize::Vector{Int64},
				stepSize::Vector{Int64},
				nBorder::Vector{Int64})

	ph = imgSize[1] - 2*nBorder[1] - patchSize[1] + 1
	pw = imgSize[2] - 2*nBorder[2] - patchSize[2] + 1

	I = [nBorder[1] .+ 1:stepSize[1]:ph...]
	J = [nBorder[2] .+ 1:stepSize[2]:pw...]

	# Make sure there is a patch touching the lower and right borders
	if maximum(I) < nBorder[1] + ph
		I = [I; nBorder[1] + ph]
	end
	if maximum(J) < nBorder[2] + pw
		J = [J; nBorder[2] + pw]
	end

	return (I,J)
end