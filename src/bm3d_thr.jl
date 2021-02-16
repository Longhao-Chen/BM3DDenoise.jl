"""
	bm3d_thr(img::Matrix{Float64}, sigma::AbstractFloat)

1st step of BM3D denoising: hard thresholding
"""
function bm3d_thr(img::Matrix{Float64}, sigma::AbstractFloat)

	# Hard code algorithm parameters for now...
	patchSize = [8;8] 
	stepSize = [3;3]  
	nBorder = [0;0]
	searchWin = [19;19]
	nMatch = 31
	thresh3D = 2.7

	# Block matching
	@info "1st get_reference_pixels"
	(Ilist, Jlist) = get_reference_pixels([size(img,1); size(img,2)], patchSize, stepSize, nBorder)
	@info "1st match_patches"
	matchTable = match_patches(img, Ilist, Jlist, patchSize, searchWin, nMatch)
	@info "1st match_patches end"

	Wout = zeros(Float64, size(img))
	imgOut = zeros(Float64, size(img))
	G3D = zeros(Float64, nMatch+1, patchSize[1], patchSize[2])

	# Each reference block is processed to reduce memory usage
	for J = 1:length(Jlist)
		for I = 1:length(Ilist)
			form_group!(G3D, img, matchTable, Ilist, Jlist, patchSize, (I, J))

			# Filter 3D groups by hard thresholding 
			HardThresholding!(G3D, sigma * thresh3D)

			T = nnz(G3D)
			W = T > 0 ? 1.0 / T : 1.0
			@strided G3D .*= W

			invert_group!(imgOut, G3D, matchTable, Ilist, Jlist, patchSize, (I, J)) 
			group_to_image!(Wout, W, matchTable, Ilist, Jlist, patchSize, (I, J))

		end
	end
	return imgOut ./ Wout

end

"""
	nnz(data::AbstractArray{Float64})

Returns the number of non-zero elements in the array
"""
function nnz(data::AbstractArray{Float64})
	sum(data .!= 0.)
end