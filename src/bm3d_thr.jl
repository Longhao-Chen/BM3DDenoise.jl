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
	(Ilist, Jlist) = get_reference_pixels([size(img,1); size(img,2)], patchSize, stepSize, nBorder)
	matchTable = match_patches(img, Ilist, Jlist, patchSize, searchWin, nMatch)

	G3D = form_groups(img,matchTable,Ilist,Jlist,patchSize)

	# Filter 3D groups by hard thresholding 
	HardThresholding!(G3D, sigma * thresh3D)

	W = zeros(Float64, size(G3D))
	@inbounds @views Base.Threads.@threads for j = 1:length(Jlist)
		for i = 1:length(Ilist)
			T = nnz(G3D[:,:,:,i,j])
			W[:,:,:,i,j] .= T > 0 ? 1.0/T : 1.0
		end
	end
	@strided G3D .*= W

	imgOut = invert_groups([size(img,1); size(img,2)], G3D, matchTable, Ilist, Jlist, patchSize) 

	Wout = zeros(Float64, size(img))
	groups_to_image!(Wout, W, matchTable, Ilist, Jlist, patchSize)

	return imgOut ./ Wout

end