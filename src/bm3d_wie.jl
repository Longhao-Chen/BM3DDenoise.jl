"""
	bm3d_wie(img::Matrix{Float64}, imgBasic::Matrix{Float64}, sigma::AbstractFloat)

2nd step of BM3D
img: input noisy image
imgBasic: denoised image from first step of BM3D (hard thresholding)
sigma: known or assumed standard deviation of noise
"""
function bm3d_wie(img::Matrix{Float64}, imgBasic::Matrix{Float64}, sigma::AbstractFloat)

	# Hard code algorithm parameters for now...
	patchSize = [8;8] 
	stepSize = [3;3]  
	nBorder = [0;0]
	searchWin = [11;11]
	nMatch = 15
	thresh3D = 2.7

	# block matching step
	@info "2st get_reference_pixels"
	(Ilist,Jlist) = get_reference_pixels([size(img,1);size(img,2)],patchSize,stepSize,nBorder)
	@info "2st match_patches"
	matchTable = match_patches(imgBasic,Ilist,Jlist,patchSize,searchWin,nMatch)
	@info "2st match_patches end"

	G3D = zeros(Float64, nMatch+1, patchSize[1], patchSize[2])
	G3Dbasic = zeros(Float64, nMatch+1, patchSize[1], patchSize[2])
	WC = similar(G3Dbasic)
	Wout = zeros(Float64, size(img))
	imgOut = zeros(Float64, size(img))

	# Each reference block is processed to reduce memory usage
	for J = 1:length(Jlist)
		for I = 1:length(Ilist)
			# Compute 3D group spectrum
			form_group!(G3D, img, matchTable, Ilist, Jlist, patchSize, (I, J))
			form_group!(G3Dbasic, imgBasic, matchTable, Ilist, Jlist, patchSize, (I, J))

			# Wiener filtering of 3D groups, using basic estimate as target spectrum
			WC .= @strided G3Dbasic.^2 ./ (G3Dbasic.^2 .+ sigma^2)
			@strided G3D .*= WC

			# Weight
			T = norm(WC, 2)
			W = T > 0 ? 1.0/T : 1.0
			@strided G3D .*= W

			invert_group!(imgOut, G3D, matchTable, Ilist, Jlist, patchSize, (I, J))
			group_to_image!(Wout, W, matchTable, Ilist, Jlist, patchSize, (I, J))
		end
	end
	return imgOut ./ Wout

end