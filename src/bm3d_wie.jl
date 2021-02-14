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
	(Ilist,Jlist) = get_reference_pixels([size(img,1);size(img,2)],patchSize,stepSize,nBorder)
	matchTable = match_patches(imgBasic,Ilist,Jlist,patchSize,searchWin,nMatch)

	# Compute 3D group spectrum
	G3D = form_groups(img,matchTable,Ilist,Jlist,patchSize)
	G3Dbasic = form_groups(imgBasic,matchTable,Ilist,Jlist,patchSize)

	# Wiener filtering of 3D groups, using basic estimate as target spectrum
	WC = @strided G3Dbasic.^2 ./ (G3Dbasic.^2 .+ sigma^2) # devec?
	@strided G3D .*= WC

	# Weight groups 
	W = zeros(Float64,size(G3D))
	@inbounds @views Base.Threads.@threads for j = 1:length(Jlist)
		for i = 1:length(Ilist)
			T = norm(WC[:,:,:,i,j], 2)
			W[:,:,:,i,j] .= T > 0 ? 1.0/T : 1.0
		end
	end

	@strided G3D .*= W

	imgOut = invert_groups([size(img,1);size(img,2)], G3D, matchTable, Ilist, Jlist, patchSize) 

	Wout = zeros(Float64,size(img))
	groups_to_image!(Wout,W,matchTable,Ilist,Jlist,patchSize)

	return imgOut ./ Wout

end

"""
	nnz(data::AbstractArray{Float64})

Returns the number of non-zero elements in the array
"""
function nnz(data::AbstractArray{Float64})
	sum(data .!= 0.)
end