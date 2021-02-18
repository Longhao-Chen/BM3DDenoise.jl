"""
	bm3d_wie(img::Array{Float64}, imgBasic::Array{Float64}, sigma::AbstractFloat, config::bm3d_config)

2nd step of BM3D
img: input noisy image
imgBasic: denoised image from first step of BM3D (hard thresholding)
sigma: known or assumed standard deviation of noise
"""
function bm3d_wie(img::Array{Float64}, imgBasic::Array{Float64}, sigma::AbstractFloat, config::bm3d_config)

	# parameters
	patchSize = config.wie_patchSize
	stepSize = config.wie_stepSize
	nBorder = config.wie_nBorder
	searchWin = config.wie_searchWin
	nMatch = config.wie_nMatch
	threshSimilar = config.wie_threshSimilar
	thresh3D =config.wie_thresh3D

	# block matching step
	@info "2st get_reference_pixels"
	(Ilist,Jlist) = get_reference_pixels([size(img,1);size(img,2)],patchSize,stepSize,nBorder)
	@info "2st match_patches"
	matchTable = match_patches(imgBasic,Ilist,Jlist,patchSize,searchWin,nMatch)
	@info "2st match_patches end"

	G3D = zeros(Float64, nMatch+1, patchSize[1], patchSize[2], axes(img)[3:end]...)
	G3Dbasic = similar(G3D)
	WC = similar(G3Dbasic)
	Wout = zeros(Float64, size(img))
	imgOut = zeros(Float64, size(img))

	# 3D filtering
	@info "2st 3D filtering"
	@views valid_match = matchTable[3, :, :, :] .< threshSimilar	# Determine which matching blocks meet the requirements.
	wie_3D_filtering!(Wout, imgOut, G3D, G3Dbasic, img, imgBasic, matchTable, valid_match, WC, Ilist, Jlist, patchSize, sigma)
	
	return imgOut ./ Wout

end

"""
3D filtering
"""
function wie_3D_filtering!(Wout::AbstractArray{<:AbstractFloat, 2},
			imgOut::AbstractArray{<:AbstractFloat, 2},
			G3D::AbstractArray{<:AbstractFloat, 3},
			G3Dbasic::AbstractArray{<:AbstractFloat, 3},
			img::AbstractArray{<:AbstractFloat, 2},
			imgBasic::AbstractArray{<:AbstractFloat, 2},
			matchTable::Array{<:AbstractFloat},
			vaild_match::BitArray,
			WC::AbstractArray{<:AbstractFloat},
			Ilist::Array{Int}, Jlist::Array{Int},
			patchSize::Array{Int}, sigma::AbstractFloat)
	# Each reference block is processed to reduce memory usage
	@views @inbounds for J = 1:length(Jlist)
		for I = 1:length(Ilist)
			# Compute 3D group spectrum
			form_group!(G3D, img, matchTable, vaild_match, Ilist, Jlist, patchSize, (I, J))
			form_group!(G3Dbasic, imgBasic, matchTable, vaild_match, Ilist, Jlist, patchSize, (I, J))

			# Wiener filtering of 3D groups, using basic estimate as target spectrum
			WC .= @strided G3Dbasic.^2 ./ (G3Dbasic.^2 .+ sigma^2)
			G3D .*= WC

			# Weight
			T = norm(WC, 2)
			W = T > 0 ? 1.0/T : 1.0
			G3D .*= W

			invert_group!(imgOut, G3D, matchTable, vaild_match, Ilist, Jlist, patchSize, (I, J))
			group_to_image!(Wout, W, matchTable, vaild_match, Ilist, Jlist, patchSize, (I, J))
		end
	end
end

function wie_3D_filtering!(Wout::Array{<:AbstractFloat, 3},
			imgOut::Array{<:AbstractFloat, 3},
			G3D::Array{<:AbstractFloat, 4},
			G3Dbasic::Array{<:AbstractFloat, 4},
			img::Array{<:AbstractFloat, 3},
			imgBasic::Array{<:AbstractFloat, 3},
			matchTable::Array{<:AbstractFloat},
			vaild_match::BitArray,
			WC::Array{<:AbstractFloat},
			Ilist::Array{Int}, Jlist::Array{Int},
			patchSize::Array{Int}, sigma::AbstractFloat)
	@views @inbounds Threads.@threads for i = 1:size(img, 3)
		wie_3D_filtering!(Wout[:, :, i], imgOut[:, :, i], G3D[:, :, :, i], G3Dbasic[:, :, :, i], img[:, :, i], imgBasic[:, :, i],
			matchTable, vaild_match, WC[:, :, :, i], Ilist, Jlist, patchSize, sigma)
	end
end