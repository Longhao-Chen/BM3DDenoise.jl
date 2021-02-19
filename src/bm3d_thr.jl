"""
	bm3d_thr(img::Array{Float64}, sigma::AbstractFloat, config::bm3d_config)

1st step of BM3D denoising: hard thresholding
"""
function bm3d_thr(img::Array{Float64}, sigma::AbstractFloat, config::bm3d_config)

	# parameters
	patchSize = config.thr_patchSize
	stepSize = config.thr_stepSize
	nBorder = config.thr_nBorder
	searchWin = config.thr_searchWin
	nMatch = config.thr_nMatch
	threshSimilar = config.thr_threshSimilar
	thresh3D = config.thr_thresh3D

	# Block matching
	@info "1st get_reference_pixels"
	(Ilist, Jlist) = get_reference_pixels([size(img,1); size(img,2)], patchSize, stepSize, nBorder)
	@info "1st match_patches"
	matchTable = match_patches(img, Ilist, Jlist, patchSize, searchWin, nMatch)
	@info "1st match_patches end"

	Wout = zeros(Float64, size(img))
	imgOut = zeros(Float64, size(img))
	G3D = zeros(Float64, nMatch+1, patchSize[1], patchSize[2], axes(img)[3:end]...)

	# 3D filtering
	@info "1st 3D filtering"
	@views valid_match = matchTable[3, :, :, :] .< threshSimilar	# Determine which matching blocks meet the requirements.
	thr_3D_filtering!(Wout, imgOut, G3D, img, matchTable, valid_match, Ilist, Jlist, patchSize, thresh3D, sigma)

	return imgOut ./ Wout
end

"""
3D filtering
"""
function thr_3D_filtering!(Wout::AbstractArray{<:AbstractFloat, 2},
			imgOut::AbstractArray{<:AbstractFloat, 2},
			G3D::AbstractArray{<:AbstractFloat, 3},
			img::AbstractArray{<:AbstractFloat, 2},
			matchTable::Array{<:AbstractFloat},
			valid_match::BitArray,
			Ilist::Array{Int}, Jlist::Array{Int},
			patchSize::Array{Int}, thresh3D::AbstractFloat, sigma::AbstractFloat)
	# Each reference block is processed to reduce memory usage
	@views @inbounds for J = 1:length(Jlist)
		for I = 1:length(Ilist)
			form_group!(G3D, img, matchTable, valid_match, Ilist, Jlist, patchSize, (I, J))

			# Filter 3D groups by hard thresholding 
			HardThresholding!(G3D, sigma * thresh3D)

			T = nnz(G3D)
			W = T > 0 ? 1.0 / (T * sigma^2) : 1.0
			G3D .*= W

			invert_group!(imgOut, G3D, matchTable, valid_match, Ilist, Jlist, patchSize, (I, J))
			group_to_image!(Wout, W, matchTable, valid_match, Ilist, Jlist, patchSize, (I, J))

		end
	end
end

# For color images
function thr_3D_filtering!(Wout::Array{<:AbstractFloat, 3},
	imgOut::Array{<:AbstractFloat, 3},
	G3D::Array{<:AbstractFloat, 4},
	img::Array{<:AbstractFloat, 3},
	matchTable::Array{<:AbstractFloat},
	valid_match::BitArray,
	Ilist::Array{Int}, Jlist::Array{Int},
	patchSize::Array{Int}, thresh3D::AbstractFloat, sigma::AbstractFloat)
	@views @inbounds Threads.@threads for i = 1:size(img, 3)
		thr_3D_filtering!(Wout[:, :, i], imgOut[:, :, i], G3D[:, :, :, i], img[:, :, i], matchTable, valid_match, Ilist, Jlist, patchSize, thresh3D, sigma)
	end
end

"""
	nnz(data::AbstractArray{Float64})

Returns the number of non-zero elements in the array
"""
function nnz(data::AbstractArray{Float64})
	sum(data .!= 0.)
end