"""
	bm3d_thr(img::Array{Float64}, sigma::AbstractFloat, config::bm3d_config)

1st step of BM3D denoising: hard thresholding

img: input noisy image
sigma: known or assumed standard deviation of noise
config: bm3d_config
"""
function bm3d_thr(img::Array{Float64}, sigma::AbstractFloat, config::bm3d_config)

	# parameters
	patchSize = config.thr_patchSize
	searchStride = config.thr_searchStride
	searchWin = config.thr_searchWin
	nMatch = config.thr_nMatch
	thresh3D = config.thr_thresh3D

	# Block matching
	config.show_info && @info "1st get_reference_pixels"
	if config.show_info
		refIndex = @time get_reference_pixels(
			[size(img, 1); size(img, 2)],
			patchSize,
			searchStride,
		)
	else
		refIndex = get_reference_pixels(
			[size(img, 1); size(img, 2)],
			patchSize,
			searchStride,
		)
	end
	config.show_info && @info "1st get_reference_pixels end"
	config.show_info && @info "1st match_patches"
	if config.show_info
		matchTable = @time match_patches(
			img,
			refIndex,
			patchSize,
			searchWin,
			nMatch,
			config.thr_distFunction,
		)
	else
		matchTable = match_patches(
			img,
			refIndex,
			patchSize,
			searchWin,
			nMatch,
			config.thr_distFunction,
		)
	end
	config.show_info && @info "1st match_patches end"

	# Don't use similar(), because it need initialize to 0.0
	Wout = zeros(Float64, size(img)...)
	imgOut = zeros(Float64, size(img)...)

	# 3D filtering
	config.show_info && @info "1st 3D filtering"
	if config.show_info
		@time thr_3D_filtering!(
			Wout,
			imgOut,
			img,
			matchTable,
			refIndex,
			patchSize,
			nMatch,
			thresh3D,
			sigma,
			config,
		)
	else
		thr_3D_filtering!(
			Wout,
			imgOut,
			img,
			matchTable,
			refIndex,
			patchSize,
			nMatch,
			thresh3D,
			sigma,
			config,
		)
	end
	config.show_info && @info "1st 3D filtering end"

	return imgOut ./ Wout
end

"""
3D filtering
"""
function thr_3D_filtering!(
	Wout::AbstractArray{<:AbstractFloat,2},
	imgOut::AbstractArray{<:AbstractFloat,2},
	img::AbstractArray{<:AbstractFloat,2},
	matchTable::Array{CartesianIndex{2},3},
	refIndex::Array{CartesianIndex{2},2},
	patchSize::CartesianIndex{2},
	nMatch::Int,
	thresh3D::AbstractFloat,
	sigma::AbstractFloat,
	config::bm3d_config,
)
	# Each reference block is processed to reduce memory usage
	I_end, J_end = size(refIndex)

	# Preventing conflicts in group_to_image!
	imgLockPool = Array{ReentrantLock}(undef, I_end, J_end)
	for i in 1:length(imgLockPool)
		imgLockPool[i] = ReentrantLock()
	end

	@views @inbounds Threads.@threads for J in 1:J_end
		# Preventing conflicts in parallel computing
		G3D = zeros(Float64, nMatch + 1, patchSize[1], patchSize[2])
		for I in 1:I_end
			form_group!(
				G3D,
				img,
				matchTable,
				refIndex,
				patchSize,
				CartesianIndex(I, J),
				config.thr_group_transform!,
			)

			# Filter 3D groups by hard thresholding
			T = HardThresholding!(G3D, sigma * thresh3D)

			if T > zero(T)
				W = 1.0 / (T * sigma^2)
				G3D .*= W
			else
				W = 1.0
			end

			invert_group!(
				imgOut,
				G3D,
				matchTable,
				refIndex,
				patchSize,
				CartesianIndex(I, J),
				config.thr_group_itransform!,
				imgLockPool,
			)
			group_to_image!(
				Wout,
				W,
				matchTable,
				refIndex,
				patchSize,
				CartesianIndex(I, J),
				imgLockPool,
			)
		end
	end
end

# For color images
function thr_3D_filtering!(
	Wout::Array{<:AbstractFloat,3},
	imgOut::Array{<:AbstractFloat,3},
	img::Array{<:AbstractFloat,3},
	matchTable::Array{CartesianIndex{2},3},
	refIndex::Array{CartesianIndex{2},2},
	patchSize::CartesianIndex{2},
	nMatch::Int,
	thresh3D::AbstractFloat,
	sigma::AbstractFloat,
	config::bm3d_config,
)
	@views @inbounds for i in 1:size(img, 3)
		thr_3D_filtering!(
			Wout[:, :, i],
			imgOut[:, :, i],
			img[:, :, i],
			matchTable,
			refIndex,
			patchSize,
			nMatch,
			thresh3D,
			sigma,
			config,
		)
	end
end
