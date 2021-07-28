"""
	bm3d_wie(img::Array{Float64}, imgBasic::Array{Float64}, sigma::AbstractFloat, config::bm3d_config)

2nd step of BM3DDenoise

img: input noisy image
imgBasic: denoised image from first step of BM3D (hard thresholding)
sigma: known or assumed standard deviation of noise
config: bm3d_config
"""
function bm3d_wie(
	img::Array{Float64},
	imgBasic::Array{Float64},
	sigma::AbstractFloat,
	config::bm3d_config,
)

	# parameters
	patchSize = config.wie_patchSize
	searchStride = config.wie_searchStride
	searchWin = config.wie_searchWin
	nMatch = config.wie_nMatch

	# block matching step
	@info "2st get_reference_pixels"
	@time refIndex =
		get_reference_pixels([size(img, 1); size(img, 2)], patchSize, searchStride)
	@info "2st get_reference_pixels end"
	@info "2st match_patches"
	@time matchTable = match_patches(
		imgBasic,
		refIndex,
		patchSize,
		searchWin,
		nMatch,
		config.wie_distFunction,
	)
	@info "2st match_patches end"

	# Don't use similar(), because it need initialize to 0.0
	Wout = zeros(Float64, size(img)...)
	imgOut = zeros(Float64, size(img)...)

	# 3D filtering
	@info "2st 3D filtering"
	@time wie_3D_filtering!(
		Wout,
		imgOut,
		img,
		imgBasic,
		matchTable,
		refIndex,
		patchSize,
		nMatch,
		sigma,
		config,
	)
	@info "2st 3D filtering end"

	return imgOut ./ Wout

end

"""
3D filtering
"""
function wie_3D_filtering!(
	Wout::AbstractArray{<:AbstractFloat,2},
	imgOut::AbstractArray{<:AbstractFloat,2},
	img::AbstractArray{<:AbstractFloat,2},
	imgBasic::AbstractArray{<:AbstractFloat,2},
	matchTable::Array{CartesianIndex{2},3},
	refIndex::Array{CartesianIndex{2},2},
	patchSize::CartesianIndex{2},
	nMatch::Int,
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
		G3Dbasic = similar(G3D)
		WC = similar(G3Dbasic)
		for I in 1:I_end
			# Compute 3D group spectrum
			form_group!(
				G3D,
				img,
				matchTable,
				refIndex,
				patchSize,
				CartesianIndex(I, J),
				config.wie_group_transform!,
			)
			form_group!(
				G3Dbasic,
				imgBasic,
				matchTable,
				refIndex,
				patchSize,
				CartesianIndex(I, J),
				config.wie_group_transform!,
			)

			# Wiener filtering of 3D groups, using basic estimate as target spectrum
			WC .= G3Dbasic .^ 2 ./ (G3Dbasic .^ 2 .+ sigma^2)
			G3D .*= WC

			# Weight
			T = norm(WC, 1)
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
				config.wie_group_itransform!,
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

function wie_3D_filtering!(
	Wout::Array{<:AbstractFloat,3},
	imgOut::Array{<:AbstractFloat,3},
	img::Array{<:AbstractFloat,3},
	imgBasic::Array{<:AbstractFloat,3},
	matchTable::Array{CartesianIndex{2},3},
	refIndex::Array{CartesianIndex{2},2},
	patchSize::CartesianIndex{2},
	nMatch::Int,
	sigma::AbstractFloat,
	config::bm3d_config,
)
	@views @inbounds for i in 1:size(img, 3)
		wie_3D_filtering!(
			Wout[:, :, i],
			imgOut[:, :, i],
			img[:, :, i],
			imgBasic[:, :, i],
			matchTable,
			refIndex,
			patchSize,
			nMatch,
			sigma,
			config,
		)
	end
end
