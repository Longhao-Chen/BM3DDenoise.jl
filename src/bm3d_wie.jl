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
	thresh3D =config.wie_thresh3D

	# block matching step
	@info "2st get_reference_pixels"
	(Ilist,Jlist) = get_reference_pixels([size(img,1);size(img,2)],patchSize,stepSize,nBorder)
	@info "2st match_patches"
	matchTable = match_patches(imgBasic,Ilist,Jlist,patchSize,searchWin,nMatch)
	@info "2st match_patches end"

	Wout = zeros(Float64, size(img))
	imgOut = zeros(Float64, size(img))

	# 3D filtering
	@info "2st 3D filtering"
	wie_3D_filtering!(Wout, imgOut, img, imgBasic, matchTable, Ilist, Jlist, patchSize, searchWin, nMatch, sigma)
	
	return imgOut ./ Wout

end

"""
3D filtering
"""
function wie_3D_filtering!(Wout::AbstractArray{<:AbstractFloat, 2},
			imgOut::AbstractArray{<:AbstractFloat, 2},
			img::AbstractArray{<:AbstractFloat, 2},
			imgBasic::AbstractArray{<:AbstractFloat, 2},
			matchTable::Array{<:AbstractFloat},
			Ilist::Array{Int}, Jlist::Array{Int},
			patchSize::Array{Int}, searchWin::Array{Int},
			nMatch::Int, sigma::AbstractFloat)
	# Each reference block is processed to reduce memory usage
	I_end = length(Ilist)
	J_end = length(Jlist)
	# Preventing conflicts in parallel computing
	@views @inbounds for offset in 0:2searchWin[2] - 1
		Threads.@threads for J = 1 + offset:2searchWin[2]:J_end
			G3D = zeros(Float64, nMatch+1, patchSize[1], patchSize[2])
			G3Dbasic = similar(G3D)
			WC = similar(G3Dbasic)
			for I = 1:I_end
				# Compute 3D group spectrum
				form_group!(G3D, img, matchTable, Ilist, Jlist, patchSize, (I, J))
				form_group!(G3Dbasic, imgBasic, matchTable, Ilist, Jlist, patchSize, (I, J))

				# Wiener filtering of 3D groups, using basic estimate as target spectrum
				WC .= G3Dbasic.^2 ./ (G3Dbasic.^2 .+ sigma^2)
				G3D .*= WC

				# Weight
				T = norm(WC, 1)
				W = T > 0 ? 1.0/(T *sigma^2) : 1.0
				G3D .*= W

				invert_group!(imgOut, G3D, matchTable, Ilist, Jlist, patchSize, (I, J))
				group_to_image!(Wout, W, matchTable, Ilist, Jlist, patchSize, (I, J))
			end
		end
	end
end

function wie_3D_filtering!(Wout::Array{<:AbstractFloat, 3},
			imgOut::Array{<:AbstractFloat, 3},
			img::Array{<:AbstractFloat, 3},
			imgBasic::Array{<:AbstractFloat, 3},
			matchTable::Array{<:AbstractFloat},
			Ilist::Array{Int}, Jlist::Array{Int},
			patchSize::Array{Int}, searchWin::Array{Int},
			nMatch::Int, sigma::AbstractFloat)
	@views @inbounds Threads.@threads for i = 1:size(img, 3)
		wie_3D_filtering!(Wout[:, :, i], imgOut[:, :, i], img[:, :, i], imgBasic[:, :, i],
			matchTable, Ilist, Jlist, patchSize, searchWin, nMatch, sigma)
	end
end