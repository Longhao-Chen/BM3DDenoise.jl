module BM3D
export bm3d, bm3d_config

import FFTW
using LinearAlgebra
using Strided
using ImageCore

include("parameter.jl")
include("bm3d_thr.jl")
include("bm3d_wie.jl")
include("get_reference_pixels.jl")
include("groups.jl")
include("HardThresholding.jl")
include("match_patches.jl")

"""
	bm3d(img, σ)
	bm3d(img, σ, config)

```
img: input grayscale image(Float64 or Gray)
σ: known or assumed standard deviation of noise
config: see bm3d_config.
return: denoised image
```
"""
function bm3d(img, σ::AbstractFloat)
	config = bm3d_config()	# Use default parameters
	bm3d(img, σ, config)
end

function bm3d(img::Array{Float64}, σ::AbstractFloat, config::bm3d_config)
	print(config)
	imgBasic = bm3d_thr(img, σ, config)
	bm3d_wie(img, imgBasic, σ, config)
end

function bm3d(img::Matrix{Gray{T}}, σ::AbstractFloat, config::bm3d_config) where {T}
	img_float = bm3d(Float64.(img), σ, config)
	Gray.(img_float)
end

# Here is the color image
function bm3d(img::Matrix{YCbCr{T}}, σ::AbstractFloat, config::bm3d_config) where {T}
	# Split into 3-dimensional array, each dimension data is:[:, :, 1] - Y; [:, :, 2] - Cb; [:, :, 3] -Cr
	img_Array = permutedims(Float64.(channelview(img)),[2,3,1])
	# normalization
	@strided img_Array .-= 16.
	@strided img_Array[:, :, 1] ./= (235. - 16.)
	@strided img_Array[:, :, 2:3] ./= (240. - 16.)

	img_Array .= bm3d(img_Array, σ, config)

	# Return to the original range
	@strided img_Array[:, :, 1] .*= (235. - 16.)
	@strided img_Array[:, :, 2:3] .*= (240. - 16.)
	@strided img_Array .+= 16.

	colorview(YCbCr, permutedims(T.(img_Array),[3,1,2]))
end

function bm3d(img::Matrix{RGB{T}}, σ::AbstractFloat, config::bm3d_config) where {T}
	img_YCbCr = bm3d(YCbCr.(img), σ, config)
	RGB.(img_YCbCr)
end

end # module