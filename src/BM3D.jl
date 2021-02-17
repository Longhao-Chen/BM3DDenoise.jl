module BM3D
export bm3d, bm3d_config

import FFTW
using LinearAlgebra
using Strided
using Images

include("bm3d_thr.jl")
include("bm3d_wie.jl")
include("get_reference_pixels.jl")
include("groups.jl")
include("HardThresholding.jl")
include("match_patches.jl")
include("parameter.jl")

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

function bm3d(img::Matrix{Float64}, σ::AbstractFloat, config::bm3d_config)
	imgBasic = bm3d_thr(img, σ, config)
	bm3d_wie(img, imgBasic, σ, config)
end

function bm3d(img::Matrix{Gray{T}}, σ::AbstractFloat, config::bm3d_config) where {T}
	img_float = bm3d(Float64.(img), σ, config)
	Gray.(img_float)
end

end # module