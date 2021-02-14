module BM3D
export bm3d

import FFTW
using LinearAlgebra
using Strided
using Images

include("bm3d_thr.jl")
include("bm3d_wie.jl")
include("get_reference_pixels.jl")
include("groups.jl")
include("HardThresholding.jl")
include("kaiser_window.jl")
include("match_patches.jl")

"""
	bm3d(img, σ)

img: input grayscale image(Float64 or Gray)
σ: known or assumed standard deviation of noise
return: denoised image
"""
function bm3d(img::Matrix{Float64}, σ::AbstractFloat)
	imgBasic = bm3d_thr(img, σ)
	bm3d_wie(img, imgBasic, σ)
end

function bm3d(img::Matrix{Gray{T}}, σ::AbstractFloat) where {T}
	img_float = bm3d(Float64.(img), σ)
	Gray.(img_float)
end

end # module