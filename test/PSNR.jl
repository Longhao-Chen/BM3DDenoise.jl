import Statistics
using Strided

"""
	PSNR(img1, img2; max_value::Float64=1.)

Compute Peak signal-to-noise ratio of img1 and img2.
max_value: The maximum pixel value of the image. For floating-point data, the maximum pixel value is 1.
"""
function PSNR(img1::AbstractArray{<:AbstractFloat}, img2::AbstractArray{<:AbstractFloat}; max_value::Float64=1.)
	10 * log10(max_value^2 / MSE(img1, img2))
end

function PSNR(img1::Matrix{Gray{T}}, img2::Matrix{Gray{T}}) where {T}
	PSNR(Float64.(img1), Float64.(img2))
end

function PSNR(img1::Matrix{RGB{T}}, img2::Matrix{RGB{T}}) where {T}
	MSE_mean = Statistics.mean(
		MSE(Float64.(red.(img1)), Float64.(red.(img2))),
		MSE(Float64.(green.(img1)), Float64.(green.(img2))),
		MSE(Float64.(blue.(img1)), Float64.(blue.(img2))),
		)
	10 * log10(1 / MSE_mean)
end

"""
	MSE(I::AbstractArray, K::AbstractArray)

Calculation of mean square error.
"""
function MSE(I::AbstractArray, K::AbstractArray)
	@strided Statistics.mean(abs2.(I .- K))
end