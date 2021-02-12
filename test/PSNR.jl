import Statistics
using Strided

"""
	PSNR(img1::AbstractArray, img2::AbstractArray; max_value::Float64=1.)

Compute Peak signal-to-noise ratio of img1 and img2.
max_value: The maximum pixel value of the image. For floating-point data, the maximum pixel value is 1.
"""
function PSNR(img1::AbstractArray, img2::AbstractArray; max_value::Float64=1.)
	10 * log10(max_value^2 / MSE(img1, img2))
end

"""
	MSE(I::AbstractArray, K::AbstractArray)

Calculation of mean square error.
"""
function MSE(I::AbstractArray, K::AbstractArray)
	@strided Statistics.mean(abs2.(I .- K))
end