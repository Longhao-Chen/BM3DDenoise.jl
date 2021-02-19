"""
	HardThresholding!(data::AbstractArray, λ::AbstractFloat)

Hard Thresholding function  
λ is threshold parameter
"""
function HardThresholding!(data::AbstractArray{T, N}, λ::AbstractFloat) where {T <:Number, N}
	i_end = length(data)
	@views @inbounds for i in 1:i_end
		if abs(data[i]) <= λ
			data[i] = 0
		end
	end
end