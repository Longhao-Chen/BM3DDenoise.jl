"""
	HardThresholding!(data::AbstractArray, λ::AbstractFloat)

Hard Thresholding function  
λ is threshold parameter
"""
function HardThresholding!(data::AbstractArray{T, N}, λ::AbstractFloat) where {T <:Number, N}
	@views @inbounds for i in 1:length(data)
		if abs(data[i]) <= λ
			data[i] = 0
		end
	end
end