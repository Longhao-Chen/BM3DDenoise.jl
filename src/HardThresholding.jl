"""
	HardThresholding!(data::AbstractArray, λ::AbstractFloat)

Hard Thresholding function  
λ is threshold parameter
Return the number of non-zero elements in the array
"""
function HardThresholding!(data::AbstractArray{T, N}, λ::AbstractFloat) where {T <:Number, N}
	i_end = length(data)
	nnz = i_end
	@views @inbounds for i in 1:i_end
		if abs(data[i]) <= λ
			data[i] = zero(T)
			nnz -= 1
		end
	end
	return nnz
end