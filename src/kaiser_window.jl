"""
	kaiser_window(W::Int, H::Int, α::Float64)

Generate a kaiser window.

W - Width
H - Height
α - Parameters
"""
function kaiser_window(W::Int, H::Int, α::Float64)
	window_1D_W = zeros(Float64, W)
	window_1D_H = zeros(Float64, H)
	window_1D_W .= K.([0 : W - 1...], W, α)
	window_1D_H .= K.([0 : H - 1...], H, α)
	window_1D_H * window_1D_W'
end

function K(n::Int, k::Int, α::Float64)
	besseli0(α * sqrt(1 - (2 * n / (k - 1) - 1)^2)) / besseli0(α)
end

"""
	besseli0(x::Float64)

Modified Bessel function of the first kind of zero order.
"""
function besseli0(x::Float64)
	ax = abs(x)
	if ax < 3.75
		y = x / 3.75
		y = y * y
		ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
			+ y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))))
	else
		y = 3.75 / ax
		ans = (exp(ax) / sqrt(ax)) * (0.39894228 + y * (0.1328592e-1
			+ y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2
			+ y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1
			+ y * 0.392377e-2))))))))
	end
	return ans
end