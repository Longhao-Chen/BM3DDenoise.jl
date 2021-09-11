include("../../src/Transform/dct.jl")

@testset "Transform/dct" begin
	a=rand(10, 20, 30)
	b=copy(a)
	dct!(b)
	idct!(b)
	@test isapprox(a, b)
	a=rand(12, 22, 32)
	b=copy(a)
	dct!(b)
	idct!(b)
	@test isapprox(a, b)
end