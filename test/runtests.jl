using Test
using Images
import BM3D
include("PSNR.jl")

l1 = load("data/Lena512.png")
l2 = load("data/Lena512_noi_s100.png")

@testset "BM3D" begin
	l3 = BM3D.bm3d(l2, 0.39)
	@test PSNR(Float64.(l1), Float64.(l3)) > 23.4
end