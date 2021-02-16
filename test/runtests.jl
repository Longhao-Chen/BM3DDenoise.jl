using Test
using Images
import BM3D
include("PSNR.jl")

σ = [5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 75, 80, 90, 100]
PSNR_data = [38, 35, 34, 33, 32, 31, 30, 29, 29, 28, 27, 27, 26, 26, 25]

# download file
if !isfile(joinpath(tempdir(), "Lena512.png"))
	println("Downloading Lena512.png")
	download("http://www.cs.tut.fi/~foi/GCF-BM3D/images/Lena512.png", joinpath(tempdir(), "Lena512.png"))
end
Threads.@threads for i in σ
	if !isfile(joinpath(tempdir(), "Lena512_noi_s" * string(i) * ".png"))
		println("Downloading Lena512_noi_s" * string(i) * ".png")
		download("http://www.cs.tut.fi/~foi/GCF-BM3D/images/Lena512_noi_s" * string(i) * ".png",
			joinpath(tempdir(), "Lena512_noi_s" * string(i) * ".png"))
	end
end


Lena = load(joinpath(tempdir(), "Lena512.png"))
@testset "BM3D" begin
	for i = 1:length(σ)
		@testset "σ=$(σ[i])" begin
			img = load(joinpath(tempdir(), "Lena512_noi_s" * string(σ[i]) * ".png"))
			@test PSNR(Float64.(Lena), Float64.(BM3D.bm3d(img, σ[i] / 255))) >= PSNR_data[i]
		end
	end
end