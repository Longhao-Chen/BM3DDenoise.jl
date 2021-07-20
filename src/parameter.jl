"""
Some runtime parameters. The default value is after the `=`

## step 1:
```
bm3d_config.thr_patchSize = CartesianIndex(8, 8)	size of patches;
bm3d_config.thr_searchStride = CartesianIndex(3, 3)	In order to speed up the processing,the loop over the pixels of the imageis done with a stepp(integer) in row and column;
bm3d_config.thr_searchWin = CartesianIndex(19, 19)	search window size;
bm3d_config.thr_nMatch = 31	maximum number of similar patches kept;
bm3d_config.thr_thresh3D = 2.7	3D group hreshold;
bm3d_config.thr_distFunction = Euclidean2	A function that measures the distance between blocks;
bm3d_config.thr_transform_1D! = FFTW.dct!	A function that 1D transform in group;
bm3d_config.thr_transform_2D! = FFTW.dct!	A function that 2D transform in group;
bm3d_config.thr_itransform_1D! = FFTW.idct!	A function that 1D inverse transform in group;
bm3d_config.thr_itransform_2D! = FFTW.idct!	A function that 2D inverse transform in group;
```

## step2:
```
bm3d_config.wie_patchSize = CartesianIndex(8, 8)	size of patches;
bm3d_config.wie_searchStride = CartesianIndex(3, 3)	In order to speed up the processing,the loop over the pixels of the imageis done with a stepp(integer) in row and column;
bm3d_config.wie_searchWin = CartesianIndex(11, 11)	search window size;
bm3d_config.wie_nMatch = 15	maximum number of similar patches kept;
bm3d_config.wie_distFunction = Euclidean2	A function that measures the distance between blocks;
bm3d_config.wie_transform_1D! = FFTW.dct!	A function that 1D transform in group;
bm3d_config.wie_transform_2D! = FFTW.dct!	A function that 2D transform in group;
bm3d_config.wie_itransform_1D! = FFTW.idct!	A function that 1D inverse transform in group;
bm3d_config.wie_itransform_2D! = FFTW.idct!	A function that 2D inverse transform in group;
```
"""
mutable struct bm3d_config
	thr_patchSize::CartesianIndex{2}
	thr_searchStride::CartesianIndex{2}
	thr_searchWin::CartesianIndex{2}
	thr_nMatch::Int
	thr_thresh3D::Float64
	thr_distFunction::Function
	thr_transform_1D!::Function
	thr_transform_2D!::Function
	thr_itransform_1D!::Function
	thr_itransform_2D!::Function

	wie_patchSize::CartesianIndex{2}
	wie_searchStride::CartesianIndex{2}
	wie_searchWin::CartesianIndex{2}
	wie_nMatch::Int
	wie_distFunction::Function
	wie_transform_1D!::Function
	wie_transform_2D!::Function
	wie_itransform_1D!::Function
	wie_itransform_2D!::Function

	bm3d_config() = new(
		CartesianIndex(8, 8),
		CartesianIndex(3, 3),
		CartesianIndex(19, 19),
		31,
		2.7,
		Euclidean2,
		FFTW.dct!,
		FFTW.dct!,
		FFTW.idct!,
		FFTW.idct!,

		CartesianIndex(8, 8),
		CartesianIndex(3, 3),
		CartesianIndex(11, 11),
		15,
		Euclidean2,
		FFTW.dct!,
		FFTW.dct!,
		FFTW.idct!,
		FFTW.idct!,
	)
end

"""
Euclidean2(a, b)

(Euclidean distance)^2
"""
function Euclidean2(a::AbstractArray{T, N}, b::AbstractArray{T, N}) where {T, N}
	norm(a .- b)^2
end

function print(config::bm3d_config)
	@info "========================================"
	@info "|       bm3d config                    |"
	@info "========================================"
	@info "| thr_patchSize: $(config.thr_patchSize)"
	@info "| thr_searchStride: $(config.thr_searchStride)"
	@info "| thr_searchWin: $(config.thr_searchWin)"
	@info "| thr_nMatch: $(config.thr_nMatch)"
	@info "| thr_thresh3D: $(config.thr_thresh3D)"
	@info "| thr_distFunction: $(config.thr_distFunction)"
	@info "| thr_transform_1D!: $(config.thr_transform_1D!)"
	@info "| thr_transform_2D!: $(config.thr_transform_2D!)"
	@info "| thr_itransform_1D!: $(config.thr_itransform_1D!)"
	@info "| thr_itransform_2D!: $(config.thr_itransform_2D!)"
	@info "| wie_patchSize: $(config.wie_patchSize)"
	@info "| wie_searchStride: $(config.wie_searchStride)"
	@info "| wie_searchWin: $(config.wie_searchWin)"
	@info "| wie_nMatch: $(config.wie_nMatch)"
	@info "| wie_distFunction: $(config.wie_distFunction)"
	@info "| wie_transform_1D!: $(config.wie_transform_1D!)"
	@info "| wie_transform_2D!: $(config.wie_transform_2D!)"
	@info "| wie_itransform_1D!: $(config.wie_itransform_1D!)"
	@info "| wie_itransform_2D!: $(config.wie_itransform_2D!)"
	@info "========================================"
end
