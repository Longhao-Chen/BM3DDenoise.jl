include("Transform/dct.jl")
"""
Some runtime parameters. The default value is after the `=`

## step 1:
```julia
bm3d_config.thr_patchSize = CartesianIndex(8, 8)	size of patches;
bm3d_config.thr_searchStride = CartesianIndex(3, 3)	In order to speed up the processing,the loop over the pixels of the imageis done with a stepp(integer) in row and column;
bm3d_config.thr_searchWin = CartesianIndex(19, 19)	search window size;
bm3d_config.thr_nMatch = 31	maximum number of similar patches kept;
bm3d_config.thr_thresh3D = 2.7	3D group hreshold;
bm3d_config.thr_distFunction = Euclidean2	A function that measures the distance between blocks;
bm3d_config.thr_group_transform! = dct!	A function that 3D transform in group
bm3d_config.thr_group_itransform! = idct!	A function that 3D inverse transform in group;
```

## step2:
```julia
bm3d_config.wie_patchSize = CartesianIndex(8, 8)	size of patches;
bm3d_config.wie_searchStride = CartesianIndex(3, 3)	In order to speed up the processing,the loop over the pixels of the imageis done with a stepp(integer) in row and column;
bm3d_config.wie_searchWin = CartesianIndex(11, 11)	search window size;
bm3d_config.wie_nMatch = 15	maximum number of similar patches kept;
bm3d_config.wie_distFunction = Euclidean2	A function that measures the distance between blocks;
bm3d_config.wie_group_transform! = dct!	A function that 3D transform in group;
bm3d_config.wie_group_itransform! = idct!	A function that 3D inverse transform in group;
```

## others
bm3d_config.show_info = false	show some running message

note:

The 3D group data struct:
```
G3D[k, :, :]	# k-th matching block.
size(G3D)	# (nMatch, patchSize[1], patchSize[2])
```
"""
mutable struct bm3d_config
	show_info::Bool

	thr_patchSize::CartesianIndex{2}
	thr_searchStride::CartesianIndex{2}
	thr_searchWin::CartesianIndex{2}
	thr_nMatch::Int
	thr_thresh3D::Float64
	thr_distFunction::Function
	thr_group_transform!::Function
	thr_group_itransform!::Function

	wie_patchSize::CartesianIndex{2}
	wie_searchStride::CartesianIndex{2}
	wie_searchWin::CartesianIndex{2}
	wie_nMatch::Int
	wie_distFunction::Function
	wie_group_transform!::Function
	wie_group_itransform!::Function

	bm3d_config() = new(
		false,
		CartesianIndex(8, 8),
		CartesianIndex(2, 2),
		CartesianIndex(24, 24),
		31,
		2.6,
		Euclidean2,
		dct!,
		idct!,

		CartesianIndex(8, 8),
		CartesianIndex(2, 2),
		CartesianIndex(19, 19),
		15,
		Euclidean2,
		dct!,
		idct!,
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
	@info "| thr_group_transform!: $(config.thr_group_transform!)"
	@info "| thr_group_itransform!: $(config.thr_group_itransform!)"
	@info "| wie_patchSize: $(config.wie_patchSize)"
	@info "| wie_searchStride: $(config.wie_searchStride)"
	@info "| wie_searchWin: $(config.wie_searchWin)"
	@info "| wie_nMatch: $(config.wie_nMatch)"
	@info "| wie_distFunction: $(config.wie_distFunction)"
	@info "| wie_group_transform!: $(config.wie_group_transform!)"
	@info "| wie_group_itransform!: $(config.wie_group_itransform!)"
	@info "========================================"
end
