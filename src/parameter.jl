"""
Some runtime parameters. The default value is after the =.

## step 1:
```
bm3d_config.thr_patchSize = [8; 8]	size of patches;
bm3d_config.thr_stepSize = [3; 3]	In order to speed up the processing,the loop over the pixels of the imageis done with a stepp(integer) in row and column;
bm3d_config.thr_nBorder = [0; 0]	number of border pixels on each side in x and y;
bm3d_config.thr_searchWin = [19; 19]	search window size;
bm3d_config.thr_nMatch = 31	maximum number of similar patches kept;
bm3d_config.thr_thresh3D = 2.7	3D group hreshold.
```

## step2:
```
bm3d_config.wie_patchSize = [8; 8]	size of patches;
bm3d_config.wie_stepSize = [3; 3]	In order to speed up the processing,the loop over the pixels of the imageis done with a stepp(integer) in row and column;
bm3d_config.wie_nBorder = [0; 0]	number of border pixels on each side in x and y;
bm3d_config.wie_searchWin = [11; 11]	search window size;
bm3d_config.wie_nMatch = 15	maximum number of similar patches kept;
bm3d_config.wie_thresh3D = 2.7	3D group hreshold.
```
"""
mutable struct bm3d_config
	thr_patchSize::Array{Int, 1}
	thr_stepSize::Array{Int, 1}
	thr_nBorder::Array{Int, 1}
	thr_searchWin::Array{Int, 1}
	thr_nMatch::Int
	thr_thresh3D::Float64

	wie_patchSize::Array{Int, 1}
	wie_stepSize::Array{Int, 1}
	wie_nBorder::Array{Int, 1}
	wie_searchWin::Array{Int, 1}
	wie_nMatch::Int
	wie_thresh3D::Float64

	bm3d_config() = new(
		[8; 8],
		[3; 3],
		[0; 0],
		[19; 19],
		31,
		2.7,

		[8; 8],
		[3; 3],
		[0; 0],
		[11; 11],
		15,
		2.7
	)
end