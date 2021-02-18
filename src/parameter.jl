"""
Some runtime parameters. The default value is after the =.

## step 1:
```
bm3d_config.thr_patchSize = [8; 8]	size of patches;
bm3d_config.thr_stepSize = [3; 3]	In order to speed up the processing,the loop over the pixels of the imageis done with a stepp(integer) in row and column;
bm3d_config.thr_nBorder = [0; 0]	number of border pixels on each side in x and y;
bm3d_config.thr_searchWin = [39; 39]	search window size;
bm3d_config.thr_nMatch = 32	maximum number of similar patches kept;
bm3d_config.thr_threshSimilar = 2500/255^2	maximum thresholds for the distance between two similar patches;
bm3d_config.thr_thresh3D = 2.7	3D group threshold.
```

## step2:
```
bm3d_config.wie_patchSize = [8; 8]	size of patches;
bm3d_config.wie_stepSize = [3; 3]	In order to speed up the processing,the loop over the pixels of the imageis done with a stepp(integer) in row and column;
bm3d_config.wie_nBorder = [0; 0]	number of border pixels on each side in x and y;
bm3d_config.wie_searchWin = [39; 39]	search window size;
bm3d_config.wie_nMatch = 32	maximum number of similar patches kept;
bm3d_config.wie_threshSimilar = 400/255^2	maximum thresholds for the distance between two similar patches;
bm3d_config.wie_thresh3D = 2.7	3D group threshold.
```
"""
mutable struct bm3d_config
	thr_patchSize::Array{Int, 1}
	thr_stepSize::Array{Int, 1}
	thr_nBorder::Array{Int, 1}
	thr_searchWin::Array{Int, 1}
	thr_nMatch::Int
	thr_threshSimilar::Float64
	thr_thresh3D::Float64

	wie_patchSize::Array{Int, 1}
	wie_stepSize::Array{Int, 1}
	wie_nBorder::Array{Int, 1}
	wie_searchWin::Array{Int, 1}
	wie_nMatch::Int
	wie_threshSimilar::Float64
	wie_thresh3D::Float64

	bm3d_config() = new(
		[8; 8],
		[3; 3],
		[0; 0],
		[39; 39],
		16,
		2500/255^2,
		2.7,

		[8; 8],
		[3; 3],
		[0; 0],
		[39; 39],
		32,
		400/255^2,
		2.7
	)
end

function print(config::bm3d_config)
	@info "========================================"
	@info "|       bm3d config                    |"
	@info "========================================"
	@info "| thr_patchSize: $(config.thr_patchSize)"
	@info "| thr_stepSize: $(config.thr_stepSize)"
	@info "| thr_nBorder: $(config.thr_nBorder)"
	@info "| thr_searchWin: $(config.thr_searchWin)"
	@info "| thr_nMatch: $(config.thr_nMatch)"
	@info "| thr_threshSimilar: $(config.thr_threshSimilar)"
	@info "| thr_thresh3D: $(config.thr_thresh3D)"
	@info "| wie_patchSize: $(config.wie_patchSize)"
	@info "| wie_stepSize: $(config.wie_stepSize)"
	@info "| wie_nBorder: $(config.wie_nBorder)"
	@info "| wie_searchWin: $(config.wie_searchWin)"
	@info "| wie_nMatch: $(config.wie_nMatch)"
	@info "| wie_threshSimilar: $(config.wie_threshSimilar)"
	@info "| wie_thresh3D: $(config.wie_thresh3D)"
	@info "========================================"
end