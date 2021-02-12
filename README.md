# BM3D.jl

An implementation of the BM3D denoising algorithm for Julia.

Some code reference: https://github.com/rcrandall/BM3D.jl

## Installation

Within Julia, use the [package manager][pkg]:
```julia
using Pkg
Pkg.add("https://github.com/Longhao-Chen/BM3D.jl")
```

## Usage

```julia
using BM3D
img = load("image path")
noise_img = Gary.(img)	# Now only supports grayscale images
denoise_img = bm3d(img, noise_variance)
```

## Todo

- [ ] Add color image support
- [ ] Performance optimization
- [ ] Add noise variance estimation

## Known issues

* Too much memory required for large pictures.