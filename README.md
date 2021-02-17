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
using Images
img = load("noise_image_path")
noise_img = Gray.(img)	# Now only supports grayscale images
denoise_img = bm3d(noise_img, noise_variance)	# noise_variance noise_variance is the variance of the noise.
```

### Example

```julia
using BM3D
using Images
noise_image_path = download("http://www.cs.tut.fi/~foi/GCF-BM3D/images/Lena512_noi_s100.png")
noise_variance = 100 / 255
img = load(noise_image_path)
noise_img = Gray.(img)	# Now only supports grayscale images
denoise_img = bm3d(noise_img, noise_variance)	# noise_variance noise_variance is the variance of the noise.
```

If you want to customize parameters, see:

```julia
?bm3d_config
```


## Todo

- [ ] Add color image support
- [ ] Performance optimization
- [ ] Add noise variance estimation