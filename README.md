# BM3D.jl

An implementation of the BM3D denoising algorithm for Julia.

## Installation

Within Julia, use the [package manager][pkg]:
```julia
using Pkg
Pkg.add("https://github.com/Longhao-Chen/BM3D.jl.git")
```

## Usage

```julia
using BM3D
using Images
img = load("noise_image_path")
denoise_img = bm3d(img, noise_variance)	# noise_variance noise_variance is the variance of the noise.
```

### Example

```julia
using BM3D
using Images
noise_image_path = download("http://www.cs.tut.fi/~foi/GCF-BM3D/images/Lena512_noi_s10.png")
noise_variance = 10 / 255
img = load(noise_image_path)
denoise_img = bm3d(img, noise_variance)	# noise_variance noise_variance is the variance of the noise.
```

If you want to customize parameters, see:

```julia
?bm3d_config
```


## Todo

- [x] Add color image support
- [ ] Performance optimization
- [ ] Add noise variance estimation

## Reference

> https://github.com/rcrandall/BM3D.jl
> https://blog.csdn.net/qq_33552519/article/details/108632146
> http://www.ipol.im/pub/art/2012/l-bm3d/article.pdf
> https://blog.csdn.net/github_28260175/article/details/81101457
> https://blog.csdn.net/qq_36955294/article/details/83443317