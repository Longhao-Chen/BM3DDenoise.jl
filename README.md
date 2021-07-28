# BM3D.jl
![](https://github.com/Longhao-Chen/BM3D.jl/workflows/Unit%20test/badge.svg)

An implementation of the BM3D(sparse 3D transform-domain collaborative filtering) denoising algorithm for Julia.

## Installation

Within Julia, use the `Pkg`:
```julia
using Pkg
Pkg.add(url = "https://github.com/Longhao-Chen/BM3D.jl.git")
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
using Downloads
noise_image_path = download("http://www.cs.tut.fi/~foi/GCF-BM3D/images/Lena512_noi_s10.png")
noise_variance = 10 / 255
img = load(noise_image_path)
denoise_img = bm3d(img, noise_variance)	# noise_variance noise_variance is the variance of the noise.
```

If you want to customize parameters, see:

```julia
?bm3d_config
```

## Noise estimator

For the image that does not know the noise variance, we can use [Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl/) to make noise estimation.

```julia
using Pkg; Pkg.add("Wavelets")
using Wavelets
using Images

noise_img = load("image_path")

# Convert to YCrCb color space for estimation
noise_img = YCbCr.(noise_img)
noise_img_y = channelview(img_y)[1,:,:]	# Extract y channel.
# normalization
noise_img_y .-= 16.
noise_img_y ./= (235. - 16.)

noise_variance = Threshold.noisest(noise_img_y)
```

# Known Issues
* 3D filtering is slow.

# Reference

> https://webpages.tuni.fi/foi/GCF-BM3D/  
> https://github.com/rcrandall/BM3D.jl  
> https://blog.csdn.net/qq_33552519/article/details/108632146  
> http://www.ipol.im/pub/art/2012/l-bm3d/article.pdf  
> https://blog.csdn.net/github_28260175/article/details/81101457  
> https://blog.csdn.net/qq_36955294/article/details/83443317
