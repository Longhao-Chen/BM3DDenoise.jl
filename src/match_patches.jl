"""
	match_patches(img::AbstractArray{Float64},
			refIndex::Array{CartesianIndex{2},2},
			patchSize::CartesianIndex{2},
			searchWindow::CartesianIndex{2},
			nMatch::Int64,
			distF::Function)

Full-search block matching algorithm for BM3D
"""
function match_patches(
	img::AbstractArray{Float64,2},
	refIndex::Array{CartesianIndex{2},2},
	patchSize::CartesianIndex{2},
	searchWindow::CartesianIndex{2},
	nMatch::Int64,
	distF::Function,
)

	refIndexH, refIndexW = size(refIndex)

	matches = Array{CartesianIndex{2}}(undef, refIndexH, refIndexW, nMatch)

	dist = Array{Float64}(undef, refIndexH, refIndexW, nMatch + 1)# dist[i, j, end] is the largest item in current dist[i, j, 1:end-1]
	dist .= typemax(Float64)

	maxDistMatch = Array{UInt}(undef, refIndexH, refIndexW)# maxDistMatch[i, j] is the position in matches[i, j, :] corresponding to the largest distance block
	maxDistMatch .= 1

	lockPool = Array{ReentrantLock}(undef, refIndexH, refIndexW)
	for i in 1:length(lockPool)
		lockPool[i] = ReentrantLock()
	end


	@inbounds Threads.@threads for j in 1:refIndexW
		for i in 1:refIndexH
			p = refIndex[i, j]
			iE = min(refIndexH, i + searchWindow[1])
			jE = min(refIndexW, j + searchWindow[2])
			imgRef = @view img[p:(p + patchSize - CartesianIndex(1, 1))]

			# In order to skip the current block and save unnecessary judgment, it is split into three parts.
			for sj in (j + 1):jE, si in (i + 1):iE
				@views weight = distF(
					imgRef,
					img[refIndex[si, sj]:(refIndex[
						si,
						sj,
					] + patchSize - CartesianIndex(1, 1))],
				)
				lock(lockPool[i, j]) do
					if weight < dist[i, j, nMatch + 1]
						maxDistImg = maxDistMatch[i, j]
						matches[i, j, maxDistImg] = CartesianIndex(si, sj)
						dist[i, j, maxDistImg] = weight
						dist[i, j, nMatch + 1], maxDistMatch[i, j] =
							findmax(dist[i, j, 1:nMatch])
					end
				end

				lock(lockPool[si, sj]) do
					if weight < dist[si, sj, nMatch + 1]
						maxDistImg = maxDistMatch[si, sj]
						matches[si, sj, maxDistImg] = CartesianIndex(i, j)
						dist[si, sj, maxDistImg] = weight
						dist[si, sj, nMatch + 1],
						maxDistMatch[si, sj] =
							findmax(dist[si, sj, 1:nMatch])
					end
				end
			end
			for si in (i + 1):iE
				weight = distF(
					imgRef,
					img[refIndex[si, j]:(refIndex[si, j] + patchSize - CartesianIndex(1, 1))],
				)
				lock(lockPool[i, j]) do
					if weight < dist[i, j, nMatch + 1]
						maxDistImg = maxDistMatch[i, j]
						matches[i, j, maxDistImg] = CartesianIndex(si, j)
						dist[i, j, maxDistImg] = weight
						dist[i, j, nMatch + 1], maxDistMatch[i, j] =
							findmax(dist[i, j, 1:nMatch])
					end
				end

				lock(lockPool[si, j]) do
					if weight < dist[si, j, nMatch + 1]
						maxDistImg = maxDistMatch[si, j]
						matches[si, j, maxDistImg] = CartesianIndex(i, j)
						dist[si, j, maxDistImg] = weight
						dist[si, j, nMatch + 1],
						maxDistMatch[si, j] =
							findmax(dist[si, j, 1:nMatch])
					end
				end
			end
			for sj in (j + 1):jE
				weight = distF(
					imgRef,
					img[refIndex[i, sj]:(refIndex[i, sj] + patchSize - CartesianIndex(1, 1))],
				)
				lock(lockPool[i, j]) do
					if weight < dist[i, j, nMatch + 1]
						maxDistImg = maxDistMatch[i, j]
						matches[i, j, maxDistImg] = CartesianIndex(i, sj)
						dist[i, j, maxDistImg] = weight
						dist[i, j, nMatch + 1], maxDistMatch[i, j] =
							findmax(dist[i, j, 1:nMatch])
					end
				end

				lock(lockPool[i, sj]) do
					if weight < dist[i, sj, nMatch + 1]
						maxDistImg = maxDistMatch[i, sj]
						matches[i, sj, maxDistImg] = CartesianIndex(i, j)
						dist[i, sj, maxDistImg] = weight
						dist[i, sj, nMatch + 1],
						maxDistMatch[i, sj] =
							findmax(dist[i, sj, 1:nMatch])
					end
				end
			end
		end
	end
	matches
end

# For color images
function match_patches(
	img::Array{Float64,3},
	refIndex::Array{CartesianIndex{2},2},
	patchSize::CartesianIndex{2},
	searchWindow::CartesianIndex{2},
	nMatch::Int64,
	distF::Function,
)
	match_patches(view(img, :, :, 1), refIndex, patchSize, searchWindow, nMatch, distF)
end
