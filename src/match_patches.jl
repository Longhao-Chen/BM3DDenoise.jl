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
	for i in length(lockPool)
		lockPool[i] = ReentrantLock()
	end


	@inbounds Threads.@threads for j = 1:refIndexW
		for i = 1:refIndexH
			p = @view refIndex[i, j]
			iE = min(refIndexH, i + searchWindow[1])
			jE = min(refIndexW, j + searchWindow[2])
			imgRef = @view img[p:p+patchSize]

			# In order to skip the current block and save unnecessary judgment, it is split into three parts.
			for sj = (j+1):jE, si = (i+1):iE
				@views weight = distF(
					imgRef,
					img[refIndex[si, sj]:(refIndex[si, sj]+patchSize)],
				)
				lock(lockPool[i, j]) do
					if weight < dist[i, j, nMatch+1]
						maxDistImg = maxDistMatch[i, j]
						matches[i, j, maxDistImg] = refIndex[si, sj]
						dist[i, j, maxDistImg] = weight
						dist[i, j, nMatch+1], maxDistMatch[i, j] =
							findmax(dist[i, j, 1:nMatch])
					end
				end

				lock(lockPool[si, sj]) do
					if weight < dist[si, sj, nMatch+1]
						maxDistImg = maxDistMatch[si, sj]
						matches[si, sj, maxDistImg] = p
						dist[si, sj, maxDistImg] = weight
						dist[si, sj, nMatch+1],
						maxDistMatch[si, sj] =
							findmax(dist[si, sj, 1:nMatch])
					end
				end
			end
			for si = (i+1):iE
				weight = distF(
					imgRef,
					img[refIndex[si, j]:(refIndex[si, j]+patchSize)],
				)
				lock(lockPool[i, j]) do
					if weight < dist[i, j, nMatch+1]
						maxDistImg = maxDistMatch[i, j]
						matches[i, j, maxDistImg] = refIndex[si, j]
						dist[i, j, maxDistImg] = weight
						dist[i, j, nMatch+1], maxDistMatch[i, j] =
							findmax(dist[i, j, 1:nMatch])
					end
				end

				lock(lockPool[si, j]) do
					if weight < dist[si, j, nMatch+1]
						maxDistImg = maxDistMatch[si, sj]
						matches[si, j, maxDistImg] = p
						dist[si, j, maxDistImg] = weight
						dist[si, j, nMatch+1], maxDistMatch[si, j] =
							findmax(dist[si, j, 1:nMatch])
					end
				end
			end
			for sj = (j+1):jE
				weight = distF(
					imgRef,
					img[refIndex[i, sj]:(refIndex[i, sj]+patchSize)],
				)
				lock(lockPool[i, j]) do
					if weight < dist[i, j, nMatch+1]
						maxDistImg = maxDistMatch[i, j]
						matches[i, j, maxDistImg] = refIndex[i, sj]
						dist[i, j, maxDistImg] = weight
						dist[i, j, nMatch+1], maxDistMatch[i, j] =
							findmax(dist[i, j, 1:nMatch])
					end
				end

				lock(lockPool[si, j]) do
					if weight < dist[i, sj, nMatch+1]
						maxDistImg = maxDistMatch[i, sj]
						matches[i, sj, maxDistImg] = p
						dist[i, sj, maxDistImg] = weight
						dist[i, sj, nMatch+1], maxDistMatch[i, sj] =
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
