"""
	match_patches(img::Matrix{Float64},
				Ilist::Vector{Int64},
				Jlist::Vector{Int64},
				patchSize::Vector{Int64},
				searchWin::Vector{Int64},
				nMatch::Int64)

Full-search block matching algorithm for BM3D
"""
function match_patches(img::Matrix{Float64},
			Ilist::Vector{Int64},
			Jlist::Vector{Int64},
			patchSize::Vector{Int64},
			searchWin::Vector{Int64},
			nMatch::Int64)

	# Dimensions of patch table
	N1 = length(Ilist)
	N2 = length(Jlist)

	# Allocate 
	matchTable = zeros(Float64,(3,nMatch,N1,N2))
	matchMaxTable = ones(Float64,(2,N1,N2))
	matchTable[3,:,:,:] .= typemax(Float64)
	matchMaxTable[2,:,:] .= typemax(Float64)

	@inbounds @views Base.Threads.@threads for j1 = 1:N2
		for i1 = 1:N1
			for j2 = j1:minimum([N2;j1 + searchWin[2]])
			# Lower bound on columns
			LB = (j1 == j2) ? (i1+1) : (i1 - searchWin[1])
			for i2 = maximum([1,LB]):minimum([i1+searchWin[1];N1])

				d2 = norm(img[Ilist[i1]:(patchSize[1]-1+Ilist[i1]), Jlist[j1]:(patchSize[2]-1+Jlist[j1])]
					.- img[Ilist[i2]:(patchSize[1]-1+Ilist[i2]), Jlist[j2]:(patchSize[2]-1+Jlist[j2])])^2
				d2 /= prod(patchSize)

				# Check current maximum for patch (i1,j1)
				if (d2 < matchMaxTable[2,i1,j1])
					kmatch = Int(matchMaxTable[1,i1,j1])
					matchTable[1, kmatch, i1, j1] = i2 - i1
					matchTable[2, kmatch, i1, j1] = j2 - j1
					matchTable[3, kmatch, i1, j1] = d2

					(tmp2,tmp1) = findmax(matchTable[3,:,i1,j1])
					matchMaxTable[1,i1,j1] = tmp1
					matchMaxTable[2,i1,j1] = tmp2
				end

				# Check current maximum for patch (i2,j2)
				if (d2 < matchMaxTable[2,i2,j2])
					kmatch = Int(matchMaxTable[1,i2,j2])
					matchTable[1, kmatch, i2, j2] = i1 - i2
					matchTable[2, kmatch, i2, j2] = j1 - j2
					matchTable[3, kmatch, i2, j2] = d2

					(tmp2,tmp1) = findmax(matchTable[3,:,i2,j2])
					matchMaxTable[1,i2,j2] = tmp1
					matchMaxTable[2,i2,j2] = tmp2
				end

				end
			end
		end
	end

	return matchTable[1:2, :, :, :]
end