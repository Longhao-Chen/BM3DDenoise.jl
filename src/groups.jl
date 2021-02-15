"""
	form_groups(img::Matrix{Float64},
		matchTable::Array{Float64,4},
		Ilist::Vector{Int64},
		Jlist::Vector{Int64},
		patchSize::Vector{Int64})

Forward BM3D groupings (full transform... inefficient!)
"""
function form_groups(img::Matrix{Float64},
			matchTable::Array{Float64, 4},
			Ilist::Vector{Int64},
			Jlist::Vector{Int64},
			patchSize::Vector{Int64})

	(t,Nmatch,N1,N2) = size(matchTable)

	G3D = zeros(Float64, Nmatch+1, patchSize[1], patchSize[2], N1, N2)

	# Form table of 3D groups
	image_to_groups!(img, G3D, matchTable, Ilist, Jlist, patchSize)

	# Apply 3D DCT on groups. To prevent OutOfMemoryError, we will transform it separately
	@views @inbounds for n1 in 1:N1
		for n2 in 1:N2
			FFTW.dct!(G3D[:, :, :, n1, n2])
		end
	end

	return G3D

end

"""
	invert_groups(imgSize::Vector{Int64},
					G3D::Array{Float64,5},
					matchTable::Array{Float64,4},
					Ilist::Vector{Int64},
					Jlist::Vector{Int64},
					patchSize::Vector{Int64})

Inverse BM3D groupings
"""
function invert_groups(imgSize::Vector{Int64},
				G3D::Array{Float64, 5},
				matchTable::Array{Float64, 4},
				Ilist::Vector{Int64},
				Jlist::Vector{Int64},
				patchSize::Vector{Int64})

	(t, Nmatch, N1, N2) = size(matchTable)

	# Allocate image and weight table
	img = zeros(Float64, imgSize[1], imgSize[2])

	# Apply inverse 3D DCT on groups. To prevent OutOfMemoryError, we will transform it separately
	@views @inbounds for n1 in 1:N1
		for n2 in 1:N2
			FFTW.idct!(G3D[:, :, :, n1, n2])
		end
	end

	groups_to_image!(img, G3D, matchTable, Ilist, Jlist, patchSize)

	return img

end

"""
	groups_to_image!(img::Matrix{Float64},
					G3D::Array{Float64,5},
					matchTable::Array{Float64,4},
					Ilist::Vector{Int64},
					Jlist::Vector{Int64},
					patchSize::Vector{Int64})

Return filtered patches to their place in the image
"""
function groups_to_image!(img::Matrix{Float64},
				G3D::Array{Float64,5},
				matchTable::Array{Float64,4},
				Ilist::Vector{Int64},
				Jlist::Vector{Int64},
				patchSize::Vector{Int64})

	Nmatch = size(matchTable,2)

	@inbounds @views Base.Threads.@threads for j1 = 1:length(Jlist)
		for i1 = 1:length(Ilist)
			img[Ilist[i1]:(patchSize[1]-1+Ilist[i1]), Jlist[j1]:(patchSize[2]-1+Jlist[j1])] .+= G3D[1, 1:patchSize[1], 1:patchSize[2], i1, j1]

			for k = 1:Nmatch
				i2 = i1 + Int(matchTable[1, k, i1, j1])
				j2 = j1 + Int(matchTable[2, k, i1, j1])
				img[Ilist[i2]:(patchSize[1]-1+Ilist[i2]), Jlist[j2]:(patchSize[2]-1+Jlist[j2])] .+= G3D[k + 1, 1:patchSize[1], 1:patchSize[2], i1, j1]
			end

		end
	end
end

function image_to_groups!(img::Matrix{Float64},
				G3D::Array{Float64,5},
				matchTable::Array{Float64,4},
				Ilist::Vector{Int64},
				Jlist::Vector{Int64},
				patchSize::Vector{Int64})

	Nmatch = size(matchTable,2)

	@inbounds @views Base.Threads.@threads for j1 = 1:length(Jlist)
		for i1 = 1:length(Ilist)
			G3D[1, 1:patchSize[1], 1:patchSize[2], i1, j1] .= img[Ilist[i1]:(patchSize[1]-1+Ilist[i1]), Jlist[j1]:(patchSize[2]-1+Jlist[j1])]

			for k = 1:Nmatch
				i2 = i1 + Int(matchTable[1,k,i1,j1])
				j2 = j1 + Int(matchTable[2,k,i1,j1])
				G3D[k + 1, 1:patchSize[1], 1:patchSize[2], i1, j1] .= img[Ilist[i2]:(patchSize[1]-1+Ilist[i2]), Jlist[j2]:(patchSize[2]-1+Jlist[j2])]
			end

		end
	end
end