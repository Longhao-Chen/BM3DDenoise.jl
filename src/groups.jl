"""
	form_group!(G3D::Array{Float64, 3},
					img::Matrix{Float64},
					matchTable::Array{Float64, 4},
					vaild_match::BitArray,
					Ilist::Vector{Int64},
					Jlist::Vector{Int64},
					patchSize::Vector{Int64},
					position::Tuple{Int, Int})

Forward BM3D groupings
"""
function form_group!(G3D::AbstractArray{Float64, 3},
			img::AbstractArray{Float64, 2},
			matchTable::Array{Float64, 4},
			vaild_match::BitArray,
			Ilist::Vector{Int64},
			Jlist::Vector{Int64},
			patchSize::Vector{Int64},
			position::Tuple{Int, Int})

	# Form table of 3D groups
	image_to_group!(img, G3D, matchTable, vaild_match, Ilist, Jlist, patchSize, position)

	# Apply 3D DCT on groups.
	FFTW.dct!(G3D)

end

"""
	invert_group!(img::AbstractArray{Float64, 2},
				G3D::AbstractArray{Float64, 3},
				matchTable::Array{Float64, 4},
				vaild_match::BitArray,
				Ilist::Vector{Int64},
				Jlist::Vector{Int64},
				patchSize::Vector{Int64},
				position::Tuple{Int, Int})

Inverse BM3D groupings
"""
function invert_group!(img::AbstractArray{Float64, 2},
				G3D::AbstractArray{Float64, 3},
				matchTable::Array{Float64, 4},
				vaild_match::BitArray,
				Ilist::Vector{Int64},
				Jlist::Vector{Int64},
				patchSize::Vector{Int64},
				position::Tuple{Int, Int})

	(t, Nmatch) = size(matchTable)

	# Apply inverse 3D DCT on groups
	FFTW.idct!(G3D)

	group_to_image!(img, G3D, matchTable, vaild_match, Ilist, Jlist, patchSize, position)

end

"""
	group_to_image!(img::AbstractArray{Float64, 2},
				G3D::AbstractArray{Float64, 3},
				matchTable::Array{Float64, 4},
				vaild_match::BitArray,
				Ilist::Vector{Int64},
				Jlist::Vector{Int64},
				patchSize::Vector{Int64},
				position::Tuple{Int, Int})

group to image
"""
function group_to_image!(img::AbstractArray{Float64, 2},
				G3D::AbstractArray{Float64, 3},
				matchTable::Array{Float64, 4},
				vaild_match::BitArray,
				Ilist::Vector{Int64},
				Jlist::Vector{Int64},
				patchSize::Vector{Int64},
				position::Tuple{Int, Int})

	Nmatch = size(matchTable,2)
	(i1, j1) = position

	@strided img[Ilist[i1]:(patchSize[1] - 1 + Ilist[i1]), Jlist[j1]:(patchSize[2] - 1 + Jlist[j1])] .+= G3D[1, 1:patchSize[1], 1:patchSize[2]]

	@inbounds @views for k = 1:Nmatch
		if vaild_match[k, i1, j1]
			i2 = i1 + Int(matchTable[1, k, i1, j1])
			j2 = j1 + Int(matchTable[2, k, i1, j1])
			img[Ilist[i2]:(patchSize[1] - 1 + Ilist[i2]), Jlist[j2]:(patchSize[2] - 1 + Jlist[j2])] .+= G3D[k + 1, 1:patchSize[1], 1:patchSize[2]]
		else
			continue
		end
	end
end

"""
	group_to_image!(img::AbstractArray{Float64, 2},
				W::Float64,
				matchTable::Array{Float64, 4},
				vaild_match::BitArray,
				Ilist::Vector{Int64},
				Jlist::Vector{Int64},
				patchSize::Vector{Int64},
				position::Tuple{Int, Int})
"""
function group_to_image!(img::AbstractArray{Float64, 2},
				W::Float64,
				matchTable::Array{Float64, 4},
				vaild_match::BitArray,
				Ilist::Vector{Int64},
				Jlist::Vector{Int64},
				patchSize::Vector{Int64},
				position::Tuple{Int, Int})

	Nmatch = size(matchTable,2)
	(i1, j1) = position

	@strided img[Ilist[i1]:(patchSize[1] - 1 + Ilist[i1]), Jlist[j1]:(patchSize[2] - 1 + Jlist[j1])] .+= W

	@inbounds @views for k = 1:Nmatch
		if vaild_match[k, i1, j1]
			i2 = i1 + Int(matchTable[1, k, i1, j1])
			j2 = j1 + Int(matchTable[2, k, i1, j1])
			img[Ilist[i2]:(patchSize[1] - 1 + Ilist[i2]), Jlist[j2]:(patchSize[2] - 1 + Jlist[j2])] .+= W
		else
			continue
		end
	end
end

"""
	image_to_group!(img::AbstractArray{Float64, 2},
					G3D::AbstractArray{Float64, 3},
					matchTable::Array{Float64, 4},
					vaild_match::BitArray,
					Ilist::Vector{Int64},
					Jlist::Vector{Int64},
					patchSize::Vector{Int64},
					position::Tuple{Int, Int})

"""
function image_to_group!(img::AbstractArray{Float64, 2},
				G3D::AbstractArray{Float64, 3},
				matchTable::Array{Float64, 4},
				vaild_match::BitArray,
				Ilist::Vector{Int64},
				Jlist::Vector{Int64},
				patchSize::Vector{Int64},
				position::Tuple{Int, Int})

	Nmatch = size(matchTable,2)
	i1, j1 = position

	G3D[1, 1:patchSize[1], 1:patchSize[2]] .= img[Ilist[i1]:(patchSize[1] - 1 + Ilist[i1]), Jlist[j1]:(patchSize[2] - 1 + Jlist[j1])]

	@views @inbounds for k = 1:Nmatch
		if vaild_match[k, i1, j1]
			i2 = i1 + Int(matchTable[1, k, i1, j1])
			j2 = j1 + Int(matchTable[2, k, i1, j1])
			G3D[k + 1, 1:patchSize[1], 1:patchSize[2]] .= img[Ilist[i2]:(patchSize[1] - 1 + Ilist[i2]), Jlist[j2]:(patchSize[2] - 1 + Jlist[j2])]
		else
			# To prevent the impact on nnz() and norm().
			G3D[k + 1, 1:patchSize[1], 1:patchSize[2]] .= 0
		end
	end
end