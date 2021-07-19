"""
	form_group!(G3D::AbstractArray{Float64, 3},
			img::AbstractArray{Float64, 2},
			matchTable::Array{CartesianIndex{2}, 3},
			refIndex::Array{CartesianIndex{2},2},
			patchSize::CartesianIndex{2},
			position::CartesianIndex{2},
			transform_1D!::Function,
			transform_2D!::Function)

Forward BM3D groupings
"""
function form_group!(
	G3D::AbstractArray{Float64,3},
	img::AbstractArray{Float64,2},
	matchTable::Array{CartesianIndex{2},3},
	refIndex::Array{CartesianIndex{2},2},
	patchSize::CartesianIndex{2},
	position::CartesianIndex{2},
	transform_1D!::Function,
	transform_2D!::Function,
)

	# Form table of 3D groups
	image_to_group!(img, G3D, matchTable, refIndex, patchSize, position)

	# Apply 2D transform
	@inbounds @views Threads.@threads for k in 1:size(G3D, 1)
		transform_2D!(G3D[k, :, :])
	end

	# Apply 1D transform
	@inbounds @views Threads.@threads for j in 1:patchSize[2]
		for i in 1:patchSize[1]
			transform_1D!(G3D[:, i, j])
		end
	end
end

"""
	invert_group!(img::AbstractArray{Float64, 2},
				G3D::AbstractArray{Float64, 3},
				matchTable::Array{CartesianIndex{2},3},
				refIndex::Array{CartesianIndex{2},2},
				patchSize::CartesianIndex{2},
				position::CartesianIndex{2},
				itransform_1D!::Function,
				itransform_2D!::Function)

Inverse BM3D groupings
"""
function invert_group!(
	img::AbstractArray{Float64,2},
	G3D::AbstractArray{Float64,3},
	matchTable::Array{CartesianIndex{2},3},
	refIndex::Array{CartesianIndex{2},2},
	patchSize::CartesianIndex{2},
	position::CartesianIndex{2},
	itransform_1D!::Function,
	itransform_2D!::Function,
)

	# Apply inverse 1D transform
	@inbounds @views Threads.@threads for j in 1:patchSize[2]
		for i in 1:patchSize[1]
			itransform_1D!(G3D[:, i, j])
		end
	end

	# Apply inverse 2D transform
	@inbounds @views Threads.@threads for k in 1:size(G3D, 1)
		itransform_2D!(G3D[k, :, :])
	end

	group_to_image!(img, G3D, matchTable, refIndex, patchSize, position)

end

"""
	group_to_image!(img::AbstractArray{Float64, 2},
			G3D::AbstractArray{Float64, 3},
			matchTable::Array{CartesianIndex{2},3},
			refIndex::Array{CartesianIndex{2},2},
			patchSize::CartesianIndex{2},
			position::CartesianIndex{2})

group to image
"""
function group_to_image!(
	img::AbstractArray{Float64,2},
	G3D::AbstractArray{Float64,3},
	matchTable::Array{CartesianIndex{2},3},
	refIndex::Array{CartesianIndex{2},2},
	patchSize::CartesianIndex{2},
	position::CartesianIndex{2},
)

	Nmatch = size(matchTable, 3)

	@views img[refIndex[position]:(patchSize - CartesianIndex(
		1,
		1,
	) + refIndex[position])] .+= G3D[1, 1:patchSize[1], 1:patchSize[2]]

	@inbounds @views for k in 1:Nmatch
		position2 = position + matchTable[position, k]
		img[refIndex[position2]:(patchSize - CartesianIndex(
			1,
			1,
		) + refIndex[position2])] .+= G3D[k + 1, 1:patchSize[1], 1:patchSize[2]]
	end
end

"""
	group_to_image!(img::AbstractArray{Float64, 2},
			W::Float64,
			matchTable::Array{CartesianIndex{2},3},
			refIndex::Array{CartesianIndex{2},2},
			patchSize::CartesianIndex{2},
			position::CartesianIndex{2})
"""
function group_to_image!(
	img::AbstractArray{Float64,2},
	W::Float64,
	matchTable::Array{CartesianIndex{2},3},
	refIndex::Array{CartesianIndex{2},2},
	patchSize::CartesianIndex{2},
	position::CartesianIndex{2},
)

	nMatch = size(matchTable, 3)

	@views img[refIndex[position]:(patchSize - CartesianIndex(
		1,
		1,
	) + refIndex[position])] .+= W

	@inbounds @views for k in 1:nMatch
		position2 = position + matchTable[position, k]
		img[refIndex[position2]:(patchSize - CartesianIndex(
			1,
			1,
		) + refIndex[position2])] .+= W
	end
end

"""
	image_to_group!(img::AbstractArray{Float64, 2},
			G3D::AbstractArray{Float64, 3},
			matchTable::Array{CartesianIndex{2},3},
			refIndex::Array{CartesianIndex{2},2},
			patchSize::CartesianIndex{2},
			position::CartesianIndex{2})

"""
function image_to_group!(
	img::AbstractArray{Float64,2},
	G3D::AbstractArray{Float64,3},
	matchTable::Array{CartesianIndex{2},3},
	refIndex::Array{CartesianIndex{2},2},
	patchSize::CartesianIndex{2},
	position::CartesianIndex{2},
)

	nMatch = size(matchTable, 3)

	@views G3D[1, CartesianIndex(1, 1):patchSize] .=
		img[refIndex[position]:(patchSize - CartesianIndex(
			1,
			1,
		) + refIndex[position])]

	@views @inbounds for k in 1:nMatch
		position2 = position + matchTable[position, k]
		G3D[k + 1, CartesianIndex(1, 1) + patchSize] .=
			img[refIndex[position2]:(patchSize - CartesianIndex(
				1,
				1,
			) + refIndex[position2])]
	end
end
