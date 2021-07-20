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
	@inbounds @views for k in 1:size(G3D, 1)
		transform_2D!(G3D[k, :, :])
	end

	# Apply 1D transform
	@inbounds @views for j in 1:patchSize[2]
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
				itransform_2D!::Function,
				imgLockPool::Array{ReentrantLock})

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
	imgLockPool::Array{ReentrantLock},
)

	# Apply inverse 1D transform
	@inbounds @views for j in 1:patchSize[2]
		for i in 1:patchSize[1]
			itransform_1D!(G3D[:, i, j])
		end
	end

	# Apply inverse 2D transform
	@inbounds @views for k in 1:size(G3D, 1)
		itransform_2D!(G3D[k, :, :])
	end

	group_to_image!(img, G3D, matchTable, refIndex, patchSize, position, imgLockPool)

end

"""
	group_to_image!(img::AbstractArray{Float64, 2},
			G3D::AbstractArray{Float64, 3},
			matchTable::Array{CartesianIndex{2},3},
			refIndex::Array{CartesianIndex{2},2},
			patchSize::CartesianIndex{2},
			position::CartesianIndex{2},
			imgLockPool::Array{ReentrantLock})

group to image
"""
function group_to_image!(
	img::AbstractArray{Float64,2},
	G3D::AbstractArray{Float64,3},
	matchTable::Array{CartesianIndex{2},3},
	refIndex::Array{CartesianIndex{2},2},
	patchSize::CartesianIndex{2},
	position::CartesianIndex{2},
	imgLockPool::Array{ReentrantLock},
)

	nMatch = size(matchTable, 3)

	lock(imgLockPool[position]) do
		@views img[refIndex[position]:(patchSize - CartesianIndex(
			1,
			1,
		) + refIndex[position])] .+= G3D[1, :, :]
	end

	@inbounds @views for k in 1:nMatch
		position2 = matchTable[position, k]
		lock(imgLockPool[position2]) do
			img[refIndex[position2]:(patchSize - CartesianIndex(
				1,
				1,
			) + refIndex[position2])] .+= G3D[k + 1, :, :]
		end
	end
end

"""
	group_to_image!(img::AbstractArray{Float64, 2},
			W::Float64,
			matchTable::Array{CartesianIndex{2},3},
			refIndex::Array{CartesianIndex{2},2},
			patchSize::CartesianIndex{2},
			position::CartesianIndex{2},
			imgLockPool::Array{ReentrantLock})
"""
function group_to_image!(
	img::AbstractArray{Float64,2},
	W::Float64,
	matchTable::Array{CartesianIndex{2},3},
	refIndex::Array{CartesianIndex{2},2},
	patchSize::CartesianIndex{2},
	position::CartesianIndex{2},
	imgLockPool::Array{ReentrantLock},
)

	nMatch = size(matchTable, 3)

	lock(imgLockPool[position]) do
		@views img[refIndex[position]:(patchSize - CartesianIndex(
			1,
			1,
		) + refIndex[position])] .+= W
	end

	@inbounds @views for k in 1:nMatch
		position2 = matchTable[position, k]
		lock(imgLockPool[position2]) do
			img[refIndex[position2]:(patchSize - CartesianIndex(
				1,
				1,
			) + refIndex[position2])] .+= W
		end
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

	@views G3D[1, :, :] .=
		img[refIndex[position]:(patchSize - CartesianIndex(
			1,
			1,
		) + refIndex[position])]

	@views @inbounds for k in 1:nMatch
		position2 = matchTable[position, k]
		G3D[k + 1, :, :] .=
			img[refIndex[position2]:(patchSize - CartesianIndex(
				1,
				1,
			) + refIndex[position2])]
	end
end
