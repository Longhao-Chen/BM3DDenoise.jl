"""
	form_group!(G3D::AbstractArray{Float64, 3},
			img::AbstractArray{Float64, 2},
			matchTable::Array{CartesianIndex{2}, 3},
			refIndex::Array{CartesianIndex{2},2},
			patchSize::CartesianIndex{2},
			position::CartesianIndex{2},
			group_transform!::Function)

Forward BM3DDenoise groupings
"""
function form_group!(
	G3D::AbstractArray{Float64,3},
	img::AbstractArray{Float64,2},
	matchTable::Array{CartesianIndex{2},3},
	refIndex::Array{CartesianIndex{2},2},
	patchSize::CartesianIndex{2},
	position::CartesianIndex{2},
	group_transform!::Function,
)

	# Form table of 3D groups
	image_to_group!(img, G3D, matchTable, refIndex, patchSize, position)

	# Apply group transform
	group_transform!(G3D)
end

"""
	invert_group!(img::AbstractArray{Float64, 2},
				G3D::AbstractArray{Float64, 3},
				matchTable::Array{CartesianIndex{2},3},
				refIndex::Array{CartesianIndex{2},2},
				patchSize::CartesianIndex{2},
				position::CartesianIndex{2},
				group_itransform!::Function,
				imgLockPool::Array{ReentrantLock})

Inverse BM3DDenoise groupings
"""
function invert_group!(
	img::AbstractArray{Float64,2},
	G3D::AbstractArray{Float64,3},
	matchTable::Array{CartesianIndex{2},3},
	refIndex::Array{CartesianIndex{2},2},
	patchSize::CartesianIndex{2},
	position::CartesianIndex{2},
	group_itransform!::Function,
	imgLockPool::Array{ReentrantLock},
)

	# Apply inverse group transform
	group_itransform!(G3D)

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
