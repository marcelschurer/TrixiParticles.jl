struct FaceNeighborhoodSearch{NDIMS, ELTYPE, PB}
    cell_list         :: Dict{NTuple{NDIMS, Int}, Vector{Int}}
    cell_size         :: NTuple{NDIMS, ELTYPE}
    neighbor_iterator :: Dict{NTuple{NDIMS, Int}, Vector{Int}}
    empty_vector      :: Vector{Int} # Just an empty vector (used in `eachneighbor`)
    n_cells           :: NTuple{NDIMS, Int}
    periodic_box      :: Nothing
    search_radius     :: ELTYPE

    function FaceNeighborhoodSearch{NDIMS}(search_radius) where {NDIMS}
        ELTYPE = eltype(search_radius)

        cell_list = Dict{NTuple{NDIMS, Int}, Vector{Int}}()
        neighbor_iterator = Dict{NTuple{NDIMS, Int}, Vector{Int}}()

        cell_size = ntuple(dim -> search_radius, NDIMS)
        empty_vector = Int[]

        n_cells = ntuple(_ -> -1, Val(NDIMS))

        new{NDIMS, ELTYPE, Nothing}(cell_list, cell_size, neighbor_iterator, empty_vector,
                                    n_cells, nothing, search_radius)
    end
end

@inline Base.ndims(::FaceNeighborhoodSearch{NDIMS}) where {NDIMS} = NDIMS

@inline function eachneighbor(coords, neighborhood_search::FaceNeighborhoodSearch)
    (; neighbor_iterator, empty_vector) = neighborhood_search
    cell = PointNeighbors.cell_coords(coords, neighborhood_search)

    haskey(neighbor_iterator, cell) && return neighbor_iterator[cell]

    return empty_vector
end

function faces_in_cell(cell, neighborhood_search)
    (; cell_list, empty_vector) = neighborhood_search

    haskey(cell_list, cell) && return cell_list[cell]

    return empty_vector
end

function initialize!(neighborhood_search::FaceNeighborhoodSearch, geometry;
                     pad=ntuple(_ -> 1, ndims(geometry)))
    (; cell_list, neighbor_iterator, search_radius) = neighborhood_search

    empty!(cell_list)

    # Fill cells with intersecting faces
    for face in eachface(geometry)

        # Check if any face intersects a cell in the face-embedding cell grid
        for cell in cell_grid(face, geometry, neighborhood_search)
            if cell_intersection(face, geometry, cell, neighborhood_search)
                if haskey(cell_list, cell) && !(face in cell_list[cell])
                    # Add face to corresponding cell
                    append!(cell_list[cell], face)
                else
                    # Create cell
                    cell_list[cell] = [face]
                end
            end
        end
    end

    empty!(neighbor_iterator)

    min_cell = PointNeighbors.cell_coords(geometry.min_corner, neighborhood_search) .- pad
    max_cell = PointNeighbors.cell_coords(geometry.max_corner, neighborhood_search) .+ pad

    # Merge all lists of faces in the neighboring cells into one iterator
    face_ids = Int[]
    for cell_runner in meshgrid(min_cell, max_cell)
        resize!(face_ids, 0)
        for neighbor in PointNeighbors.neighboring_cells(cell_runner, neighborhood_search,
                                                         search_radius)
            append!(face_ids, faces_in_cell(Tuple(neighbor), neighborhood_search))
        end

        unique!(face_ids)

        if isempty(face_ids)
            continue
        end

        neighbor_iterator[cell_runner] = copy(face_ids)
    end

    return neighborhood_search
end

function cell_intersection(face, geometry, cell,
                           neighborhood_search::FaceNeighborhoodSearch{NDIMS}) where {NDIMS}
    (; cell_size) = neighborhood_search

    vertice_list = face_vertices(face, geometry)

    # Check if one of the vertices is inside cell
    for v in vertice_list
        cell == PointNeighbors.cell_coords(v, neighborhood_search) && return true
    end

    # Check if line segments intersect cell
    min_corner = SVector(cell .* cell_size...)
    max_corner = min_corner + SVector(cell_size...)

    ray_direction = vertice_list[2] - vertice_list[1]
    ray_origin = vertice_list[1]

    ray_intersection(min_corner, max_corner, ray_origin, ray_direction) && return true

    if NDIMS == 3
        ray_direction = vertice_list[2] - vertice_list[3]
        ray_origin = vertice_list[3]

        ray_intersection(min_corner, max_corner, ray_origin, ray_direction) && return true

        ray_direction = vertice_list[3] - vertice_list[1]
        ray_origin = vertice_list[1]

        ray_intersection(min_corner, max_corner, ray_origin, ray_direction) && return true

        # For 3D,  Check if triangle plane intersects cell (for very large triangles)
        normal = face_normal(face, geometry)

        return triangle_plane_intersection(ray_origin, normal, min_corner, max_corner,
                                           cell_size)
    end

    return false
end

# See https://tavianator.com/2022/ray_box_boundary.html
function ray_intersection(min_corner, max_corner, ray_origin, ray_direction;
                          pad=sqrt(eps()))
    NDIMS = length(ray_origin)

    inv_dir = SVector(ntuple(@inline(dim->1 / ray_direction[dim]), NDIMS))

    tmin = zero(eltype(ray_direction))
    tmax = Inf
    @inbounds for dim in 1:NDIMS
        # `pad` is to handle rays on the boundary
        t1 = (min_corner[dim] - pad - ray_origin[dim]) * inv_dir[dim]
        t2 = (max_corner[dim] + pad - ray_origin[dim]) * inv_dir[dim]

        tmin = min(max(t1, tmin), max(t2, tmin))
        tmax = max(min(t1, tmax), min(t2, tmax))
    end

    return tmin <= tmax
end

# Check if each cell vertex is located on the same side of the plane.
# Otherwise the plane intersects the cell.
function triangle_plane_intersection(point, plane_normal, min_corner, max_corner, cell_size)
    dirx = SVector(cell_size[1], zero(eltype(point)), zero(eltype(point)))
    diry = SVector(zero(eltype(point)), cell_size[2], zero(eltype(point)))
    dirz = SVector(zero(eltype(point)), zero(eltype(point)), cell_size[3])

    cell_vertex_1 = min_corner
    cell_vertex_2 = min_corner + dirx

    pos_diff_1 = cell_vertex_1 - point
    pos_diff_2 = cell_vertex_2 - point

    # Corners: bottom north-west and bottom north-east
    dot1 = sign(dot(pos_diff_1, plane_normal))
    dot2 = sign(dot(pos_diff_2, plane_normal))
    !(dot1 == dot2) && return true

    cell_vertex_3 = min_corner + diry
    pos_diff_3 = cell_vertex_3 - point

    # Corners: bottom north-east and top north-west
    dot3 = sign(dot(pos_diff_3, plane_normal))
    !(dot2 == dot3) && return true

    cell_vertex_4 = min_corner + dirz
    pos_diff_4 = cell_vertex_4 - point

    # Corners: top north-west and bottom south-west
    dot4 = sign(dot(pos_diff_4, plane_normal))
    !(dot3 == dot4) && return true

    cell_vertex_5 = max_corner
    pos_diff_5 = cell_vertex_5 - point

    # Corners: bottom south-west and top south-east
    dot5 = sign(dot(pos_diff_5, plane_normal))
    !(dot4 == dot5) && return true

    cell_vertex_6 = max_corner - dirx
    pos_diff_6 = cell_vertex_6 - point

    # Corners: top south-east and top south-west
    dot6 = sign(dot(pos_diff_6, plane_normal))
    !(dot5 == dot6) && return true

    cell_vertex_7 = max_corner - diry
    pos_diff_7 = cell_vertex_7 - point

    # Corners: top south-west and bottom south-east
    dot7 = sign(dot(pos_diff_7, plane_normal))
    !(dot6 == dot7) && return true

    cell_vertex_8 = max_corner - dirz
    pos_diff_8 = cell_vertex_8 - point

    # Corners: bottom south-east and top north-east
    dot8 = sign(dot(pos_diff_8, plane_normal))
    !(dot7 == dot8) && return true

    # All edge vertices are on one side of the plane
    return false
end

@inline function cell_grid(face, geometry,
                           neighborhood_search::FaceNeighborhoodSearch{NDIMS}) where {NDIMS}
    vertice_list = face_vertices(face, geometry)

    # Compute the cell coordinates for each vertex
    cells = [PointNeighbors.cell_coords(v, neighborhood_search) for v in vertice_list]

    # Compute the element-wise minimum and maximum cell coordinates across all vertices
    mins = reduce((a, b) -> min.(a, b), cells)
    maxs = reduce((a, b) -> max.(a, b), cells)

    return meshgrid(mins, maxs)
end

@inline function meshgrid(min_corner, max_corner; increment=1)
    min_ = collect(min_corner)
    max_ = collect(max_corner)

    ranges = ntuple(dim -> (min_[dim]:increment:max_[dim]), length(min_corner))

    return Iterators.product(ranges...)
end