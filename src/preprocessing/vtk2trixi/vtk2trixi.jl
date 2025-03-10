"""
    vtk2trixi(file::String)

Convert data from VTK-file to InitialCondition

# Arguments
- `file`: Name of the file to be loaded.

# Example
```jldoctest; output = false
# Create a rectangular shape
rectangular = RectangularShape(0.1, (10, 10), (0, 0), density=1.5, velocity=(1.0, -2.0),
                               pressure=1000.0)

# Write the `InitialCondition` to a vtk file
trixi2vtk(rectangular; filename="rectangular", output_directory="out")

# Read the vtk file and convert it to 'InitialCondition'
ic = vtk2trixi(joinpath("out", "rectangular.vtu")

# output
┌──────────────────────────────────────────────────────────────────────────────────────────────────┐
│ InitialCondition{Float64}                                                                        │
│ ═════════════════════════                                                                        │
│ #dimensions: ……………………………………………… 2                                                                │
│ #particles: ………………………………………………… 100                                                              │
│ particle spacing: ………………………………… 0.1                                                              │
└──────────────────────────────────────────────────────────────────────────────────────────────────┘
"""
function vtk2trixi(file)
    vtk_file = ReadVTK.VTKFile(file)

    # Retrieve data fields (e.g., pressure, velocity, ...)
    point_data = ReadVTK.get_point_data(vtk_file)
    field_data = ReadVTK.get_field_data(vtk_file)

    # Retrieve fields
    ndims = first(ReadVTK.get_data(field_data["ndims"]))
    particle_spacing = first(ReadVTK.get_data(field_data["particle_spacing"]))

    coordinates = ReadVTK.get_points(vtk_file)[1:ndims, :]

    fields = ["velocity", "density", "pressure", "mass"]
    results = Dict{String, Array{Float64}}()

    for field in fields
        found = false
        for k in ReadVTK.keys(point_data)
            if !isnothing(match(Regex("$field"), k))
                results[field] = ReadVTK.get_data(point_data[k])
                found = true
                break
            end
        end
        if !found
            # Set fields to zero if not found
            if field in ["density", "pressure", "mass"]
                results[field] = zeros(size(coordinates, 2))
            else
                results[field] = zero(coordinates)
            end
            @info "No '$field' field found in VTK file. Will be set to zero."
        end
    end

    return InitialCondition(
                            ; coordinates, particle_spacing,
                            velocity=results["velocity"],
                            mass=results["mass"],
                            density=results["density"],
                            pressure=results["pressure"])
end