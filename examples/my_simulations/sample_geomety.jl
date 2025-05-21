using TrixiParticles
using OrdinaryDiffEq
using LinearAlgebra

# Create Dictionaries
files = Dict()
geometry = Dict()

filenames = [
    "intake_manifold",
    "intake_manifold_sum",
    "inlet",
    "outlet_left",
    "outlet_right"
]
for filename in filenames
    files[filename] = pkgdir(TrixiParticles, "examples", "preprocessing", "data",
                             filename * ".stl")
end

# ==========================================================================================
# ==== Packing parameters
tlsph = true

# ==========================================================================================
# ==== Resolution
particle_spacing = 0.002

# The following depends on the sampling of the particles. In this case `boundary_thickness`
# means literally the thickness of the boundary packed with boundary particles and *not*
# how many rows of boundary particles will be sampled.
boundary_thickness = 8 * particle_spacing

# ==========================================================================================
# ==== Load complex geometry
density = 1000.0
pressure = 1.0e5
state_equation = nothing

for filename in filenames
    geometry[filename] = load_geometry(files[filename])
end

point_in_geometry_algorithm = WindingNumberJacobson(;
                                                    geometry=geometry["intake_manifold_sum"],
                                                    winding_number_factor=0.4,
                                                    hierarchical_winding=true)
# Returns `InitialCondition`
shape_sampled = ComplexShape(geometry["intake_manifold_sum"]; particle_spacing, density,
                             point_in_geometry_algorithm=point_in_geometry_algorithm)

signed_distance_field = SignedDistanceField(geometry["intake_manifold_sum"],
                                            particle_spacing;
                                            use_for_boundary_packing=true,
                                            max_signed_distance=boundary_thickness)
# Returns `InitialCondition`
boundary_sampled = sample_boundary(signed_distance_field; boundary_density=density,
                                   boundary_thickness, tlsph=tlsph)

# ==========================================================================================
# ==== Packing

# Large `background_pressure` can cause high accelerations. That is, the adaptive
# time-stepsize will be adjusted properly. We found that the following order of
# `background_pressure` result in appropriate stepsizes.
background_pressure = 1e6 * particle_spacing^ndims(geometry["intake_manifold_sum"])

packing_system = ParticlePackingSystem(shape_sampled;
                                       signed_distance_field, tlsph=tlsph,
                                       background_pressure)

boundary_system = ParticlePackingSystem(boundary_sampled;
                                        is_boundary=true, signed_distance_field,
                                        tlsph=tlsph, boundary_compress_factor=0.8,
                                        background_pressure)

# ==========================================================================================
# ==== Simulation
semi = Semidiscretization(packing_system, boundary_system)

# Use a high `tspan` to guarantee that the simulation runs at least for `maxiters`
tspan = (0, 10.0)
ode = semidiscretize(semi, tspan)

# Use this callback to stop the simulation when it is sufficiently close to a steady state
steady_state = SteadyStateReachedCallback(; interval=1, interval_size=10,
                                          abstol=1.0e-5, reltol=1.0e-3)

info_callback = InfoCallback(interval=50)

save_intervals = false
saving_callback = save_intervals ?
                  SolutionSavingCallback(interval=10, prefix="", ekin=kinetic_energy) :
                  nothing

callbacks = CallbackSet(UpdateCallback(), saving_callback, info_callback, steady_state)

sol = solve(ode, RDPK3SpFSAL35();
            save_everystep=false, maxiters=1000, callback=callbacks, dtmax=1e-2)

ic_packed = InitialCondition(sol, packing_system, semi)
ic_boundary = InitialCondition(sol, boundary_system, semi)

ic_inlet = intersect(ic_packed, geometry["inlet"])
ic_outlet_left = intersect(ic_packed, geometry["outlet_left"])
ic_outlet_right = intersect(ic_packed, geometry["outlet_right"])
ic_intake_manifold = intersect(ic_packed, geometry["intake_manifold"])

trixi2vtk(ic_inlet, filename="ic_inlet", output_directory="sampled_geometries")
trixi2vtk(ic_outlet_left, filename="ic_outlet_left", output_directory="sampled_geometries")
trixi2vtk(ic_outlet_right, filename="ic_outlet_right",
          output_directory="sampled_geometries")
trixi2vtk(ic_intake_manifold, filename="ic_intake_manifold",
          output_directory="sampled_geometries")
trixi2vtk(ic_boundary, filename="ic_boundary", output_directory="sampled_geometries")
