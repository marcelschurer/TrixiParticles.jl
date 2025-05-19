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

# trixi2vtk(ic_inlet, filename="ic_inlet")
# trixi2vtk(ic_outlet_left, filename="ic_outlet_left")
# trixi2vtk(ic_outlet_right, filename="ic_outlet_right")
# trixi2vtk(ic_packed, filename="ic_intake_manifold")
# trixi2vtk(ic_boundary, filename="ic_boundary")

# ==========================================================================================
# ==== Fluid
t_span_sim = [0.0, 0.5]

reynolds_number = 100
const prescribed_velocity = 2.0

smoothing_length = 1.5 * particle_spacing
smoothing_kernel = WendlandC2Kernel{3}()

fluid_density_calculator = ContinuityDensity()

kinematic_viscosity = prescribed_velocity * 50.0 / reynolds_number

viscosity = ViscosityAdami(nu=kinematic_viscosity)

n_buffer_particles = 40000
sound_speed = 10*prescribed_velocity

fluid_system = EntropicallyDampedSPHSystem(ic_intake_manifold, smoothing_kernel, smoothing_length,
                                           sound_speed, viscosity=viscosity,
                                           density_calculator=fluid_density_calculator,
                                           buffer_size=n_buffer_particles)

# ==========================================================================================
# ==== Open Boundary
open_boundary_model = BoundaryModelLastiwka()

A = [-0.010606, 0.11675, -0.002136]
B = [0.072195, 0.092376, -0.001068]
C = [0.023133, 0.10682, 0.042448]
flow_direction = normalize(cross(B .- A, C .- A))
plane_in = (A, B, C)

function velocity_function3d(pos, t)
    # Use this for a time-dependent inflow velocity
    # return SVector(0.5prescribed_velocity * sin(2pi * t) + prescribed_velocity, 0)

    return SVector(prescribed_velocity, 0.0, 0.0)
end

inflow = BoundaryZone(; plane=plane_in, plane_normal=flow_direction,
                      density=density, particle_spacing=particle_spacing,
                      open_boundary_layers=6, initial_condition=ic_inlet,
                      boundary_type=InFlow())

reference_velocity_in = velocity_function3d
reference_pressure_in = pressure
reference_density_in = density
open_boundary_in = OpenBoundarySPHSystem(inflow; fluid_system,
                                         boundary_model=open_boundary_model,
                                         buffer_size=n_buffer_particles,
                                         reference_density=reference_density_in,
                                         reference_pressure=reference_pressure_in,
                                         reference_velocity=reference_velocity_in)

# ==========================================================================================
# ==== Boundary
viscosity_boundary = ViscosityAdami(nu=1e-4)
boundary_model = BoundaryModelDummyParticles(ic_boundary.density, ic_boundary.mass,
                                             AdamiPressureExtrapolation(),
                                             state_equation=state_equation,
                                             viscosity=viscosity_boundary,
                                             smoothing_kernel, smoothing_length)

boundary_system = BoundarySPHSystem(ic_boundary, boundary_model)

# ==========================================================================================
# ==== Simulation
semi = Semidiscretization(fluid_system, boundary_system, open_boundary_in)

ode = semidiscretize(semi, t_span_sim)

info_callback = InfoCallback(interval=100)
saving_callback = SolutionSavingCallback(dt=0.02, prefix="")

extra_callback = nothing

callbacks = CallbackSet(info_callback, saving_callback, UpdateCallback(), extra_callback)

sol = solve(ode, RDPK3SpFSAL35(),
            abstol=1e-5, # Default abstol is 1e-6 (may need to be tuned to prevent boundary penetration)
            reltol=1e-3, # Default reltol is 1e-3 (may need to be tuned to prevent boundary penetration)
            dtmax=1e-2, # Limit stepsize to prevent crashing
            save_everystep=false, callback=callbacks);
