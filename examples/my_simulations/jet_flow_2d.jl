using TrixiParticles
using OrdinaryDiffEq

# ==========================================================================================
# ==== Resolution
particle_spacing = 0.02

# ==========================================================================================
# ==== Experiment Setup
gravity = 9.81
tspan = (0.0, 1)

boundary_layers = 3
wall_size = (60, 40) # wall bevor jet (length, height)
jet_number_of_stages = 3
jet_stage_size = 2

fluid_size = (30, wall_size[2] - 2 * boundary_layers)
velocity = (2.0, 0.0)
acceleration = (0.0, 0.0)

fluid_density = 1000.0
sound_speed = 10.0

# ==========================================================================================
# ==== Boundary geometry and initial fluid particle positions

state_equation = StateEquationCole(; sound_speed, reference_density=fluid_density,
                                   exponent=7, clip_negative_pressure=false)

wall_bottom = RectangularShape(particle_spacing, (wall_size[1], boundary_layers),
                               (0.0, 0.0), density=fluid_density)

wall_top = RectangularShape(particle_spacing, (wall_size[1], boundary_layers),
                            (0.0, (wall_size[2] - boundary_layers) * particle_spacing),
                            density=fluid_density)

jet_bottom = Vector{Any}(undef, jet_number_of_stages)
jet_top = Vector{Any}(undef, jet_number_of_stages)
for i in 1:jet_number_of_stages
    jet_bottom[i] = RectangularShape(particle_spacing, (jet_stage_size, boundary_layers),
                                     ((wall_size[1] - jet_stage_size + i * jet_stage_size) *
                                      particle_spacing,
                                      i * particle_spacing),
                                     density=fluid_density)
    jet_top[i] = RectangularShape(particle_spacing, (jet_stage_size, boundary_layers),
                                  ((wall_size[1] - jet_stage_size + i * jet_stage_size) *
                                   particle_spacing,
                                   (wall_size[2] - boundary_layers -
                                    i) * particle_spacing),
                                  density=fluid_density)
end

wall_after_jet_bottom = RectangularShape(particle_spacing, (wall_size[1], boundary_layers),
                                         ((wall_size[1] +
                                           jet_number_of_stages * jet_stage_size) *
                                          particle_spacing,
                                          jet_number_of_stages *
                                          particle_spacing),
                                         density=fluid_density)
wall_after_jet_top = RectangularShape(particle_spacing, (wall_size[1], boundary_layers),
                                      ((wall_size[1] +
                                        jet_number_of_stages * jet_stage_size) *
                                       particle_spacing,
                                       (wall_size[2] - boundary_layers -
                                        jet_number_of_stages) *
                                       particle_spacing),
                                      density=fluid_density)

channel = union(wall_bottom, wall_top, jet_bottom..., jet_top..., wall_after_jet_bottom,
                wall_after_jet_top)

fluid = RectangularShape(particle_spacing, fluid_size,
                         (0.0,
                          boundary_layers * particle_spacing),
                         density=fluid_density, velocity=velocity)

# ==========================================================================================
# ==== Fluid
smoothing_length = 1.2 * particle_spacing
smoothing_kernel = SchoenbergCubicSplineKernel{2}()

viscosity = ArtificialViscosityMonaghan(alpha=0.02, beta=0.0)

fluid_density_calculator = ContinuityDensity()

fluid_system = WeaklyCompressibleSPHSystem(fluid,
                                           fluid_density_calculator,
                                           state_equation,
                                           smoothing_kernel,
                                           smoothing_length,
                                           viscosity=viscosity,
                                           acceleration=acceleration, #-gravity,
                                           source_terms=nothing)

# ==========================================================================================
# ==== Boundary

# This is to set another boundary density calculation with `trixi_include`
boundary_density_calculator = AdamiPressureExtrapolation()

# This is to set wall viscosity with `trixi_include`
viscosity_wall = nothing
boundary_model = BoundaryModelDummyParticles(channel.density,
                                             channel.mass,
                                             boundary_density_calculator, smoothing_kernel,
                                             smoothing_length; viscosity=nothing,
                                             state_equation=state_equation,
                                             correction=nothing)

boundary_system = BoundarySPHSystem(channel, boundary_model,
                                    movement=nothing)

# ==========================================================================================
# ==== Simulation
semi = Semidiscretization(boundary_system, fluid_system)
ode = semidiscretize(semi, tspan)

info_callback = InfoCallback(interval=50)
saving_callback = SolutionSavingCallback(dt=0.02, prefix="",
                                         output_directory="out_jet_flow_2d")
fluid_flow_callback = FluidFlowCallback(dt=0.02)

# This is to easily add a new callback with `trixi_include`
extra_callback = fluid_flow_callback

callbacks = CallbackSet(info_callback, saving_callback, extra_callback)

# Use a Runge-Kutta method with automatic (error based) time step size control
sol = solve(ode, RDPK3SpFSAL35(), save_everystep=false, callback=callbacks);
