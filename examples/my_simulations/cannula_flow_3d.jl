using TrixiParticles
using OrdinaryDiffEq

# ==========================================================================================
# ==== Resolution
particle_spacing = 0.02

# ==========================================================================================
# ==== Experiment Setup
gravity = 9.81
tspan = (0.0, 1.0)

boundary_layers = 3
wall_size = (40.0, 40.0, 100.0) # wall bevor jet (length, height)

fluid_size = (30.0, 30.0, 0.0)

velocity = (0.0, 0.0, 0.0)
acceleration = (0.0, 0.0, 0.0)

fluid_density = 1000.0
sound_speed = 10.0

# ==========================================================================================
# ==== Boundary geometry and initial fluid particle positions

state_equation = StateEquationCole(; sound_speed, reference_density=fluid_density,
                                   exponent=7, clip_negative_pressure=false)

cannula = RectangularShape(particle_spacing, wall_size,
                           (0.0, 0.0, 0.0), density=fluid_density)

fluid = RectangularShape(particle_spacing, fluid_size,
                         (0.0, 0.0, 0.0),
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
                                           acceleration=acceleration,
                                           source_terms=nothing)

# ==========================================================================================
# ==== Boundary

# This is to set another boundary density calculation with `trixi_include`
boundary_density_calculator = AdamiPressureExtrapolation()

# This is to set wall viscosity with `trixi_include`
viscosity_wall = nothing
boundary_model = BoundaryModelDummyParticles(cannula.density,
                                             cannula.mass,
                                             boundary_density_calculator, smoothing_kernel,
                                             smoothing_length; viscosity=nothing,
                                             state_equation=state_equation,
                                             correction=nothing)

boundary_system = BoundarySPHSystem(cannula, boundary_model,
                                    movement=nothing)

# ==========================================================================================
# ==== Simulation
semi = Semidiscretization(boundary_system, fluid_system)
ode = semidiscretize(semi, tspan)

info_callback = InfoCallback(interval=50)
saving_callback = SolutionSavingCallback(dt=0.02, prefix="",
                                         output_directory="out_my_simulation")
fluid_flow_callback = FluidFlowCallback(dt=0.02)

# This is to easily add a new callback with `trixi_include`
extra_callback = fluid_flow_callback

callbacks = CallbackSet(info_callback, saving_callback, extra_callback)

# Use a Runge-Kutta method with automatic (error based) time step size control
sol = solve(ode, RDPK3SpFSAL35(), save_everystep=false, callback=callbacks);
