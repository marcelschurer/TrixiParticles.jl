# 2D channel flow simulation with open boundaries.

using TrixiParticles
using OrdinaryDiffEq

# ==========================================================================================
# ==== Resolution
particle_spacing = 0.1 # following https://doi.org/10.1016/j.cma.2020.113119

# Make sure that the kernel support of fluid particles at a boundary is always fully sampled
boundary_layers = 4

# Make sure that the kernel support of fluid particles at an open boundary is always
# fully sampled.
# Note: Due to the dynamics at the inlets and outlets of open boundaries,
# it is recommended to use `open_boundary_layers > boundary_layers`
open_boundary_layers = 6 # following https://doi.org/10.1016/j.cma.2020.113119

# ==========================================================================================
# ==== Experiment Setup
tspan = (0.0, 10.0)
flow_direction = [1.0, 0.0]
const prescribed_velocity = 1.0 # following https://doi.org/10.1016/j.cma.2020.113119 
sound_speed = 10 * prescribed_velocity # following https://doi.org/10.1016/j.cma.2020.113119 

# length of domain and initial fluid particle positions
L = 10 # following http://dx.doi.org/10.1017/S0022112083002839 (200mm/500mm)
h = 4.9 # Stepsize following http://dx.doi.org/10.1017/S0022112083002839
reynolds_number = 389 # following https://doi.org/10.1016/j.cma.2020.113119
fluid_density = 1225.0 # following https://doi.org/10.1016/j.cma.2020.113119

# For this particular example, it is necessary to have a background pressure.
# Otherwise the suction at the outflow is to big and the simulation becomes unstable.
pressure = 1.0 # 1.0 following https://doi.org/10.1016/j.cma.2020.113119

state_equation = StateEquationCole(; sound_speed, reference_density=fluid_density,
                                   exponent=7, background_pressure=pressure)

fluid_size_inlet = (Int(floor(0.25 * L / particle_spacing)),
                    Int(floor(5.2 / particle_spacing)))
fluid_size_outlet = (Int(floor(0.75 * L / particle_spacing)),
                     Int(floor((h + 5.2) / particle_spacing)))

fluid_inlet = RectangularShape(particle_spacing,
                               fluid_size_inlet, (0.0, 0.0),
                               pressure=pressure,
                               density=fluid_density)

fluid_outlet = RectangularShape(particle_spacing,
                                fluid_size_outlet, (0.25 * L, -h), pressure=pressure,
                                density=fluid_density)

ic_fluid = union(fluid_inlet, fluid_outlet)

boundary_top = RectangularShape(particle_spacing,
                                ((fluid_size_inlet[1] + fluid_size_outlet[1] +
                                  2 * open_boundary_layers),
                                 boundary_layers),
                                (-open_boundary_layers * particle_spacing,
                                 fluid_size_inlet[2] * particle_spacing),
                                pressure=pressure,
                                density=fluid_density)

boundary_bottom_left = RectangularShape(particle_spacing,
                                        (fluid_size_inlet[1] + open_boundary_layers,
                                         boundary_layers),
                                        (-open_boundary_layers * particle_spacing,
                                         -boundary_layers * particle_spacing),
                                        pressure=pressure,
                                        density=fluid_density)

boundary_bottom_right = RectangularShape(particle_spacing,
                                         (fluid_size_outlet[1] + open_boundary_layers,
                                          boundary_layers),
                                         (fluid_size_inlet[1] * particle_spacing,
                                          -(fluid_size_outlet[2] - fluid_size_inlet[2] +
                                            boundary_layers) * particle_spacing),
                                         pressure=pressure,
                                         density=fluid_density)

boundary_left = RectangularShape(particle_spacing,
                                 (boundary_layers,
                                  fluid_size_outlet[2] - fluid_size_inlet[2]),
                                 ((fluid_size_inlet[1] - boundary_layers) *
                                  particle_spacing,
                                  -(fluid_size_outlet[2] - fluid_size_inlet[2] +
                                    boundary_layers) *
                                  particle_spacing),
                                 pressure=pressure,
                                 density=fluid_density)

ic_boundary = union(boundary_top, boundary_bottom_left, boundary_bottom_right,
                    boundary_left)

# ==========================================================================================
# ==== Fluid
smoothing_length = 1.4 * particle_spacing
smoothing_kernel = SchoenbergQuinticSplineKernel{2}() # SchoenbergQuinticSplineKernel/WendlandC2Kernel following https://doi.org/10.1016/j.cma.2020.113119

fluid_density_calculator = ContinuityDensity()

kinematic_viscosity = 4 * h / (3 * reynolds_number * prescribed_velocity) # following https://doi.org/10.1016/j.cma.2020.113119

n_buffer_particles = open_boundary_layers * fluid_size_outlet[2] * fluid_size_inlet[2]

viscosity = ViscosityAdami(nu=kinematic_viscosity)

fluid_system = EntropicallyDampedSPHSystem(ic_fluid, smoothing_kernel,
                                           smoothing_length,
                                           sound_speed, viscosity=viscosity,
                                           density_calculator=fluid_density_calculator,
                                           buffer_size=n_buffer_particles)

#Alternatively the WCSPH scheme can be used
# alpha = 8 * kinematic_viscosity / (smoothing_length * sound_speed)
# viscosity = ArtificialViscosityMonaghan(; alpha, beta=0.0)
# density_diffusion = DensityDiffusionAntuono(ic_fluid; delta=0.1)

# fluid_system = WeaklyCompressibleSPHSystem(ic_fluid, fluid_density_calculator,
#                                            state_equation, smoothing_kernel,
#                                            smoothing_length, viscosity=viscosity,
#                                            buffer_size=n_buffer_particles,
#                                            density_diffusion=density_diffusion)

# ==========================================================================================
# ==== Open Boundary
function velocity_function_inlet(pos, t)
    # Use this for a time-dependent inflow velocity
    # return SVector(0.5prescribed_velocity * sin(2pi * t) + prescribed_velocity, 0)

    return SVector(prescribed_velocity, 0.0)
end

function velocity_function_outlet(pos, t)
    # Use this for a time-dependent inflow velocity
    #return SVector(0.5prescribed_velocity * sin(2pi * t) + prescribed_velocity, 0)
    return SVector(0.15(1 + tanh(5(t - 0.9))), 0)
    # return SVector((5.2 / (4.9 + 5.2)) * prescribed_velocity, 0.0)
    #return SVector((0.0, 0.0))
end

inflow = InFlow(; plane=([0.0, 0.0], [0.0, fluid_size_inlet[2] * particle_spacing]),
                flow_direction,
                open_boundary_layers, density=fluid_density, particle_spacing)

open_boundary_in = OpenBoundarySPHSystem(inflow; sound_speed, fluid_system,
                                         buffer_size=n_buffer_particles,
                                         reference_pressure=pressure,
                                         reference_velocity=velocity_function_inlet)

outflow = OutFlow(;
                  plane=([
                             (fluid_size_inlet[1] + fluid_size_outlet[1]) *
                             particle_spacing,
                             (fluid_size_inlet[2] - fluid_size_outlet[2]) *
                             particle_spacing
                         ],
                         [
                             (fluid_size_inlet[1] + fluid_size_outlet[1]) *
                             particle_spacing,
                             fluid_size_inlet[2] * particle_spacing
                         ]),
                  flow_direction, open_boundary_layers, density=fluid_density,
                  particle_spacing)

open_boundary_out = OpenBoundarySPHSystem(outflow; sound_speed, fluid_system,
                                          buffer_size=n_buffer_particles,
                                          reference_pressure=pressure,
                                          reference_velocity=velocity_function_outlet)

# ==========================================================================================
# ==== Boundary

boundary_model = BoundaryModelDummyParticles(ic_boundary.density,
                                             ic_boundary.mass,
                                             AdamiPressureExtrapolation(),
                                             state_equation=state_equation,
                                             viscosity=viscosity,
                                             smoothing_kernel, smoothing_length)

boundary_system = BoundarySPHSystem(ic_boundary, boundary_model)

# ==========================================================================================
# ==== Simulation
semi = Semidiscretization(fluid_system, open_boundary_in, open_boundary_out,
                          boundary_system)

ode = semidiscretize(semi, tspan)

info_callback = InfoCallback(interval=50)
saving_callback = SolutionSavingCallback(dt=0.02, prefix="",
                                         output_directory="out_my_simulation")

callbacks = CallbackSet(info_callback, saving_callback, UpdateCallback())

sol = solve(ode, RDPK3SpFSAL35(),
            abstol=1e-5, # Default abstol is 1e-6 (may need to be tuned to prevent boundary penetration)
            reltol=1e-3, # Default reltol is 1e-3 (may need to be tuned to prevent boundary penetration)
            dtmax=1e-2, # Limit stepsize to prevent crashing
            save_everystep=false, callback=callbacks);
