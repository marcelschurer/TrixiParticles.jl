using TrixiParticles
using OrdinaryDiffEq

# ==========================================================================================
# ==== Resolution
particle_spacing = 0.05

# ==========================================================================================
# ==== Experiment Setup
tspan = (0.0, 2.0)
flow_direction = [0.0, 0.0, 1.0]
structure_size = (40, 40, 60) # size to taper

boundary_layers = 4
open_boundary_layers = 6
tapering_layers = 4 # number of tapering layers / stages
tapering_stage_size = 3 # length of one tapering stage

n_buffer_particles = 100000
const prescribed_velocity = 1.0
fluid_density = 1000.0
pressure = 1000.0
sound_speed = 10.0 * prescribed_velocity
reynolds_number = 100

structure = Vector{Any}(undef, tapering_layers)
bore = Vector{Any}(undef, tapering_layers)
fluid = Vector{Any}(undef, tapering_layers)

# ==========================================================================================
# ==== Boundary geometry and initial fluid particle positions

state_equation = StateEquationCole(; sound_speed, reference_density=fluid_density,
                                   exponent=7, background_pressure=pressure)

for i in 1:tapering_layers
    if i == 1
        structure[i] = RectangularShape(particle_spacing, structure_size,
                                        (0.0, 0.0, 0.0),
                                        density=fluid_density, pressure=pressure)
        bore[i] = RectangularShape(particle_spacing,
                                   (structure_size[1] - 2 * boundary_layers,
                                    structure_size[2] - 2 * boundary_layers,
                                    structure_size[3]),
                                   (boundary_layers * particle_spacing,
                                    boundary_layers * particle_spacing,
                                    0.0),
                                   density=fluid_density, pressure=pressure)
        fluid[i] = RectangularShape(particle_spacing,
                                    (structure_size[1] - 2 * boundary_layers,
                                     structure_size[2] - 2 * boundary_layers,
                                     structure_size[3] - open_boundary_layers),
                                    (boundary_layers * particle_spacing,
                                     boundary_layers * particle_spacing,
                                     open_boundary_layers * particle_spacing),
                                    density=fluid_density, pressure=pressure)
    elseif i < tapering_layers
        structure[i] = RectangularShape(particle_spacing,
                                        (structure_size[1] - 2 * (i - 1),
                                         structure_size[2] - 2 * (i - 1),
                                         tapering_stage_size),
                                        ((i - 1) * particle_spacing,
                                         (i - 1) * particle_spacing,
                                         (structure_size[3] + tapering_stage_size * (i - 2)) *
                                         particle_spacing),
                                        density=fluid_density, pressure=pressure)
        bore[i] = RectangularShape(particle_spacing,
                                   (structure_size[1] - 2 * (boundary_layers + (i - 1)),
                                    structure_size[2] - 2 * (boundary_layers + (i - 1)),
                                    tapering_stage_size),
                                   ((boundary_layers + (i - 1)) * particle_spacing,
                                    (boundary_layers + (i - 1)) * particle_spacing,
                                    (structure_size[3] + tapering_stage_size * (i - 2)) *
                                    particle_spacing),
                                   density=fluid_density, pressure=pressure)
        fluid[i] = bore[i]
    else
        structure[i] = RectangularShape(particle_spacing,
                                        (structure_size[1] - 2 * (i - 1),
                                         structure_size[2] - 2 * (i - 1),
                                         structure_size[3]),
                                        ((i - 1) * particle_spacing,
                                         (i - 1) * particle_spacing,
                                         (structure_size[3] + tapering_stage_size * (i - 2)) *
                                         particle_spacing),
                                        density=fluid_density, pressure=pressure)
        bore[i] = RectangularShape(particle_spacing,
                                   (structure_size[1] - 2 * (boundary_layers + (i - 1)),
                                    structure_size[2] - 2 * (boundary_layers + (i - 1)),
                                    structure_size[3]),
                                   ((boundary_layers + (i - 1)) * particle_spacing,
                                    (boundary_layers + (i - 1)) * particle_spacing,
                                    (structure_size[3] + tapering_stage_size * (i - 2)) *
                                    particle_spacing),
                                   density=fluid_density, pressure=pressure)
        fluid[i] = RectangularShape(particle_spacing,
                                    (structure_size[1] - 2 * (boundary_layers + (i - 1)),
                                     structure_size[2] - 2 * (boundary_layers + (i - 1)),
                                     structure_size[3] - open_boundary_layers),
                                    ((boundary_layers + (i - 1)) * particle_spacing,
                                     (boundary_layers + (i - 1)) * particle_spacing,
                                     (structure_size[3] + tapering_stage_size * (i - 2)) *
                                     particle_spacing),
                                    density=fluid_density, pressure=pressure)
    end
end

initial_structure = setdiff(union(structure...), union(bore...))
initial_fluid = union(fluid...)

# ==========================================================================================
# ==== Fluid
smoothing_length = 3.0 * particle_spacing
smoothing_kernel = WendlandC2Kernel{3}()

fluid_density_calculator = ContinuityDensity()
kinematic_viscosity = prescribed_velocity * 0.4 / reynolds_number
viscosity = ViscosityAdami(nu=kinematic_viscosity)

fluid_system = EntropicallyDampedSPHSystem(initial_fluid, smoothing_kernel,
                                           smoothing_length,
                                           sound_speed, viscosity=viscosity,
                                           density_calculator=fluid_density_calculator,
                                           buffer_size=n_buffer_particles)

# ==========================================================================================
# ==== Open Boundary
function velocity_function(pos, t)
    # Use this for a time-dependent inflow velocity
    # return SVector(0.5prescribed_velocity * sin(2pi * t) + prescribed_velocity, 0)

    return SVector(0.0, 0.0, prescribed_velocity)
end

inflow = InFlow(;
                plane=([
                           (boundary_layers) * particle_spacing,
                           (boundary_layers) * particle_spacing,
                           (open_boundary_layers) * particle_spacing
                       ],
                       [
                           (structure_size[1] - boundary_layers) * particle_spacing,
                           (boundary_layers) * particle_spacing,
                           (open_boundary_layers) * particle_spacing
                       ],
                       [
                           (boundary_layers) * particle_spacing,
                           (structure_size[2] - boundary_layers) * particle_spacing,
                           (open_boundary_layers) * particle_spacing
                       ]),
                flow_direction,
                open_boundary_layers, density=fluid_density, particle_spacing)

open_boundary_in = OpenBoundarySPHSystem(inflow; sound_speed, fluid_system,
                                         buffer_size=n_buffer_particles,
                                         reference_pressure=pressure,
                                         reference_velocity=velocity_function)

outflow = OutFlow(;
                  plane=([
                             (boundary_layers) * particle_spacing,
                             (boundary_layers) * particle_spacing,
                             (2 * structure_size[3] - open_boundary_layers +
                              tapering_stage_size * (tapering_layers - 2)) *
                             particle_spacing
                         ],
                         [
                             (structure_size[1] - boundary_layers) * particle_spacing,
                             (boundary_layers) * particle_spacing,
                             (2 * structure_size[3] - open_boundary_layers +
                              tapering_stage_size * (tapering_layers - 2)) *
                             particle_spacing
                         ],
                         [
                             (boundary_layers) * particle_spacing,
                             (structure_size[2] - boundary_layers) * particle_spacing,
                             (2 * structure_size[3] - open_boundary_layers +
                              tapering_stage_size * (tapering_layers - 2)) *
                             particle_spacing
                         ]),
                  flow_direction,
                  open_boundary_layers, density=fluid_density, particle_spacing)

open_boundary_out = OpenBoundarySPHSystem(outflow; sound_speed, fluid_system,
                                          buffer_size=n_buffer_particles,
                                          reference_pressure=pressure,
                                          reference_velocity=velocity_function)

# ==========================================================================================
# ==== Boundary

boundary_model = BoundaryModelDummyParticles(initial_structure.density,
                                             initial_structure.mass,
                                             AdamiPressureExtrapolation(),
                                             state_equation=state_equation,
                                             #viscosity=ViscosityAdami(nu=1e-4),
                                             smoothing_kernel, smoothing_length)

boundary_system = BoundarySPHSystem(initial_structure, boundary_model)

# ==========================================================================================
# ==== Simulation
semi = Semidiscretization(boundary_system, open_boundary_in, open_boundary_out,
                          fluid_system)
# semi = Semidiscretization(boundary_system, fluid_system)
ode = semidiscretize(semi, tspan)

info_callback = InfoCallback(interval=100)
saving_callback = SolutionSavingCallback(dt=0.02, prefix="",
                                         output_directory="out_channel_flow_3d")

callbacks = CallbackSet(info_callback, saving_callback, UpdateCallback())

sol = solve(ode, RDPK3SpFSAL35(),
            abstol=1e-5, # Default abstol is 1e-6 (may need to be tuned to prevent boundary penetration)
            reltol=1e-3, # Default reltol is 1e-3 (may need to be tuned to prevent boundary penetration)
            dtmax=1e-2, # Limit stepsize to prevent crashing
            save_everystep=false, callback=callbacks);
