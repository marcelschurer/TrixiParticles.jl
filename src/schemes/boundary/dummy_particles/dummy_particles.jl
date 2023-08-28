@doc raw"""
    BoundaryModelDummyParticles(initial_density, hydrodynamic_mass, state_equation,
                                density_calculator, smoothing_kernel, smoothing_length)

Boundaries modeled as dummy particles, which are treated like fluid particles,
but their positions and velocities are not evolved in time. Since the force towards the fluid
should not change with the material density when used with a `TotalLagrangianSPHSystem`, the
dummy particles need to have a mass corresponding to the fluid's rest density, which we call
"hydrodynamic mass", as opposed to mass corresponding to the material density of a
`TotalLagrangianSPHSystem`.

Here, `initial_density` and `hydrodynamic_mass` are vectors that contains the initial density
and the hydrodynamic mass respectively for each boundary particle.
Note that when used with [`SummationDensity`](@ref) (see below), this is only used to determine
the element type and the number of boundary particles.

To establish a relationship between density and pressure, a `state_equation` has to be passed,
which should be the same as for the adjacent fluid systems.
To sum over neighboring particles, a `smoothing_kernel` and `smoothing_length` needs to be passed.
This should be the same as for the adjacent fluid system with the largest smoothing length.

In the literature, this kind of boundary particles is referred to as
"dummy particles" (Adami et al., 2012 and Valizadeh & Monaghan, 2015),
"frozen fluid particles" (Akinci et al., 2012) or "dynamic boundaries (Crespo et al., 2007).
The key detail of this boundary condition and the only difference between the boundary models
in these references is the way the density and pressure of boundary particles is computed.

Since boundary particles are treated like fluid particles, the force
on fluid particle ``a`` due to boundary particle ``b`` is given by
```math
f_{ab} = m_a m_b \left( \frac{p_a}{\rho_a^2} + \frac{p_b}{\rho_b^2} \right) \nabla_{r_a} W(\Vert r_a - r_b \Vert, h).
```
The quantities to be defined here are the density ``\rho_b`` and pressure ``p_b``
of the boundary particle ``b``.

We provide three options to compute the boundary density and pressure, determined by the `density_calculator`:
1. With [`SummationDensity`](@ref), the density is calculated by summation over the neighboring particles,
   and the pressure is computed from the density with the state equation.
2. With [`ContinuityDensity`](@ref), the density is integrated from the continuity equation,
   and the pressure is computed from the density with the state equation.
   Note that this causes a gap between fluid and boundary where the boundary is initialized
   without any contact to the fluid. This is due to overestimation of the boundary density
   as soon as the fluid comes in contact with boundary particles that initially did not have
   contact to the fluid.
   Therefore, in dam break simulations, there is a visible "step", even though the boundary is supposed to be flat.
   See also [dual.sphysics.org/faq/#Q_13](https://dual.sphysics.org/faq/#Q_13).
3. With [`AdamiPressureExtrapolation`](@ref), the pressure is extrapolated from the pressure of the
   fluid according to (Adami et al., 2012), and the density is obtained by applying the inverse of the state equation.

## References:
- S. Adami, X. Y. Hu, N. A. Adams.
  "A generalized wall boundary condition for smoothed particle hydrodynamics".
  In: Journal of Computational Physics 231, 21 (2012), pages 7057–7075.
  [doi: 10.1016/J.JCP.2012.05.005](https://doi.org/10.1016/J.JCP.2012.05.005)
- Alireza Valizadeh, Joseph J. Monaghan.
  "A study of solid wall models for weakly compressible SPH".
  In: Journal of Computational Physics 300 (2015), pages 5–19.
  [doi: 10.1016/J.JCP.2015.07.033](https://doi.org/10.1016/J.JCP.2015.07.033)
- Nadir Akinci, Markus Ihmsen, Gizem Akinci, Barbara Solenthaler, Matthias Teschner.
  "Versatile rigid-fluid coupling for incompressible SPH".
  ACM Transactions on Graphics 31, 4 (2012), pages 1–8.
  [doi: 10.1145/2185520.2185558](https://doi.org/10.1145/2185520.2185558)
- A. J. C. Crespo, M. Gómez-Gesteira, R. A. Dalrymple.
  "Boundary conditions generated by dynamic particles in SPH methods"
  In: Computers, Materials and Continua 5 (2007), pages 173-184.
  [doi: 10.3970/cmc.2007.005.173](https://doi.org/10.3970/cmc.2007.005.173)
"""
struct BoundaryModelDummyParticles{ELTYPE <: Real, SE, DC, K, V, C}
    pressure           :: Vector{ELTYPE}
    hydrodynamic_mass  :: Vector{ELTYPE}
    state_equation     :: SE
    density_calculator :: DC
    smoothing_kernel   :: K
    smoothing_length   :: ELTYPE
    viscosity          :: V
    cache              :: C

    function BoundaryModelDummyParticles(initial_density, hydrodynamic_mass, state_equation,
                                         density_calculator, smoothing_kernel,
                                         smoothing_length; viscosity=NoViscosity())
        pressure = similar(initial_density)

        n_particles = length(initial_density)

        cache = (; create_cache(viscosity, n_particles, ndims(smoothing_kernel))...,
                 create_cache(initial_density, density_calculator)...)

        new{eltype(initial_density), typeof(state_equation),
            typeof(density_calculator), typeof(smoothing_kernel), typeof(viscosity),
            typeof(cache)}(pressure, hydrodynamic_mass, state_equation, density_calculator,
                           smoothing_kernel, smoothing_length, viscosity, cache)
    end
end

function Base.show(io::IO, model::BoundaryModelDummyParticles)
    @nospecialize model # reduce precompilation time

    print(io, "BoundaryModelDummyParticles(")
    print(io, model.density_calculator |> typeof |> nameof)
    print(io, ", ")
    print(io, model.viscosity |> typeof |> nameof)
    print(io, ")")
end

@inline function boundary_particle_impact(particle, boundary_particle,
                                          boundary_model::BoundaryModelDummyParticles,
                                          v_particle_system, v_boundary_system,
                                          particle_system, boundary_system,
                                          pos_diff, distance, m_b)
    rho_a = particle_density(v_particle_system, particle_system, particle)
    rho_b = particle_density(v_boundary_system, boundary_system, boundary_particle)

    grad_kernel = smoothing_kernel_grad(particle_system, pos_diff, distance)

    return -m_b *
           (particle_system.pressure[particle] / rho_a^2 +
            boundary_model.pressure[boundary_particle] / rho_b^2) *
           grad_kernel
end

@doc raw"""
    AdamiPressureExtrapolation()

The pressure of the boundary particles is obtained by extrapolating the pressure of the fluid
according to (Adami et al., 2012).
The pressure of a boundary particle ``b`` is given by
```math
p_b = \frac{\sum_f (p_f + \rho_f (\bm{g} - \bm{a}_b) \cdot \bm{r}_{bf}) W(\Vert r_{bf} \Vert, h)}{\sum_f W(\Vert r_{bf} \Vert, h)},
```
where the sum is over all fluid particles, ``\rho_f`` and ``p_f`` denote the density and pressure of fluid particle ``f``, respectively,
``r_{bf} = r_b - r_f`` denotes the difference of the coordinates of particles ``b`` and ``f``,
``\bm{g}`` denotes the gravitational acceleration acting on the fluid, and ``\bm{a}_b`` denotes the acceleration of the boundary particle ``b``.

## References:
- S. Adami, X. Y. Hu, N. A. Adams.
  "A generalized wall boundary condition for smoothed particle hydrodynamics".
  In: Journal of Computational Physics 231, 21 (2012), pages 7057–7075.
  [doi: 10.1016/J.JCP.2012.05.005](https://doi.org/10.1016/J.JCP.2012.05.005)
"""
struct AdamiPressureExtrapolation end

function create_cache(initial_density, ::SummationDensity)
    density = similar(initial_density)

    return (; density)
end

function create_cache(initial_density, ::ContinuityDensity)
    return (; initial_density)
end

function create_cache(initial_density, ::AdamiPressureExtrapolation)
    density = similar(initial_density)
    volume = similar(initial_density)

    return (; density, volume)
end

function create_cache(viscosity::NoViscosity, n_particles, n_dims)
    return (;)
end

function create_cache(viscosity::ViscosityAdami, n_particles, n_dims)
    ELTYPE = eltype(viscosity.nu)

    wall_velocity = zeros(ELTYPE, n_dims, n_particles)

    return (; wall_velocity)
end

function reset_cache!(cache, viscosity)
    (; density, volume) = cache

    set_zero!(density)
    set_zero!(volume)

    return cache
end

function reset_cache!(cache, viscosity::ViscosityAdami)
    (; density, volume, wall_velocity) = cache

    set_zero!(density)
    set_zero!(volume)
    set_zero!(wall_velocity)

    return cache
end

@inline function particle_density(v, model::BoundaryModelDummyParticles, system, particle)
    return particle_density(v, model.density_calculator, model, particle)
end

# Note that the other density calculators are dispatched in `density_calculators.jl`
@inline function particle_density(v, ::AdamiPressureExtrapolation, boundary_model, particle)
    (; cache) = boundary_model

    return cache.density[particle]
end

@inline function update_density!(boundary_model::BoundaryModelDummyParticles,
                                 system, system_index, v, u, v_ode, u_ode, semi)
    (; density_calculator) = boundary_model

    compute_density!(boundary_model, density_calculator, system, system_index, v, u, v_ode,
                     u_ode, semi)

    return boundary_model
end

function compute_density!(boundary_model,
                          ::Union{ContinuityDensity, AdamiPressureExtrapolation},
                          system, system_index, v, u, v_ode, u_ode, semi)
    # No density update for `ContinuityDensity`.
    # For `AdamiPressureExtrapolation`, the density is updated in `compute_pressure!`.
    return boundary_model
end

@inline function update_pressure!(boundary_model::BoundaryModelDummyParticles,
                                  system, system_index, v, u, v_ode, u_ode, semi)
    (; density_calculator) = boundary_model

    compute_pressure!(boundary_model, density_calculator, system, system_index, v, u, v_ode,
                      u_ode, semi)

    return boundary_model
end

function compute_density!(boundary_model, ::SummationDensity,
                          system, system_index, v, u, v_ode, u_ode, semi)
    (; cache) = boundary_model
    (; density) = cache # Density is in the cache for SummationDensity

    summation_density!(system, system_index, semi, u, u_ode, density,
                       particles=eachparticle(system))
end

function compute_pressure!(boundary_model, ::Union{SummationDensity, ContinuityDensity},
                           system, system_index, v, u, v_ode, u_ode, semi)
    (; state_equation, pressure) = boundary_model

    for particle in eachparticle(system)
        pressure[particle] = state_equation(particle_density(v, boundary_model, particle))
    end

    return boundary_model
end

function compute_pressure!(boundary_model, ::AdamiPressureExtrapolation,
                           system, system_index, v, u, v_ode, u_ode, semi)
    (; systems, neighborhood_searches) = semi
    (; pressure, state_equation, cache, viscosity) = boundary_model
    (; volume, density) = cache

    set_zero!(pressure)

    # Set `volume` and `density` to zero.
    # For `ViscosityAdami` the `wall_velocity` is also set to zero.
    reset_cache!(cache, viscosity)

    system_coords = current_coordinates(u, system)

    # Use all other systems for the pressure extrapolation
    @trixi_timeit timer() "compute boundary pressure" foreach_enumerate(systems) do (neighbor_system_index,
                                                                                     neighbor_system)
        v_neighbor_system = wrap_v(v_ode, neighbor_system_index,
                                   neighbor_system, semi)
        u_neighbor_system = wrap_u(u_ode, neighbor_system_index,
                                   neighbor_system, semi)

        neighborhood_search = neighborhood_searches[system_index][neighbor_system_index]

        neighbor_coords = current_coordinates(u_neighbor_system, neighbor_system)

        adami_pressure_extrapolation!(boundary_model, system, neighbor_system,
                                      system_coords, neighbor_coords,
                                      v_neighbor_system, neighborhood_search)
    end

    for particle in eachparticle(system)

        # The summation is only over fluid particles, thus the volume stays zero when a boundary
        # particle isn't surrounded by fluid particles.
        # Check the volume to avoid NaNs in pressure and velocity.
        if volume[particle] > eps()
            pressure[particle] /= volume[particle]

            # To impose no-slip condition
            compute_wall_velocity!(viscosity, system, system_coords, particle)
        end

        density[particle] = inverse_state_equation(state_equation, pressure[particle])
    end
end

@inline function adami_pressure_extrapolation!(boundary_model, system,
                                               neighbor_system::WeaklyCompressibleSPHSystem,
                                               system_coords, neighbor_coords,
                                               v_neighbor_system, neighborhood_search)
    (; pressure, cache, viscosity) = boundary_model

    # Loop over all pairs of particles and neighbors within the kernel cutoff.
    for_particle_neighbor(system, neighbor_system,
                          system_coords, neighbor_coords,
                          neighborhood_search;
                          particles=eachparticle(system)) do particle, neighbor,
                                                             pos_diff, distance
        density_neighbor = particle_density(v_neighbor_system, neighbor_system, neighbor)

        resulting_acc = neighbor_system.acceleration -
                        current_acceleration(system, particle)

        kernel_weight = smoothing_kernel(boundary_model, distance)

        pressure[particle] += (neighbor_system.pressure[neighbor] +
                               dot(resulting_acc, density_neighbor * pos_diff)) *
                              kernel_weight

        cache.volume[particle] += kernel_weight

        compute_smoothed_velocity!(cache, viscosity, neighbor_system, v_neighbor_system,
                                   kernel_weight, particle, neighbor)
    end

    for particle in eachparticle(system)
        # Limit pressure to be non-negative to avoid negative pressures at free surfaces
        pressure[particle] = max(pressure[particle], 0.0)
    end
end

@inline function adami_pressure_extrapolation!(boundary_model, system, neighbor_system,
                                               system_coords, neighbor_coords,
                                               v_neighbor_system, neighborhood_search)
    return boundary_model
end

function compute_smoothed_velocity!(cache, viscosity, neighbor_system, v_neighbor_system,
                                    kernel_weight, particle, neighbor)
    return cache
end

function compute_smoothed_velocity!(cache, viscosity::ViscosityAdami,
                                    neighbor_system, v_neighbor_system, kernel_weight,
                                    particle, neighbor)
    v_b = current_velocity(v_neighbor_system, neighbor_system, neighbor)

    for dim in 1:ndims(neighbor_system)
        cache.wall_velocity[dim, particle] += kernel_weight * v_b[dim]
    end

    return cache
end

@inline function compute_wall_velocity!(viscosity, system, system_coords, particle)
    return viscosity
end

@inline function compute_wall_velocity!(viscosity::ViscosityAdami, system,
                                        system_coords, particle)
    (; boundary_model) = system
    (; cache) = boundary_model
    (; volume, wall_velocity) = cache

    # Prescribed velocity of the boundary particle.
    # This velocity is zero when not using moving boundaries.
    v_boundary = current_velocity(system_coords, system, particle)

    for dim in 1:ndims(system)
        # The second term is the precalculated smoothed velocity field of the fluid.
        wall_velocity[dim, particle] = 2 * v_boundary[dim] -
                                       wall_velocity[dim, particle] / volume[particle]
    end
    return viscosity
end