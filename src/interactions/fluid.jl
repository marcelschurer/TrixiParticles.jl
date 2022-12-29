# Fluid-fluid interaction
function interact!(du, u_particle_container, u_neighbor_container, neighborhood_search,
                   particle_container::FluidParticleContainer,
                   neighbor_container::FluidParticleContainer)
    @unpack density_calculator, smoothing_kernel, smoothing_length = particle_container

    @threaded for particle in each_moving_particle(particle_container)
        particle_coords = get_current_coords(particle, u_particle_container, particle_container)
        for neighbor in eachneighbor(particle_coords, neighborhood_search)
            neighbor_coords = get_current_coords(neighbor, u_neighbor_container, neighbor_container)

            pos_diff = particle_coords - neighbor_coords
            distance = norm(pos_diff)

            if sqrt(eps()) < distance <= compact_support(smoothing_kernel, smoothing_length)
                calc_dv!(du, u_particle_container, u_neighbor_container,
                         particle, neighbor, pos_diff, distance,
                         particle_container, neighbor_container)

                continuity_equation!(du, density_calculator,
                                     u_particle_container, u_neighbor_container,
                                     particle, neighbor, pos_diff, distance,
                                     particle_container, neighbor_container)
            end
        end
    end

    return du
end


@inline function calc_dv!(du, u_particle_container, u_neighbor_container,
                          particle, neighbor, pos_diff, distance,
                          particle_container, neighbor_container)
    @unpack smoothing_kernel, smoothing_length, state_equation, viscosity = particle_container

    density_particle = get_particle_density(particle, u_particle_container, particle_container)
    density_neighbor = get_particle_density(neighbor, u_neighbor_container, neighbor_container)

    # Viscosity
    v_diff = get_particle_vel(particle, u_particle_container, particle_container) -
        get_particle_vel(neighbor, u_neighbor_container, neighbor_container)
    density_mean = (density_particle + density_neighbor) / 2
    pi_ab = viscosity(state_equation.sound_speed, v_diff, pos_diff,
                      distance, density_mean, smoothing_length)

    grad_kernel = kernel_deriv(smoothing_kernel, distance, smoothing_length) * pos_diff / distance
    m_b = neighbor_container.mass[neighbor]
    dv_pressure = -m_b * (particle_container.pressure[particle] / density_particle^2 +
                          neighbor_container.pressure[neighbor] / density_neighbor^2) * grad_kernel
    dv_viscosity = m_b * pi_ab * grad_kernel

    dv = dv_pressure + dv_viscosity

    for i in 1:ndims(particle_container)
        du[ndims(particle_container) + i, particle] += dv[i]
    end

    return du
end


@inline function continuity_equation!(du, density_calculator::ContinuityDensity,
                                      u_particle_container, u_neighbor_container,
                                      particle, neighbor, pos_diff, distance,
                                      particle_container::FluidParticleContainer,
                                      neighbor_container)
    @unpack smoothing_kernel, smoothing_length = particle_container

    vdiff = get_particle_vel(particle, u_particle_container, particle_container) -
            get_particle_vel(neighbor, u_neighbor_container, neighbor_container)

    du[2 * ndims(particle_container) + 1, particle] += sum(neighbor_container.mass[neighbor] * vdiff *
                                             kernel_deriv(smoothing_kernel, distance, smoothing_length) .*
                                             pos_diff) / distance

    return du
end

@inline function continuity_equation!(du, density_calculator::ContinuityDensity,
                                      u_particle_container, u_neighbor_container,
                                      particle, neighbor, pos_diff, distance,
                                      particle_container::FluidParticleContainer,
                                      neighbor_container::BoundaryParticleContainer)
    @unpack smoothing_kernel, smoothing_length = particle_container

    vdiff = get_particle_vel(particle, u_particle_container, particle_container)

    du[2 * ndims(particle_container) + 1, particle] += sum(neighbor_container.mass[neighbor] * vdiff *
                                             kernel_deriv(smoothing_kernel, distance, smoothing_length) .*
                                             pos_diff) / distance

    return du
end


# Fluid-boundary and fluid-solid interaction
function interact!(du, u_particle_container, u_neighbor_container, neighborhood_search,
                   particle_container::FluidParticleContainer,
                   neighbor_container::Union{BoundaryParticleContainer, SolidParticleContainer})
    @unpack density_calculator, state_equation, viscosity, smoothing_kernel, smoothing_length = particle_container
    @unpack sound_speed = state_equation

    @threaded for particle in each_moving_particle(particle_container)
        density_a = get_particle_density(particle, u_particle_container, particle_container)
        v_a = get_particle_vel(particle, u_particle_container, particle_container)

        particle_coords = get_current_coords(particle, u_particle_container, particle_container)
        for neighbor in eachneighbor(particle_coords, neighborhood_search)
            neighbor_coords = get_current_coords(neighbor, u_neighbor_container, neighbor_container)

            pos_diff = particle_coords - neighbor_coords
            distance = norm(pos_diff)

            if sqrt(eps()) < distance <= compact_support(smoothing_kernel, smoothing_length)
                m_b = neighbor_container.mass[neighbor]

                continuity_equation!(du, density_calculator,
                                     u_particle_container, u_neighbor_container,
                                     particle, neighbor, pos_diff, distance,
                                     particle_container, neighbor_container)

                pi_ab = viscosity(sound_speed, v_a, pos_diff, distance, density_a, smoothing_length)
                dv_viscosity = m_b * pi_ab * kernel_deriv(smoothing_kernel, distance, smoothing_length) * pos_diff / distance

                dv_boundary = boundary_particle_impact(particle, particle_container, neighbor_container,
                                                       pos_diff, distance, density_a, m_b)

                dv = dv_boundary + dv_viscosity

                for i in 1:ndims(particle_container)
                    du[ndims(particle_container) + i, particle] += dv[i]
                end
            end
        end
    end

    return du
end