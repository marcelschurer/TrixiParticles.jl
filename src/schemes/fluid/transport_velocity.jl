struct TransportVelocityAdami{ELTYPE}
    background_pressure::ELTYPE
    function TransportVelocityAdami(background_pressure)
        new{typeof(background_pressure)}(background_pressure)
    end
end

@inline update_transport_velocity!(system, v_ode, semi) = system

@inline function update_transport_velocity!(system::FluidSystem, v_ode, semi)
    update_transport_velocity!(system, v_ode, semi, system.transport_velocity)
end

@inline update_transport_velocity!(system, v_ode, semi, ::Nothing) = system

@inline function update_transport_velocity!(system, v_ode, semi, ::TransportVelocityAdami)
    v = wrap_v(v_ode, system, semi)

    for particle in each_moving_particle(system)
        for i in 1:ndims(system)
            v[ndims(system) + i, particle] = v[i, particle]
        end
    end

    return system
end

# Add momentum velocity.
@inline function add_velocity!(du, v, particle, system, ::Nothing)
    for i in 1:ndims(system)
        du[i, particle] = v[i, particle]
    end

    return du
end

# Add advection velocity.
@inline function add_velocity!(du, v, particle, system, ::TransportVelocityAdami)
    for i in 1:ndims(system)
        du[i, particle] = v[ndims(system) + i, particle]
    end

    return du
end

@inline function advection_velocity(v, system, particle)
    return SVector(ntuple(@inline(dim->v[ndims(system) + dim, particle]), ndims(system)))
end

@inline function momentum_convection(system, neighbor_system,
                                     v_particle_system, v_neighbor_system, rho_a, rho_b,
                                     particle, neighbor, grad_kernel, volume_term)
    return SVector(ntuple(_ -> 0.0, Val(ndims(system))))
end

@inline function momentum_convection(system::FluidSystem,
                                     neighbor_system::FluidSystem,
                                     v_particle_system, v_neighbor_system, rho_a, rho_b,
                                     particle, neighbor, grad_kernel, volume_term)
    momentum_convection(system, neighbor_system, system.transport_velocity,
                        v_particle_system, v_neighbor_system, rho_a, rho_b,
                        particle, neighbor, grad_kernel, volume_term)
end

@inline function momentum_convection(system, neighbor_system, ::Nothing,
                                     v_particle_system, v_neighbor_system, rho_a, rho_b,
                                     particle, neighbor, grad_kernel, volume_term)
    return SVector(ntuple(_ -> 0.0, Val(ndims(system))))
end

@inline function momentum_convection(system, neighbor_system, ::TransportVelocityAdami,
                                     v_particle_system, v_neighbor_system, rho_a, rho_b,
                                     particle, neighbor, grad_kernel, volume_term)
    momentum_velocity_a = current_velocity(v_particle_system, system, particle)
    advection_velocity_a = advection_velocity(v_particle_system, system, particle)

    momentum_velocity_b = current_velocity(v_neighbor_system, neighbor_system, neighbor)
    advection_velocity_b = advection_velocity(v_neighbor_system, neighbor_system, neighbor)

    A_a = rho_a * momentum_velocity_a * (advection_velocity_a - momentum_velocity_a)'
    A_b = rho_b * momentum_velocity_b * (advection_velocity_b - momentum_velocity_b)'

    return volume_term * (0.5 * (A_a + A_b)) * grad_kernel
end

@inline function transport_velocity!(dv, system::FluidSystem, volume_term,
                                     grad_kernel, particle)
    transport_velocity!(dv, system, system.transport_velocity, volume_term, grad_kernel,
                        particle)
end

@inline transport_velocity!(dv, system, volume_term, grad_kernel, particle) = dv

@inline transport_velocity!(dv, system, ::Nothing, volume_term, grad_kernel, particle) = dv

@inline function transport_velocity!(dv, system, ::TransportVelocityAdami, volume_term,
                                     grad_kernel, particle)
    (; transport_velocity) = system
    (; background_pressure) = transport_velocity
    n_dims = ndims(system)

    for dim in 1:n_dims
        dv[n_dims + dim, particle] -= volume_term * background_pressure * grad_kernel[dim]
    end

    return dv
end
