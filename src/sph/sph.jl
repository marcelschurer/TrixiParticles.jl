struct SummationDensity end

struct ContinuityDensity end

struct SPHSemidiscretization{NDIMS, ELTYPE<:Real, DC, SE, K, V, BC, C}
    density_calculator  ::DC
    state_equation      ::SE
    smoothing_kernel    ::K
    smoothing_length    ::ELTYPE
    viscosity           ::V
    boundary_conditions ::BC
    gravity             ::SVector{NDIMS, ELTYPE}
    cache               ::C

    function SPHSemidiscretization{NDIMS}(particle_masses,
                                          density_calculator, state_equation,
                                          smoothing_kernel, smoothing_length;
                                          viscosity=NoViscosity(),
                                          boundary_conditions=nothing,
                                          gravity=ntuple(_ -> 0.0, Val(NDIMS))) where NDIMS
        ELTYPE = eltype(particle_masses)
        nparticles = length(particle_masses)

        # Make boundary_conditions a tuple
        boundary_conditions_ = digest_boundary_conditions(boundary_conditions)

        # Make gravity an SVector
        gravity_ = SVector(gravity...)

        cache = (; create_cache(particle_masses, density_calculator, ELTYPE, nparticles)...)

        return new{NDIMS, ELTYPE, typeof(density_calculator), typeof(state_equation),
                   typeof(smoothing_kernel), typeof(viscosity),
                   typeof(boundary_conditions_), typeof(cache)}(
            density_calculator, state_equation, smoothing_kernel, smoothing_length,
            viscosity, boundary_conditions_,  gravity_, cache)
    end
end


function create_cache(mass, density_calculator, eltype, nparticles)
    pressure = Vector{eltype}(undef, nparticles)

    return (; mass, pressure, create_cache(density_calculator, eltype, nparticles)...)
end

function create_cache(::SummationDensity, eltype, nparticles)
    density = Vector{eltype}(undef, nparticles)

    return (; density)
end

function create_cache(::ContinuityDensity, eltype, nparticles)
    return (; )
end


function semidiscretize(semi::SPHSemidiscretization{NDIMS, ELTYPE, SummationDensity},
                        particle_coordinates, particle_velocities, tspan) where {NDIMS, ELTYPE}

    u0 = Array{eltype(particle_coordinates), 2}(undef, 2 * ndims(semi), nparticles(semi))

    for particle in eachparticle(semi)
        # Set particle coordinates
        for dim in 1:ndims(semi)
            u0[dim, particle] = particle_coordinates[dim, particle]
        end

        # Set particle velocities
        for dim in 1:ndims(semi)
            u0[dim + ndims(semi), particle] = particle_velocities[dim, particle]
        end
    end

    # Compute quantities like density and pressure
    compute_quantities!(u0, semi)

    return ODEProblem(rhs!, u0, tspan, semi)
end


function semidiscretize(semi::SPHSemidiscretization{NDIMS, ELTYPE, ContinuityDensity},
                        particle_coordinates, particle_velocities, particle_densities, tspan) where {NDIMS, ELTYPE}

    u0 = Array{eltype(particle_coordinates), 2}(undef, 2 * ndims(semi) + 1, nparticles(semi))

    for particle in eachparticle(semi)
        # Set particle coordinates
        for dim in 1:ndims(semi)
            u0[dim, particle] = particle_coordinates[dim, particle]
        end

        # Set particle velocities
        for dim in 1:ndims(semi)
            u0[dim + ndims(semi), particle] = particle_velocities[dim, particle]
        end

        # Set particle densities
        u0[2 * ndims(semi) + 1, particle] = particle_densities[particle]
    end

    # Compute quantities like pressure
    compute_quantities!(u0, semi)

    return ODEProblem(rhs!, u0, tspan, semi)
end


function compute_quantities!(u, semi::SPHSemidiscretization{NDIMS, ELTYPE, SummationDensity}) where {NDIMS, ELTYPE}
    @unpack smoothing_kernel, smoothing_length, state_equation, cache = semi
    @unpack mass, density, pressure = cache

    for particle in eachparticle(semi)
        @pixie_timeit timer() "Compute density" begin
            density[particle] = sum(eachparticle(semi)) do neighbor
                distance = norm(get_particle_coords(u, semi, particle) -
                                get_particle_coords(u, semi, neighbor))

                if distance > compact_support(smoothing_kernel, smoothing_length)
                    return 0.0
                end

                return mass[neighbor] * kernel(smoothing_kernel, distance, smoothing_length)
            end
        end

        pressure[particle] = state_equation(density[particle])
    end
end

function compute_quantities!(u, semi::SPHSemidiscretization{NDIMS, ELTYPE, ContinuityDensity}) where {NDIMS, ELTYPE}
    @unpack density_calculator, state_equation, cache = semi
    @unpack pressure = cache

    for particle in eachparticle(semi)
        pressure[particle] = state_equation(get_particle_density(u, cache, density_calculator, particle))
    end
end


function rhs!(du, u, semi, t)
    @unpack smoothing_kernel, smoothing_length,
            density_calculator, state_equation, viscosity,
            boundary_conditions, gravity, cache = semi

    @pixie_timeit timer() "rhs!" begin
        @unpack mass, pressure = cache

        compute_quantities!(u, semi)

        # u[1:3] = coordinates
        # u[4:6] = velocity
        for particle in eachparticle(semi)
            # dr = v
            for i in 1:ndims(semi)
                du[i, particle] = u[i + ndims(semi), particle]
            end

            # dv (constant smoothing length, Price (31))
            r1 = get_particle_coords(u, semi, particle)
            @pixie_timeit timer() "Compute dv" begin
                dv = sum(eachparticle(semi)) do neighbor
                    m = mass[neighbor]
                    r2 = get_particle_coords(u, semi, neighbor)

                    pos_diff = r1 - r2
                    distance = norm(pos_diff)

                    if eps() < distance <= compact_support(smoothing_kernel, smoothing_length)
                        density_particle = get_particle_density(u, cache, density_calculator, particle)
                        density_neighbor = get_particle_density(u, cache, density_calculator, particle)

                        # Viscosity
                        v_diff = get_particle_vel(u, semi, particle) - get_particle_vel(u, semi, neighbor)
                        density_diff = (density_particle + density_neighbor) / 2
                        pi_ij = viscosity(state_equation.sound_speed, v_diff, pos_diff,
                                          distance, density_diff, smoothing_length)

                        result = -m * (pressure[particle] / density_particle^2 +
                                      pressure[neighbor] / density_neighbor^2 + pi_ij) *
                            kernel_deriv(smoothing_kernel, distance, smoothing_length) * pos_diff / distance
                    else
                        # Don't compute pressure and density terms, just return zero
                        result = zeros(SVector{ndims(semi), eltype(pressure)})
                    end

                    return result
                end
            end

            for i in 1:ndims(semi)
                # Gravity
                du[i + ndims(semi), particle] = dv[i] + gravity[i]
            end
        end

        continuity_equation!(du, u, semi)

        # Boundary conditions
        @pixie_timeit timer() "Boundary conditions" begin
            for boundary_condition in boundary_conditions
                calc_boundary_condition!(du, u, boundary_condition, semi)
            end
        end
    end

    return du
end


function continuity_equation!(du, u, semi::SPHSemidiscretization{NDIMS, ELTYPE, ContinuityDensity}) where {NDIMS, ELTYPE}
    @unpack smoothing_kernel, smoothing_length, cache = semi
    @unpack mass = cache

    for particle in eachparticle(semi)
        r1 = get_particle_coords(u, semi, particle)
        @pixie_timeit timer() "Compute drho" begin
            du[2 * ndims(semi) + 1, particle] = sum(eachparticle(semi)) do neighbor
                m = mass[neighbor]
                r2 = get_particle_coords(u, semi, neighbor)

                diff = r1 - r2
                distance = norm(diff)

                if eps() < distance <= compact_support(smoothing_kernel, smoothing_length)
                    vdiff = get_particle_vel(u, semi, particle) -
                            get_particle_vel(u, semi, neighbor)

                    result = sum(m * vdiff * kernel_deriv(smoothing_kernel, distance, smoothing_length) .* diff) / distance
                else
                    # Don't compute pressure and density terms, just return zero
                    result = 0.0
                end

                return result
            end
        end
    end

    return du
end

function continuity_equation!(du, u, ::SPHSemidiscretization{NDIMS, ELTYPE, SummationDensity}) where {NDIMS, ELTYPE}
    return du
end


function calc_boundary_condition!(du, u, boundary_condition::BoundaryConditionMonaghanKajtar, semi)
    @unpack smoothing_kernel, smoothing_length,
            density_calculator, state_equation, viscosity, cache = semi
    @unpack K, coordinates, mass, spacing = boundary_condition

    for particle in eachparticle(semi)
        dv = sum(eachparticle(boundary_condition)) do boundary_particle
            pos_diff = get_particle_coords(u, semi, particle) -
                       get_boundary_coords(boundary_condition, semi, boundary_particle)
            distance = norm(pos_diff)

            if eps() < distance <= compact_support(smoothing_kernel, smoothing_length)
                # Viscosity
                v_diff = get_particle_vel(u, semi, particle)
                pi_ij = viscosity(state_equation.sound_speed, v_diff, pos_diff, distance,
                                  get_particle_density(u, cache, density_calculator, particle), smoothing_length)

                m_b = mass[boundary_particle]

                return K * spacing[boundary_particle] * pos_diff / distance^2 *
                    kernel(smoothing_kernel, distance, smoothing_length) * 2 * m_b / (cache.mass[particle] + m_b) -
                    kernel_deriv(smoothing_kernel, distance, smoothing_length) * m_b * pi_ij * pos_diff / distance
            else
                return zeros(SVector{ndims(semi), eltype(cache.mass)})
            end
        end

        for i in 1:ndims(semi)
            du[i + ndims(semi), particle] += dv[i]
        end
    end

    return du
end


@inline function get_particle_coords(u, semi, particle)
    return SVector(ntuple(@inline(dim -> u[dim, particle]), Val(ndims(semi))))
end

@inline function get_particle_vel(u, semi, particle)
    return SVector(ntuple(@inline(dim -> u[dim + ndims(semi), particle]), Val(ndims(semi))))
end


@inline function get_particle_density(u, cache, ::SummationDensity, particle)
    return cache.density[particle]
end

@inline function get_particle_density(u, cache, ::ContinuityDensity, particle)
    return u[end, particle]
end


@inline function get_boundary_coords(boundary_container, semi, particle)
    @unpack coordinates = boundary_container
    SVector(ntuple(@inline(dim -> coordinates[dim, particle]), Val(ndims(semi))))
end


# This can be used both for Semidiscretization or boundary container types
@inline eachparticle(container) = Base.OneTo(nparticles(container))
@inline nparticles(semi) = length(semi.cache.mass)
@inline Base.ndims(::SPHSemidiscretization{NDIMS}) where NDIMS = NDIMS
