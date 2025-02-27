"""
    write2json!(meta_data, system; output_directory="out", prefix="",
                system_name=vtkname(system_))

Write simulation metadata to a JSON file.

# Arguments
- `meta_data`:      Dictionary where metadata will be stored.
- `system`:         The simulation system whose metadata should be recorded.

# Keywords
- `output_directory="out"`: Output directory path.
- `prefix=""`:              Prefix for output files.
- `system_name=vtkname(system_)`: Name of the system, used for naming the output file.
- `git_hash=compute_git_hash()`: Git hash of the current repository.

"""
function write2json!(system; output_directory="out", prefix="",
                     system_name=vtkname(system_), git_hash=compute_git_hash())
    meta_data = Dict{String, Any}(
        "solver_version" => git_hash,
        "julia_version" => string(VERSION)
    )

    get_meta_data!(meta_data, system)

    # handle "_" on optional prefix strings
    add_opt_str_pre(str) = (str === "" ? "" : "$(str)_")

    # Write meta data to json file
    json_file = joinpath(output_directory,
                         add_opt_str_pre(prefix) * "$(system_name)_metadata.json")

    open(json_file, "w") do file
        JSON.print(file, meta_data, 2)
    end
end

function get_meta_data!(meta_data, system::FluidSystem)
    meta_data["acceleration"] = system.acceleration
    meta_data["viscosity"] = type2string(system.viscosity)
    get_meta_data!(meta_data, system.viscosity)
    meta_data["smoothing_kernel"] = type2string(system.smoothing_kernel)
    meta_data["smoothing_length"] = system.smoothing_length
    meta_data["density_calculator"] = type2string(system.density_calculator)

    if system isa WeaklyCompressibleSPHSystem
        meta_data["state_equation"] = type2string(system.state_equation)
        meta_data["state_equation_rho0"] = system.state_equation.reference_density
        meta_data["state_equation_pa"] = system.state_equation.background_pressure
        meta_data["state_equation_c"] = system.state_equation.sound_speed
        meta_data["solver"] = "WCSPH"

        meta_data["correction_method"] = type2string(system.correction)
        if system.correction isa AkinciFreeSurfaceCorrection
            meta_data["correction_rho0"] = system.correction.rho0
        end
        if system.state_equation isa StateEquationCole
            meta_data["state_equation_exponent"] = system.state_equation.exponent
        end
        if system.state_equation isa StateEquationIdealGas
            meta_data["state_equation_gamma"] = system.state_equation.gamma
        end
    else
        meta_data["solver"] = "EDAC"
        meta_data["sound_speed"] = system.sound_speed
        meta_data["background_pressure_TVF"] = system.transport_velocity isa Nothing ? "-" :
                                               system.transport_velocity.background_pressure
    end

    return meta_data
end

get_meta_data!(meta_data, viscosity::Nothing) = meta_data

function get_meta_data!(meta_data, viscosity::Union{ViscosityAdami, ViscosityMorris})
    meta_data["viscosity_nu"] = viscosity.nu
    meta_data["viscosity_epsilon"] = viscosity.epsilon
end

function get_meta_data!(meta_data, viscosity::ArtificialViscosityMonaghan)
    meta_data["viscosity_alpha"] = viscosity.alpha
    meta_data["viscosity_beta"] = viscosity.beta
    meta_data["viscosity_epsilon"] = viscosity.epsilon
end

function get_meta_data!(meta_data, system::TotalLagrangianSPHSystem)
    meta_data["young_modulus"] = system.young_modulus
    meta_data["poisson_ratio"] = system.poisson_ratio
    meta_data["lame_lambda"] = system.lame_lambda
    meta_data["lame_mu"] = system.lame_mu
    meta_data["smoothing_kernel"] = type2string(system.smoothing_kernel)
    meta_data["smoothing_length"] = system.smoothing_length

    get_meta_data!(meta_data, system.boundary_model, system)
end

function get_meta_data!(meta_data, system::OpenBoundarySPHSystem)
    meta_data["boundary_zone"] = type2string(system.boundary_zone)
    meta_data["width"] = round(system.boundary_zone.zone_width, digits=3)
    meta_data["flow_direction"] = system.flow_direction
    meta_data["velocity_function"] = type2string(system.reference_velocity)
    meta_data["pressure_function"] = type2string(system.reference_pressure)
    meta_data["density_function"] = type2string(system.reference_density)
end

function get_meta_data!(meta_data, system::BoundarySPHSystem)
    get_meta_data!(meta_data, system.boundary_model, system)
end

function get_meta_data!(meta_data, model, system)
    return meta_data
end

function get_meta_data!(meta_data, model::BoundaryModelMonaghanKajtar, system)
    meta_data["boundary_model"] = "BoundaryModelMonaghanKajtar"
    meta_data["boundary_spacing_ratio"] = model.beta
    meta_data["boundary_K"] = model.K
end

function get_meta_data!(meta_data, model::BoundaryModelDummyParticles, system)
    meta_data["boundary_model"] = "BoundaryModelDummyParticles"
    meta_data["smoothing_kernel"] = type2string(model.smoothing_kernel)
    meta_data["smoothing_length"] = model.smoothing_length
    meta_data["density_calculator"] = type2string(model.density_calculator)
    meta_data["state_equation"] = type2string(model.state_equation)
    meta_data["viscosity_model"] = type2string(model.viscosity)
end

function get_meta_data!(meta_data, system::BoundaryDEMSystem)
    return meta_data
end