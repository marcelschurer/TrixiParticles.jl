function write_simulation_script(callback)
    starting_file = joinpath(examples_dir(), "test_file", "moving_wall_2d.jl")
    #TODO: Pfad muss gespeichert werden. Eigener Callback? oder in SolutionSavingCallback übergeben?
    saving_file = joinpath(pkgdir(TrixiParticles), callback.output_directory,
                           "simulation_script.jl")
    # TODO: Make filename settable

    cp(starting_file, saving_file; force=true)

    return
end

function restart_simulation(input_directory, output_directory, script_name, start_iter, end_time)
    starting_file = joinpath(input_directory, script_name * ".jl")

    trixi_include(@__MODULE__, starting_file, sol=nothing)
    ode = Main.ode
    #@autoinfiltrate
    #semi = ode.p

    # Simulationsdaten des benötigten Zeitpunkts mit `vtk2trixi()`laden.
    fluid = vtk2trixi(joinpath(input_directory, "1_fluid_1_" * string(start_iter) * ".vtu"))

    # Erzeugen der neuen "Start"-Daten
    v0_ode = vec(TrixiParticles.vcat(fluid.velocity, fluid.density'))
    u0_ode = vec(fluid.coordinates)

    ode.u0.x[1] .= v0_ode
    ode.u0.x[2] .= u0_ode

    start_time = fluid.time
    ode.tspan = (start_time, end_time)

    saving_callback = SolutionSavingCallback(dt=0.02, prefix="2",
                                             output_directory=output_directory,
                                             save_initial_solution=true,
                                             save_script=false)

    # "simulation_script.jl" mit veränderten Werten ausführen.
    trixi_include(starting_file, ode=ode, saving_callback=saving_callback)

    return
end
