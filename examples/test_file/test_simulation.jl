using TrixiParticles
using OrdinaryDiffEq

# Simulation restarten
input_directory = joinpath(pkgdir(TrixiParticles), "out")
script_name = "simulation_script"

output_directory = input_directory
start_iter = 50
end_time = 2.0

restart_simulation(input_directory, output_directory, script_name, start_iter, end_time)
