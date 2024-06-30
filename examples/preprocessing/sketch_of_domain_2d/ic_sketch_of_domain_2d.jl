using TrixiParticles

particle_spacing = 0.05
fluid_density = 1000.0
pressure = 1000.0
domain_inlet_size = (1.0, 0.4)
domain_outlet_size = (1.0, 2.0)
boundary_layers = 3

file_boundary_out = joinpath("examples", "preprocessing", "sketch_of_domain_2d",
                             "boundary_out.asc")
shape_boundary_out = load_shape(file_boundary_out)
ic_boundary_out = ComplexShape(shape_boundary_out; particle_spacing, density=fluid_density)

# file_boundary_in = joinpath("examples", "preprocessing", "sketch_of_domain_2d",
#                             "boundary_in.asc")
# shape_boundary_in = load_shape(file_boundary_in)
# ic_boundary_in = ComplexShape(shape_boundary_in; particle_spacing, density=fluid_density)

# ic_boundary = setdiff(ic_boundary_out, ic_boundary_in)
# ic_boundary = ic_boundary_out
# ic_boundary = ic_boundary_in

domain = ic_boundary_out

trixi2vtk(domain)
