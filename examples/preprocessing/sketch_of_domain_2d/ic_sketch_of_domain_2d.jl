using TrixiParticles

particle_spacing = 0.05
fluid_density = 1000.0
pressure = 1000.0
domain_inlet_size = (1.0, 0.4)
domain_outlet_size = (1.0, 2.0)
boundary_layers = 3
boundary_size_inlet = (domain_inlet_size[1] + 2 * particle_spacing * open_boundary_layers,
                       domain_inlet_size[2])

# file_boundary_out = joinpath("examples", "preprocessing", "sketch_of_domain_2d",
#                              "boundary_out.asc")
# shape_boundary_out = load_shape(file_boundary_out)
# ic_boundary_out = ComplexShape(shape_boundary_out; particle_spacing, density=fluid_density)

# file_boundary_in = joinpath("examples", "preprocessing", "sketch_of_domain_2d",
#                             "boundary_in.asc")
# shape_boundary_in = load_shape(file_boundary_in)
# ic_boundary_in = ComplexShape(shape_boundary_in; particle_spacing, density=fluid_density)

# ic_boundary = setdiff(ic_boundary_out, ic_boundary_in)
# ic_boundary = ic_boundary_out
# ic_boundary = ic_boundary_in

domain_inlet = RectangularTank(particle_spacing, domain_inlet_size, boundary_size,
                               fluid_density,
                               pressure=pressure, n_layers=boundary_layers,
                               faces=(false, false, true, true))

domain_outlet = RectangularTank(particle_spacing, domain_outlet_size, boundary_size,
                                fluid_density,
                                pressure=pressure, n_layers=boundary_layers,
                                min_coordinates=[
                                    domain_inlet_size[1] * 2,
                                    domain_inlet_size[2] * 2
                                ],
                                faces=(false, false, true, true))

domain = union(domain_inlet.boundary, domain_outlet.boundary)

trixi2vtk(domain)
