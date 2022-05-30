import FESTIM as F
import properties
from decay_output import DecayXDMF

id_W = 6  # volume W
id_Cu = 7  # volume Cu
id_CuCrZr = 8  # volume CuCrZr
id_W_top = 9
id_coolant = 10
id_poloidal_gap = 11
id_toroidal_gap = 12
id_bottom = 13
id_top_pipe = 14

my_model = F.Simulation()


my_model.mesh = F.MeshFromXDMF("mesh/mesh_cells.xdmf", "mesh/mesh_facets.xdmf")


tungsten = F.Material(
    id=id_W,
    D_0=4.1e-7,
    E_D=0.39,
    S_0=1.87e24,
    E_S=1.04,
    thermal_cond=properties.thermal_cond_W,
    heat_capacity=properties.rhoCp_W,
    rho=1,
)

copper = F.Material(
    id=id_Cu,
    D_0=6.6e-7,
    E_D=0.387,
    S_0=3.14e24,
    E_S=0.572,
    thermal_cond=properties.thermal_cond_Cu,
    heat_capacity=properties.rhoCp_Cu,
    rho=1,
)

cucrzr = F.Material(
    id=id_CuCrZr,
    D_0=3.92e-7,
    E_D=0.418,
    S_0=4.28e23,
    E_S=0.387,
    thermal_cond=properties.thermal_cond_CuCrZr,
    heat_capacity=properties.rhoCp_CuCrZr,
    rho=1,
)

my_model.materials = F.Materials([tungsten, copper, cucrzr])

# traps
my_model.traps = F.Traps(
    [
        F.Trap(
            k_0=8.96e-17,
            E_k=0.39,
            p_0=1e13,
            E_p=0.87,
            density=1.3e-3 * properties.atom_density_W,
            materials=tungsten,
        ),
        F.Trap(
            k_0=[8.96e-17, 6e-17, 1.2e-16],
            E_k=[0.39, 0.39, 0.42],
            p_0=[1e13, 8e13, 8e13],
            E_p=[1.0, 0.5, 0.85],
            density=[
                4e-4 * properties.atom_density_W,
                5e-5 * properties.atom_density_Cu,
                5e-5 * properties.atom_density_CuCrZr,
            ],
            materials=[tungsten, copper, cucrzr],
        ),
    ]
)

# temperature
my_model.T = F.HeatTransferProblem(transient=False)


# boundary conditions
heat_flux_top = F.FluxBC(surfaces=id_W_top, value=10e6, field="T")
convective_heat_flux_coolant = F.ConvectiveFlux(
    h_coeff=7e04, T_ext=323, surfaces=id_coolant
)

heat_transfer_bcs = [heat_flux_top, convective_heat_flux_coolant]

instantaneous_recombination_poloidal = F.DirichletBC(value=0, surfaces=id_poloidal_gap)
instantaneous_recombination_toroidal = F.DirichletBC(value=0, surfaces=id_toroidal_gap)
instantaneous_recombination_bottom = F.DirichletBC(value=0, surfaces=id_bottom)
instantaneous_recombination_top_pipe = F.DirichletBC(value=0, surfaces=id_top_pipe)

recombination_flux_coolant = F.RecombinationFlux(
    Kr_0=2.9e-14, E_Kr=1.92, order=2, surfaces=id_coolant
)
h_implantation_top = F.ImplantationDirichlet(
    surfaces=id_W_top, phi=1.61e22, R_p=9.52e-10, D_0=4.1e-7, E_D=0.39
)

h_transport_bcs = [
    h_implantation_top,
    recombination_flux_coolant,
    instantaneous_recombination_poloidal,
    instantaneous_recombination_toroidal,
    instantaneous_recombination_bottom,
    instantaneous_recombination_top_pipe,
]


my_model.boundary_conditions = heat_transfer_bcs + h_transport_bcs

my_model.settings = F.Settings(
    absolute_tolerance=1e4,
    relative_tolerance=1e-5,
    maximum_iterations=15,
    traps_element_type="DG",
    chemical_pot=True,
    transient=False,
    linear_solver="mumps",
)

folder = "results"
my_model.exports = F.Exports(
    [
        F.XDMFExport("T", folder=folder),
        F.XDMFExport("solute", folder=folder),
        F.XDMFExport("retention", folder=folder),
        DecayXDMF("decay", folder=folder)
    ]
)

my_model.initialise()


# modify formulation
# dc/dt = ..... - decay_constant*c
tritium_half_life = 12.4  # years
tritium_half_life *= 364.25  # days
tritium_half_life *= 24  # hours
tritium_half_life *= 3600  # seconds

decay_constant = 0.69 / tritium_half_life  # ln(2)/half_life  in s-1

decay_mobile = (
    decay_constant
    * my_model.h_transport_problem.mobile.solution
    * my_model.h_transport_problem.mobile.test_function
    * my_model.mesh.dx
)
decay_traps = [
    decay_constant * trap.solution * trap.test_function * my_model.mesh.dx
    for trap in my_model.traps.traps
]

my_model.h_transport_problem.F += decay_mobile + sum(decay_traps)


my_model.run()
