""" Operating conditions """

"""Electrochemistry settings"""
# if current_control = True, provided current_density value or array is used
# as operating point
# if current_control = False, provided average_cell_voltage value or array is
# used as operating point (not working as well, so current_control should be
# used at the moment)
current_control = True

# target current density [A/m^2]
current_density = 10000.0
# current_density = \
#     (445.5,	891.1,	1336.6,	1782.1,	2227.7,	2673.2,	3118.7,	3564.3,	4009.8,
#      4455.3, 4900.9, 5346.4, 5569.2, 6683.0, 7128.5, 7574.1, 8019.6, 8465.1,
#      8910.7, 9356.2, 9801.7, 10247.3, 10692.8, 11138.3)

# average cell voltage [V] (not used, when current_control = True)
# average_cell_voltage = 0.5

# open circuit voltage [V]
open_circuit_voltage = 0.96

# cathode stoichiometry
stoichiometry_cathode = 2.0

# anode stoichiometry
stoichiometry_anode = 1.5

# reaction stoichiometry
cathode_reaction_stoich = [-1.0, 0.0, 2.0]
anode_reaction_stoich = [-2.0, 0.0, 0.0]

# reaction charge number
cathode_charge_number = 4.0
anode_charge_number = 4.0

"""Fluid settings"""
# components
cathode_species = {'O2': 'gas', 'N2': 'gas', 'H2O': 'gas-liquid'}
anode_species = {'H2': 'gas', 'N2': 'gas', 'H2O': 'gas-liquid'}

# inlet composition (molar fractions)
cathode_inlet_composition = [0.21, 0.79, 0.0]
anode_inlet_composition = [0.667, 0.333, 0.0]


""""Thermal Settings"""
# temperature value used for most boundary conditions (see below and edit
# as required) [K]
temperature = 433.15

# set the total mass flow of the coolant for the stack [kg/s]
# coolant_mass_flow = 1.e-2
# if coolant_mass_flow is not provided, set the desired coolant temperature
# difference [K] (mass flow overrides the temperature difference, if provided)
coolant_temperature_difference = 5.0

# coolant pinch point temperature difference [K]
temp_pinch = 10.0

# reactant inlet temperatures [K]
temp_cathode_in = temperature  # 443.15
temp_anode_in = temperature  # 443.15

# coolant inlet temperature [K]
temp_coolant_in = temperature - coolant_temperature_difference - temp_pinch

# environment temperature [K]
temp_environment = 298.15

# heat boundary conditions at the endplates [W]
endplates_heat_power = 0.0

# convection coefficient between the stack walls and the environment [W/(m^2K)]
convection_coefficient_environment = 10.0

# initial temperature [K]
temp_initial = temperature  # 443.15


"""Pressure settings"""
# pressure at the outlet of the manifolds [Pa]
p_manifold_cathode_out = 101325.0
# pressure at the outlet of the anode manifold [Pa]
p_manifold_anode_out = 101325.0




