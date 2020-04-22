""" This file contains the geometry data"""

"""Electrochemistry"""
# target current density [A/m^2]
current_control = True
# current_density = 10000.0
# current_density = (2.23E+02, 1.11E+03, 2.23E+03,
#                    3.34E+03, 4.46E+03, 5.57E+03, 6.68E+03)

# current_density = (445.5, 891.1, 1336.6, 2227.7, 3118.7, 4009.8, 4900.9, 5569.2,
#                    6683.0, 7574.1, 8465.1, 9356.2, 10247.3, 11138.3,
#                    12000.0, 13000.0)
current_density = 12000.0
# current_density = \
#     (445.5,	891.1,	1336.6,	1782.1,	2227.7,	2673.2,	3118.7,	3564.3,	4009.8,
#      4455.3, 4900.9, 5346.4, 5569.2, 6683.0, 7128.5, 7574.1, 8019.6, 8465.1,
#      8910.7, 9356.2, 9801.7, 10247.3, 10692.8, 11138.3)

average_cell_voltage = 0.5
# open circuit voltage [V]
open_circuit_voltage = 0.96
# cathode stoichiometry
stoichiometry_cathode = 2.5
# anode stoichiometry
stoichiometry_anode = 1.34
# reaction stoichiometry
cathode_reaction_stoich = [-1.0, 0.0, 2.0]
anode_reaction_stoich = [-2.0, 0.0, 0.0]
# anode_reaction_stoich = cathode_reaction_stoich
# reaction charge number
cathode_charge_number = 2.0
anode_charge_number = 4.0

""""Thermal Settings"""
temperature = 433.15
# total mass flow of the coolant per channel [kg/s]
# coolant_mass_flow = 1.e-2
coolant_temperature_difference = 5.0
# air inlet temperature [K]
temp_air_in = temperature  # 443.15
# anode gas inlet temperature [K]
temp_anode_gas_in = temperature  # 443.15
# coolant inlet temperature [K]
temp_coolant_in = temperature - coolant_temperature_difference - 5.0  # 443.15
# environment temperature [K]
temp_environment = 298.15
# heat power of the endplates [W]
endplates_heat_power = 0.
# convection coefficient between the stack walls and the environment [W/(m^2K)]
convection_coefficient_stack_environment = 0.0
# initial temperature [K]
temp_initial = temperature  # 443.15

"""Humidification"""
# # cathode inlet gas relative humidity
# inlet_humidity_cathode = 0.
# # anode inlet gas relative humidity
# inlet_humidity_anode = 0.

"""Fluid Mechanic Settings"""
# pressure at the outlet of the cathode manifold [Pa]
p_manifold_cathode_out = 101325.0
# pressure at the outlet of the anode manifold [Pa]
p_manifold_anode_out = 101325.0

"""Fluid settings"""
# species names
# cathode_species = {'O2': 'gas', 'N2': 'gas', 'H2O': 'gas-liquid'}
cathode_species = {'O2': 'gas', 'N2': 'gas', 'H2O': 'gas-liquid'}
anode_species = {'H2': 'gas', 'N2': 'gas', 'H2O': 'gas-liquid'}
# anode_species = cathode_species

# inlet composition (molar fractions)
cathode_inlet_composition = [0.21, 0.79, 0.0]
anode_inlet_composition = [0.667, 0.333, 0.0]
# anode_inlet_composition = [1.0, 0.0, 0.0]
# anode_inlet_composition = cathode_inlet_composition


