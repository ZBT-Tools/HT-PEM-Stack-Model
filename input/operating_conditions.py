import numpy as np


""" This file contains the geometry data"""

"""Electrochemistry"""
# target current density [A/m^2]
target_current_density = 6000.
# open circuit voltage [V]
open_circuit_voltage = 1.00
# cathode stoichiometry
stoichiometry_cathode = 2.
# anode stoichiometry
stoichiometry_anode = 2.
# pem-type (True = HT-PEM, False = NT-PEM)
is_ht_pem = False

""""Thermal Settings"""
# mass flow of the coolant per channel [kg/s]
mass_flow_coolant = 1.e-2
# air inlet temperature [K]
temp_air_in = 340.
# anode gas inlet temperature [K]
temp_anode_gas_in = 340.
# coolant inlet temperature [K]
temp_coolant_in = 340.
# environment temperature [K]
temp_environment = 298.15
# heat power of the endplates [W]
endplates_heat_power = 50.
# convection coefficient between the stack walls and the environment [W/(m^2K)]
convection_coefficient_stack_environment = 10.0
# initial temperature [K]
temp_initial = 350.

"""Humidification"""
# cathode inlet gas relative humidity
inlet_humidity_cathode = 0.
# anode inlet gas relative humidity
inlet_humidity_anode = 0.

"""Fluid Mechanic Settings"""
# pressure at the outlet of the cathode manifold [Pa]
p_manifold_cathode_out = 1.e5
# pressure at the outlet of the anode manifold [Pa]
p_manifold_anode_out = 1.e5

"""Species settings"""
# species names
cathode_species = {'O2': 'gas', 'N2': 'gas', 'H2O': 'gas-liquid'}
anode_species = {'H2': 'gas', 'N2': 'gas', 'H2O': 'gas-liquid'}
# molar masses (kg/mol)
cathode_molar_mass = [0.032, 0.028, 0.018]
anode_molar_mass = [0.002, 0.028, 0.018]
# inlet composition (molar fractions)
cathode_inlet_composition = [0.21, 0.79, 0.0]
anode_inlet_composition = [0.5, 0.5, 0.0]
# reaction stoichiometry
cathode_reaction_stoich = [-1.0, 0.0, 2.0]
anode_reaction_stoich = [-2.0, 0.0, 0.0]
# reaction charge number
cathode_charge_number = 4.0
anode_charge_number = 4.0
