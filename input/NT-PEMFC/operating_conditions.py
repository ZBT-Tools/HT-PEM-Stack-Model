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
stoichiometry_anode = 5.
# pem-type (True = HT-PEM, False = NT-PEM)
pem_type = False

""""Thermal Settings"""
# air inlet temperature [K]
temp_air_in = 298.15
# anode gas inlet temperature [K]
temp_anode_gas_in = 298.15
# coolant inlet temperature [K]
temp_coolant_in = 298.15
# environment temperature [K]
temp_environment = 298.15
# heat power of the endplates [W]
endplates_heat_power = 0.
# initial temperature [K]
temp_initial = 298.15

"""Fluid Mechanic Settings"""
# pressure at the outlet of the cathode manifold [Pa]
p_manifold_cathode_out = 1.e5
# pressure at the outlet of the anode manifold [Pa]
p_manifold_anode_out = 1.e5
# mass flow of the coolant per channel [kg/s]
mass_flow_coolant = 1.e-4

"""Cell Settings"""
# number of gas channels
gas_channel_number = 10
# number of coolant channels
coolant_channel_number = 10

"""Stack Settings"""
# number of the pemfc
cell_number = 40
# coolant channel configuration (no cooling channel at the endplates = False)
cooling_bc = True
