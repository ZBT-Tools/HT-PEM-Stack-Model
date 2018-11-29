""" This file contains the geometry data"""

"""Electrochemistry"""
# target current density [A/m^2]
target_current_density = [10., 100., 1000., 3000., 4000., 5000., 6000., 7000.]
# open circuit voltage [V]
open_circuit_voltage = 1.0
# cathode stoichiometry
stoichiometry_cathode = 1.6
# anode stoichiometry
stoichiometry_anode = 1.6
# pem-type (True = HT-PEM, False = NT-PEM)
pem_type = True

""""Thermal Settings"""
# air inlet temperature [K]
temp_air_in = 410.15
# anode gas inlet temperature [K]
temp_anode_gas_in = 410.15
# coolant inlet temperature [K]
temp_coolant_in = 410.15
# environment temperature [K]
temp_environment = 298.15
# heat power of the endplates [W]
endplates_heat_power = 0.e0
# initial temperature [K]
temp_initial = 433.15

"""Fluid Mechanic Settings"""
# pressure at the outlet of the cathode manifold [Pa]
p_manifold_cathode_out = 1.e5
# pressure at the outlet of the anode manifold [Pa]
p_manifold_anode_out = 1.e5
# mass flow of the coolant per channel [kg/s]
mass_flow_coolant = 1.e-6

"""Stack Settings"""
# number of the pemfc
cell_number = 5
# coolant channel configuration (no cooling channel at the endplates = False)
cooling_bc = True
