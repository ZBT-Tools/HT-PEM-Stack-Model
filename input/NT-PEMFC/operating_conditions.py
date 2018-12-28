import numpy as np


""" This file contains the geometry data"""

"""Electrochemistry"""
# target current density [A/m^2]
a = [3.17782859, 2.998470663, 2.798749363, 2.598393607,
                         2.396849473, 2.197724053, 1.997807213, 1.796569687,
                         1.5972518, 1.396593663, 1.196849567, 0.997384197,
                         0.796733187, 0.597545647, 0.49775168, 0.396710477,
                         0.29834263, 0.1982076, 0.148961347, 0.098618743,
                         0.04888981, 0.01861823, 0.000440543]
target_current_density = np.array(a)*1.e4 #np.linspace(1., 31778.2859, 100)
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

"""Stack Settings"""
# number of the pemfc
cell_number = 3
# coolant channel configuration (no cooling channel at the endplates = False)
cooling_bc = True
