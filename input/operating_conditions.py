""" This file contains the geometry data"""

"""Electrochemistry"""
# target current density [A/m^2]
target_current_density = [10.,100.,500.,1000.,2000.,3000.,4000.,5000.,6000.,7000.]
# open circuit volatage [V]
open_circuit_voltage = 1.0
# cathode stoichiometry
stoichiometry_cathode = 1.6
# anode stoichiometry
stoichiometry_anode = 1.6
# pem-type (True = HT-PEM, False = NT-PEM)
pem_type = True

""""Thermal Settings"""
# inlet temperature cathode gas [K]
t_cathode_in = 410.15
# inlet temperature anode gas [K]
t_anode_in = 410.15
# inlet temperature coolant [K]
t_coolant_in = 410.15
# environment temperature [K]
t_environment = 298.15
# heat power of the endplates [W]
heat_power_endplates = 0.e0
# initial temperature [K]
t_initial = 433.15

"""Fluid Mechanic Settings"""
# outlet pressure cathode gas [Pa]
p_cathode_basic = 1.e5
# outlet pressure anode gas [Pa]
p_anode_basic = 1.e5
# mass flow of the coolant per channel [kg/s]
mass_flow_coolant = 1.e-6

"""Stack Settings"""
# number of cells
cell_number = 5
# coolant channel configuration (no cooling channel at the endplates = False)
cooling_bc = True
