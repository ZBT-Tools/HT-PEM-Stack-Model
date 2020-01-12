""" This file contains the geometry settings"""

"""Gas Channel Geometry"""
# channel length [m]
channel_length = 0.651
# channel width [m]
channel_width = 1.e-3
# channel height [m]
channel_height = 1.e-3
# rack width [m]
#rib_width = 1.e-3
# number of channels
gas_channel_number = 10.
# channel bends [n]
channel_bends = 48
# bend pressure loss coefficient of the channel bends
bend_pressure_loss_coefficient = 0.1
# flow direction in cathode channel along x-axis
cathode_flow_direction = 1
# flow direction in anode channel along x-axis
anode_flow_direction = -1

"""Coolant Channel Geometry"""
# channel length [m]
coolant_channel_length = 0.651
# height of the coolant channel [m]
coolant_channel_height = 1.e-3
# width of the coolant channel [m]
coolant_channel_width = 1.e-3
# number of coolant channels
coolant_channel_number = 10.
# channel bends [n]
coolant_channel_bends = 48.
# bend pressure loss coefficient of the channel bends
coolant_bend_pressure_loss_coefficient = 0.1


""""Cell Geometry """
# thickness of the membrane [m]
membrane_thickness = 50.e-6
# thickness of the catalyst layer [m]
catalyst_layer_thickness = 10.e-6
# catalyst layer porosity (ionomer considered as porous volume) [-]
catalyst_layer_porosity = 0.5
# thickness of the gas diffusion layer [m]
gas_diffusion_layer_thickness = 250.e-6
# gas diffusion layer porosity [-]
gas_diffusion_layer_porosity = 0.8
# thickness of the bipolar plate [m]
bipolar_plate_thickness = 1.11275e-3
# length of the cell, a side of the active area [m]
cell_length = 100.e-3
# height of the cell, b side of the active area [m]
cell_width = 100.e-3


"""Manifold Geometry"""
# Configuration: U- or Z-shape
anode_manifold_configuration = 'U'
cathode_manifold_configuration = 'U'
# manifold height [m]
manifold_height = 10.5e-3
# manifold width [m]
manifold_width = 5.5e-3
# geometrical pressure loss coefficient of the manifold header
manifold_pressure_loss_coefficient = 0.1

"""Stack Settings"""
# number of the pemfc
cell_number = 3
# coolant channel configuration (no cooling channel at the endplates = False)
cooling_bc = True
