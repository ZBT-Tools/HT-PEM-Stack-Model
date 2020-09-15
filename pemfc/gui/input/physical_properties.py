# Thermal conductivities
thermal_conductivity_bpp = \
    {'label': 'Bipolar Plate Thermal Conductivity:', 'value': [100, 100],
     'sim_name': [['anode', 'thermal_conductivity_bpp'],
                  ['cathode', 'thermal_conductivity_bpp']],
     'dtype': 'float', 'dimensions': 'W/(m K)', 'type': 'EntrySet'}

thermal_conductivity_gde = \
    {'label': 'Gas Diffusion Electrode Thermal Conductivity:',
     'value': [[10.0, 10.0], [10.0, 10.0]],
     'sim_name': [['anode', 'thermal_conductivity_gde'],
                  ['cathode', 'thermal_conductivity_gde']],
     'dtype': 'float', 'dimensions': 'W/(m K)', 'type': 'EntrySet'}

electrical_conductivity_bpp = \
    {'label': 'Bipolar Plate Electrical Conductivity:', 'value': [6e4, 6e4],
     'sim_name': [['anode', 'electrical_conductivity_bpp'],
                  ['cathode', 'electrical_conductivity_bpp']],
     'dtype': 'float', 'dimensions': 'S/m', 'type': 'EntrySet'}
electrical_conductivity_gde = \
    {'label': 'Gas Diffusion Electrode Electrical Conductivity:',
     'value': [500.0, 500.0],
     'sim_name': [['anode', 'electrical_conductivity_gde'],
                  ['cathode', 'electrical_conductivity_gde']],
     'dtype': 'float', 'dimensions': 'S/m', 'type': 'EntrySet'}

porosity_gdl = \
    {'label': 'Gas Diffusion Layer Porosity:', 'number': 2, 'value': 0.8,
     'sim_name': [['anode', 'porosity_gdl'], ['cathode', 'porosity_gdl']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}
porosity_cl = \
    {'label': 'Catalyst Layer Porosity:', 'number': 2, 'value': 0.5,
     'sim_name': [['anode', 'porosity_cl'], ['cathode', 'porosity_cl']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}

