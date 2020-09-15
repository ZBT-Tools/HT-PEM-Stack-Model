# Thermal conductivities
thermal_conductivity_bpp = \
    {'label': 'Bipolar Plate Thermal Conductivity:',
     'value': [[100.0, 100.0], [100.0, 100.0]],
     'sim_name': [['anode', 'thermal_conductivity_bpp', [0, 1]],
                  ['cathode', 'thermal_conductivity_bpp', [2, 3]]],
     'dtype': 'float', 'dimensions': 'W/(m K)', 'type': 'EntrySet'}

thermal_conductivity_gde = \
    {'label': 'Gas Diffusion Electrode Thermal Conductivity:',
     'value': [[10.0, 10.0], [10.0, 10.0]],
     'sim_name': [['anode', 'thermal_conductivity_gde', [0, 1]],
                  ['cathode', 'thermal_conductivity_gde', [2, 3]]],
     'dtype': 'float', 'dimensions': 'W/(m K)', 'type': 'EntrySet'}

electrical_conductivity_bpp = \
    {'label': 'Bipolar Plate Electrical Conductivity:',
     'number': 4, 'value': 6e4,
     'sim_name': [['anode', 'electrical_conductivity_bpp', [0, 1]],
                  ['cathode', 'electrical_conductivity_bpp', [2, 3]]],
     'dtype': 'float', 'dimensions': 'S/m', 'type': 'EntrySet'}

electrical_conductivity_gde = \
    {'label': 'Gas Diffusion Electrode Electrical Conductivity:',
     'number': 4, 'value': 500.0,
     'sim_name': [['anode', 'electrical_conductivity_gde', [0, 1]],
                  ['cathode', 'electrical_conductivity_gde', [2, 3]]],
     'dtype': 'float', 'dimensions': 'S/m', 'type': 'EntrySet'}

porosity_gdl = \
    {'label': 'Gas Diffusion Layer Porosity:', 'number': 2, 'value': 0.8,
     'sim_name': [['anode', 'porosity_gdl'], ['cathode', 'porosity_gdl']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}
porosity_cl = \
    {'label': 'Catalyst Layer Porosity:', 'number': 2, 'value': 0.5,
     'sim_name': [['anode', 'porosity_cl'], ['cathode', 'porosity_cl']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}

frame_dict = \
    {'title': 'Physical Properties', 'show_title': False, 'font': 'Arial 10 bold',
     'sticky': 'WEN',
     'widget_dicts': [thermal_conductivity_bpp,
                      thermal_conductivity_gde,
                      electrical_conductivity_bpp,
                      electrical_conductivity_gde,
                      porosity_gdl,
                      porosity_cl],
     'highlightbackground': 'grey', 'highlightthickness': 1}

tab_dict = {'title': 'Physical Properties', 'show_title': False,
            'sub_frame_dicts': [frame_dict]}
