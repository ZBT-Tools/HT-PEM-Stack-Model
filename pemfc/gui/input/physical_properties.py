# Thermal conductivities
anode_label = {'label': 'Anode', 'row': 1, 'column': 1, 'columnspan': 2, 'pady': 0,
               'type': 'Label', 'sticky': 'WENS'}
cathode_label = {'label': 'Cathode', 'row': 1, 'column': 3, 'columnspan': 2, 'pady': 0,
                 'type': 'Label', 'sticky': 'WENS'}
in_plane_label_1 = {'label': 'ip', 'row': 2, 'column': 1, 'pady': 0,
                  'type': 'Label', 'sticky': 'WENS'}
through_plane_label_1 = {'label': 'tp', 'row': 2, 'column': 2, 'pady': 0,
                         'type': 'Label', 'sticky': 'WENS'}
in_plane_label_2 = {'label': 'ip', 'row': 2, 'column': 3, 'pady': 0,
                    'type': 'Label', 'sticky': 'WENS'}
through_plane_label_2 = {'label': 'tp', 'row': 2, 'column': 4, 'pady': 0,
                         'type': 'Label', 'sticky': 'WENS'}
electrical_conductivity_bpp = \
    {'label': 'Bipolar Plate Electrical Conductivity:',
     'number': 4, 'value': 6e4, 'width': 5,
     'sim_name': [['anode', 'electrical_conductivity_bpp', [0, 1]],
                  ['cathode', 'electrical_conductivity_bpp', [2, 3]]],
     'dtype': 'float', 'dimensions': 'S/m', 'type': 'EntrySet'}

electrical_conductivity_gde = \
    {'label': 'Gas Diffusion Electrode Electrical Conductivity:',
     'number': 4, 'value': 500.0, 'width': 5,
     'sim_name': [['anode', 'electrical_conductivity_gde', [0, 1]],
                  ['cathode', 'electrical_conductivity_gde', [2, 3]]],
     'dtype': 'float', 'dimensions': 'S/m', 'type': 'EntrySet'}

thermal_conductivity_bpp = \
    {'label': 'Bipolar Plate Thermal Conductivity:',
     'value': [[100.0, 100.0], [100.0, 100.0]], 'width': 5,
     'sim_name': [['anode', 'thermal_conductivity_bpp', [0, 1]],
                  ['cathode', 'thermal_conductivity_bpp', [2, 3]]],
     'dtype': 'float', 'dimensions': 'W/(m K)', 'type': 'EntrySet'}

thermal_conductivity_gde = \
    {'label': 'Gas Diffusion Electrode Thermal Conductivity:',
     'value': [[1.0, 1.0], [1.0, 1.0]], 'width': 5,
     'sim_name': [['anode', 'thermal_conductivity_gde', [0, 1]],
                  ['cathode', 'thermal_conductivity_gde', [2, 3]]],
     'dtype': 'float', 'dimensions': 'W/(m K)', 'type': 'EntrySet'}

empty_row = {'label': ' ',  'font': 'Arial 1',  # 'row': 1, 'column': 1,
             'type': 'Label', 'sticky': 'WENS'}

thermal_conductivity_mem = \
    {'label': 'Membrane Thermal Conductivity:',
     'value': [[0.26, 0.26]], 'width': 5,
     'sim_name': [['membrane', 'thermal_conductivity', [0, 1]]],
     'dtype': 'float', 'dimensions': 'W/(m K)', 'type': 'EntrySet'}

# empty_row = {'label': ' ',  'font': 'Arial 1',  # 'row': 1, 'column': 1,
#              'type': 'Label', 'sticky': 'WENS'}

porosity_gdl = \
    {'label': 'Gas Diffusion Layer Porosity:', 'number': 2, 'value': 0.8,
     'width': 5, #'column': 1, 'columnspan': 2,
     'sim_name': [['anode', 'porosity_gdl'], ['cathode', 'porosity_gdl']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}
porosity_cl = \
    {'label': 'Catalyst Layer Porosity:', 'number': 2, 'value': 0.5,
     'width': 5, #'column': 1,
     'sim_name': [['anode', 'porosity_cl'], ['cathode', 'porosity_cl']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}

porous_frame_dict = \
    {'title': 'Porous Layers', 'show_title': True,
     'font': 'Arial 10 bold', 'sticky': 'WEN',
     'widget_dicts': [anode_label,
                      cathode_label,
                      through_plane_label_1, in_plane_label_1,
                      through_plane_label_2, in_plane_label_2,
                      thermal_conductivity_bpp,
                      thermal_conductivity_gde,
                      electrical_conductivity_bpp,
                      electrical_conductivity_gde,
                      empty_row,
                      porosity_gdl, porosity_cl],
     'highlightbackground': 'grey', 'highlightthickness': 1}

membrane_model = \
    {'label': 'Membrane Model:', 'number': 1,
     'sim_name': ['membrane', 'type'],
     'value': ['Constant', 'Springer', 'Kvesic'],
     'type': 'ComboboxSet',
     'command': {'function': 'show_connected_widgets',
                 'args': [[[[[3, 0]], [[4, 0], [5, 0]]],
                           [[[4, 0]], [[3, 0], [5, 0]]],
                           [[[5, 0]], [[3, 0], [4, 0]]]]]}
     }

mem_loss = \
    {'label': 'Calculate Membrane Loss:', 'value': True,
     'sim_name': ['membrane', 'calc_loss'],
     'dtype': 'bool', 'type': 'CheckButtonSet'}

mem_ionic_conductivity = \
    {'label': 'Ionic Conductivity:', 'value': 5.0,
     'sim_name': ['membrane', 'ionic_conductivity'],
     'dtype': 'float', 'dimensions': 'S/m', 'type': 'EntrySet'}

constant_frame = \
    {'title': 'Constant Ionic Conductivity',
     'widget_dicts': [mem_ionic_conductivity],
     'sticky': 'WEN', 'columnspan': 2}

mem_constant_resistance = \
    {'label': 'Constant Resistance Coefficient:', 'value': 4.3e-5,
     'sim_name': ['membrane', 'basic_resistance'],
     'dtype': 'float', 'dimensions': 'Ohm-m²', 'type': 'EntrySet'}

mem_temp_coefficient = \
    {'label': 'Linear Temperature Coefficient:', 'value': 7e-8,
     'sim_name': ['membrane', 'temperature_coefficient'],
     'dtype': 'float', 'dimensions': 'Ohm-m²/K', 'type': 'EntrySet'}
kvesic_frame = \
    {'title': 'Kvesic Ionic Conductivity',
     'widget_dicts': [mem_constant_resistance, mem_temp_coefficient],
     'sticky': 'WEN', 'columnspan': 2}

mem_acid_group_conc = \
    {'label': 'Acid Group Concentration:', 'value': 1.2e3,
     'sim_name': ['membrane', 'acid_group_concentration'],
     'dtype': 'float', 'dimensions': 'mol/m³', 'type': 'EntrySet'}
mem_vapour_transport_coefficient = \
    {'label': 'Vapour Transport Coefficient:', 'value': 6.2e-6,
     'sim_name': ['membrane', 'vapour_transport_coefficient'],
     'dtype': 'float', 'dimensions': 'm/s', 'type': 'EntrySet'}
springer_frame = \
    {'title': 'Springer Ionic Conductivity',
     'widget_dicts': [mem_acid_group_conc,
                      mem_vapour_transport_coefficient],
     'sticky': 'WEN', 'columnspan': 2}

membrane_frame_dict = \
    {'title': 'Membrane Settings', 'show_title': True,
     'font': 'Arial 10 bold', 'sticky': 'WEN',
     'widget_dicts': [mem_loss,
                      membrane_model,
                      constant_frame,
                      springer_frame,
                      kvesic_frame],
     'highlightbackground': 'grey', 'highlightthickness': 1}


frame_dict = \
    {'title': 'Physical Properties', 'show_title': False,
     'font': 'Arial 10 bold', 'sticky': 'WEN',
     'sub_frame_dicts': [porous_frame_dict,
                         membrane_frame_dict],
     'highlightbackground': 'grey', 'highlightthickness': 1}

tab_dict = {'title': 'Physical Properties', 'show_title': False,
            'sub_frame_dicts': [frame_dict]}
