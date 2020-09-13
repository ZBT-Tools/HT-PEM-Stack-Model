"""Geometry"""
# Cell Settings
#

channel_length = \
    {'label': 'Channel Length:', 'number': 2, 'value': [0.4, 0.4],
     'sim_name': [['anode', 'channel', 'length'],
                  ['cathode', 'channel', 'length']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}
channel_width = \
    {'label': 'Channel Width:', 'number': 2, 'value': [1e-3, 1e-3],
     'sim_name': [['anode', 'channel', 'width'],
                  ['cathode', 'channel', 'width']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}
channel_height = \
    {'label': 'Channel Height:', 'number': 2, 'value': [1e-3, 1e-3],
     'sim_name': [['anode', 'channel', 'height'],
                  ['cathode', 'channel', 'height']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}
channel_bends = \
    {'label': 'Number of Channel Bends:', 'value': [48, 48],
     'sim_name': [['anode', 'channel', 'bend_number'],
                  ['cathode', 'channel', 'bend_number']],
     'dtype': 'int', 'dimensions': '-', 'type': 'EntrySet'}
bend_pressure_loss_coefficient = \
    {'label': 'Pressure Loss Coefficient for Channel Bend:', 'value': 0.5,
     'sim_name': [['anode', 'channel', 'bend_friction_factor'],
                  ['cathode', 'channel', 'bend_friction_factor']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}
channel_flow_direction = \
    {'label': 'Channel Flow Direction (1 or -1):', 'value': [1, -1],
     'sim_name': [['anode', 'channel', 'flow_direction'],
                  ['cathode', 'channel', 'flow_direction']],
     'dtype': 'int', 'dimensions': '-', 'type': 'EntrySet'}
thickness_bpp = \
    {'label': 'Bipolar Plate Thickness:', 'number': 2, 'value': 2e-3,
     'sim_name': [['anode', 'thickness_bpp'], ['cathode', 'thickness_bpp']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}
thickness_gdl = \
    {'label': 'Gas Diffusion Layer Thickness:', 'number': 2, 'value': 2e-4,
     'sim_name': [['anode', 'thickness_gdl'], ['cathode', 'thickness_gdl']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}
thickness_cl = \
    {'label': 'Catalyst Layer Thickness:', 'number': 2, 'value': 1e-5,
     'sim_name': [['anode', 'thickness_cl'], ['cathode', 'thickness_cl']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}
porosity_gdl = \
    {'label': 'Gas Diffusion Layer Porosity:', 'number': 2, 'value': 0.8,
     'sim_name': [['anode', 'porosity_gdl'], ['cathode', 'porosity_gdl']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}
porosity_cl = \
    {'label': 'Catalyst Layer Porosity:', 'number': 2, 'value': 0.5,
     'sim_name': [['anode', 'porosity_cl'], ['cathode', 'porosity_cl']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}
thickness_mem = \
    {'label': 'Membrane Thickness:',  'value': 1.5e-5,
     'sim_name': ['membrane', 'thickness'],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}
cool_circuit = \
    {'label': 'Activate Cooling', 'value': True,
     'sim_name': ['stack', 'cool_flow'], 'type': 'CheckButtonSet'}
cool_channel_length = \
    {'label': 'Coolant Channel Length:', 'value': 0.4,
     'sim_name': ['coolant channel', 'length'],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}
cool_channel_height = \
    {'label': 'Coolant Channel Height:', 'value': 1e-3,
     'sim_name': ['coolant channel', 'height'],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}
cool_channel_width = \
    {'label': 'Coolant Channel Height:', 'value': 1e-3,
     'sim_name': ['coolant channel', 'width'],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}
cool_channel_number = {'label': 'Coolant Channel number:', 'type': 'EntrySet'}
cool_channel_bends = {'label': 'Number of Channel Bends:', 'type': 'EntrySet'}
cool_bend_pressure_loss_coefficient = \
    {'label': 'Pressure Loss Coefficient for Coolant Channel Bend:',
     'dimensions': '-', 'type': 'EntrySet'}
cool_flow_end_cells = \
    {'label': 'Activate Cooling Flow at End Plates', 'type': 'CheckButtonSet'}

anode_label_channel = {'label': 'Anode', 'row': 0, 'column': 1,
                       'type': 'Label', 'sticky': 'WENS'}
cathode_label_channel = \
    {'label': 'Cathode', 'row': 0, 'column': 2,
     'type': 'Label', 'sticky': 'WENS'}
empty_row = {'label': ' ', 'font': 'Arial 1',  # 'row': 1, 'column': 1,
             'type': 'Label', 'sticky': 'WENS'}

"""Manifold Settings"""
calc_distribution = \
    {'label': 'Activate Distribution Calculation', 'number': 3,
     'sim_name': [['anode', 'flow circuit', 'calc_distribution'],
                  ['cathode', 'flow circuit', 'calc_distribution'],
                  ['coolant flow circuit', 'calc_distribution']],
     'type': 'CheckButtonSet'}

anode_label_manifold = {'label': 'Anode', 'row': 1, 'column': 1,
                        'type': 'Label', 'sticky': 'WENS'}
cathode_label_manifold = \
    {'label': 'Cathode', 'row': 1, 'column': 2,
     'type': 'Label', 'sticky': 'WENS'}
cooling_label_manifold = \
    {'label': 'Cathode', 'row': 1, 'column': 3,
     'type': 'Label', 'sticky': 'WENS'}

manifold_configuration = \
    {'label': 'Manifold Flow Configuration', 'number': 3,
     'sim_name': [['anode', 'flow circuit', 'shape'],
                  ['cathode', 'flow circuit', 'shape'],
                  ['coolant flow circuit', 'shape']],
     'options': ['U', 'Z'], 'type': 'ComboboxSet'}

channel_frame_dict = \
    {'title': 'Channel Settings', #'grid_location': (1, 1),
     'font': 'Arial 10 bold',
     'widget_set_dicts': [anode_label_channel,
                          cathode_label_channel,
                          channel_length, channel_width,
                          channel_height,
                          channel_bends,
                          channel_flow_direction,
                          empty_row,
                          bend_pressure_loss_coefficient],
     'highlightbackground': 'grey', 'highlightthickness': 1,
     'sticky': 'WENS'}

cool_frame_dict = \
    {'title': 'Cooling Settings', 'font': 'Arial 10 bold',
     'widget_set_dicts': [cool_circuit, cool_channel_length],
     'highlightbackground': 'grey', 'highlightthickness': 1, 'sticky': 'WENS'}

# Cell Settings (sub frame)
cell_number = {'label': 'Cell Number:', 'number': 1, 'value': 10,
               'sim_name': ['stack', 'cell_number'], 'dtype': 'int',
               # 'grid_location': (1, 0),
               'type': 'EntrySet'}
cell_length = {'label': 'Cell Length:', 'dimensions': 'm', 'value': 0.1,
               'sim_name': ['cell', 'length'], 'dtype': 'float',
               # 'grid_location': (2, 0),
               'type': 'EntrySet'}
cell_width = {'label': 'Cell Width:', 'number': 1, 'value': 0.1,
              'sim_name': ['cell', 'width'],
              # 'grid_location': (3, 0),
              'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}
anode_label_cell = {'label': 'Anode', 'row': 5, 'column': 1,
                    'type': 'Label', 'sticky': 'WENS'}
cathode_label_cell = \
    {'label': 'Cathode', 'row': 5, 'column': 2,
     'type': 'Label', 'sticky': 'WENS'}
empty_row = {'label': ' ',  'font': 'Arial 1',  # 'row': 1, 'column': 1,
             'type': 'Label', 'sticky': 'WENS'}
cell_frame_sub_dict = \
    {'title': 'Cell Lengths', 'show_title': False,
     'grid_location': (1, 0), 'font': 'Arial 10 bold',
     'widget_set_dicts': [cell_number, cell_length,
                          cell_width,
                          anode_label_cell,
                          cathode_label_cell,
                          thickness_bpp,
                          thickness_gdl,
                          thickness_cl,
                          empty_row,
                          thickness_mem],
     #'highlightbackground': 'grey', 'highlightthickness': 1,
     'sticky': 'WESN'}

cell_frame_dict = \
    {'title': 'Cell Settings', 'show_title': True, 'font': 'Arial 10 bold',
     'sub_frame_dicts': [cell_frame_sub_dict, channel_frame_dict],
     # 'widget_set_dicts': [cell_number, cell_length, cell_width],
     'highlightbackground': 'grey', 'highlightthickness': 1}

manifold_frame_dict = \
    {'title': 'Manifold Settings', 'show_title': True, 'font': 'Arial 10 bold',
     'widget_set_dicts': [calc_distribution, manifold_configuration],
     'highlightbackground': 'grey', 'highlightthickness': 1}


output_dir_button_dict = \
    {'label': 'Open', 'type': 'OpenDirectoryButton'}
output_dir = \
    {'label': 'Output Directory:', 'button_dict': output_dir_button_dict,
     'sim_name': ['output', 'directory'], 'width': 20,
     'dtype': 'string', 'type': 'EntryButtonSet', 'sticky': 'W'}

run_button_dict = {'label': 'Run Simulation', 'type': 'RunButton',
                   'columnspan': 3}

output_frame_dict = \
    {'title': 'Output Settings', 'show_title': True, 'font': 'Arial 10 bold',
     'widget_set_dicts': [output_dir, run_button_dict],
     # 'button_dicts': [run_button_dict],
     'sticky': 'WEN',
     'highlightbackground': 'grey', 'highlightthickness': 1}

geometry_frame_dict = \
    {'title': 'Geometry', 'show_title': False, 'font': 'Arial 10 bold',
     'sub_frame_dicts': [cell_frame_dict, cool_frame_dict, manifold_frame_dict],
     'highlightbackground': 'grey', 'highlightthickness': 1}

simulation_frame_dict = \
    {'title': 'Simulation', 'show_title': False, 'font': 'Arial 10 bold',
     'sub_frame_dicts': [output_frame_dict],
     'highlightbackground': 'grey', 'highlightthickness': 1}


main_frame_dicts = [geometry_frame_dict, simulation_frame_dict]







