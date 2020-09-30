# Cell Settings

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

thickness_mem = \
    {'label': 'Membrane Thickness:',  'value': 1.5e-5,
     'sim_name': ['membrane', 'thickness'],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}

empty_row = {'label': ' ',  'font': 'Arial 1',  # 'row': 1, 'column': 1,
             'type': 'Label', 'sticky': 'WENS'}

cell_frame_sub_dict = \
    {'title': 'Cell Lengths', 'show_title': False,
     'grid_location': (1, 0), 'font': 'Arial 10 bold',
     'widget_dicts': [cell_number, cell_length,
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

# Channel Settings
channel_number = \
    {'label': 'Channel Number:', 'number': 2, 'value': [10, 10],
     'sim_name': [['anode', 'channel', 'number'],
                  ['cathode', 'channel', 'number']],
     'dtype': 'int', 'dimensions': 'm', 'type': 'EntrySet'}

channel_shape = \
    {'label': 'Shape of Cross-Section:', 'number': 1,
     'sim_name': ['membrane', 'type'],
     'value': ['rectangular', 'trapezoidal', 'triangular'],
     'type': 'ComboboxSet', 'sticky': 'WN',
     'command': {'function': 'show_connected_widgets',
                 'args': [[[[[1, 0]], [[2, 0], [3, 0]]],
                           [[[2, 0]], [[1, 0], [3, 0]]],
                           [[[1, 0]], [[2, 0], [3, 0]]]]]}}


channel_width = \
    {'label': 'Width:', 'number': 2, 'value': [1e-3, 1e-3],
     'sim_name': [['anode', 'channel', 'width'],
                  ['cathode', 'channel', 'width']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}

channel_rib_width = \
    {'label': 'Width of Ribs:', 'number': 2, 'value': [1e-3, 1e-3],
     'sim_name': [['anode', 'channel', 'rib_width'],
                  ['cathode', 'channel', 'width']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}

channel_height = \
    {'label': 'Height:', 'number': 2, 'value': [1e-3, 1e-3],
     'sim_name': [['anode', 'channel', 'height'],
                  ['cathode', 'channel', 'height']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}

channel_base_width = \
    {'label': 'Base Width:', 'number': 2, 'value': [1e-3, 1e-3],
     'sim_name': [['anode', 'channel', 'base_width'],
                  ['cathode', 'channel', 'base_width']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}

rectangular_frame = \
    {'widget_dicts': [channel_width,
                      channel_rib_width,
                      channel_height],
     'sticky': 'WEN', 'columnspan': 4, 'padx': 0.0, 'pady': 0.0}

trapezoidal_frame = \
    {'widget_dicts': [channel_width,
                      channel_rib_width,
                      channel_height,
                      channel_base_width],
     'sticky': 'WEN', 'columnspan': 4, 'padx': 0.0, 'pady': 0.0}

channel_shape_frame = \
    {'widget_dicts': [channel_shape,
                      rectangular_frame,
                      trapezoidal_frame,
                      rectangular_frame],
     'sticky': 'WEN', 'columnspan': 4, 'padx': 0.0, 'pady': 0.0}

channel_length = \
    {'label': 'Channel Length:', 'number': 2, 'value': [0.1, 0.1],
     'sim_name': [['anode', 'channel', 'length'],
                  ['cathode', 'channel', 'length']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}

channel_bends = \
    {'label': 'Number of Channel Bends:', 'value': [48, 48],
     'sim_name': [['anode', 'channel', 'bend_number'],
                  ['cathode', 'channel', 'bend_number']],
     'dtype': 'int', 'type': 'EntrySet'}
bend_pressure_loss_coefficient = \
    {'label': 'Pressure Loss Coefficient for Channel Bend:', 'value': 0.5,
     'sim_name': [['anode', 'channel', 'bend_friction_factor'],
                  ['cathode', 'channel', 'bend_friction_factor']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}
channel_flow_direction = \
    {'label': 'Channel Flow Direction (1 or -1):', 'value': [-1, 1],
     'sim_name': [['anode', 'channel', 'flow_direction'],
                  ['cathode', 'channel', 'flow_direction']],
     'dtype': 'int', 'dimensions': '-', 'type': 'EntrySet'}
anode_label_channel = {'label': 'Anode', 'row': 1, 'column': 1,
                       'type': 'Label', 'sticky': 'WENS'}
cathode_label_channel = \
    {'label': 'Cathode', 'row': 1, 'column': 2,
     'type': 'Label', 'sticky': 'WENS'}
channel_first_frame = \
    {'widget_dicts': [anode_label_channel,
                      cathode_label_channel,
                      channel_number],
     'sticky': 'WEN', 'columnspan': 4, 'padx': 0.0, 'pady': 0.0}

channel_last_frame = \
    {'widget_dicts': [channel_length,
                      channel_bends,
                      bend_pressure_loss_coefficient,
                      channel_flow_direction],
     'sticky': 'WEN', 'columnspan': 4, 'padx': 0.0, 'pady': 0.0}
channel_frame_dict = \
    {'title': 'Channel Settings', #'grid_location': (1, 1),
     'font': 'Arial 10 bold',
     'widget_dicts': [channel_first_frame,
                      channel_shape_frame,
                      channel_last_frame],
     'highlightbackground': 'grey', 'highlightthickness': 1,
     'sticky': 'WENS'}

cell_frame_dict = \
    {'title': 'Cell Settings', 'show_title': False, 'font': 'Arial 10 bold',
     'sticky': 'WEN',
     'sub_frame_dicts': [cell_frame_sub_dict, channel_frame_dict],
     # 'widget_dicts': [cell_number, cell_length, cell_width],
     'highlightbackground': 'grey', 'highlightthickness': 1}

tab_dict = {'title': 'Cells', 'show_title': False,
            'sub_frame_dicts': [cell_frame_dict]}
