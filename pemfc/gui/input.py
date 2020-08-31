# entry set dictionaries for geometry frame
cell_number = {'label': 'Cell Number:', 'number': 1,
               # 'grid_location': (1, 0),
               'type': 'EntrySet'}
cell_length = {'label': 'Cell Length:', 'dimensions': 'm',
               # 'grid_location': (2, 0),
               'type': 'EntrySet'}
cell_width = {'label': 'Cell Width:', 'number': 1, 'dimensions': 'm',
              # 'grid_location': (3, 0),
              'type': 'EntrySet'}
channel_length = {'label': 'Channel Length:', 'number': 2,
                  'dimensions': 'm', 'type': 'EntrySet'}
channel_width = {'label': 'Channel Width:', 'number': 2,
                 'dimensions': 'm', 'type': 'EntrySet'}
channel_height = {'label': 'Channel Height:', 'number': 2,
                  'dimensions': 'm', 'type': 'EntrySet'}
channel_bends = {'label': 'Number of Channel Bends:', 'number': 2,
                 'dimensions': '-', 'type': 'EntrySet'}
bend_pressure_loss_coefficient = \
    {'label': 'Pressure Loss Coefficient for Channel Bend:', 'dimensions': '-',
     'type': 'EntrySet'}
channel_flow_direction = {'label': 'Channel Flow Direction (1 or -1):',
                          'number': 2, 'dimensions': '-', 'type': 'EntrySet'}
bpp_thickness = {'label': 'Bipolar Plate Thickness:', 'number': 2,
                 'dimensions': 'm', 'type': 'EntrySet'}
gdl_thickness = {'label': 'Gas Diffusion Layer Thickness:', 'number': 2,
                 'dimensions': 'm', 'type': 'EntrySet'}
cl_thickness = {'label': 'Catalyst Layer Thickness:', 'number': 2,
                'dimensions': 'm', 'type': 'EntrySet'}
gdl_porosity = {'label': 'Gas Diffusion Layer Porosity:', 'number': 2,
                'dimensions': '-', 'type': 'EntrySet'}
cl_porosity = {'label': 'Catalyst Layer Porosity:', 'number': 2,
               'dimensions': '-', 'type': 'EntrySet'}
mem_thickness = {'label': 'Bipolar Plate Thickness:', 'dimensions': 'm',
                 'type': 'EntrySet'}
cool_circuit = {'label': 'Activate Cooling', 'type': 'CheckButtonSet'}
cool_channel_length = {'label': 'Coolant Channel Length:', 'dimensions': 'm',
                       'type': 'EntrySet'}
cool_channel_height = {'label': 'Coolant Channel Height:', 'dimensions': 'm',
                       'type': 'EntrySet'}
cool_channel_width = {'label': 'Coolant Channel Height:', 'dimensions': 'm',
                      'type': 'EntrySet'}
cool_channel_number = {'label': 'Coolant Channel number:', 'type': 'EntrySet'}
cool_channel_bends = {'label': 'Number of Channel Bends:', 'type': 'EntrySet'}
cool_bend_pressure_loss_coefficient = \
    {'label': 'Pressure Loss Coefficient for Coolant Channel Bend:',
     'dimensions': '-', 'type': 'EntrySet'}
cool_flow_end_cells = \
    {'label': 'Activate Cooling Flow at End Plates', 'type': 'CheckButtonSet'}

anode_label_cell = {'label': 'Anode', 'row': 1, 'column': 1,
                    'type': 'Label', 'sticky': 'WENS'}
cathode_label_cell = \
    {'label': 'Cathode', 'row': 1, 'column': 2,
     'type': 'Label', 'sticky': 'WENS'}

"""Manifold Settings"""
calc_distribution = \
    {'label': 'Activate Distribution Calculation', 'number': 3,
     'type': 'CheckButtonSet'}

anode_label_manifold = {'label': 'Anode', 'row': 1, 'column': 1,
                        'type': 'Label', 'sticky': 'WENS'}
cathode_label_manifold = \
    {'label': 'Cathode', 'row': 1, 'column': 2,
     'type': 'Label', 'sticky': 'WENS'}
cooling_label_manifold = \
    {'label': 'Cathode', 'row': 1, 'column': 3,
     'type': 'Label', 'sticky': 'WENS'}

channel_frame_dict = \
    {'title': 'Channel Settings', 'grid_location': (1, 1),
     'widget_set_dicts': [anode_label_cell,
                          cathode_label_cell,
                          channel_length, channel_width,
                          channel_height,
                          channel_bends, bend_pressure_loss_coefficient,
                          channel_flow_direction],
     'highlightbackground': 'grey', 'highlightthickness': 1}
     # 'sticky': 'WENS'}
cool_frame_dict = \
    {'title': 'Cooling Settings', 'font': 'Arial 10 bold',
     'widget_set_dicts': [cool_circuit, cool_channel_length],
     'highlightbackground': 'grey', 'highlightthickness': 1, 'sticky': 'WENS'}
cell_frame_sub_dict = \
    {'title': 'Cell Lengths', 'grid_location': (1, 0),
     'widget_set_dicts': [cell_number, cell_length, cell_width],
     'highlightbackground': 'grey', 'highlightthickness': 1, 'sticky': 'WEN'}

cell_frame_dict = \
    {'title': 'Cell Settings', 'show_title': True,
     'sub_frame_dicts': [cell_frame_sub_dict, channel_frame_dict],
     'highlightbackground': 'grey', 'highlightthickness': 1}

manifold_frame_dict = \
    {'title': 'Manifold Settings', 'show_title': True,
     'widget_set_dicts': [calc_distribution],
     'highlightbackground': 'grey', 'highlightthickness': 1}

geometry_frame_dict = \
    {'title': 'Geometry', 'show_title': False,
     'sub_frame_dicts': [cell_frame_dict, cool_frame_dict, manifold_frame_dict],
     'highlightbackground': 'grey', 'highlightthickness': 1}

main_frame_dicts = [geometry_frame_dict]

