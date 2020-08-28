# entry set dictionaries for geometry frame
cell_number = {'label': 'Cell Number:', 'type': 'EntrySet'}
cell_length = {'label': 'Cell Length:', 'dimensions': 'm', 'type': 'EntrySet'}
cell_width = {'label': 'Cell Width:', 'dimensions': 'm', 'type': 'EntrySet'}
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
cool_circuit = {'label': 'Activate Cooling Circuit', 'type': 'CheckButtonSet'}
cool_channel_length = {'label': 'Coolant Channel Length:', 'dimensions': 'm',
                       'type': 'EntrySet'}
anode_label = {'label': 'Anode', 'row': 1, 'column': 1,
               'type': 'Label', 'sticky': 'WENS'}
cathode_label = \
    {'label': 'Cathode', 'row': 1, 'column': 2,
     'type': 'Label', 'sticky': 'WENS'}

channel_frame_dict = \
    {'title': 'Channel Settings',
     'widget_set_dicts': [anode_label,
                          cathode_label,
                          channel_length, channel_width,
                          channel_height,
                          channel_bends, bend_pressure_loss_coefficient,
                          channel_flow_direction],
     # 'highlightbackground': 'grey', 'highlightthickness': 1,
     'sticky': 'WENS'}
cool_frame_dict = \
    {'title': 'Cooling Settings', 'font': 'Arial 10 bold',
     'widget_set_dicts': [cool_circuit, cool_channel_length],
     'highlightbackground': 'grey', 'highlightthickness': 1} #, 'sticky': 'WENS'}

cell_frame_dict = {'title': 'Cell Settings',
                   'sub_frame_dicts': [channel_frame_dict],
                   'widget_set_dicts': [cell_number, cell_length, cell_width],
                   'highlightbackground': 'grey', 'highlightthickness': 1}
geometry_frame_dict = {'title': 'Geometry', 'show_title': False,
                       'sub_frame_dicts': [cell_frame_dict, cool_frame_dict]}
main_frame_dicts = [geometry_frame_dict]

