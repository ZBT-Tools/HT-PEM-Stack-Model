# Cooling Settings
cool_circuit = \
    {'label': 'Activate Cooling:', 'value': True, 'sticky': ['NW', 'NWE'],
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
cool_channel_number = \
    {'label': 'Coolant Channel Number:', 'value': 2, 'type': 'EntrySet'}
cool_channel_bends = \
    {'label': 'Number of Coolant Channel Bends:', 'value': 0, 'type':
        'EntrySet'}
cool_bend_pressure_loss_coefficient = \
    {'label': 'Pressure Loss Coefficient for Coolant Channel Bend:',
     'value': 0.5, 'dimensions': '-', 'type': 'EntrySet'}
cool_flow_end_cells = \
    {'label': 'Activate Cooling Flow at End Plates:', 'value': False,
     'sticky': ['NW', 'NWE'], 'type': 'CheckButtonSet'}
cool_frame_dict = \
    {'title': 'Cooling Settings', 'show_title': False, 'font': 'Arial 10 bold',
     'sticky': 'WEN',
     'widget_dicts': [cool_circuit, cool_channel_number,
                      cool_channel_length, cool_channel_height,
                      cool_channel_width, cool_channel_bends,
                      cool_bend_pressure_loss_coefficient,
                      cool_flow_end_cells],
     'highlightbackground': 'grey', 'highlightthickness': 1}

tab_dict = {'title': 'Cooling', 'show_title': False,
            'sub_frame_dicts': [cool_frame_dict]}

# geometry_frame_dict = \
#     {'title': 'Geometry', 'show_title': False, 'font': 'Arial 10 bold',
#      'sub_frame_dicts': [cool_frame_dict, manifold_frame_dict,
#                          cell_frame_dict],
#      'highlightbackground': 'grey', 'highlightthickness': 1}

# main_frame_dicts = [geometry_frame_dict, simulation_frame_dict]







