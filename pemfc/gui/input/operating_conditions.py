# # Electrochemistry Settings
# current_control = \
#     {'label': 'Current Control:', 'value': True, 'sticky': ['NW', 'NWE'],
#      'sim_name': ['simulation', 'current_control'], 'type': 'CheckButtonSet',
#      'command': {'function': 'set_visibility',
#                  'args': [[[1, 0], [1, 1], [1, 2]],
#                           [[2, 0], [2, 1], [2, 2]]]}}
#
# current_density = \
#     {'label': 'Current Density:', 'value': 0.0,
#      'sim_name': ['simulation', 'current_density'],
#      'dtype': 'float', 'dimensions': 'A/m²', 'type': 'EntrySet'}
#
# electrochem_frame_dict = \
#     {'title': 'Electrochemistry Settings', 'show_title': True,
#      'widget_dicts': [current_control,
#                       current_density]}
#
#
# tab_dict = \
#     {'title': 'Operating Conditions', 'show_title': False,
#      'sub_frame_dicts': [electrochem_frame_dict]}
