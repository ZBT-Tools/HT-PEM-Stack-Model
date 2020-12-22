# Electrochemistry Settings
# current_control = \
#     {'label': 'Current Control:', 'value': True, 'sticky': ['NW', 'NWE'],
#      'sim_name': ['simulation', 'current_control'], 'type': 'CheckButtonSet',
#      # 'command': {'function': 'set_visibility',
#      #             'args': [[[1, 0], [1, 1], [1, 2]],
#      #                      [[2, 0], [2, 1], [2, 2]]]}
#      }

current_control = \
    {'label': 'Operation Control:', 'number': 1,
     'sim_name': ['simulation', 'operation_control'],
     'value': ['Current', 'Voltage'],
     'type': 'ComboboxSet',
     'command': {'function': 'show_connected_widgets',
                 'args': [[[[[2, 0], [2, 1], [2, 2]],
                            [[3, 0], [3, 1], [3, 2]]],
                           [[[3, 0], [3, 1], [3, 2]],
                            [[2, 0], [2, 1], [2, 2]]]]]}
     }

current_density = \
    {'label': 'Current Density:', 'value': 0.0,
     'sim_name': ['simulation', 'current_density'],
     'dtype': 'float', 'dimensions': 'A/m²', 'type': 'EntrySet'}

average_cell_voltage = \
    {'label': 'Cell Voltage:', 'value': 0.5,
     'sim_name': ['simulation', 'average_cell_voltage'],
     'dtype': 'float', 'dimensions': 'V', 'type': 'EntrySet'}

open_circuit_voltage = \
    {'label': 'Open Circuit Voltage:', 'value': 1.0,
     'sim_name': ['cell', 'open_circuit_voltage'],
     'dtype': 'float', 'dimensions': 'V', 'type': 'EntrySet'}

electrochem_frame_dict = \
    {'title': 'Electrochemistry Settings', 'show_title': True, 'sticky': 'NWE',
     'widget_dicts': [current_control,
                      current_density,
                      average_cell_voltage,
                      open_circuit_voltage]}

tab_dict = \
    {'title': 'Operating Conditions', 'show_title': False,
     'sub_frame_dicts': [electrochem_frame_dict]}
