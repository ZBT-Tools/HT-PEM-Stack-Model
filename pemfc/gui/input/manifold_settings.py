# Manifold Settings
calc_distribution = \
    {'label': 'Activate Distribution Calculation:', 'number': 3,
     'sim_name': [['anode', 'flow_circuit', 'calc_distribution'],
                  ['cathode', 'flow_circuit', 'calc_distribution'],
                  ['coolant_flow_circuit', 'calc_distribution']],
     'sticky': ['NW', 'NWE'], 'type': 'CheckButtonSet',
     'command': {'function': 'set_status',
                 'args': [[[3, 1], [4, 1], [5, 1], [6, 1], [7, 1], [8, 1],
                           [9, 1], [10, 1], [11, 1], [12, 1], [13, 1]],
                          [[3, 2], [4, 2], [5, 2], [6, 2], [7, 2], [8, 2],
                           [9, 2], [10, 2], [11, 2], [12, 2], [13, 2]],
                          [[3, 3], [4, 3], [5, 3], [6, 3], [7, 3], [8, 3],
                           [9, 3], [10, 3], [11, 3], [12, 3], [13, 3]]]}}

anode_label_manifold = {'label': 'Anode', 'row': 1, 'column': 1,
                        'type': 'Label', 'sticky': 'WENS'}
cathode_label_manifold = \
    {'label': 'Cathode', 'row': 1, 'column': 2,
     'type': 'Label', 'sticky': 'WENS'}
cooling_label_manifold = \
    {'label': 'Cooling', 'row': 1, 'column': 3,
     'type': 'Label', 'sticky': 'WENS'}

manifold_configuration = \
    {'label': 'Manifold Flow Configuration:', 'number': 3,
     'sim_name': [['anode', 'flow_circuit', 'shape'],
                  ['cathode', 'flow_circuit', 'shape'],
                  ['coolant_flow_circuit', 'shape']],
     'value': ['U', 'Z'], 'type': 'ComboboxSet'}


inlet_manifold_cross_section = \
    {'label': 'Inlet Manifold Cross-Section:', 'number': 3,
     'sim_name': [['anode', 'flow_circuit', 'inlet_manifold',
                   'cross_sectional_shape'],
                  ['cathode', 'flow_circuit', 'inlet_manifold',
                   'cross_sectional_shape'],
                  ['coolant_flow_circuit', 'inlet_manifold',
                   'cross_sectional_shape']],
     'value': ['circular', 'rectangular'], 'type': 'ComboboxSet',
     'command': {'function': 'show_connected_widgets',
                 'args': [[[[[6, 1]], [[8, 1], [10, 1]]],
                           [[[8, 1], [10, 1]], [[6, 1]]]],
                          [[[[6, 2], [7, 2], [8, 2]], [[9, 2]]],
                           [[[6, 2], [7, 2], [8, 2], [9, 2]], []]],
                          [[[[6, 3], [7, 3], [8, 3]], [[9, 3]]],
                           [[[6, 3], [7, 3], [8, 3], [9, 3]], []]]]
                          }
     }

outlet_manifold_cross_section = \
    {'label': 'Outlet Manifold Cross-Section:', 'number': 3,
     'sim_name': [['anode', 'flow_circuit', 'outlet_manifold',
                   'cross_sectional_shape'],
                  ['cathode', 'flow_circuit', 'outlet_manifold',
                   'cross_sectional_shape'],
                  ['coolant_flow_circuit', 'outlet_manifold',
                   'cross_sectional_shape']],
     'value': ['circular', 'rectangular'], 'type': 'ComboboxSet'}

inlet_manifold_diameter = \
    {'label': 'Inlet Manifold Diameter:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow_circuit', 'inlet_manifold', 'diameter'],
                  ['cathode', 'flow_circuit', 'inlet_manifold', 'diameter'],
                  ['coolant_flow_circuit', 'inlet_manifold', 'diameter']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': 'm',
     'type': 'EntrySet'}

outlet_manifold_diameter = \
    {'label': 'Outlet Manifold Diameter:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow_circuit', 'outlet_manifold', 'diameter'],
                  ['cathode', 'flow_circuit', 'outlet_manifold', 'diameter'],
                  ['coolant_flow_circuit', 'outlet_manifold', 'diameter']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': 'm',
     'type': 'EntrySet'}

inlet_manifold_height = \
    {'label': 'Inlet Manifold Height:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow_circuit', 'inlet_manifold', 'height'],
                  ['cathode', 'flow_circuit', 'inlet_manifold', 'height'],
                  ['coolant_flow_circuit', 'inlet_manifold', 'height']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': 'm',
     'type': 'EntrySet'}

outlet_manifold_height = \
    {'label': 'Outlet Manifold Height:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow_circuit', 'outlet_manifold', 'height'],
                  ['cathode', 'flow_circuit', 'outlet_manifold', 'height'],
                  ['coolant_flow_circuit', 'outlet_manifold', 'height']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': 'm',
     'type': 'EntrySet'}

inlet_manifold_width = \
    {'label': 'Inlet Manifold Width:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow_circuit', 'inlet_manifold', 'width'],
                  ['cathode', 'flow_circuit', 'inlet_manifold', 'width'],
                  ['coolant_flow_circuit', 'inlet_manifold', 'width']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': 'm',
     'type': 'EntrySet'}

outlet_manifold_width = \
    {'label': 'Outlet Manifold Width:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow_circuit', 'outlet_manifold', 'width'],
                  ['cathode', 'flow_circuit', 'outlet_manifold', 'width'],
                  ['coolant_flow_circuit', 'outlet_manifold', 'width']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': 'm',
     'type': 'EntrySet'}

inlet_pressure_loss_coeff = \
    {'label': 'Inlet Manifold T-Junction Loss Coefficient:',
     'value': [0.4, 0.4, 0.4],
     'sim_name': [['anode', 'flow_circuit', 'inlet_manifold',
                   'constant_friction_factor'],
                  ['cathode', 'flow_circuit', 'inlet_manifold',
                   'constant_friction_factor'],
                  ['coolant_flow_circuit', 'inlet_manifold',
                   'constant_friction_factor']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': '-',
     'type': 'EntrySet'}

outlet_pressure_loss_coeff = \
    {'label': 'Outlet Manifold T-Junction Loss Coefficient:',
     'value': [0.4, 0.4, 0.4],
     'sim_name': [['anode', 'flow_circuit', 'outlet_manifold',
                   'constant_friction_factor'],
                  ['cathode', 'flow_circuit', 'outlet_manifold',
                   'constant_friction_factor'],
                  ['coolant_flow_circuit', 'outlet_manifold',
                   'constant_friction_factor']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': '-',
     'type': 'EntrySet'}

manifold_frame_dict = \
    {'title': 'Manifold Settings', 'show_title': False, 'font': 'Arial 10 bold',
     'sticky': 'WEN',
     'widget_dicts': [anode_label_manifold, cathode_label_manifold,
                      cooling_label_manifold, calc_distribution,
                      manifold_configuration,
                      inlet_manifold_cross_section,
                      outlet_manifold_cross_section,
                      inlet_manifold_diameter, outlet_manifold_diameter,
                      inlet_manifold_width, outlet_manifold_width,
                      inlet_manifold_height, outlet_manifold_height,
                      inlet_pressure_loss_coeff, outlet_pressure_loss_coeff],
     'highlightbackground': 'grey', 'highlightthickness': 1}

tab_dict = {'title': 'Manifolds', 'show_title': False,
            'sub_frame_dicts': [manifold_frame_dict]}
