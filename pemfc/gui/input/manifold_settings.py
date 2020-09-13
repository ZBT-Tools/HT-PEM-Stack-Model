# Manifold Settings
calc_distribution = \
    {'label': 'Activate Distribution Calculation:', 'number': 3,
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
    {'label': 'Cooling', 'row': 1, 'column': 3,
     'type': 'Label', 'sticky': 'WENS'}

manifold_configuration = \
    {'label': 'Manifold Flow Configuration:', 'number': 3,
     'sim_name': [['anode', 'flow circuit', 'shape'],
                  ['cathode', 'flow circuit', 'shape'],
                  ['coolant flow circuit', 'shape']],
     'options': ['U', 'Z'], 'type': 'ComboboxSet'}

inlet_manifold_cross_section = \
    {'label': 'Inlet Manifold Cross-Section:', 'number': 3,
     'sim_name': [['anode', 'flow circuit', 'inlet manifold',
                   'cross_sectional_shape'],
                  ['cathode', 'flow circuit', 'inlet manifold',
                   'cross_sectional_shape'],
                  ['coolant flow circuit', 'inlet manifold',
                   'cross_sectional_shape']],
     'options': ['circular', 'rectangular'], 'type': 'ComboboxSet'}

outlet_manifold_cross_section = \
    {'label': 'Outlet Manifold Cross-Section:', 'number': 3,
     'sim_name': [['anode', 'flow circuit', 'outlet manifold',
                   'cross_sectional_shape'],
                  ['cathode', 'flow circuit', 'outlet manifold',
                   'cross_sectional_shape'],
                  ['coolant flow circuit', 'outlet manifold',
                   'cross_sectional_shape']],
     'options': ['circular', 'rectangular'], 'type': 'ComboboxSet'}

inlet_manifold_diameter = \
    {'label': 'Inlet Manifold Diameter:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow circuit', 'inlet manifold', 'diameter'],
                  ['cathode', 'flow circuit', 'inlet manifold', 'diameter'],
                  ['coolant flow circuit', 'inlet manifold', 'diameter']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}

outlet_manifold_diameter = \
    {'label': 'Outlet Manifold Diameter:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow circuit', 'outlet manifold', 'diameter'],
                  ['cathode', 'flow circuit', 'outlet manifold', 'diameter'],
                  ['coolant flow circuit', 'outlet manifold', 'diameter']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}

inlet_manifold_height = \
    {'label': 'Inlet Manifold Height:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow circuit', 'inlet manifold', 'height'],
                  ['cathode', 'flow circuit', 'inlet manifold', 'height'],
                  ['coolant flow circuit', 'inlet manifold', 'height']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}

outlet_manifold_height = \
    {'label': 'Outlet Manifold Height:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow circuit', 'outlet manifold', 'height'],
                  ['cathode', 'flow circuit', 'outlet manifold', 'height'],
                  ['coolant flow circuit', 'outlet manifold', 'height']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}

inlet_manifold_width = \
    {'label': 'Inlet Manifold Width:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow circuit', 'inlet manifold', 'width'],
                  ['cathode', 'flow circuit', 'inlet manifold', 'width'],
                  ['coolant flow circuit', 'inlet manifold', 'width']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}

outlet_manifold_width = \
    {'label': 'Outlet Manifold Width:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow circuit', 'outlet manifold', 'width'],
                  ['cathode', 'flow circuit', 'outlet manifold', 'width'],
                  ['coolant flow circuit', 'outlet manifold', 'width']],
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}

inlet_pressure_loss_coeff = \
    {'label': 'Inlet Manifold T-Junction Loss Coefficient:',
     'value': [0.4, 0.4, 0.4],
     'sim_name': [['anode', 'flow circuit', 'inlet manifold',
                   'constant_friction_factor'],
                  ['cathode', 'flow circuit', 'inlet manifold',
                   'constant_friction_factor'],
                  ['coolant flow circuit', 'inlet manifold',
                   'constant_friction_factor']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}

outlet_pressure_loss_coeff = \
    {'label': 'Outlet Manifold T-Junction Loss Coefficient:',
     'value': [0.4, 0.4, 0.4],
     'sim_name': [['anode', 'flow circuit', 'outlet manifold',
                   'constant_friction_factor'],
                  ['cathode', 'flow circuit', 'outlet manifold',
                   'constant_friction_factor'],
                  ['coolant flow circuit', 'outlet manifold',
                   'constant_friction_factor']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}

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
