load_settings_button_dict = \
    {'label': 'Load Settings', 'type': 'OpenDirectoryButton'}

output_dir_button_dict = \
    {'label': 'Open', 'type': 'OpenDirectoryButton'}
output_dir = \
    {'label': 'Output Directory:', 'button_dict': output_dir_button_dict,
     'sim_name': ['output', 'directory'], 'width': 20,
     'dtype': 'string', 'type': 'EntryButtonSet', 'sticky': 'W'}


output_frame_dict = \
    {'title': 'Output Settings', 'show_title': True, 'font': 'Arial 10 bold',
     'widget_dicts': [output_dir],
     'sticky': 'WEN',
     'highlightbackground': 'grey', 'highlightthickness': 1}

run_button_dict = {'label': 'Run Simulation', 'type': 'RunButton',
                   'columnspan': 3, 'width': 20, 'sticky': 'WNE'}

empty_row = {'label': ' ',  'font': 'Arial 1',  'height': 100,
             'type': 'Label', 'sticky': 'WENS'}
run_button_frame_dict = \
    {'title': 'Run Simulation', 'show_title': False,
     'widget_dicts': [run_button_dict, empty_row], 'sticky': ''}

tab_dict = \
    {'title': 'Simulation', 'show_title': False,
     'sub_frame_dicts': [output_frame_dict, run_button_frame_dict]}
