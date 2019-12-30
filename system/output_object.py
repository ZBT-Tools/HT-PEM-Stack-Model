import numpy as np


class OutputObject:
    def __init__(self):
        self.print_data_1d = {}
        self.print_data_2d = {}
        self.print_data = [self.print_data_1d, self.print_data_2d]

    def add_print_data(self, data_array, name, units='-', sub_names=None):
        if len(np.shape(data_array)) == 2:
            if sub_names is None:
                sub_names = [str(i+1) for i in range(len(data_array))]
            self.print_data_2d[name] = \
                {sub_names[i]:
                 {'value': data_array[i], 'units': str(units)}
                 for i in range(len(sub_names))}
        elif len(np.shape(data_array)) == 1:
            self.print_data_1d[name] = \
                {'value': data_array, 'units': str(units)}

        else:
            raise ValueError('argument data_array must be 1- or 2-dimensional')
