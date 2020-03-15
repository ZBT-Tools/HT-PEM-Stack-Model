import numpy as np
import string
import weakref
from copy import deepcopy
import collections


class OutputObject:

    PRINT_HIERARCHY = 3
    CLUSTER_NAMES = [['Cell', 'Flow Circuit']]
    _instances = set()

    def __init__(self, name):
        assert isinstance(name, str)
        self._name = name
        self.active = True

        self.print_data_1d = {}
        self.print_data_2d = {}
        self.print_data = [self.print_data_1d, self.print_data_2d]
        self._instances.add(weakref.ref(self))

    def _get_name(self):
        return self._name

    def _set_name(self, name):
        self._name = name

    @property
    def name(self):
        return self._get_name()

    @name.setter
    def name(self, name):
        self._set_name(name)

    @classmethod
    def getinstances(cls):
        dead = set()
        for ref in cls._instances:
            obj = ref()
            if obj is not None:
                yield obj
            else:
                dead.add(ref)
        cls._instances -= dead

    def copy(self):
        copy = deepcopy(self)
        self._instances.add(weakref.ref(copy))
        # self._names.append(self.name)
        # self._instances_strong_ref.add(self)
        return copy

    def add_print_data(self, data_array, name, units='-', sub_names=None):
        if data_array.ndim == 2:
            if sub_names is None:
                sub_names = [str(i+1) for i in range(len(data_array))]
            self.print_data_2d[name] = \
                {sub_names[i]:
                 {'value': data_array[i], 'units': str(units)}
                 for i in range(len(sub_names))}
        elif data_array.ndim == 1:
            self.print_data_1d[name] = \
                {'value': data_array, 'units': str(units)}
        else:
            raise ValueError('argument data_array must be 1- or 2-dimensional')

    def add_print_variables(self, print_variables):
        for i, name in enumerate(print_variables['names']):
            attr = eval('self.' + name)
            description = string.capwords(name.replace('_', ' '))
            units = print_variables['units'][i]
            sub_names = eval(print_variables['sub_names'][i])
            self.add_print_data(attr, description,
                                units=units, sub_names=sub_names)

    @staticmethod
    def combine_print_variables(dict_a, dict_b):
        if dict_b is not None:
            for key in dict_a.keys():
                dict_a[key] += dict_b[key]
        return dict_a

    @classmethod
    def make_name_list(cls):
        name_list = []
        for obj in cls.getinstances():
            obj.name_list = obj.name.split(': ')
            name_list.append(obj.name_list)
        return name_list

    @classmethod
    def cluster_objects(cls):
        cluster = []


