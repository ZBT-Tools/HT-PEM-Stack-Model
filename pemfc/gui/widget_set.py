# global imports
import tkinter as tk
from tkinter import Grid, ttk
from abc import ABC, abstractmethod

# local imports
from pemfc.src import global_functions as gf
from . import base
from . import entry_value


class WidgetSetFactory:

    def create_set(self, frame, **kwargs):
        widget_type = kwargs.pop('type', None)
        if widget_type == 'EntrySet':
            return self.create_entry_set(frame, **kwargs)
        elif widget_type == 'CheckButtonSet':
            return MultiCheckButtonSet(frame, **kwargs)
        elif widget_type == 'Label':
            return Label(frame, **kwargs)
        elif widget_type == 'OptionMenuSet':
            return OptionMenuSet(frame, **kwargs)
        elif widget_type == 'ComboboxSet':
            return ComboboxSet(frame, **kwargs)
        else:
            raise NotImplementedError('type of WidgetSet not implemented')

    @staticmethod
    def create_entry_set(frame, **kwargs):
        dimensions = kwargs.get('dimensions', None)
        if dimensions is None:
            return MultiEntrySet(frame, **kwargs)
        else:
            return DimensionedEntrySet(frame, **kwargs)


class WidgetSet(base.Base, ABC):

    def __init__(self, frame, label, **kwargs):

        self.columns = 1
        self.name = label.lower().strip(':')
        super().__init__(self.name, **kwargs)
        self.frame = frame
        kwargs = self.remove_dict_entries(
            kwargs, ['padx', 'pady', 'row', 'column', 'grid_location',
                     'sticky', 'sim_name', 'dtype'])
        kwargs['text'] = label
        self.label = tk.Label(frame, **kwargs)
        # self.label.grid(row=self.row, column=self.column, padx=self.padx,
        #                 pady=self.pady, sticky=kwargs.pop('sticky', 'W'))

    def set_grid(self, widget=None, **kwargs):
        row = kwargs.pop('row', self.row)
        column = kwargs.pop('column', self.column)
        if widget is None:
            widget = self.label
        self._set_grid(widget, row=row, column=column, **kwargs)
        return row, column

    def _get_values(self):
        if self.sim_name is None:
            return {'gui_name': self.label.cget('text')}
        else:
            return {'sim_name': self.sim_name,
                    'gui_name': self.label.cget('text')}

    @abstractmethod
    def get_values(self):
        return self._get_values()


class Label(WidgetSet):
    def __init__(self, frame, label, **kwargs):
        super().__init__(frame, label, **kwargs)

    def get_values(self, values=None):
        return super().get_values()


class MultiWidgetSet(WidgetSet, ABC):

    def __init__(self, frame, label, **kwargs):
        super().__init__(frame, label, **kwargs)
        self.dtype = kwargs.pop('dtype', None)
        self.entry_value_factory = entry_value.EntryValueFactory()
        self.widgets = []

    def set_grid(self, widgets=None, **kwargs):
        row = kwargs.pop('row', self.row)
        column = kwargs.pop('column', self.column)
        row, column = super().set_grid(row=row, column=column, **kwargs)
        if widgets is None:
            widgets = self.widgets
        for i, widget in enumerate(widgets):
            column += 1
            super().set_grid(widget=widget, row=row, column=column, **kwargs)
        return row, column

    def get_tk_values(self, tk_objects):
        values = super().get_values()
        if len(tk_objects) > 1:
            values['value'] = []
            for item in tk_objects:
                values['value'].append(
                    self.entry_value_factory.create(item.get(),
                                                    self.dtype, self))
        else:
            values['value'] = \
                self.entry_value_factory.create(tk_objects[0].get(),
                                                self.dtype, self)
        return values

    def get_values(self):
        return self.get_tk_values(self.widgets)


class MultiEntrySet(MultiWidgetSet):

    def __init__(self, frame, label, number=1, value=None, **kwargs):
        super().__init__(frame, label, **kwargs)
        self.dtype = kwargs.pop('dtype', 'float')
        kwargs = self.remove_dict_entries(kwargs, ['grid_location', 'sim_name'])
        if value is not None:
            value = gf.ensure_list(value, length=number)
            number = len(value)
        for i in range(number):
            self.columns += 1
            entry = tk.Entry(frame, justify='right', **kwargs)
            # entry.grid(row=self.row, column=self.column + 1 + i,
            #            padx=self.padx, pady=self.pady)
            entry.delete(0, -1)
            if value is not None:
                entry.insert(0, value[i])
            self.widgets.append(entry)


class DimensionedEntrySet(MultiEntrySet):
    def __init__(self, frame, label, number=1, dimensions='-',
                 value=None, **kwargs):
        super().__init__(frame, label, number=number, value=value, **kwargs)
        kwargs['text'] = dimensions
        kwargs = \
            self.remove_dict_entries(kwargs,
                                     ['grid_location', 'sim_name', 'dtype'])
        self.dimensions = tk.Label(frame, **kwargs)
        self.columns += 1
        # self.dimensions.grid(row=self.row, column=self.column + number + 1,
        #                      padx=kwargs.get('padx', self.PADX),
        #                      pady=kwargs.get('pady', self.PADY))

    def set_grid(self, **kwargs):
        row = kwargs.pop('row', self.row)
        column = kwargs.pop('column', self.column)
        row, column = super().set_grid(row=row, column=column, **kwargs)
        column += 1
        self._set_grid(self.dimensions, row=row, column=column, **kwargs)
        return row, column


class MultiCheckButtonSet(MultiWidgetSet):
    def __init__(self, frame, label, number=1, value=None, **kwargs):
        super().__init__(frame, label, **kwargs)
        self.dtype = kwargs.pop('dtype', 'boolean')
        kwargs = self.remove_dict_entries(kwargs, ['grid_location', 'sim_name'])
        if value is not None:
            value = gf.ensure_list(value, length=number)
        self.check_vars = []
        for i in range(number):
            self.columns += 1
            check_var = tk.BooleanVar()
            self.check_vars.append(check_var)
            check_button = \
                tk.Checkbutton(frame, variable=check_var, onvalue=True,
                               offvalue=False, **kwargs)
            # check_button.grid(row=self.row, column=self.column + 1 + i,
            #                   padx=self.padx, pady=self.pady)
            if value is not None and value[i] is True:
                check_button.select()
            self.widgets.append(check_button)

    def get_values(self):
        return super().get_tk_values(self.check_vars)


class OptionMenuSet(MultiWidgetSet):
    def __init__(self, frame, label, number=1, **kwargs):
        options = kwargs.pop('options', [])
        super().__init__(frame, label, **kwargs)
        self.dtype = kwargs.pop('dtype', 'string')
        self.option_vars = []
        for i in range(number):
            self.columns += 1
            option_var = tk.StringVar()
            self.option_vars.append(option_var)
            option_menu = tk.OptionMenu(self.frame, option_var, *options)
            self.widgets.append(option_menu)


class ComboboxSet(MultiWidgetSet):
    def __init__(self, frame, label, number=1, **kwargs):
        options = kwargs.pop('options', [])
        super().__init__(frame, label, **kwargs)
        self.dtype = kwargs.pop('dtype', 'string')
        for i in range(number):
            self.columns += 1
            option_menu = ttk.Combobox(self.frame, values=options)
            option_menu.current(0)
            self.widgets.append(option_menu)


