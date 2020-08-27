# global imports
import tkinter as tk
from abc import ABC

# local imports
from pemfc.src import global_functions as gf


class WidgetSetFactory:

    def create_set(self, frame, **kwargs):
        type = kwargs.pop('type', None)
        if type == 'EntrySet':
            return self.create_entry_set(frame, **kwargs)
        elif type == 'CheckButtonSet':
            return MultiCheckButtonSet(frame, **kwargs)
        elif type == 'Label':
            return Label(frame, **kwargs)

    @staticmethod
    def create_entry_set(frame, **kwargs):
        dimensions = kwargs.get('dimensions', None)
        if dimensions is None:
            return MultiEntrySet(frame, **kwargs)
        else:
            return DimensionedEntrySet(frame, **kwargs)


class WidgetSet(ABC):

    PADX = 1
    PADY = 1

    def __init__(self, frame, label, **kwargs):
        self.padx = kwargs.pop('padx', self.PADX)
        self.pady = kwargs.pop('pady', self.PADY)
        self.column = kwargs.pop('column', 0)
        self.row = kwargs.pop('row', frame.grid_size()[1] + 1)
        kwargs['text'] = label
        sticky = kwargs.pop('sticky', 'W')
        self.label = tk.Label(frame, **kwargs)
        self.label.grid(row=self.row, column=self.column, padx=self.padx,
                        pady=self.pady, sticky=sticky)


class Label(WidgetSet):
    def __init__(self, frame, label, **kwargs):
        super().__init__(frame, label, **kwargs)


class MultiEntrySet(WidgetSet):

    def __init__(self, frame, label, number=1, value=None, **kwargs):
        super().__init__(frame, label, **kwargs)
        if value is not None:
            value = gf.ensure_list(value, length=number)
        self.entries = []
        for i in range(number):
            entry = tk.Entry(frame, justify='right')
            entry.grid(row=self.row, column=self.column + 1 + i,
                       padx=self.padx, pady=self.pady)
            entry.delete(0, -1)
            if value is not None:
                entry.insert(0, value[i])
            self.entries.append(entry)


class DimensionedEntrySet(MultiEntrySet):
    def __init__(self, frame, label, number=1, dimensions='-',
                 value=None, **kwargs):
        super().__init__(frame, label, number=number, value=value, **kwargs)
        kwargs['text'] = dimensions
        self.dimensions = tk.Label(frame, **kwargs)
        self.dimensions.grid(row=self.row, column=self.column + number + 1,
                             padx=kwargs.pop('padx', self.PADX),
                             pady=kwargs.pop('pady', self.PADY))


class MultiCheckButtonSet(WidgetSet):
    def __init__(self, frame, label, number=1, value=None, **kwargs):
        super().__init__(frame, label, **kwargs)
        if value is not None:
            value = gf.ensure_list(value, length=number)
        self.check_buttons = []
        self.check_vars = []
        for i in range(number):
            check_var = tk.BooleanVar()
            self.check_vars.append(check_var)
            check_button = \
                tk.Checkbutton(frame, variable=check_var, onvalue=True,
                               offvalue=False, **kwargs)
            check_button.grid(row=self.row, column=self.column + 1 + i,
                              padx=self.padx, pady=self.pady)
            if value is not None and value[i] is True:
                check_button.select()
            self.check_buttons.append(check_button)
