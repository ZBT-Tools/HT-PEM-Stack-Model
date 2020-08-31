# global imports
import tkinter as tk
from tkinter import Grid
from abc import ABC, abstractmethod

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
        grid_location = kwargs.pop('grid_location', (None, None))
        self.row = grid_location[0]
        self.column = grid_location[1]
        row = kwargs.pop('row', None)
        column = kwargs.pop('column', None)
        if row is not None:
            self.row = row
        if column is not None:
            self.column = column
        self.name = label.lower()
        kwargs['text'] = label
        self.sticky = kwargs.pop('sticky', 'NW')
        self.columns = 1
        self.frame = frame
        self.label = tk.Label(frame, **kwargs)
        # self.label.grid(row=self.row, column=self.column, padx=self.padx,
        #                 pady=self.pady, sticky=kwargs.pop('sticky', 'W'))

    def _set_grid(self, widget, **kwargs):
        # Grid.rowconfigure(self.frame, row, weight=1)
        # Grid.columnconfigure(self.frame, column, weight=1)
        # self.frame.rowconfigure(row, weight=1)
        # self.frame.columnconfigure(column, weight=1)
        row = kwargs.pop('row', self.row)
        column = kwargs.pop('column', self.column)
        widget.grid(row=row, column=column,
                    padx=kwargs.get('padx', self.PADX),
                    pady=kwargs.get('pady', self.PADY),
                    sticky=kwargs.pop('sticky', self.sticky), **kwargs)
        return row, column

    def set_grid(self, widget=None, **kwargs):
        row = kwargs.pop('row', self.row)
        column = kwargs.pop('column', self.column)
        if widget is None:
            widget = self.label
        self._set_grid(widget, row=row, column=column, **kwargs)
        return row, column


class Label(WidgetSet):
    def __init__(self, frame, label, **kwargs):
        super().__init__(frame, label, **kwargs)


class MultiWidgetSet(WidgetSet):

    def __init__(self, frame, label, number=1, **kwargs):
        super().__init__(frame, label, **kwargs)
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


class MultiEntrySet(MultiWidgetSet):

    def __init__(self, frame, label, number=1, value=None, **kwargs):
        super().__init__(frame, label, number=number, **kwargs)
        kwargs.pop('grid_location', None)
        if value is not None:
            value = gf.ensure_list(value, length=number)
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
        kwargs.pop('grid_location', None)
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

