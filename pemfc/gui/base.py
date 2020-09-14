# global imports
from abc import ABC


class Base(ABC):

    PADX = 1
    PADY = 1

    def __init__(self, name, **kwargs):
        self.name = name
        self.sim_name = kwargs.pop('sim_name', None)
        self.padx = kwargs.pop('padx', self.PADX)
        self.pady = kwargs.pop('pady', self.PADY)
        self.rowspan = kwargs.pop('rowspan', None)
        self.columnspan = kwargs.pop('columnspan', None)
        grid_location = kwargs.pop('grid_location', (None, None))
        self.row = grid_location[0]
        self.column = grid_location[1]
        row = kwargs.pop('row', None)
        column = kwargs.pop('column', None)
        if row is not None:
            self.row = row
        if column is not None:
            self.column = column
        self.sticky = kwargs.pop('sticky', 'NW')

    def _set_grid(self, widget, **kwargs):
        # Grid.rowconfigure(self.frame, row, weight=1)
        # Grid.columnconfigure(self.frame, column, weight=1)
        # self.frame.rowconfigure(row, weight=1)
        # self.frame.columnconfigure(column, weight=1)
        row = kwargs.pop('row', self.row)
        column = kwargs.pop('column', self.column)
        widget.grid(row=row, column=column,
                    padx=kwargs.pop('padx', self.PADX),
                    pady=kwargs.pop('pady', self.PADY),
                    rowspan=kwargs.pop('rowspan', self.rowspan),
                    columnspan=kwargs.pop('columnspan', self.columnspan),
                    sticky=kwargs.pop('sticky', self.sticky), **kwargs)
        return row, column

    @staticmethod
    def remove_dict_entries(dictionary, entries):
        for entry in entries:
            dictionary.pop(entry, None)
        return dictionary

