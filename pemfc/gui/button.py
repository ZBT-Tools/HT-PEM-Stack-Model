# global imports
import tkinter as tk
from tkinter import filedialog

# local imports
from . import base


class ButtonFactory:

    @staticmethod
    def create(frame, **kwargs):
        button_type = kwargs.pop('type', None)
        if button_type == 'Button':
            return Button(frame, **kwargs)
        elif button_type == 'RunButton':
            return RunButton(frame, **kwargs)
        elif button_type == 'OpenDirectoryButton':
            return OpenDirectoryButton(frame, kwargs.pop('entry', None),
                                       **kwargs)
        else:
            raise NotImplementedError('type of WidgetSet not implemented')


class Button(base.Base):

    PADX = 1
    PADY = 1

    def __init__(self, frame, **kwargs):

        label = kwargs.pop('label', '')
        self.name = label.lower()
        super().__init__(self.name, **kwargs)
        self.frame = frame
        kwargs = self.remove_dict_entries(
            kwargs, ['row', 'column', 'grid_location', 'columnspan',
                     'rowspan', 'sticky', 'sim_name', 'dtype'])
        self.button = tk.Button(self.frame, text=label, command=self.command,
                                **kwargs)

    def set_grid(self, **kwargs):
        row = kwargs.pop('row', self.row)
        column = kwargs.pop('column', self.column)
        self._set_grid(self.button, row=row, column=column, **kwargs)
        return row, column

    def command(self, *args):
        pass

    def _get_values(self):
        if self.sim_name is None:
            return {'gui_name': self.button.cget('text')}
        else:
            return {'sim_name': self.sim_name,
                    'gui_name': self.button.cget('text')}

    def get_values(self):
        return self._get_values()


class RunButton(Button):
    def __init__(self, frame, **kwargs):
        super().__init__(frame, **kwargs)


class OpenDirectoryButton(Button):
    def __init__(self, frame, entry, **kwargs):
        self.entry = entry
        self.directory = kwargs.pop('directory', None)
        super().__init__(frame, **kwargs)

    def command(self):
        directory = filedialog.askdirectory()
        try:
            self.entry.insert(0, directory)
        except AttributeError:
            raise AttributeError('member entry hast not been referenced '
                                 'correctly to tk.Entry object')
        print(directory)
        return directory


