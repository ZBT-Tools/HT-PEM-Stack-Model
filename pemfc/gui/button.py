# global imports
import tkinter as tk

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
        else:
            raise NotImplementedError('type of WidgetSet not implemented')


class Button(base.Base):

    PADX = 1
    PADY = 1

    def __init__(self, frame, **kwargs):

        label = kwargs.get('label', '')
        kwargs['text'] = label
        self.name = label.lower()
        super().__init__(self.name, **kwargs)
        self.frame = frame
        self.button = tk.Button(self.frame, text=label, command=self.command)

    def set_grid(self, **kwargs):
        row = kwargs.pop('row', self.row)
        column = kwargs.pop('column', self.column)
        self._set_grid(self.button, row=row, column=column, **kwargs)
        return row, column

    def command(self, *args):
        pass


class RunButton(Button):
    def __init__(self, frame, **kwargs):
        super().__init__(frame, **kwargs)

