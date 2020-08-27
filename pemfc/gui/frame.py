# global imports
import tkinter as tk
from tkinter import ttk

# local imports
from . import widget_set as ws


class FrameFactory:
    @staticmethod
    def create_frame(master, sub_frame_dicts: list = None,
                     widget_set_dicts: list = None, **kwargs):
        if sub_frame_dicts is None:
            return BaseFrame(master, widget_set_dicts=widget_set_dicts,
                             **kwargs)
        else:
            return MainFrame(master, sub_frame_dicts=sub_frame_dicts,
                             widget_set_dicts=widget_set_dicts, **kwargs)


class BaseFrame(tk.Frame):

    PADX = 2
    PADY = 2

    def __init__(self, master, widget_set_dicts: list = None, **kwargs):
        title = kwargs.pop('title', None)
        self.name = None
        if title is not None:
            self.name = title.lower()
        super().__init__(master, name=self.name, **kwargs)
        # self.grid()
        self.widget_factory = ws.WidgetSetFactory()
        # self.widget_sets = []
        if widget_set_dicts is not None:
            self.widget_sets = [self.widget_factory.create_set(self, **ws_dict)
                                for ws_dict in widget_set_dicts]
            # for ws_dict in widget_set_dicts:
            #     self.widget_sets.append(self.widget_factory.create_set(self, **ws_dict))
            #     print(title, self.grid_size(), 'BaseFrame')
        self.title = None
        print(title, self.grid_size(), 'BaseFrame')
        print('master', master)
        if title is not None:
            self.title = self.set_title(title)

    def set_title(self, text, row=0, column=0, font=None, bold=True,
                  **kwargs):
        if text is not None:
            columns = self.grid_size()[0]
            if font is None:
                if bold:
                    font = 'Arial 10 bold'
                else:
                    font = 'Arial 10'
            if columns == 0:
                columnspan = 1
            else:
                columnspan = columns
            title = tk.Label(self, text=text, font=font)
            title.grid(row=row, column=column, columnspan=columnspan,
                       padx=kwargs.get('padx', self.PADX),
                       pady=kwargs.get('pady', self.PADY), **kwargs)
            return title

    def set_grid(self, grid_list=None, **kwargs):
        if grid_list is None:
            for i, widget_set in enumerate(self.widget_sets):
                kwargs.pop('row', None)
                kwargs.pop('column', None)
                row = i + 1
                column = 0
                widget_set.set_grid(row=row, column=column, **kwargs)


class MainFrame(BaseFrame):
    def __init__(self, master, sub_frame_dicts: list = None,
                 widget_set_dicts: list = None, **kwargs):
        title = kwargs.pop('title')
        show_title = kwargs.pop('show_title', False)
        super().__init__(master, widget_set_dicts=widget_set_dicts, **kwargs)
        self.frame_factory = FrameFactory()
        self.sub_frames = []
        if sub_frame_dicts is not None:
            for i, frame_dict in enumerate(sub_frame_dicts):
                sub_frame = self.frame_factory.create_frame(self, **frame_dict)
                # sub_frame.grid(row=i)
                self.sub_frames.append(sub_frame)
        print(title, self.grid_size(), 'MainFrame')
        print('master', master)
        self.set_grid()
        if show_title:
            self.set_title(title)

    def set_grid(self, grid_list=None, **kwargs):
        if grid_list is None:
            for i, sub_frame in enumerate(self.sub_frames):
                kwargs.pop('row', None)
                kwargs.pop('column', None)
                row = i + 1
                column = 0
                sub_frame.set_grid(row=row, column=column, **kwargs)





