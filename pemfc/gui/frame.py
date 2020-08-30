# global imports
import tkinter as tk
from tkinter import ttk
import numpy as np

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
        show_title = kwargs.pop('show_title', True)
        title = kwargs.pop('title', None)
        font = kwargs.pop('font', None)
        self.grid_location = kwargs.pop('grid_location', (None, None))
        if isinstance(master, ttk.Notebook):
            self.notebook_tab = True
        else:
            self.notebook_tab = False
        self.title = None
        if title is None:
            self.name = None
        else:
            self.name = title.lower()
        self.sticky = kwargs.pop('sticky', 'WENS')
        super().__init__(master, name=self.name, **kwargs)
        if title is not None and show_title:
            self.title = self.set_title(title, font=font)
        # self.grid(sticky=kwargs.pop('sticky', self.sticky),
        #           row=kwargs.pop('row', self.grid_location[0]),
        #           column=kwargs.pop('column', self.grid_location[1]),
        #           padx=kwargs.get('padx', self.PADX),
        #           pady=kwargs.get('pady', self.PADY))
        self.widget_factory = ws.WidgetSetFactory()
        self.widget_sets = []
        if widget_set_dicts is not None:
            self.widget_sets = [self.widget_factory.create_set(self, **ws_dict)
                                for ws_dict in widget_set_dicts]
            self.columns = np.max([widget_set.columns for widget_set
                                   in self.widget_sets])

    def set_title(self, text, font=None, **kwargs):
        title = ws.Label(self, label=text, font=font, sticky='WENS', **kwargs)
        return title

    def set_grid(self, tk_objects=None, grid_list=None, **kwargs):
        row = kwargs.pop('row', 0)
        column = kwargs.pop('column', 0)
        if self.grid_location[0] is not None:
            row = self.grid_location[0]
        if self.grid_location[1] is not None:
            column = self.grid_location[1]
        self.rowconfigure(row, weight=1)
        self.columnconfigure(column, weight=1)
        if not self.notebook_tab:
            self.grid(sticky=kwargs.pop('sticky', self.sticky),
                      row=row, column=column,
                      padx=kwargs.get('padx', self.PADX),
                      pady=kwargs.get('pady', self.PADY), **kwargs)

        if tk_objects is None:
            tk_objects = self.widget_sets
        if grid_list is None:
            kwargs.pop('row', None)
            kwargs.pop('column', None)
            for i, tk_object in enumerate(tk_objects):
                if hasattr(tk_object, 'row') and tk_object.row is not None:
                    row = tk_object.row
                else:
                    row = i + 1
                if hasattr(tk_object, 'column') and tk_object.row is not None:
                    column = tk_object.column
                else:
                    column = 0
                tk_object.set_grid(row=row, column=column, **kwargs)
            if self.title is not None:
                self.title.set_grid(row=0, column=0,
                                    columnspan=self.grid_size()[0])
        # for i in range(self.grid_size()[1]):
        #     for j in range(self.grid_size()[0]):
        #         self.rowconfigure(i, weight=1)
        #         self.columnconfigure(j, weight=1)


class MainFrame(BaseFrame):
    def __init__(self, master, sub_frame_dicts: list = None,
                 widget_set_dicts: list = None, **kwargs):
        super().__init__(master, widget_set_dicts=widget_set_dicts, **kwargs)
        self.frame_factory = FrameFactory()
        self.sub_frames = []
        if sub_frame_dicts is not None:
            for i, frame_dict in enumerate(sub_frame_dicts):
                sub_frame = self.frame_factory.create_frame(self, **frame_dict)
                self.sub_frames.append(sub_frame)

    def set_grid(self, grid_list=None, **kwargs):
        super().set_grid(tk_objects=self.sub_frames, grid_list=grid_list,
                         **kwargs)
        # if grid_list is None:
        #     kwargs.pop('row', None)
        #     kwargs.pop('column', None)
        #     for i, sub_frame in enumerate(self.sub_frames):
        #         row = i + 1
        #         column = 0
        #         sub_frame.set_grid(row=row, column=column, **kwargs)
        #     if self.title is not None:
        #         self.title.set_grid(row=0, column=0,
        #                             columnspan=self.grid_size()[0])





