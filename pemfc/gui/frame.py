# global imports
import tkinter as tk
from tkinter import Grid
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
        self.title = None
        if title is None:
            self.name = None
        else:
            self.name = title.lower()
        self.sticky = kwargs.pop('sticky', 'WENS')
        super().__init__(master, name=self.name, **kwargs)
        if title is not None and show_title:
            self.title = self.set_title(title, font=font)
        self.grid(sticky=self.sticky)
        self.widget_factory = ws.WidgetSetFactory()
        self.widget_sets = []
        if widget_set_dicts is not None:
            self.widget_sets = [self.widget_factory.create_set(self, **ws_dict)
                                for ws_dict in widget_set_dicts]
            self.columns = np.max([widget_set.columns for widget_set
                                   in self.widget_sets])

    def set_title(self, text, row=0, column=0, font=None, bold=True,
                  **kwargs):
        # if text is not None:
        #     columns = self.grid_size()[0]
        #     if font is None:
        #         if bold:
        #             font = 'Arial 10 bold'
        #         else:
        #             font = 'Arial 10'
        #     if columns == 0:
        #         columnspan = 1
        #     else:
        #         columnspan = columns
        title = ws.Label(self, label=text, font=font, sticky='WENS', **kwargs)
        # title.grid(row=row, column=column, columnspan=columnspan,
        #            padx=kwargs.get('padx', self.PADX),
        #            pady=kwargs.get('pady', self.PADY), **kwargs)
        return title

    def set_grid(self, grid_list=None, **kwargs):
        # columns, rows = self.grid_size()
        # for i in range(rows):
        #     self.rowconfigure(i, weight=1)
        #     for j in range(columns):
        #         self.columnconfigure(j, weight=1)



        # self.grid(sticky=self.sticky)
        self.rowconfigure(kwargs.get('row', 0), weight=1)
        self.columnconfigure(kwargs.get('column', 0), weight=1)

        if grid_list is None:
            for i, widget_set in enumerate(self.widget_sets):
                kwargs.pop('row', None)
                kwargs.pop('column', None)
                row = None
                column = None
                if widget_set.row is None:
                    row = i + 1
                # else:
                #     row = widget_set.row
                if widget_set.column is None:
                    column = 0
                # else:
                #     column = widget_set.column
                Grid.rowconfigure(self, 0, weight=1)
                Grid.columnconfigure(self, 0, weight=1)
                widget_set.set_grid(row=row, column=column, **kwargs)
            if self.title is not None:
                self.title.set_grid(row=0, column=0,
                                    columnspan=self.grid_size()[0])


class MainFrame(BaseFrame):
    def __init__(self, master, sub_frame_dicts: list = None,
                 widget_set_dicts: list = None, **kwargs):
        super().__init__(master, widget_set_dicts=widget_set_dicts, **kwargs)
        self.frame_factory = FrameFactory()
        self.sub_frames = []
        if sub_frame_dicts is not None:
            for i, frame_dict in enumerate(sub_frame_dicts):
                sub_frame = self.frame_factory.create_frame(self, **frame_dict)
                # sub_frame.grid(row=i)
                self.sub_frames.append(sub_frame)

    def set_grid(self, grid_list=None, **kwargs):
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)
        super().set_grid()
        if grid_list is None:
            for i, sub_frame in enumerate(self.sub_frames):
                kwargs.pop('row', None)
                kwargs.pop('column', None)
                row = i + 1
                column = 0
                sub_frame.set_grid(row=row, column=column, **kwargs)
            if self.title is not None:
                row = 0
                column = 0
                Grid.rowconfigure(self, row, weight=1)
                Grid.columnconfigure(self, column, weight=1)
                self.title.set_grid(row=row, column=column,
                                    columnspan=self.grid_size()[0])





