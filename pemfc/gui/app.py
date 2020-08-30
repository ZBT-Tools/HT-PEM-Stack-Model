# global imports
import tkinter as tk
from tkinter import ttk
import os

# local imports
from pemfc.gui import frame
from pemfc.gui import widget_set as ws
from pemfc.gui import input


class NotebookApp:
    def __init__(self, master, main_frame_dicts=None, **kwargs):
        self.notebook = ttk.Notebook(master)

        # Make tabs and corresponding frames
        self.frame_factory = frame.FrameFactory()
        self.frames = []
        if main_frame_dicts is not None:
            for i, main_frame_dict in enumerate(main_frame_dicts):
                main_frame = self.frame_factory.create_frame(self.notebook,
                                                             **main_frame_dict)
                self.frames.append(main_frame)
                self.notebook.add(main_frame, text=main_frame_dict['title'])

        self.notebook.select(self.frames[0])
        self.notebook.enable_traversal()
        self.set_grid()

    def set_grid(self, grid_list=None, **kwargs):
        self.notebook.rowconfigure(0, weight=1)
        self.notebook.columnconfigure(0, weight=1)
        # self.notebook.grid(sticky='WENS', **kwargs)

        for fr in self.frames:
            fr.set_grid(grid_list=grid_list, **kwargs)
        self.notebook.grid(sticky='WENS', **kwargs)


if __name__ == "__main__":
    root = tk.Tk()
    root.title('PEMFC Model')
    gui_dir = os.path.dirname(os.path.abspath(__file__))
    root.iconbitmap(os.path.join(gui_dir, 'logo-zbt.ico'))

    # root.resizable(False, False)
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)

    base_app = NotebookApp(root, main_frame_dicts=input.main_frame_dicts)
    # base_app.set_grid()
    root.mainloop()
