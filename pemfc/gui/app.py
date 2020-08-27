# global imports
import tkinter as tk
from tkinter import ttk

# local imports
from pemfc.gui import frame
from pemfc.gui import widget_set as ws
from pemfc.gui import input


class NotebookApp:
    def __init__(self, master, main_frame_dicts=None, **kwargs):
        notebook = ttk.Notebook(master)
        notebook.grid()

        # Make tabs and corresponding frames
        self.frame_factory = frame.FrameFactory()
        self.frames = []
        if main_frame_dicts is not None:
            for i, main_frame_dict in enumerate(main_frame_dicts):
                main_frame = self.frame_factory.create_frame(notebook,
                                                             **main_frame_dict)
                self.frames.append(main_frame)
                notebook.add(main_frame, text=main_frame_dict['title'])

        notebook.select(self.frames[0])
        notebook.enable_traversal()

        # # self.label = tk.Label(self.frames[-1], text=title)
        # # self.label.pack(padx=10, pady=10)
        # self.entries = []
        # for frame in self.frames:
        #     for name in ['pdc1', 'pdc2', 'pdc3', 'pdc4']:
        #         entry = tk.Entry(frame)
        #         entry.grid()
        #         self.entries.append((name, entry))
        #
        #     tk.Button(
        #         frame,
        #         text='test',
        #         command=lambda: print(self.collect_entries())
        #     ).grid()
        #     dim_entry = ws.DimensionedEntrySet(frame, label='test',
        #                                        number=3, dimensions='m')

    # def collect_entries(self):
    #     return {name: entry.get() for name, entry in self.entries}


if __name__ == "__main__":
    root = tk.Tk()
    # make frames
    # frame_names = ['Geometry', 'Physical Properties',
    #                'Operating Conditions', 'Output', 'Simulation']

    # base_app = NotebookApp(root, frame_names=frame_names)
    base_app = NotebookApp(root, main_frame_dicts=input.main_frame_dicts)

    root.mainloop()
