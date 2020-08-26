import tkinter as tk
from tkinter import ttk


class MultiEntrySet:

    PADX = 1
    PADY = 1

    def __init__(self, frame, label, number=1, **kwargs):
        padx = kwargs.pop('padx', self.PADX)
        pady = kwargs.pop('pady', self.PADY)
        self.column = kwargs.pop('column', 0)
        self.row = kwargs.get('column', None)
        if self.row is None:
            self.row = frame.grid_size()[1] + 1
        kwargs['text'] = label
        sticky = kwargs.pop('sticky', 'W')
        self.label = tk.Label(frame, **kwargs)
        self.label.grid(row=self.row, column=self.column, padx=padx, pady=pady,
                        sticky=sticky)
        self.entries = []
        for i in range(number):
            entry = tk.Entry(frame)
            entry.grid(row=self.row, column=self.column + 1 + i,
                       padx=padx, pady=pady)
            self.entries.append(tk.Entry(frame))


class DimensionedEntrySet(MultiEntrySet):
    def __init__(self, frame, label, number=1, dimensions='-', **kwargs):
        super().__init__(frame, label, number=number, **kwargs)
        kwargs['text'] = dimensions
        self.dimensions = tk.Label(frame, **kwargs)
        self.dimensions.grid(row=self.row, column=self.column + number + 1,
                             padx=kwargs.pop('padx', self.PADX),
                             pady=kwargs.pop('pady', self.PADY))


def entry_set_factory(frame, **kwargs):
    dimensions = kwargs.get('dimensions', None)
    if dimensions is None:
        return MultiEntrySet(frame, **kwargs)
    else:
        return DimensionedEntrySet(frame, **kwargs)


class BaseFrame(tk.Frame):
    PADX = 2
    PADY = 2

    def __init__(self, master, entry_set_dicts: list = None, **kwargs):
        title = kwargs.pop('title', None)
        self.name = title
        super().__init__(master, name=self.name, **kwargs)
        #self.frame = tk.Frame(master, name=self.name, **kwargs)
        self.grid()
        self.entry_sets = []
        self.entry_set_dicts = entry_set_dicts
        self.create_entry_sets()

        if title is not None:
            self.label = tk.Label(self, text=title)
            self.label.grid(row=0, column=0, padx=self.PADX, pady=self.PADY)

    def create_entry_sets(self):
        for entry_dict in self.entry_set_dicts:
            entry_set = entry_set_factory(self, **entry_dict)
            self.entry_sets.append(entry_set)


class MainFrame(tk.Frame):
    def __init__(self, master, n_subframes=0, **kwargs):
        title = kwargs.get('title', None)
        self.name = title
        super().__init__(master, name=self.name, **kwargs)
        # self.grid()
        if title is not None:
            self.label = tk.Label(self, text=title)
            self.label.pack(padx=10, pady=10)


class NotebookApp:
    def __init__(self, master, frame_input, frame_names, **kwargs):
        notebook = ttk.Notebook(master)
        notebook.grid()

        # Make tabs and corresponding frames
        self.frames = []
        names = frame_names
        for i, name in enumerate(names):
            frame = BaseFrame(notebook, entry_set_dicts=frame_input[i])
            self.frames.append(frame)
            notebook.add(frame, text=name)

        notebook.select(self.frames[0])
        notebook.enable_traversal()

        # self.label = tk.Label(self.frames[-1], text=title)
        # self.label.pack(padx=10, pady=10)
        self.entries = []
        for frame in self.frames:
            for name in ['pdc1', 'pdc2', 'pdc3', 'pdc4']:
                entry = tk.Entry(frame)
                entry.grid()
                self.entries.append((name, entry))

            tk.Button(
                frame,
                text='test',
                command=lambda: print(self.collect_entries())
            ).grid()
            dim_entry = DimensionedEntrySet(frame, label='test',
                                            number=3, dimensions='m')

    def collect_entries(self):
        return {name: entry.get() for name, entry in self.entries}


if __name__ == "__main__":
    root = tk.Tk()
    # make frames
    frame_names = ['Geometry', 'Physical Properties',
                   'Operating Conditions', 'Output', 'Simulation']
    frame_dicts = [[{'label': 'Electrical Conductivity', 'number': 3, 'dimensions': 'm/s'},
                    {'label': 'test_10', 'number': 5, 'dimensions': 'm/s'}],
                   [{'label': 'test_2', 'number': 2, 'dimensions': 'K'}],
                   [{'label': 'test_0', 'number': 2, 'dimensions': '°C'}],
                   [{'label': 'test_4', 'number': 1, 'dimensions': 's'}],
                   [{'label': 'test_4', 'number': 1, 'dimensions': 's'}]]

    # base_app = NotebookApp(root, frame_names=frame_names)
    base_app = NotebookApp(root, frame_input=frame_dicts,
                           frame_names=frame_names)

    root.mainloop()
