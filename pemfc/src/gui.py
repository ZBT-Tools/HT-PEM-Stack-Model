import tkinter as tk
from tkinter import ttk


class EntrySet:
    def __init__(self, frame, **kwargs):
        column = kwargs.pop('column', 0)
        row = kwargs.get('column', None)
        if row is None:
            row = frame.grid_size()[1] + 1
        self.label = tk.Label(frame, text=kwargs.pop('text'), **kwargs)
        self.label.grid(row=row, column=column)
        self.entry = tk.Entry(frame)
        self.entry.grid(row=row, column=column + 1)


class BaseFrame:
    def __init__(self, master, **kwargs):
        title = kwargs.get('title', None)
        self.frame = tk.Frame(master, **kwargs)
        self.frame.pack(side="top", fill="both", expand=True)
        if title is not None:
            self.label = tk.Label(self.frame, text=title)
            self.label.pack(padx=10, pady=10)


class NotebookApp:
    def __init__(self, master, frames=None, **kwargs):
        frame_names = kwargs.get('tab_names', None)
        notebook = ttk.Notebook(master)
        notebook.pack()

        # Make tabs and corresponding frames
        if frames is None:
            self.frames = []
        for name in frame_names:
            frame = tk.Frame(notebook)
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
                entry.pack()
                self.entries.append((name, entry))

            tk.Button(
                frame,
                text='test',
                command=lambda: print(self.collect_entries())
            ).pack()

    def collect_entries(self):
        return {name: entry.get() for name, entry in self.entries}


if __name__ == "__main__":
    root = tk.Tk()
    tab_names = ['Geometry', 'Physical Properties',
                 'Operating Conditions', 'Output', 'Simulation']
    base_app = NotebookApp(root, tab_names=tab_names)
    root.mainloop()
