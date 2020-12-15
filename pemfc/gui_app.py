# global imports
import tkinter as tk
from tkinter import ttk
import os
import sys
import json

# local imports
from pemfc.gui import frame
from pemfc.gui import input
from pemfc.gui import data_transfer
# from pemfc.gui import icon
from pemfc.src import stack_simulation


class NotebookApp:
    def __init__(self, master, main_frame_dicts=None, **kwargs):
        self.notebook = ttk.Notebook(master)
        self.master = master
        self.registry = {}
        # self.configure_gui()
        # Make tabs and corresponding frames
        self.frame_factory = frame.FrameFactory()
        self.frames = []
        if main_frame_dicts is not None:
            for i, main_frame_dict in enumerate(main_frame_dicts):
                main_frame = self.frame_factory.create_frame(self.notebook,
                                                             **main_frame_dict)
                self.frames.append(main_frame)
                self.notebook.add(main_frame, text=main_frame_dict['title'],
                                  sticky='WENS')

        # configure load settings button
        self.load_settings_button = self.frames[-1].sub_frames[0].widgets[0]
        self.load_settings_button.button.configure(command=self.load_settings)
        self.load_settings_button.filetypes = [('JSON Files', ['.json'])]
        self.load_settings_button.title = 'Please select settings file.'

        # configure save settings button
        self.save_settings_button = self.frames[-1].sub_frames[0].widgets[1]
        self.save_settings_button.button.configure(command=self.save_settings)

        # configure run button
        self.run_button = self.frames[-1].sub_frames[-1].widgets[0]
        self.run_button.button.configure(command=self.run)

        self.notebook.select(self.frames[0])
        self.notebook.enable_traversal()
        self.set_grid()
        self.call_commands()

    def set_grid(self, grid_list=None, **kwargs):
        self.notebook.rowconfigure(0, weight=1)
        self.notebook.columnconfigure(0, weight=1)
        # self.notebook.grid(sticky='WENS', **kwargs)
        for fr in self.frames:
            fr.set_grid(grid_list=grid_list, **kwargs)
        self.notebook.grid(sticky='WEN', **kwargs)

    def get_values(self, get_object=False):
        # return [fr.get_values() for fr in self.frames]
        return {fr.name[-1]: fr.get_values(get_object=get_object)
                for fr in self.frames}

    def call_commands(self):
        for item in self.frames:
            item.call_commands()

    # def configure_gui(self):
    #     self.master.tk.call('wm', 'iconphoto', self.master._w,
    #                         tk.PhotoImage(data=icon.icon()))
    #     self.master.title("Example")
    #     self.master.minsize(250, 50)

    def get_settings(self):
        values = self.get_values()
        return data_transfer.gui_to_sim_transfer(values, data_transfer.sim_dict)

    def save_settings(self):
        settings_dict = self.get_settings()
        self.save_settings_button.file_content = \
            json.dumps(settings_dict, indent=2)
        self.save_settings_button.command()
        return settings_dict

    def load_settings(self):
        widgets_registry = self.get_values(get_object=True)
        settings_dict = json.load(self.load_settings_button.command())
        data_transfer.sim_to_gui_transfer(settings_dict, widgets_registry)

    def run(self):
        settings_dict = self.get_settings()
        data_transfer.save_settings(settings_dict)
        avg_icd = stack_simulation.main()
        # window = tk.Toplevel(self.master)
        # testresult = tk.Label(window, text=str(avg_icd), borderwidth=1)
        # testresult.pack(ipadx=1)
        # print(values)


def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS',
                        os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)


if __name__ == "__main__":
    root = tk.Tk()
    root.title('PEMFC Model')
    if getattr(sys, 'frozen', False):
        # frozen
        gui_dir = os.path.dirname(sys.executable)
    else:
        # unfrozen
        gui_dir = os.path.dirname(os.path.abspath(__file__))
    # gui_dir = os.path.dirname(os.path.abspath(__file__))
    root.iconbitmap(os.path.join(gui_dir, 'logo-zbt.ico'))

    # root.resizable(False, False)
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)

    base_app = NotebookApp(root, main_frame_dicts=input.main_frame_dicts)

    root.mainloop()
