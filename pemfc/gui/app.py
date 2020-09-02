# global imports
import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfilename
import os

# local imports
from pemfc.gui import button
from pemfc.gui import frame
from pemfc.gui import input
from pemfc.gui import data_transfer
from pemfc import stack_simulation
from pemfc.src import global_functions as gf


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
        # self.buttons = []
        # if 'button_dicts' in kwargs:
        #     button_dicts = kwargs['button_dicts']
        #     button_dicts = gf.ensure_list(button_dicts)
        #     for button_dict in button_dicts:
        #         button_factory = button.ButtonFactory()
        #
        #         button_widget = \
        #             button_factory.create(self.frames[-1], **button_dict)
        #         self.buttons.append(button_widget)
        #
        #     self.frames[-1].add_widget(self.buttons[0])
        #     self.buttons[0].button.configure(command=self.run)
        #     self.frames[-1].add_widget(self.buttons[1])
        self.frames[-1].sub_frames[-1].buttons[-1].button.configure(command=self.run)

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

    def get_values(self):
        # return [fr.get_values() for fr in self.frames]
        return {fr.name: fr.get_values() for fr in self.frames}

    def run(self):
        values = self.get_values()
        data_transfer.transfer(values, data_transfer.sim_dict)
        stack_simulation.main()
        # print(values)


if __name__ == "__main__":
    root = tk.Tk()
    root.title('PEMFC Model')
    gui_dir = os.path.dirname(os.path.abspath(__file__))
    root.iconbitmap(os.path.join(gui_dir, 'logo-zbt.ico'))

    # root.resizable(False, False)
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)

    base_app = NotebookApp(root, main_frame_dicts=input.main_frame_dicts,
                           button_dicts=[input.run_button_dict,
                                         input.output_dir_button_dict])

    root.mainloop()
