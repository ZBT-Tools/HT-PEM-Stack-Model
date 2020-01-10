import numpy as np
import matplotlib.pyplot as plt
import os
import errno
import system.interpolation as ip
import input.geometry as geom
import shutil
import data.global_parameters as g_par
import system.global_functions as g_func
from itertools import cycle, islice


class Output:

    def __init__(self, dict_output):
        self.save_csv = dict_output['save_csv']
        # switch to save the csv data
        self.save_plot = dict_output['save_plot']
        # switch to save the plot data
        self.show_loss = dict_output['show_loss']
        # switch to show the single voltage losses in the u-i-graph
        self.delimiter = ','
        self.csv_format = '%.9e'
        # object of the class Stack
        self.output_dir = os.path.join(os.path.dirname(__file__), 'output')

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            # shutil.rmtree(self.output_dir, ignore_errors=True)

    def save(self, folder_name, fc_stack):
        if self.save_plot:
            self.output_plots(folder_name, fc_stack)
        if self.save_csv:
            self.output_csv(folder_name, fc_stack)

    @staticmethod
    def clean_directory(directory):
        for file in os.listdir(directory):
            file_path = os.path.join(directory, file)
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path, ignore_errors=True)
            except Exception as e:
                print(e)

    @staticmethod
    def set_ax_properties(ax, **kwargs):
        fontsize = kwargs.get('fontsize', 16)
        if 'xlabel' in kwargs:
            ax.set_xlabel(kwargs['xlabel'], fontsize=fontsize)
        if 'ylabel' in kwargs:
            ax.set_ylabel(kwargs['ylabel'], fontsize=fontsize)
        if 'margins' in kwargs:
            ax.margins(x=kwargs['margins'][0], y=kwargs['margins'][1])
        if 'xlim' in kwargs:
            ax.set_xlim(kwargs['xlim'])
        if 'ylim' in kwargs:
            ax.set_ylim(kwargs['ylim'])
        if 'xticks' in kwargs:
            ax.set_xticks(kwargs['xticks'])
        if 'yticks' in kwargs:
            ax.set_yticks(kwargs['yticks'])
        if 'labels' in kwargs:
            ax.legend(kwargs['labels'], fontsize=fontsize)
        if 'title' in kwargs:
            ax.set_title(kwargs['title'], fontsize=fontsize)
        return ax

    def plot_lines(self, ax, x, y, colormap='', **kwargs):
        x = np.asarray(x)
        y = np.asarray(y)
        ny = len(y)

        if x.ndim != y.ndim:
            if x.ndim == 1:
                x = np.tile(x, (ny, 1))
            else:
                raise ValueError('Outer dimension of x is not one and not '
                                 'equal to outer dimension of y')
        if y.ndim == 1:
            ax.plot(x, y, marker=kwargs.get('marker', '.'),
                    markersize=kwargs.get('markersize', 5.0),
                    fillstyle=kwargs.get('fillstyle', 'full'),
                    linewidth=kwargs.get('linewidth', 1.0),
                    linestyle=kwargs.get('linestyle', '-'),
                    color=kwargs.get('color', 'k'))
        else:
            if colormap:
                cmap = plt.get_cmap(colormap)
                colors = cmap(np.linspace(0.0, 1.0, ny))
            else:
                colors = \
                    kwargs.get('color',
                               list(islice(cycle(['k', 'b', 'r', 'g', 'y']),
                                           ny)))
            linestyles = \
                list(islice(cycle(kwargs.get('linestyle', ['-'])), ny))
            markers = \
                list(islice(cycle(kwargs.get('marker', ['.'])), ny))
            fillstyles = \
                list(islice(cycle(kwargs.get('fillstyle', ['full'])), ny))
            for i in range(ny):
                ax.plot(x[i], y[i], marker=markers[i],
                        markersize=kwargs.get('markersize', 5.0),
                        fillstyle=fillstyles[i],
                        linewidth=kwargs.get('linewidth', 1.0),
                        linestyle=linestyles[i],
                        color=colors[i])
        ax.grid(True)
        ax.use_sticky_edges = False
        ax.autoscale()
        ax.set_xscale(kwargs.get('xscale', 'linear'))
        ax.set_yscale(kwargs.get('yscale', 'linear'))
        ax = self.set_ax_properties(ax, **kwargs)
        return ax

    def create_figure(self, filepath, x_array, y_array, xlabels, ylabels,
                      xlims=None, ylims=None, xticks=None, yticks=None,
                      titles=None, legend=None, rows=1, cols=1,
                      **kwargs):
        nplots = rows*cols

        def check_dims(variable, correct_single_dim=False):
            if isinstance(variable, str):
                variable = [variable]
            if not (isinstance(variable, (list, tuple, np.ndarray))):
                raise TypeError('variable must be provided '
                                'as tuple, list or numpy array')
            if len(variable) != nplots:
                if correct_single_dim:
                    if nplots == 1:
                        variable = [variable]
                    else:
                        raise ValueError('variable must be sequence with '
                                         'length equivalent to number of plots')
                else:
                    raise ValueError('variable must be sequence with '
                                     'length equivalent to number of plots')
            return variable

        if rows > 2:
            figsize = kwargs.get('figsize', (6.4, 4.8*float(rows)/2.0))
        else:
            figsize = kwargs.get('figsize', (6.4, 4.8))
        fig = plt.figure(dpi=kwargs.get('dpi', 150), figsize=figsize)

        x_array = np.asarray(x_array)
        y_array = check_dims(np.asarray(y_array), correct_single_dim=True)

        if len(x_array) != nplots:
            if x_array.ndim == 1:
                x_array = np.tile(x_array, (nplots, 1))
            else:
                raise ValueError('Dimension of x-array is not one and does not '
                                 'match number of plot')
        fontsize = kwargs.get('fontsize', 16)
        xlabels = check_dims(xlabels)
        ylabels = check_dims(ylabels)

        for i in range(nplots):
            ax = fig.add_subplot(rows, cols, i+1)
            ax = self.plot_lines(ax, x_array[i], y_array[i],
                                 xlabel=xlabels[i], ylabel=ylabels[i], **kwargs)
            if titles is not None:
                titles = check_dims(titles)
                ax.set_title(titles[i], fontsize=fontsize)
            if legend is not None:
                legend = check_dims(legend, correct_single_dim=True)
                ax.legend(legend[i])
            if xlims is not None:
                xlims = check_dims(xlims, correct_single_dim=True)
                ax.set_xlim(xlims[i])
            if ylims is not None:
                xlims = check_dims(xlims, correct_single_dim=True)
                ax.set_ylim(ylims[i])
            if xticks is not None:
                xticks = check_dims(xticks, correct_single_dim=True)
                ax.set_xticks(xticks[i])
            if yticks is not None:
                xlims = check_dims(yticks, correct_single_dim=True)
                ax.set_yticks(yticks[i])
        plt.tight_layout()
        if filepath:
            fig.savefig(filepath, format=kwargs.get('fileformat', 'png'))
        return fig

    @staticmethod
    def plot(y_values, y_label, x_label, y_scale, colors,
             title, xlim_low, xlim_up, labels, path):
        if labels:
            for i in range(len(y_values)):
                plt.plot(y_values[i], color=colors[i],
                         marker='.', label=labels[i])
        else:
            for i in range(len(y_values)):
                plt.plot(y_values[i], color=colors[i], marker='.')

        plt.xlabel(x_label, fontsize=16)
        plt.ylabel(y_label, fontsize=16)
        plt.yscale(y_scale)
        plt.tick_params(labelsize=14)
        plt.xlim(xlim_low, xlim_up)
        plt.tight_layout()
        plt.grid()
        if labels:
            plt.legend()
        plt.savefig(os.path.join(path, title + '.png'))
        plt.close()

    @staticmethod
    def x_plot(path, x, y, x_label, y_label, x_scale='linear',
               y_scale='linear', xlim=None, ylim=None, title=None, labels=None):
        if labels:
            for i in range(len(y)):
                plt.plot(x, y[i],
                         color=plt.cm.coolwarm(i / len(y)),
                         marker='.', label=labels[i])
        else:
            for i in range(len(y)):
                plt.plot(x, y[i],
                         color=plt.cm.coolwarm(i / len(y)), marker='.')

        plt.xlabel(x_label, fontsize=16)
        plt.ylabel(y_label, fontsize=16)
        plt.yscale(y_scale)
        plt.tick_params(labelsize=14)
        plt.autoscale(tight=True, axis='both', enable=True)
        if xlim:
            plt.xlim(xlim[0], xlim[1])
        if ylim:
            plt.ylim(ylim[0], ylim[1])
        plt.tight_layout()
        plt.grid()
        if labels:
            plt.legend()
        plt.savefig(os.path.join(path, title + '.png'))
        plt.close()

    def write_array_to_csv(self, file_path, array,
                           header='', separator_lines=None):
        with open(file_path, 'w') as file:
            if header:
                file.write('# ' + header + '\n')
            if separator_lines:
                for i in range(len(separator_lines)):
                    a = array[i]
                    if a.ndim == 1:
                        a = a.reshape(1, a.shape[0])
                    file.write(separator_lines[i])
                    np.savetxt(file, a,
                               delimiter=self.delimiter, fmt=self.csv_format)
            else:
                np.savetxt(file, array,
                           delimiter=self.delimiter, fmt=self.csv_format)
        return file

    def write_stack_array_to_csv(self, file_path, stack_array,
                                 header='', separator_lines=None):
        if not separator_lines:
            n = len(stack_array)
            separator_lines = ['# Cell ' + str(i) + '\n' for i in range(n)]
        return self.write_array_to_csv(file_path, stack_array,
                                       header, separator_lines)

    # def plot_cell_var(self, path, cells, y_var, y_label,
    #                   title, x_lim, x_var=None, y_lim=None,
    #                   x_label='Channel Location $[m]$', y_scale='linear'):
    #     """
    #     Creates plots by given input values
    #     """
    #     for i, cell in enumerate(cells):
    #         plt.plot(x_var, eval('cell.' + y_var),
    #                  color=plt.cm.coolwarm(i / len(cells)), marker='.')
    #     plt.xlabel(x_label, fontsize=16)
    #     plt.ylabel(y_label, fontsize=16)
    #     plt.yscale(y_scale)
    #     plt.tick_params(labelsize=14)
    #     plt.autoscale(tight=True, axis='both', enable=True)
    #     plt.xlim(x_lim[0], x_lim[1])
    #     if y_lim:
    #         plt.ylim(y_lim[0], y_lim[1])
    #     plt.tight_layout()
    #     plt.grid()
    #     plt.savefig(os.path.join(path, title + '.png'))
    #     plt.close()
    #
    # def plot_cell_var_2(self, path, cells, y_var, y_label,
    #                   title, x_lim, x_var=None, y_lim=None,
    #                   x_label='Channel Location $[m]$', y_scale='linear'):
    #     """
    #     Creates plots by given input values
    #     """
    #     for i, cell in enumerate(cells):
    #         plt.plot(x_var, eval('cell.' + y_var),
    #                  color=plt.cm.coolwarm(i / len(cells)), marker='.')
    #     plt.xlabel(x_label, fontsize=16)
    #     plt.ylabel(y_label, fontsize=16)
    #     plt.yscale(y_scale)
    #     plt.tick_params(labelsize=14)
    #     plt.autoscale(tight=True, axis='both', enable=True)
    #     plt.xlim(x_lim[0], x_lim[1])
    #     if y_lim:
    #         plt.ylim(y_lim[0], y_lim[1])
    #     plt.tight_layout()
    #     plt.grid()
    #     plt.savefig(os.path.join(path, title + '.png'))
    #     plt.close()

    def output_plots(self, folder_name, fc_stack):
        """
        Coordinates the plot sequence
        """
        path = os.path.join(self.output_dir, folder_name, 'plots')

        if not os.path.exists(path):
            os.makedirs(path)
        else:
            self.clean_directory(path)
        length = fc_stack.cells[0].cathode.channel.length
        x_lims = [0.0, length]
        x_node = fc_stack.cells[0].cathode.channel.x
        x_ele = ip.interpolate_1d(x_node)
        # self.plot([self.mdf_criteria_process, self.i_ca_criteria_process,
        #            self.temp_criteria_process], 'ERR', 'Iteration', 'log',
        #           ['k', 'r', 'b'], 'Convergence', 0.,
        #           len(self.temp_criteria_process),
        #           ['Flow Distribution', 'Current Density', 'Temperature'], path)
        # self.x_plot(fc_stack.i_cd, x_ele, 'Current Density $[A/m²]$',
        #             'Channel Location $[m]$', 'linear', 'Current Density',
        #             False, x_lims, path)
        # self.create_figure(x_ele, fc_stack.i_cd,
        #                    'Channel Location $[m]$',
        #                    'Current Density $[A/m²]$', colormap='coolwarm',
        #                    filepath=os.path.join(path, 'current_density.png'))
        n_cells = fc_stack.n_cells
        cell_range = np.asarray(range(n_cells))
        if n_cells > 1:
            self.plot([fc_stack.manifold[0].cell_stoi,
                       fc_stack.manifold[1].cell_stoi],
                      'Stoichiometry', 'Cell Number', 'linear', ['k', 'r'],
                      'Stoichimetry Distribution', 0., n_cells - 1,
                      ['Cathode', 'Anode'], path)
            file_path = os.path.join(path,'stoichimetry_distribution.png')
            self.create_figure(file_path, cell_range,
                               [fc_stack.manifold[0].cell_stoi,
                                fc_stack.manifold[1].cell_stoi],
                               'Cell Number', 'Stoichiometry',
                               legend=['Cathode', 'Anode'], xticks=cell_range)
            self.plot([fc_stack.manifold[0].cell_stoi / 2.5],
                      'Flow Distribution', 'Cell Number', 'linear', ['k'],
                      'Distribution', 0., n_cells - 1, ['Cathode'], path)
        # self.plot_cell_var(path, fc_stack.cells, 'v',
        #                    'Cell Voltage $[V]$', 'Channel Location $[m]$',
        #                    'Cell Voltage', x_lims, x_ele, [0.0, 1.0])
        # self.x_plot(fc_stack.temp_sys.temp_cool, x_node,
        #             'Coolant Temperature [K]', 'Channel Location $[m]$',
        #             'linear', 'Coolant Temperature', False,
        #             x_lims, path)
        # self.plot_cell_var(path, fc_stack.cells, 'temp[-1]',
        #                    'Anode BPP - GDE Temperature $[K]$',
        #                    'Channel Location $[m]$',
        #                    'Anode Plate - GDE Temperature', x_lims, x_ele)
        # self.plot_cell_var(path, fc_stack.cells, 'temp[-2]',
        #                    'Anode GDE - MEM Temperature $[K]$',
        #                    'Channel Location $[m]$',
        #                    'Anode GDE - Membrane Temperature', x_lims, x_ele)
        # self.plot_cell_var(path, fc_stack.cells, 'temp[2]',
        #                    'Cathode GDE - MEM Temperature $[K]$',
        #                    'Channel Location $[m]$', 'Cathode GDL Temperature',
        #                    x_lims, x_ele)
        # self.plot_cell_var(path, fc_stack.cells, 'cathode.temp_fluid',
        #                    'Cathode Fluid Temperature $[K]$',
        #                    'Channel Location $[m]$',
        #                    'Cathode_Channel_Temperature', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'temp[1]',
        #                    'Cathode BPP-GDE Temperature $[K]$',
        #                    'Channel Location $[m]$',
        #                    'Cathode GDE - Plate Temperature', x_lims, x_ele)
        # self.plot_cell_var(path, fc_stack.cells, 'anode.temp_fluid',
        #                    'Anode Fluid Temperature $[K]$',
        #                    'Channel Location $[m]$',
        #                    'Anode_Channel_Temperature', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'temp[0]',
        #                    'BPP - BPP Temperature $[K]$',
        #                    'Channel Location $[m]$',
        #                    'Coolant Plate Temperature', x_lims, x_ele)
        # self.plot_cell_var(path, fc_stack.cells, 'cathode.mol_flow[0] * 1.e3',
        #                    'Cathode Oxygen Molar Flow $[mmol/s]$',
        #                    'Channel Location $[m]$',
        #                    'Cathode Oxygen Molar Flow', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'cathode.mol_flow[1] * 1.e3',
        #                    'Cathode Water Molar Flow $[mmol/s]$',
        #                    'Channel Location $[m]$',
        #                    'Cathode Water Molar Flow', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'cathode.mol_flow[2] * 1.e3',
        #                    'Cathode Nitrogen Molar Flow $[mmol/s]$',
        #                    'Channel Location $[m]$',
        #                    'Cathode Nitrogen Molar Flow', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'anode.mol_flow[0] * 1.e3',
        #                    'Anode Hydrogen Molar Flow $[mmol/s]$',
        #                    'Channel Location $[m]$',
        #                    'Anode Hydrogen Molar Flow', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'anode.mol_flow[1] * 1.e3',
        #                    'Anode Water Molar Flow $[mmol/s]$',
        #                    'Channel Location $[m]$',
        #                    'Anode Water Molar Flow', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'anode.mol_flow[2] * 1.e3',
        #                    'Anode Nitrogen Molar Flow $[mmol/s]$',
        #                    'Channel Location $[m]$',
        #                    'Anode Nitrogen Molar Flow', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'cathode.mol_fraction[0]',
        #                    'Oxygen Molar Fraction',
        #                    'Channel Location $[m]$',
        #                    'Oxygen Molar Fraction', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'cathode.mol_fraction[1]',
        #                    'Cathode Gas Water Molar Fraction',
        #                    'Channel Location $[m]$',
        #                    'Water Molar Fraction Cathode',
        #                    x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'cathode.mol_fraction[2]',
        #                    'Cathode Nitrogen Molar Fraction',
        #                    'Channel Location $[m]$',
        #                    'Nitrogen Molar Fraction Cathode', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'anode.mol_fraction[0]',
        #                    'Hydrogen Molar Fraction',
        #                    'Channel Location $[m]$',
        #                    'Hydrogen_Molar_Fraction_Anode', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'anode.mol_fraction[1]',
        #                    'Anode Gas Water Molar Fraction',
        #                    'Channel Location $[m]$',
        #                    'Water_Molar_Fraction_Anode', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'anode.mol_fraction[2]',
        #                    'Anode Nitrogen Molar Fraction',
        #                    'Channel Location $[m]$',
        #                    'Nitrogen_Molar_Fraction_Anode', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells,
        #                    'cathode.mol_flow_liq[1] * 1.e3',
        #                    'Cathode Liquid Water Flow $[mmol/s]$',
        #                    'Channel Location $[m]$',
        #                    'Liquid Water Flow Cathode', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'cathode.cond_rate * 1.e3',
        #                    'Cathode Water Condensation Rate $[mmol/s]$',
        #                    'Channel Location $[m]$',
        #                    'Water Condensation Rate Cathode', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'cathode.humidity',
        #                    'Cathode Relative Humidity $[-]$',
        #                    'Channel Location $[m]$',
        #                    'Relative Humidity Cathode', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells,
        #                    'cathode.mass_flow_gas_total * 1.e6',
        #                    'Cathode Channel Gas Massflow $[mg/s]$',
        #                    'Channel Location $[m]$',
        #                    'Cathode_Channel_Gas_Massflow', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells,
        #                    'cathode.mass_flow_total * 1.e6',
        #                    'Cathode Channel Fluid Massflow $[mg/s]$',
        #                    'Channel Location $[m]$',
        #                    'Cathode_Channel_Fluid_Massflow', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'cathode.g_fluid * 1.e3',
        #                    'Cathode Capacity Flow $[mW/K]$',
        #                    'Channel Location $[m]$',
        #                    'Cathode Capacity Flow', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'cathode.mass_flow[0] * 1.e6',
        #                    'Oxygen Massflow $[mg/s]$', 'Channel Location $[m]$',
        #                    'Oxygen_massflow', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells,
        #                    'cathode.mass_flow_gas[1] * 1.e6',
        #                    'Cathode H2O Vapour Massflow $[mg/s]$',
        #                    'Channel Location $[m]$',
        #                    'Cathode H2O Vapour Massflow', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'anode.mass_flow[0] * 1.e6',
        #                    'Hydrogen Mass Flow $[mg/s]$',
        #                    'Channel Location $[m]$',
        #                    'Hydrogen Mass Flow', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'cathode.cp_fluid',
        #                    'Cathode Heat Capacity $[J/(kgK)]$',
        #                    'Channel Location $[m]$',
        #                    'Cathode Heat Capacity', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'cathode.p',
        #                    'Cathode Channel Pressure $[Pa]$',
        #                    'Channel Location $[m]$',
        #                    'Cathode Channel Pressure', x_lims, x_node)
        # self.plot_cell_var(path, fc_stack.cells, 'anode.p',
        #                    'Anode Channel Pressure $[Pa]$',
        #                    'Channel Location $[m]$',
        #                    'Anode Channel Pressure', x_lims, x_node)
        # Z-Axis-Temperature Plot
        x_vec_z = np.array([0.,
                           geom.bipolar_plate_thickness,
                           geom.gas_diffusion_layer_thickness,
                           geom.membrane_thickness,
                           geom.gas_diffusion_layer_thickness])
        x_vec_e = np.array([geom.bipolar_plate_thickness,
                            geom.bipolar_plate_thickness,
                            geom.gas_diffusion_layer_thickness,
                            geom.membrane_thickness,
                            geom.gas_diffusion_layer_thickness])
        x_vec_l = np.array([geom.bipolar_plate_thickness,
                            geom.bipolar_plate_thickness,
                            geom.gas_diffusion_layer_thickness,
                            geom.membrane_thickness,
                            geom.gas_diffusion_layer_thickness,
                            geom.bipolar_plate_thickness])
        x = []
        n_cells = fc_stack.n_cells
        n_ele = g_par.dict_case['elements']
        for j in range(n_cells):
            if j is 0:
                x.append(x_vec_z)
            elif 0 < j < n_cells - 1:
                x.append(x_vec_e)
            else:
                x.append(x_vec_l)
        x = np.cumsum(np.block(x))
        t = fc_stack.temp_sys.temp_layer
        for i in range(n_ele):
            t_vec = []
            for j in range(n_cells):
                if j is not n_cells - 1:
                    t_vec.append(np.array([t[j][0, i], t[j][1, i],
                                           t[j][2, i], t[j][3, i], t[j][4, i]]))
                else:
                    t_vec.append(np.array([t[j][0, i], t[j][1, i], t[j][2, i],
                                           t[j][3, i], t[j][4, i], t[j][5, i]]))
            plt.plot(x, np.block(t_vec),
                     color=plt.cm.coolwarm((i + 1.e-20) / float(n_ele)))
        plt.xlim(0, x[-1])
        plt.xlabel('Stack Location $[m]$', fontsize=16)
        plt.ylabel('Temperature $[K]$', fontsize=16)
        plt.tick_params(labelsize=14)
        plt.autoscale(tight=True, axis='both', enable=True)
        plt.tight_layout()
        plt.savefig(os.path.join(path, 'z-Cut-Temperature.png'), format='png')
        plt.close()

        # for j in range(n_cells):
        #     print(np.average(fc_stack.i_cd[j, :]))

    def output_csv(self, folder_name, fc_stack):
        csv_path = os.path.join(self.output_dir, folder_name, 'csv_data')
        plot_path = os.path.join(self.output_dir, folder_name, 'plots')
        if not os.path.exists(csv_path):
            os.makedirs(csv_path)
        else:
            self.clean_directory(csv_path)
        if not os.path.exists(plot_path):
            os.makedirs(plot_path)
        # else:
        #    self.clean_directory(plot_path)

        cells = fc_stack.cells
        x_node = cells[0].cathode.channel.x
        x_label = 'Channel Location $[m]$'
        x_ele = ip.interpolate_1d(x_node)

        def write_data(name, array, colormap='coolwarm'):
            if self.save_csv:
                file_name = name.replace(' ', '_') + '.csv'
                file_path = os.path.join(csv_path, file_name)
                header = name + ', ' + 'Units: ' + content['units'].strip('$')
                self.write_stack_array_to_csv(file_path, array,
                                              header=header)
            if self.save_plot:
                file_name = name.replace(' ', '_') + '.png'
                file_path = os.path.join(plot_path, file_name)
                if array.shape[-1] == len(x_node):
                    x = x_node
                else:
                    x = x_ele
                y_label = name + ' ' + '$[' + content['units'] + ']$'
                self.create_figure(file_path, x, array, x_label, y_label,
                                   colormap=colormap)

        # Write cell values
        for name, content in cells[0].print_data[0].items():
            value = content['value']
            var_array = \
                g_func.construct_empty_stack_array(value, fc_stack.n_cells)
            for i, cell in enumerate(cells):
                var_array[i] = cell.print_data[0][name]['value']
            write_data(name, var_array)

        for base_name, sub_dict in cells[0].print_data[1].items():
            for sub_name, content in sub_dict.items():
                value = content['value']
                var_array = \
                    g_func.construct_empty_stack_array(value, fc_stack.n_cells)
                for i, cell in enumerate(cells):
                    var_array[i] = \
                        cell.print_data[1][base_name][sub_name]['value']
                write_data(sub_name + ' ' + base_name, var_array)


        # Write half cell values
        for i in range(len(cells[0].half_cells)):
            electrode_name = cells[0].half_cells[i].name
            for name, content in cells[0].half_cells[i].print_data[0].items():
                value = content['value']
                var_array = \
                    g_func.construct_empty_stack_array(value, fc_stack.n_cells)
                for j, cell in enumerate(cells):
                    var_array[j] = \
                        cell.half_cells[i].print_data[0][name]['value']
                write_data(electrode_name + ' ' + name, var_array)

            for base_name, sub_dict \
                    in cells[0].half_cells[i].print_data[1].items():
                for sub_name, content in sub_dict.items():
                    value = content['value']
                    var_array = \
                        g_func.construct_empty_stack_array(value,
                                                           fc_stack.n_cells)
                    for j, cell in enumerate(cells):
                        var_array[j] = \
                            cell.half_cells[i].print_data[1][base_name][sub_name]['value']
                    write_data(electrode_name + ' ' + sub_name + ' ' + base_name,
                               var_array)
                # file_name = electrode_name + '_' \
                #     + sub_name.replace(' ', '_') + '.csv'
                # file_path = os.path.join(csv_path, file_name)
                # header = electrode_name + ' ' + sub_name + ', ' + 'Units: ' \
                #     + content['units'].strip('$')
                # self.write_stack_array_to_csv(file_path, var_array,
                #                               header=header)

        # for i, item in enumerate(fc_stack.cells):
        #     self.mol_flow[0, i] = item.cathode.mol_flow[0]
        #     self.mol_flow[1, i] = item.cathode.mol_flow[1]
        #     self.mol_flow[2, i] = item.cathode.mol_flow[2]
        #     self.mol_flow[3, i] = item.anode.mol_flow[0]
        #     self.mol_flow[4, i] = item.anode.mol_flow[1]
        #     self.mol_flow[5, i] = item.anode.mol_flow[2]
        #     self.gas_con[0, i] = item.cathode.gas_con[0]
        #     self.gas_con[1, i] = item.cathode.gas_con[1]
        #     self.gas_con[2, i] = item.cathode.gas_con[2]
        #     self.gas_con[3, i] = item.anode.gas_con[0]
        #     self.gas_con[4, i] = item.anode.gas_con[1]
        #     self.gas_con[5, i] = item.anode.gas_con[2]
        #     self.m_f[0, i] = item.cathode.mass_f[0]
        #     self.m_f[1, i] = item.cathode.mass_f[1]
        #     self.m_f[2, i] = item.cathode.mass_f[2]
        #     self.m_f[3, i] = item.anode.mass_f[0]
        #     self.m_f[4, i] = item.anode.mass_f[1]
        #     self.m_f[5, i] = item.anode.mass_f[2]
        #     self.mol_f[0, i] = item.cathode.mol_f[0]
        #     self.mol_f[1, i] = item.cathode.mol_f[1]
        #     self.mol_f[2, i] = item.cathode.mol_f[2]
        #     self.mol_f[3, i] = item.anode.mol_f[0]
        #     self.mol_f[4, i] = item.anode.mol_f[1]
        #     self.mol_f[5, i] = item.anode.mol_f[2]
        #     self.v_loss[i] = item.v_loss
        #     self.v_cell[i] = item.v
        #     self.cp[0, i] = item.cathode.cp[0]
        #     self.cp[1, i] = item.cathode.cp[1]
        #     self.cp[2, i] = item.cathode.cp[2]
        #     self.cp[3, i] = item.anode.cp[0]
        #     self.cp[4, i] = item.anode.cp[1]
        #     self.cp[5, i] = item.anode.cp[2]
        #     self.lambda_gas[0, i] = item.cathode.lambda_gas[0]
        #     self.lambda_gas[1, i] = item.cathode.lambda_gas[1]
        #     self.lambda_gas[2, i] = item.cathode.lambda_gas[2]
        #     self.lambda_gas[3, i] = item.anode.lambda_gas[0]
        #     self.lambda_gas[4, i] = item.anode.lambda_gas[1]
        #     self.lambda_gas[5, i] = item.anode.lambda_gas[2]
        #     self.visc[0, i] = item.cathode.visc[0]
        #     self.visc[1, i] = item.cathode.visc[1]
        #     self.visc[2, i] = item.cathode.visc[2]
        #     self.visc[3, i] = item.anode.visc[0]
        #     self.visc[4, i] = item.anode.visc[1]
        #     self.visc[5, i] = item.anode.visc[2]
        #     self.r_gas_cat[i] = item.cathode.r_gas
        #     self.r_gas_ano[i] = item.anode.r_gas
        #     self.cp_gas_cat[i] = item.cathode.cp_gas
        #     self.cp_gas_ano[i] = item.anode.cp_gas
        #     self.visc_gas_cat[i] = item.cathode.visc_gas
        #     self.visc_gas_ano[i] = item.anode.visc_gas
        #     self.lambda_gas_cat[i] = item.cathode.lambda_gas
        #     self.lambda_gas_ano[i] = item.anode.lambda_gas
        #     self.rho_gas_cat[i] = item.cathode.rho_gas
        #     self.rho_gas_cat[i] = item.anode.rho_gas
        #     self.u_gas_cat[i] = item.cathode.u
        #     self.u_gas_ano[i] = item.anode.u
        #     self.p_cat[i] = item.cathode.p
        #     self.p_ano[i] = item.anode.p
        #     self.ht_coef_cat[i] = item.cathode.ht_coef
        #     self.ht_coef_ano[i] = item.anode.ht_coef
        #     self.cp_fluid_cat[i] = item.cathode.cp_fluid
        #     self.cp_fluid_ano[i] = item.anode.cp_fluid
        #     self.m_flow_fluid_cat[i] = item.cathode.m_flow_fluid
        #     self.m_flow_fluid_ano[i] = item.anode.m_flow_fluid
        #     self.temp_fluid_cat[i] = item.cathode.temp_fluid
        #     self.temp_fluid_ano[i] = item.anode.temp_fluid
        #     self.stoi_cat[i] = item.cathode.stoi
        #     self.stoi_ano[i] = item.anode.stoi
        #     for folder_name in range(5):
        #         self.temp_layer.append(
        #             fc_stack.temp_cpl_stack.temp_layer[i][folder_name, :])
        # np.savetxt(path + 'Temperature Layer.csv',
        #            self.temp_layer, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Coolant Temperature.csv',
        #            fc_stack.temp_cpl_stack.temp_cool,
        #            delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Current Density.csv', fc_stack.i_ca,
        #            delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Cathode Gas Temperature.csv',
        #            self.temp_fluid_cat, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Anode Gas Temperature.csv',
        #            self.temp_fluid_ano, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Cathode Channel Average Velocity.csv',
        #            self.u_gas_cat, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Anode Channel Average Velocity.csv',
        #            self.u_gas_ano, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Cathode Channel Pressure.csv',
        #            self.p_cat, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Anode Channel Pressure.csv',
        #            self.p_ano, delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Oxygen Molar Flow.csv',
        #            self.mol_flow[0], delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Cathode Water Molar Flow.csv',
        #            self.mol_flow[1], delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Cathode Nitrogen Molar Flow.csv',
        #            self.mol_flow[2], delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Hydrogen Molar Flow.csv',
        #            self.mol_flow[3], delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Anode Water Molar Flow.csv',
        #            self.mol_flow[4], delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Oxygen Molar Concentration.csv',
        #            self.gas_con[0], delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Cathode Water Molar Concentration.csv',
        #            self.gas_con[1], delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Nitrogen Molar Concentration.csv',
        #            self.gas_con[2], delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Hydrogen Molar Concentration.csv',
        #            self.gas_con[3], delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Anode Water Molar Concentration.csv',
        #            self.gas_con[4], delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Oxygen Molar Fraction.csv',
        #            self.mol_f[0], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Cathode Water Molar Fraction.csv',
        #            self.mol_f[1], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Nitrogen Molar Fraction.csv',
        #            self.mol_f[2], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Hydrogen Molar Fraction.csv',
        #            self.mol_f[3], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Anode Water Molar Fraction.csv',
        #            self.mol_f[4], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Oxygen Mass Fraction.csv',
        #            self.mol_f[0], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Cathode Mass Fraction.csv',
        #            self.m_f[1], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Nitrogen Mass Fraction.csv',
        #            self.m_f[2], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Hydrogen Mass Fraction.csv',
        #            self.m_f[3], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Anode Water Mass Fraction.csv',
        #            self.m_f[4], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Cathode Activation Loss.csv',
        #            self.act_loss_cat, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Anode Activation Loss.csv',
        #            self.act_loss_ano, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Cathode Layer Diffusion Loss.csv',
        #            self.cl_diff_loss_cat, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Anode Layer Diffusion Loss.csv',
        #            self.cl_diff_loss_ano, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Cathode GDL Diffusion Loss.csv',
        #            self.gdl_diff_loss_cat, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Anode GDL Diffusion ´Loss.csv',
        #            self.gdl_diff_loss_ano, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Membrane Conductivity Loss.csv',
        #            self.mem_loss, delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Voltage Loss.csv', self.v_loss,
        #            delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Cell Voltage.csv', self.v_cell,
        #            delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Oxygen Heat Capacity.csv',
        #            self.cp[0], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Cathode Gas Water Heat Capacity.csv',
        #            self.cp[1], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Nitrogen Heat Capacity.csv',
        #            self.cp[2], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Hydrogen Heat Capacity.csv',
        #            self.cp[3], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Anode Gas Water Heat Capacity.csv',
        #            self.cp[4], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Oxygen Dynamic Viscosity.csv',
        #            self.visc[0], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Cathode Gas Water Dynamic Viscosity.csv',
        #            self.visc[1], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Nitrogen Dynamic Viscosity.csv',
        #            self.visc[2], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Hydrogen Dynamic Viscosity.csv',
        #            self.visc[3], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Anode Gas Water Dynamic Viscosity.csv',
        #            self.visc[4], delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path + 'Oxygen Thermal Conductivity.csv',
        #            self.lambda_gas[0], delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Cathode Gas Water Thermal Conductivity.csv',
        #            self.lambda_gas[1], delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Nitrogen Thermal Conductivity.csv',
        #            self.lambda_gas[2], delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Hydrogen Thermal Conductivity.csv',
        #            self.lambda_gas[3], delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Anode Water Gas Thermal Conductivity.csv',
        #            self.lambda_gas[4], delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Cathode Channel Mixture Heat Capacity.csv',
        #            self.cp_gas_cat, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Anode Channel Mixture Heat Capacity.csv',
        #            self.cp_gas_ano, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Cathode Channel Mixture Gas Constant.csv',
        #            self.r_gas_cat, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Anode Channel Mixture Gas Constant.csv',
        #            self.r_gas_ano, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Cathode Channel Mixture Dynamic Viscosity.csv',
        #            self.visc_gas_cat,
        #            delimiter=self.delimiter, fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Anode Channel Mixture Dynamic Viscosity.csv',
        #            self.visc_gas_ano, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Cathode Channel Mixture Heat Conductivity.csv',
        #            self.lambda_gas_cat, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Anode Channel Mixture Heat Conductivity.csv',
        #            self.lambda_gas_ano, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Cathode Channel Two Phase Heat Capacity.csv',
        #            self.cp_fluid_cat, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Anode Channel Two Phase Heat Capacity.csv',
        #            self.cp_fluid_ano, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Cathode Channel Gas Phase Density.csv',
        #            self.rho_gas_cat, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Anode Channel Gas Phase Density.csv',
        #            self.rho_gas_cat, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Cathode Channel Two Phase Mass Flow.csv',
        #            self.m_flow_fluid_cat, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Anode Channel Two Phase Mass Flow.csv',
        #            self.m_flow_fluid_ano, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path + 'Air Stoichiometry Distribution.csv',
        #            self.stoi_cat, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Hydrogen Stoichiometry Distribution.csv',
        #            self.stoi_ano, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Cathode Channel Heat Convection Coefficient.csv',
        #            self.ht_coef_cat, delimiter=self.delimiter,
        #            fmt=self.csv_format)
        # np.savetxt(path
        #            + 'Anode Channel Heat Convection Coefficient.csv',
        #            self.ht_coef_ano, delimiter=self.delimiter,
        #            fmt=self.csv_format)

    def plot_polarization_curve(self, voltage_loss,
                                cell_voltages, target_current_density):
        """
        Plots the polarization curve of the given
        current densities and average stack voltages.
        """
        cd_array = np.asarray(target_current_density) * 1.e-4
        plt.plot(cd_array, cell_voltages, marker='.', color='k',
                 label='Simulation')
        if self.show_loss is True:
            plt.plot(cd_array, voltage_loss['membrane']['average'],
                     color='b', marker='.', label='Membrane Loss')
            plt.plot(cd_array, voltage_loss['activation']['anode']['average'],
                     color='g', marker='*', label='Anode Activation Loss')
            plt.plot(cd_array,
                     voltage_loss['activation']['cathode']['average'],
                     color='g', marker='+', label='Cathode Activation Loss')
            plt.plot(cd_array,
                     voltage_loss['diffusion']['CL']['anode']['average'],
                     color='y', marker='*', label='Anode CL Diff Loss')
            plt.plot(cd_array,
                     voltage_loss['diffusion']['CL']['cathode']['average'],
                     color='y', marker='+', label='Cathode CL Diff Loss')
            plt.plot(cd_array,
                     voltage_loss['diffusion']['GDL']['anode']['average'],
                     color='m', marker='*', label='Anode GDL Diff Loss')
            plt.plot(cd_array,
                     voltage_loss['diffusion']['GDL']['cathode']['average'],
                     color='m', marker='+', label='Cathode GDL Diff Loss')
        plt.ylabel('Voltage $[V]$', fontsize=16)
        plt.xlabel('Current Density $[A/cm²]$', fontsize=16)
        plt.tick_params(labelsize=14)
        plt.grid()
        plt.legend()
        plt.autoscale(tight=True, axis='both', enable=True)
        plt.ylim(0., 1.)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'polarization_curve.png'))
        plt.close()
