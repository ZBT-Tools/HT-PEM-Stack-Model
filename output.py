import numpy as np
import matplotlib.pyplot as plt
import os
import errno
import system.interpolation as ip
import input.geometry as geom


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
        self.plot_path = None
        # path where the plots of the results gets saved
        self.path_csv_data = None
        # path where the csv data of the results gets saved

    def save(self, folder_name):
        if self.save_plot:
            self.output_plots(folder_name)
        if self.save_csv:
            self.output_csv(folder_name)

    @staticmethod
    def output(y_values, y_label, x_label, y_scale, color,
               title, xlim_low, xlim_up, val_label, path):
        if val_label:
            for l in range(len(y_values)):
                plt.plot(y_values[l], color=color[l],
                         marker='.', label=val_label[l])
        else:
            for l in range(len(y_values)):
                plt.plot(y_values[l], color=color[l], marker='.')

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.yscale(y_scale)
        plt.xlim(xlim_low, xlim_up)
        plt.tight_layout()
        plt.grid()
        if val_label:
            plt.legend()
        plt.savefig(os.path.join(path + title + '.png'))
        plt.close()

    @staticmethod
    def output_x(y_values, x_values, y_label, x_label,
                 y_scale, title, val_label, lim, path):
        if val_label:
            for l in range(len(y_values)):
                plt.plot(x_values, y_values[l],
                         color=plt.cm.coolwarm(l / len(y_values)),
                         marker='.', label=val_label[l])
        else:
            for l in range(len(y_values)):
                plt.plot(x_values, y_values[l],
                         color=plt.cm.coolwarm(l / len(y_values)), marker='.')

        plt.xlabel(x_label, fontsize=16)
        plt.ylabel(y_label, fontsize=16)
        plt.yscale(y_scale)
        plt.tick_params(labelsize=14)
        plt.autoscale(tight=True, axis='both', enable=True)
        plt.xlim(lim[0], lim[1])
        plt.tight_layout()
        plt.grid()
        if val_label:
            plt.legend()
        plt.savefig(os.path.join(path + title + '.png'))
        plt.close()

    def plot_cell_var(self, y_var, y_label, x_label, title, x_lim, x_var,
                      y_lim=False, y_scale='linear'):
        """
        Creates plots by given input values
        """
        for i, item in enumerate(self.stack.cells):
            plt.plot(x_var, eval('self.stack.cells' +
                                 '['+str(i)+']'+'.' + y_var),
                     color=plt.cm.coolwarm(i / self.stack.n_cells),
                     marker='.')
        plt.xlabel(x_label, fontsize=16)
        plt.ylabel(y_label, fontsize=16)
        plt.yscale(y_scale)
        plt.tick_params(labelsize=14)
        plt.autoscale(tight=True, axis='both', enable=True)
        plt.xlim(x_lim[0], x_lim[1])
        if y_lim:
            plt.ylim(y_lim[0], y_lim[1])
        plt.tight_layout()
        plt.grid()
        plt.savefig(self.plot_path + title + '.png')
        plt.close()

    def output_plots(self, folder_name):
        """
        Coordinates the plot sequence
        """
        plot_path = os.path.join(os.path.dirname(__file__), 'output',
                                 folder_name, 'plots')
        try:
            os.makedirs(plot_path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        length = self.stack.cells[0].cathode.channel.length
        x_lims = [0.0, length]
        x_node = self.stack.cells[0].cathode.channel.x
        x_ele = ip.interpolate_1d(x_node)
        self.output([self.mdf_criteria_process, self.i_ca_criteria_process,
                     self.temp_criteria_process], 'ERR', 'Iteration', 'log',
                    ['k', 'r', 'b'], 'Convergence', 0.,
                    len(self.temp_criteria_process),
                    ['Flow Distribution', 'Current Density', 'Temperature'],
                    self.plot_path)
        self.output_x(self.stack.i_cd, x_ele, 'Current Density $[A/m²]$',
                        'Channel Location $[m]$', 'linear', 'Current Density',
                        False, x_lims, plot_path)
        n_cells = self.stack.n_cells
        if n_cells > 1:
            self.output([self.stack.manifold[0].cell_stoi,
                         self.stack.manifold[1].cell_stoi],
                        'Stoichiometry', 'Cell Number', 'linear', ['k', 'r'],
                        'Stoichimetry Distribution', 0., n_cells - 1,
                        ['Cathode', 'Anode'], plot_path)
            self.output([self.stack.manifold[0].cell_stoi/2.5],
                        'Flow Distribution', 'Cell Number', 'linear', ['k'],
                        'Distribution', 0., n_cells - 1, ['Cathode'], plot_path)
        self.plot_cell_var('v', 'Voltage $[V]$', 'Channel Location $[m]$',
                           'Cell Voltage', x_lims, x_ele, [0.52, 0.54])
        self.output_x(self.stack.temp_sys.temp_cool, x_node,
                        'Coolant Temperature [K]', 'Channel Location $[m]$',
                        'linear', 'Coolant Temperature', False,
                        x_lims, plot_path)
        self.plot_cell_var('temp[-1]', 'Anode BPP - GDE Temperature $[K]$',
                           'Channel Location $[m]$',
                           'Anode Plate - GDE Temperature', x_lims, x_ele)
        self.plot_cell_var('temp[-2]', 'Anode GDE - MEM Temperature $[K]$',
                           'Channel Location $[m]$',
                           'Anode GDE - Membrane Temperature', x_lims, x_ele)
        self.plot_cell_var('temp[2]',
                           'Cathode GDE - MEM Temperature $[K]$',
                           'Channel Location $[m]$', 'Cathode GDL Temperature',
                           x_lims, x_ele)
        self.plot_cell_var('cathode.temp_fluid',
                           'Cathode Fluid Temperature $[K]$',
                           'Channel Location $[m]$',
                           'Cathode_Channel_Temperature', x_lims, x_node)
        self.plot_cell_var('temp[1]', 'Cathode BPP-GDE Temperature $[K]$',
                           'Channel Location $[m]$',
                           'Cathode GDE - Plate Temperature', x_lims, x_ele)
        self.plot_cell_var('anode.temp_fluid', 'Anode Fluid Temperature $[K]$',
                           'Channel Location $[m]$',
                           'Anode_Channel_Temperature', x_lims, x_node)
        self.plot_cell_var('temp[0]', 'BPP - BPP Temperature $[K]$',
                           'Channel Location $[m]$',
                           'Coolant Plate Temperature', x_lims, x_ele)
        self.plot_cell_var('cathode.mol_flow[0] * 1.e3',
                           'Cathode Oxygen Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$',
                           'Cathode Oxygen Molar Flow', x_lims, x_node)
        self.plot_cell_var('cathode.mol_flow[1] * 1.e3',
                           'Cathode Water Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$',
                           'Cathode Water Molar Flow', x_lims, x_node)
        self.plot_cell_var('cathode.mol_flow[2] * 1.e3',
                           'Cathode Nitrogen Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$',
                           'Cathode Nitrogen Molar Flow', x_lims, x_node)
        self.plot_cell_var('anode.mol_flow[0] * 1.e3',
                           'Anode Hydrogen Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$',
                           'Anode Hydrogen Molar Flow', x_lims, x_node)
        self.plot_cell_var('anode.mol_flow[1] * 1.e3',
                           'Anode Water Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$',
                           'Anode Water Molar Flow', x_lims, x_node)
        self.plot_cell_var('anode.mol_flow[2] * 1.e3',
                           'Anode Nitrogen Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$',
                           'Anode Nitrogen Molar Flow', x_lims, x_node)
        self.plot_cell_var('cathode.mol_f[0]', 'Oxygen Molar Fraction',
                           'Channel Location $[m]$',
                           'Oxygen Molar Fraction', x_lims, x_node)
        self.plot_cell_var('cathode.mol_f[1]',
                           'Cathode Gas Water Molar Fraction',
                           'Channel Location $[m]$',
                           'Water Molar Fraction Cathode',
                           x_lims, x_node)
        self.plot_cell_var('cathode.mol_f[2]',
                           'Cathode Nitrogen Molar Fraction',
                           'Channel Location $[m]$',
                           'Nitrogen Molar Fraction Cathode', x_lims, x_node)
        self.plot_cell_var('anode.mol_f[0]', 'Hydrogen Molar Fraction',
                           'Channel Location $[m]$',
                           'Hydrogen_Molar_Fraction_Anode', x_lims, x_node)
        self.plot_cell_var('anode.mol_f[1]', 'Anode Gas Water  Molar Fraction',
                           'Channel Location $[m]$',
                           'Water_Molar_Fraction_Anode', x_lims, x_node)
        self.plot_cell_var('anode.mol_f[2]', 'Anode Nitrogen Molar Fraction',
                           'Channel Location $[m]$',
                           'Nitrogen_Molar_Fraction_Anode', x_lims, x_node)
        self.plot_cell_var('cathode.liq_w_flow * 1.e3',
                           'Cathode Liquid Water Flow $[mmol/s]$',
                           'Channel Location $[m]$',
                           'Liquid Water Flow Cathode', x_lims, x_node)
        self.plot_cell_var('cathode.cond_rate * 1.e3',
                           'Cathode Water Condensation Rate $[mmol/s]$',
                           'Channel Location $[m]$',
                           'Water Condensation Rate Cathode', x_lims, x_node)
        self.plot_cell_var('cathode.humidity', 'Cathode Relative Humidity',
                           'Channel Location $[m]$', 'linear',
                           'Relative Humidity Cathode', x_lims, x_node)
        self.plot_cell_var('cathode.m_flow_gas * 1.e6',
                           'Cathode Channel Gas Massflow $[mg/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode_Channel__Gas_Massflow', x_lims, x_node)
        self.plot_cell_var('cathode.m_flow_fluid * 1.e6',
                           'Cathode Channel Fluid Massflow $[mg/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode_Channel_Fluid_Massflow', x_lims, x_node)
        self.plot_cell_var('cathode.g_fluid * 1.e3',
                           'Cathode Capacity Flow $[mW/K]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode Capacity Flow', x_lims, x_node)
        self.plot_cell_var('cathode.m_flow_reac * 1.e6',
                           'Oxygen Massflow $[mg/s]$', 'Channel Location $[m]$',
                           'linear', 'Oxygen_massflow', x_lims, x_node)
        self.plot_cell_var('cathode.m_flow_vap_w * 1.e6',
                           'Cathode Vapour Massflow $[mg/s]$',
                           'Channel Location $[m]$',
                           'linear', 'Vapour Massflow', x_lims, x_node)
        self.plot_cell_var('anode.m_flow_reac * 1.e6',
                           'Hydrogen Mass Flow $[mg/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Hydrogen Mass Flow', x_lims, x_node)
        self.plot_cell_var('cathode.cp_fluid',
                           'Cathode Heat Capacity $[J/(kgK)]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode Heat Capacity', x_lims, x_node)
        self.plot_cell_var('cathode.p', 'Cathode Channel Pressure $[Pa]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode Channel Pressure', x_lims, x_node)
        self.plot_cell_var('anode.p', 'Anode Channel Pressure $[Pa]$',
                           'Channel Location $[m]$', 'linear',
                           'Anode Channel Pressure', x_lims, x_node)
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
        for i in range(self.stack.n_cells):
            if i is 0:
                x.append(x_vec_z)
            elif 0 < i < self.stack.n_cells - 1:
                x.append(x_vec_e)
            else:
                x.append(x_vec_l)
        x = np.cumsum(np.block(x))
        t = self.stack.temp_sys.temp_layer
        for w in range(g_par.dict_case['nodes'] - 1):
            t_vec = []
            for i in range(self.stack.n_cells):
                if i is not self.stack.n_cells - 1:
                    t_vec.append(np.array([t[i][0, w], t[i][1, w],
                                           t[i][2, w], t[i][3, w], t[i][4, w]]))
                else:
                    t_vec.append(np.array([t[i][0, w], t[i][1, w], t[i][2, w],
                                           t[i][3, w], t[i][4, w], t[i][5, w]]))
            plt.plot(x, np.block(t_vec),
                     color=plt.cm.coolwarm((w + 1.e-20)
                                           / float(g_par.dict_case['nodes']
                                                   - 1.)))
        plt.xlim(0, x[-1])
        plt.xlabel('Stack Location $[m]$', fontsize=16)
        plt.ylabel('Temperature $[K]$', fontsize=16)
        plt.tick_params(labelsize=14)
        plt.autoscale(tight=True, axis='both', enable=True)
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_path + 'z-Cut-Temperature' + '.png'))
        plt.close()

        for i in range(self.stack.n_cells):
            print(np.average(self.stack.i_cd[i, :]))

    def output_csv(self, folder_name):
        self.path_csv_data = os.path.join(os.path.dirname(__file__),
                                          'output', folder_name, 'csv_data')
        try:
            os.makedirs(self.path_csv_data)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        for i, item in enumerate(self.stack.cells):
            self.mol_flow[0, i] = item.cathode.mol_flow[0]
            self.mol_flow[1, i] = item.cathode.mol_flow[1]
            self.mol_flow[2, i] = item.cathode.mol_flow[2]
            self.mol_flow[3, i] = item.anode.mol_flow[0]
            self.mol_flow[4, i] = item.anode.mol_flow[1]
            self.mol_flow[5, i] = item.anode.mol_flow[2]
            self.gas_con[0, i] = item.cathode.gas_con[0]
            self.gas_con[1, i] = item.cathode.gas_con[1]
            self.gas_con[2, i] = item.cathode.gas_con[2]
            self.gas_con[3, i] = item.anode.gas_con[0]
            self.gas_con[4, i] = item.anode.gas_con[1]
            self.gas_con[5, i] = item.anode.gas_con[2]
            self.m_f[0, i] = item.cathode.mass_f[0]
            self.m_f[1, i] = item.cathode.mass_f[1]
            self.m_f[2, i] = item.cathode.mass_f[2]
            self.m_f[3, i] = item.anode.mass_f[0]
            self.m_f[4, i] = item.anode.mass_f[1]
            self.m_f[5, i] = item.anode.mass_f[2]
            self.mol_f[0, i] = item.cathode.mol_f[0]
            self.mol_f[1, i] = item.cathode.mol_f[1]
            self.mol_f[2, i] = item.cathode.mol_f[2]
            self.mol_f[3, i] = item.anode.mol_f[0]
            self.mol_f[4, i] = item.anode.mol_f[1]
            self.mol_f[5, i] = item.anode.mol_f[2]
            self.v_loss[i] = item.v_loss
            self.v_cell[i] = item.v
            self.cp[0, i] = item.cathode.cp[0]
            self.cp[1, i] = item.cathode.cp[1]
            self.cp[2, i] = item.cathode.cp[2]
            self.cp[3, i] = item.anode.cp[0]
            self.cp[4, i] = item.anode.cp[1]
            self.cp[5, i] = item.anode.cp[2]
            self.lambda_gas[0, i] = item.cathode.lambda_gas[0]
            self.lambda_gas[1, i] = item.cathode.lambda_gas[1]
            self.lambda_gas[2, i] = item.cathode.lambda_gas[2]
            self.lambda_gas[3, i] = item.anode.lambda_gas[0]
            self.lambda_gas[4, i] = item.anode.lambda_gas[1]
            self.lambda_gas[5, i] = item.anode.lambda_gas[2]
            self.visc[0, i] = item.cathode.visc[0]
            self.visc[1, i] = item.cathode.visc[1]
            self.visc[2, i] = item.cathode.visc[2]
            self.visc[3, i] = item.anode.visc[0]
            self.visc[4, i] = item.anode.visc[1]
            self.visc[5, i] = item.anode.visc[2]
            self.r_gas_cat[i] = item.cathode.r_gas
            self.r_gas_ano[i] = item.anode.r_gas
            self.cp_gas_cat[i] = item.cathode.cp_gas
            self.cp_gas_ano[i] = item.anode.cp_gas
            self.visc_gas_cat[i] = item.cathode.visc_gas
            self.visc_gas_ano[i] = item.anode.visc_gas
            self.lambda_gas_cat[i] = item.cathode.lambda_gas
            self.lambda_gas_ano[i] = item.anode.lambda_gas
            self.rho_gas_cat[i] = item.cathode.rho_gas
            self.rho_gas_cat[i] = item.anode.rho_gas
            self.u_gas_cat[i] = item.cathode.u
            self.u_gas_ano[i] = item.anode.u
            self.p_cat[i] = item.cathode.p
            self.p_ano[i] = item.anode.p
            self.ht_coef_cat[i] = item.cathode.ht_coef
            self.ht_coef_ano[i] = item.anode.ht_coef
            self.cp_fluid_cat[i] = item.cathode.cp_fluid
            self.cp_fluid_ano[i] = item.anode.cp_fluid
            self.m_flow_fluid_cat[i] = item.cathode.m_flow_fluid
            self.m_flow_fluid_ano[i] = item.anode.m_flow_fluid
            self.temp_fluid_cat[i] = item.cathode.temp_fluid
            self.temp_fluid_ano[i] = item.anode.temp_fluid
            self.stoi_cat[i] = item.cathode.stoi
            self.stoi_ano[i] = item.anode.stoi
            for folder_name in range(5):
                self.temp_layer.append(
                    self.stack.temp_cpl_stack.temp_layer[i][folder_name, :])
        np.savetxt(self.path_csv_data + 'Temperature Layer.csv',
                   self.temp_layer, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Coolant Temperature.csv',
                   self.stack.temp_cpl_stack.temp_cool,
                   delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Current Density.csv', self.stack.i_ca,
                   delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Gas Temperature.csv',
                   self.temp_fluid_cat, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Gas Temperature.csv',
                   self.temp_fluid_ano, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Channel Average Velocity.csv',
                   self.u_gas_cat, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Channel Average Velocity.csv',
                   self.u_gas_ano, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Channel Pressure.csv',
                   self.p_cat, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Channel Pressure.csv',
                   self.p_ano, delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Oxygen Molar Flow.csv',
                   self.mol_flow[0], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Water Molar Flow.csv',
                   self.mol_flow[1], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Nitrogen Molar Flow.csv',
                   self.mol_flow[2], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Hydrogen Molar Flow.csv',
                   self.mol_flow[3], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Water Molar Flow.csv',
                   self.mol_flow[4], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Oxygen Molar Concentration.csv',
                   self.gas_con[0], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Water Molar Concentration.csv',
                   self.gas_con[1], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Nitrogen Molar Concentration.csv',
                   self.gas_con[2], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Hydrogen Molar Concentration.csv',
                   self.gas_con[3], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Water Molar Concentration.csv',
                   self.gas_con[4], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Oxygen Molar Fraction.csv',
                   self.mol_f[0], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Water Molar Fraction.csv',
                   self.mol_f[1], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Nitrogen Molar Fraction.csv',
                   self.mol_f[2], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Hydrogen Molar Fraction.csv',
                   self.mol_f[3], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Water Molar Fraction.csv',
                   self.mol_f[4], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Oxygen Mass Fraction.csv',
                   self.mol_f[0], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Mass Fraction.csv',
                   self.m_f[1], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Nitrogen Mass Fraction.csv',
                   self.m_f[2], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Hydrogen Mass Fraction.csv',
                   self.m_f[3], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Water Mass Fraction.csv',
                   self.m_f[4], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Activation Loss.csv',
                   self.act_loss_cat, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Activation Loss.csv',
                   self.act_loss_ano, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Layer Diffusion Loss.csv',
                   self.cl_diff_loss_cat, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Layer Diffusion Loss.csv',
                   self.cl_diff_loss_ano, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode GDL Diffusion Loss.csv',
                   self.gdl_diff_loss_cat, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode GDL Diffusion ´Loss.csv',
                   self.gdl_diff_loss_ano, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Membrane Conductivity Loss.csv',
                   self.mem_loss, delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Voltage Loss.csv', self.v_loss,
                   delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cell Voltage.csv', self.v_cell,
                   delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Oxygen Heat Capacity.csv',
                   self.cp[0], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Gas Water Heat Capacity.csv',
                   self.cp[1], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Nitrogen Heat Capacity.csv',
                   self.cp[2], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Hydrogen Heat Capacity.csv',
                   self.cp[3], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Gas Water Heat Capacity.csv',
                   self.cp[4], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Oxygen Dynamic Viscosity.csv',
                   self.visc[0], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Gas Water Dynamic Viscosity.csv',
                   self.visc[1], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Nitrogen Dynamic Viscosity.csv',
                   self.visc[2], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Hydrogen Dynamic Viscosity.csv',
                   self.visc[3], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Gas Water Dynamic Viscosity.csv',
                   self.visc[4], delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Oxygen Thermal Conductivity.csv',
                   self.lambda_gas[0], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Gas Water Thermal Conductivity.csv',
                   self.lambda_gas[1], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Nitrogen Thermal Conductivity.csv',
                   self.lambda_gas[2], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Hydrogen Thermal Conductivity.csv',
                   self.lambda_gas[3], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Anode Water Gas Thermal Conductivity.csv',
                   self.lambda_gas[4], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Channel Mixture Heat Capacity.csv',
                   self.cp_gas_cat, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Anode Channel Mixture Heat Capacity.csv',
                   self.cp_gas_ano, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Channel Mixture Gas Constant.csv',
                   self.r_gas_cat, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Anode Channel Mixture Gas Constant.csv',
                   self.r_gas_ano, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Channel Mixture Dynamic Viscosity.csv',
                   self.visc_gas_cat,
                   delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Anode Channel Mixture Dynamic Viscosity.csv',
                   self.visc_gas_ano, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Channel Mixture Heat Conductivity.csv',
                   self.lambda_gas_cat, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Anode Channel Mixture Heat Conductivity.csv',
                   self.lambda_gas_ano, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Channel Two Phase Heat Capacity.csv',
                   self.cp_fluid_cat, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Anode Channel Two Phase Heat Capacity.csv',
                   self.cp_fluid_ano, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Channel Gas Phase Density.csv',
                   self.rho_gas_cat, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Channel Gas Phase Density.csv',
                   self.rho_gas_cat, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Channel Two Phase Mass Flow.csv',
                   self.m_flow_fluid_cat, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Channel Two Phase Mass Flow.csv',
                   self.m_flow_fluid_ano, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Air Stoichiometry Distribution.csv',
                   self.stoi_cat, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Hydrogen Stoichiometry Distribution.csv',
                   self.stoi_ano, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Channel Heat Convection Coefficient.csv',
                   self.ht_coef_cat, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Anode Channel Heat Convection Coefficient.csv',
                   self.ht_coef_ano, delimiter=self.delimiter,
                   fmt=self.csv_format)

    def plot_polarization_curve(self, voltage_loss,
                                cell_voltages, target_current_density):
        """
        Plots the polarization curve of the given
        current densities and average stack voltages.
        """
        try:
            os.makedirs(os.path.join(os.path.dirname(__file__), 'output/'))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
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
        plt.savefig(os.path.join(os.path.dirname(__file__),
                                 'output/' + 'Polarization_curve' + '.png'))
        plt.close()
