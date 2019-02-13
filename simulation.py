import data.stack_dict as st_dict
import data.channel_dict as ch_dict
import data.simulation_dict as sim
import input.operating_conditions as op_con
import system.stack as st
import numpy as np
import data.global_parameters as g_par
import system.global_functions as g_func
import input.geometry as geom
import cProfile
import matplotlib.pyplot as plt
import os
import errno
import timeit
np.set_printoptions(threshold=np.nan, linewidth=10000,
                    precision=9, suppress=True)


def do_c_profile(func):
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.print_stats('cumtime')
    return profiled_func


class Simulation:

    def __init__(self, dict_simulation):
        # Handover
        self.it_crit = dict_simulation['iteration_criteria']
        # iteration criteria
        self.max_it = dict_simulation['maximal_iteration']
        # maximal number of iterations before force termination#
        self.save_csv = dict_simulation['save_csv']
        # switch to save the csv data
        self.save_plot = dict_simulation['save_plot']
        # switch to save the plot data
        self.show_loss = dict_simulation['show_loss']
        # switch to show the single voltage losses in the u-i-graph
        cell_numb = st_dict.dict_stack['cell_numb']
        # number of stack cells
        nodes = g_par.dict_case['nodes']
        # node points of the x-grid

        """General variables"""
        self.delimiter = ','
        self.csv_format = '%.9e'
        self.stack = None
        # object of the class Stack
        self.path_plot = None
        # path where the plots of the results gets saved
        self.path_csv_data = None
        # path where the csv data of the results gets saved
        self.temp_old = None
        # defined temperature of the last iteration
        self.mdf_criteria_cat_process = []
        # array of the cathodic mdf criteria over the iterations
        self.mdf_criteria_ano_process = []
        # array of the anodic mdf criteria over the iterations
        self.i_ca_criteria_process = []
        # array of the current density criteria over the iterations
        self.temp_criteria_process = []
        # array of the temperature criteria over the iterations
        self.mdf_criteria_process = []
        # common array of the mdf criteria over the iterations
        self.i_ca_criteria = None
        # convergence criteria of the current density
        self.temp_criteria = None
        # convergence criteria of the temperature
        self.v = []
        # cell voltage
        self.mol_flow = np.full((6, cell_numb, nodes), 0.)
        # molar flow of the species in the channels
        # 0: oxygen, 1: cathode water, 2: cathode nitrogen,
        # 3: hydrogen, 4: anode water, 5: anode nitrogen
        self.gas_con = np.full((6, cell_numb, nodes), 0.)
        # molar concentration of the species in the gas mixture
        self.m_f = np.full((6, cell_numb, nodes), 0.)
        # mass fraction of the species in the gas mixture
        self.mol_f = np.full((6, cell_numb, nodes), 0.)
        # molar fraction of the species in the gas mixture
        self.act_loss_cat = np.full((cell_numb, nodes - 1), 0.)
        # cathodic activation voltage loss
        self.act_loss_ano = np.full((cell_numb, nodes - 1), 0.)
        # anodic activation voltage loss
        self.cl_diff_loss_cat = np.full((cell_numb, nodes - 1), 0.)
        # cathodic catalyst layer diffusion voltage loss
        self.cl_diff_loss_ano = np.full((cell_numb, nodes - 1), 0.)
        # anodic catalyst layer diffusion voltage loss
        self.gdl_diff_loss_cat = np.full((cell_numb, nodes - 1), 0.)
        # cathodic gas diffusion layer diffusion voltage loss
        self.gdl_diff_loss_ano = np.full((cell_numb, nodes - 1), 0.)
        # anodic gas diffusion layer diffusion voltage loss
        self.mem_loss = np.full((cell_numb, nodes - 1), 0.)
        # membrane voltage loss
        self.v_loss = np.full((cell_numb, nodes - 1), 0.)
        # voltage loss over the stack
        self.v_cell = np.full((cell_numb, nodes - 1), 0.)
        # cell voltage
        self.cp = np.full((6, cell_numb, nodes), 0.)
        # heat capacity of the species in the gas phase
        self.lambda_gas = np.full((6, cell_numb, nodes), 0.)
        # heat conductivity of the species in the gas phase
        self.visc = np.full((6, cell_numb, nodes), 0.)
        # viscosity of the species in the gas phase
        self.r_gas_cat = np.full((cell_numb, nodes), 0.)
        # gas constant of the gas phase in the cathode channels
        self.r_gas_ano = np.full((cell_numb, nodes), 0.)
        # gas constant of the gas phase in the anode channels
        self.cp_gas_cat = np.full((cell_numb, nodes), 0.)
        # heat capacitiy of the gas phase in the cathode channels
        self.cp_gas_ano = np.full((cell_numb, nodes), 0.)
        # heat capacity of the gas phase in the anode channels
        self.visc_gas_cat = np.full((cell_numb, nodes), 0.)
        # viscosity of the gas phase in the cathode channels
        self.visc_gas_ano = np.full((cell_numb, nodes), 0.)
        # viscosity of the gas phase in the anode channels
        self.lambda_gas_cat = np.full((cell_numb, nodes), 0.)
        # heat conductivity of the gas phase in the cathode channels
        self.lambda_gas_ano = np.full((cell_numb, nodes), 0.)
        # heat conductivity of the gas phase in the anode channels
        self.rho_gas_cat = np.full((cell_numb, nodes), 0.)
        # density of the gas in the cathode channels
        self.rho_gas_cat = np.full((cell_numb, nodes), 0.)
        # density of the gas in the anode channels
        self.pr_gas_cat = np.full((cell_numb, nodes), 0.)
        # prandtl number of the gas in the cathode channels
        self.pr_gas_ano = np.full((cell_numb, nodes), 0.)
        # prandtl number of the gas in the anode channels
        self.u_gas_cat = np.full((cell_numb, nodes), 0.)
        # velocity of the fluid in the cathode channels
        self.u_gas_ano = np.full((cell_numb, nodes), 0.)
        # velocity of the fluid in the anode channels
        self.re_gas_cat = np.full((cell_numb, nodes), 0.)
        # reynolds number of the fluid in the cathode channels
        self.re_gas_ano = np.full((cell_numb, nodes), 0.)
        # reynolds number of the fluid in the anode channels
        self.p_cat = np.full((cell_numb, nodes), 0.)
        # pressure in the cathode channels
        self.p_ano = np.full((cell_numb, nodes), 0.)
        # pressure in the anode channels
        self.ht_coef_cat = np.full((cell_numb, nodes), 0.)
        # heat convection coefficient in the cathode channels
        self.ht_coef_ano = np.full((cell_numb, nodes), 0.)
        # heat convection coefficient in the anode channels
        self.m_flow_fluid_cat = np.full((cell_numb, nodes), 0.)
        # mass flow of the fluid in the cathode channels
        self.m_flow_fluid_ano = np.full((cell_numb, nodes), 0.)
        # mass flow of the fluid in the anode channels
        self.cp_fluid_cat = np.full((cell_numb, nodes), 0.)
        # heat capacity of the cathode fluid
        self.cp_fluid_ano = np.full((cell_numb, nodes), 0.)
        # heat capacity of the anode fluid
        self.temp_layer = []
        # temperature layer of the stack
        self.temp_fluid_ano = np.full((cell_numb, nodes), 0.)
        # temperature of the fluid in the anode channels
        self.temp_fluid_cat = np.full((cell_numb, nodes), 0.)
        # temperature of the fluid in the cathode channels
        self.temp_cool = np.full((cell_numb, nodes), 0.)
        # temperature of the coolant
        self.stoi_cat = np.full(cell_numb, 0.)
        # inlet stoichiometry of the cathode channels
        self.stoi_ano = np.full(cell_numb, 0.)
        # inlet stoichiometry of the anode channels
        self.act_loss_ui_ano = []
        # average activation voltage loss of the anode
        self.act_loss_ui_cat = []
        # average activation voltage loss of the cathode
        self.cl_diff_loss_ui_ano = []
        # average anode catalyst layer diffusion voltage losses
        self.cl_diff_loss_ui_cat = []
        # average cathode catalyst layer diffusion voltage losses
        self.gdl_diff_loss_ui_ano = []
        # average anode gdl diffusion voltage losses
        self.gdl_diff_loss_ui_cat = []
        # average cathode gdl diffusion voltage losses
        self.mem_loss_ui = []
        # average membrane voltage losses

    # @do_c_profile
    def update(self):
        """
        This function coordinates the program sequence
        """
        for i, item in enumerate(op_con.target_current_density):
            g_par.dict_case['tar_cd'] = op_con.target_current_density[i]
            self.stack = st.Stack(st_dict.dict_stack)
            statement = True
            counter = 0
            while statement is True:
                self.save_old_value()
                self.stack.update()
                if self.stack.break_program is True:
                    break
                self.calc_convergence_criteria()
                if len(op_con.target_current_density) < 1:
                    print(counter)
                counter = counter + 1
                if ((self.i_ca_criteria < self.it_crit
                     and self.temp_criteria < self.it_crit) and counter > 10)\
                        or counter > self.max_it:
                    statement = False
            if self.stack.break_program is False:
                self.mdf_criteria_process =\
                    (np.array(self.mdf_criteria_ano_process)
                     + np.array(self.mdf_criteria_cat_process)) * .5
                self.save_voltages()
                print(item)
                if self.save_plot is True:
                    self.output_plots(str(i))
                if self.save_csv is True:
                    self.output_csv(str(i))
            else:
                op_con.target_current_density = \
                    op_con.target_current_density[0:-i]
                print(op_con.target_current_density, self.v)
                break
        if len(op_con.target_current_density) > 1:
            self.plot_polarization_curve()

    def plot_polarization_curve(self):
        """
        Plots the polarization curve of the given
        current densities and average stack voltages.
        """
        try:
            os.makedirs(os.path.join(os.path.dirname(__file__), 'output/'))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        cd_array = np.asarray(op_con.target_current_density) * 1.e-4
        plt.plot(cd_array, self.v, marker='.', color='k', label='Simulation')
        if self.show_loss is True:
            plt.plot(cd_array, self.mem_loss_ui, color='b', marker='.',
                     label='Membrane Loss')
            plt.plot(cd_array, self.act_loss_ui_ano, color='g', marker='*',
                     label='Anode Activation Loss')
            plt.plot(cd_array, self.act_loss_ui_cat, color='g', marker='+',
                     label='Cathode Activation Loss')
            plt.plot(cd_array, self.cl_diff_loss_ui_ano, color='y', marker='*',
                     label='Anode Cl Diff Loss')
            plt.plot(cd_array, self.cl_diff_loss_ui_cat, color='y', marker='+',
                     label='Cathode Cl Diff Loss')
            plt.plot(cd_array, self.gdl_diff_loss_ui_ano, color='m', marker='*',
                     label='Anode GDL Diff Loss')
            plt.plot(cd_array, self.gdl_diff_loss_ui_cat, color='m', marker='+',
                     label='Cathode GDL Diff Loss')
        plt.ylabel('Voltage $[V]$', fontsize=16)
        plt.xlabel('Current Density $[A/cm²]$', fontsize=16)
        plt.tick_params(labelsize=14)
        plt.grid()
        plt.legend()
        plt.autoscale(tight=True, axis='both', enable=True)
        plt.ylim(0., 1.)
        plt.tight_layout()
        plt.savefig(os.path.join(os.path.dirname(__file__),
                                 'output/' + 'Polarization_curve' + '.jpg'))
        plt.close()

    def save_voltages(self):
        """
        Saves the average voltage losses of the stack
        """
        for w, item in enumerate(self.stack.cells):
            self.act_loss_cat[w] = item.cathode.act_loss
            self.act_loss_ano[w] = item.anode.act_loss
            self.cl_diff_loss_cat[w] = item.cathode.cl_diff_loss
            self.cl_diff_loss_ano[w] = item.anode.cl_diff_loss
            self.gdl_diff_loss_cat[w] = item.cathode.gdl_diff_loss
            self.gdl_diff_loss_ano[w] = item.anode.gdl_diff_loss
            self.mem_loss[w] = item.mem_loss
        self.v.append(np.average(self.stack.v_cell))
        self.act_loss_ui_ano.append(np.average(self.act_loss_ano))
        self.act_loss_ui_cat.append(np.average(self.act_loss_cat))
        self.cl_diff_loss_ui_ano.append(np.average(self.cl_diff_loss_ano))
        self.cl_diff_loss_ui_cat.append(np.average(self.cl_diff_loss_cat))
        self.gdl_diff_loss_ui_ano.append(np.average(self.gdl_diff_loss_ano))
        self.gdl_diff_loss_ui_cat.append(np.average(self.gdl_diff_loss_cat))
        self.mem_loss_ui.append(np.average(self.mem_loss))

    def calc_convergence_criteria(self):
        """
        Calculates the convergence criteria according to (Koh, 2003)
        """
        self.i_ca_criteria = np.abs(sum(((self.stack.i_cd.flatten()
                                          - self.stack.i_cd_old.flatten())
                                         / self.stack.i_cd.flatten()) ** 2.))
        self.temp_criteria =\
            np.abs(np.sum(((self.temp_old
                            - self.stack.temp_cpl_stack.temp_layer[0][0, 0]))
                          / self.stack.temp_cpl_stack.temp_layer[0][0, 0]))

        self.temp_criteria_process.append(self.temp_criteria)
        self.mdf_criteria_cat_process.append(self.stack.cathode_mfd_criteria)
        self.mdf_criteria_ano_process.append(self.stack.anode_mfd_criteria)
        self.i_ca_criteria_process.append(self.i_ca_criteria)

    def save_old_value(self):
        """
        Saves an defined temperature value of the current iteration
        as the old temperature value for the next iteration.
        """
        self.temp_old = self.stack.temp_cpl_stack.temp_layer[0][0, 0]

    def plot_cell_var(self, y_var, y_label, x_label,
                      y_scale, title, x_lim, x_var, y_lim):
        """
        Creates plots by given input values
        """
        for l, item in enumerate(self.stack.cells):
            plt.plot(x_var, eval('self.stack.cells' +
                                 '['+str(l)+']'+'.' + y_var),
                     color=plt.cm.coolwarm(l / self.stack.cell_numb),
                     marker='.')
        plt.xlabel(x_label, fontsize=16)
        plt.ylabel(y_label, fontsize=16)
        plt.yscale(y_scale)
        plt.tick_params(labelsize=14)
        plt.autoscale(tight=True, axis='both', enable=True)
        plt.xlim(x_lim[0], x_lim[1])
        if y_lim is not False:
            plt.ylim(y_lim[0], y_lim[1])
        plt.tight_layout()
        plt.grid()
        plt.savefig(self.path_plot + title + '.png')
        plt.close()

    def output_plots(self, q):
        """
        Coordinates the plot sequence
        """
        self.path_plot = os.path.join(os.path.dirname(__file__),
                                      'output/' + 'case' + q + '/plots' + '/')
        try:
            os.makedirs(self.path_plot)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        x_node = np.linspace(0., ch_dict.dict_cathode_channel['channel_length'],
                             g_par.dict_case['nodes'])
        x_ele = g_func.calc_elements_1_d(x_node)
        g_func.output([self.mdf_criteria_process, self.i_ca_criteria_process,
                       self.temp_criteria_process], 'ERR', 'Iteration', 'log',
                      ['k', 'r', 'b'], 'Convergence', 0.,
                      len(self.temp_criteria_process),
                      ['Flow Distribution', 'Current Density', 'Temperature'],
                      self.path_plot)
        self.mdf_criteria_process = []
        self.mdf_criteria_ano_process = []
        self.mdf_criteria_cat_process = []
        self.temp_criteria_process = []
        self.i_ca_criteria_process = []
        g_func.output_x(self.stack.i_cd, x_ele, 'Current Density $[A/m²]$',
                        'Channel Location $[m]$', 'linear', 'Current Density',
                        False,
                        [0., ch_dict.dict_cathode_channel['channel_length']],
                        self.path_plot)
        if self.stack.cell_numb > 1:
            g_func.output([self.stack.manifold[0].cell_stoi,
                           self.stack.manifold[1].cell_stoi],
                          'Stoichiometry', 'Cell Number', 'linear', ['k', 'r'],
                          'Stoichimetry Distribution', 0.,
                          self.stack.cell_numb - 1, ['Cathode', 'Anode'],
                          self.path_plot)
            g_func.output([self.stack.manifold[0].cell_stoi/2.5],
                          'Flow Distribution', 'Cell Number', 'linear', ['k'],
                          'Distribution', 0.,
                          self.stack.cell_numb - 1, ['Cathode'],
                          self.path_plot)
        self.plot_cell_var('v', 'Voltage $[V]$', 'Channel Location $[m]$',
                           'linear', 'Cell Voltage',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_ele, [0.52, 0.54])
        g_func.output_x(self.stack.temp_cpl_stack.temp_cool, x_node,
                        'Coolant Temperature [K]', 'Channel Location $[m]$',
                        'linear', 'Coolant Temperature', False,
                        [0., ch_dict.dict_cathode_channel['channel_length']],
                        self.path_plot)
        self.plot_cell_var('temp[-1]',
                           'Anode BPP - GDE Temperature $[K]$',
                           'Channel Location $[m]$', 'linear',
                           'Anode Plate - GDE Temperature',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_ele, False)
        self.plot_cell_var('temp[-2]', 'Anode GDE - MEM Temperature $[K]$',
                           'Channel Location $[m]$', 'linear',
                           'Anode GDE - Membrane Temperature',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_ele, False)
        self.plot_cell_var('temp[2]',
                           'Cathode GDE - MEM Temperature $[K]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode GDL Temperature',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_ele, False)
        self.plot_cell_var('cathode.temp_fluid',
                           'Cathode Fluid Temperature $[K]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode_Channel_Temperature',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('temp[1]', 'Cathode BPP-GDE Temperature $[K]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode GDE - Plate Temperature',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_ele, False)
        self.plot_cell_var('anode.temp_fluid', 'Anode Fluid Temperature $[K]$',
                           'Channel Location $[m]$', 'linear',
                           'Anode_Channel_Temperature',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('temp[0]',
                           'BPP - BPP Temperature $[K]$',
                           'Channel Location $[m]$', 'linear',
                           'Coolant Plate Temperature',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_ele, False)
        self.plot_cell_var('cathode.mol_flow[0] * 1.e3',
                           'Cathode Oxygen Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode Oxygen Molar Flow',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('cathode.mol_flow[1] * 1.e3',
                           'Cathode Water Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode Water Molar Flow',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('cathode.mol_flow[2] * 1.e3',
                           'Cathode Nitrogen Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode Nitrogen Molar Flow',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('anode.mol_flow[0] * 1.e3',
                           'Anode Hydrogen Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Anode Hydrogen Molar Flow',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('anode.mol_flow[1] * 1.e3',
                           'Anode Water Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Anode Water Molar Flow',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('anode.mol_flow[2] * 1.e3',
                           'Anode Nitrogen Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Anode Nitrogen Molar Flow',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('cathode.mol_f[0]', 'Oxygen  Molar Fraction',
                           'Channel Location $[m]$', 'linear',
                           'Oxygen_Molar_Fraction',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('cathode.mol_f[1]',
                           'Cathode Gas Water Molar Fraction',
                           'Channel Location $[m]$', 'linear',
                           'Water Molar Fraction Cathode',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('cathode.mol_f[2]',
                           'Cathode Nitrogen Molar Fraction',
                           'Channel Location $[m]$', 'linear',
                           'Nitrogen_Molar_Fraction_Cathode',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('anode.mol_f[0]', 'Hydrogen Molar Fraction',
                           'Channel Location $[m]$', 'linear',
                           'Hydrogen_Molar_Fraction_Anode',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('anode.mol_f[1]', 'Anode Gas Water  Molar Fraction',
                           'Channel Location $[m]$', 'linear',
                           'Water_Molar_Fraction_Anode',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('anode.mol_f[2]', 'Anode Nitrogen Molar Fraction',
                           'Channel Location $[m]$', 'linear',
                           'Nitrogen_Molar_Fraction_Anode',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('cathode.liq_w_flow * 1.e3',
                           'Cathode Liquid Water Flow $[mmol/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Liquid Water Flow Cathode',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('cathode.cond_rate * 1.e3',
                           'Cathode Water Condensation Rate $[mmol/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Water Condensation Rate Cathode',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('cathode.humidity', 'Cathode Relative Humidity',
                           'Channel Location $[m]$', 'linear',
                           'Relative Humidity Cathode',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('cathode.m_flow_gas * 1.e6',
                           'Cathode Channel Gas Massflow $[mg/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode_Channel__Gas_Massflow',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('cathode.m_flow_fluid * 1.e6',
                           'Cathode Channel Fluid Massflow $[mg/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode_Channel_Fluid_Massflow',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('cathode.g_fluid * 1.e3',
                           'Cathode Capacity Flow $[mW/K]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode Capacity Flow',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('cathode.m_flow_reac * 1.e6',
                           'Oxygen Massflow $[mg/s]$', 'Channel Location $[m]$',
                           'linear', 'Oxygen_massflow',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('cathode.m_flow_vap_w * 1.e6',
                           'Cathode Vapour Massflow $[mg/s]$',
                           'Channel Location $[m]$',
                           'linear', 'Vapour Massflow',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('anode.m_flow_reac * 1.e6',
                           'Hydrogen Massflow $[mg/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Hydrogen_massflow',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('cathode.cp_fluid',
                           'Cathode Heat Capacity $[J/(kgK)]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode Heat Capacity',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('cathode.p', 'Cathode Channel Pressure $[Pa]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode Channel Pressure',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
        self.plot_cell_var('anode.p', 'Anode Channel Pressure $[Pa]$',
                           'Channel Location $[m]$', 'linear',
                           'Anode Channel Pressure',
                           [0., ch_dict.dict_cathode_channel['channel_length']],
                           x_node, False)
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
        for l in range(self.stack.cell_numb):
            if l is 0:
                x.append(x_vec_z)
            elif 0 < l < self.stack.cell_numb - 1:
                x.append(x_vec_e)
            else:
                x.append(x_vec_l)
        x = np.cumsum(np.block(x))
        t = self.stack.temp_cpl_stack.temp_layer
        for w in range(g_par.dict_case['nodes'] - 1):
            t_vec = []
            for l in range(self.stack.cell_numb):
                if l is not self.stack.cell_numb - 1:
                    t_vec.append(np.array([t[l][0, w], t[l][1, w],
                                           t[l][2, w], t[l][3, w], t[l][4, w]]))
                else:
                    t_vec.append(np.array([t[l][0, w], t[l][1, w], t[l][2, w],
                                           t[l][3, w], t[l][4, w], t[l][5, w]]))
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
        plt.savefig(os.path.join(self.path_plot+'Z-Cut-Temperature' + '.png'))
        plt.close()

        for q in range(self.stack.cell_numb):
            print(np.average(self.stack.i_cd[q, :]))

    def output_csv(self, q):
        self.path_csv_data = os.path.join(os.path.dirname(__file__),
                                          'output/' + 'case' + q
                                          + '/csv_data' + '/')
        try:
            os.makedirs(self.path_csv_data)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        for w, item in enumerate(self.stack.cells):
            self.mol_flow[0, w] = item.cathode.mol_flow[0]
            self.mol_flow[1, w] = item.cathode.mol_flow[1]
            self.mol_flow[2, w] = item.cathode.mol_flow[2]
            self.mol_flow[3, w] = item.anode.mol_flow[0]
            self.mol_flow[4, w] = item.anode.mol_flow[1]
            self.mol_flow[5, w] = item.anode.mol_flow[2]
            self.gas_con[0, w] = item.cathode.gas_con[0]
            self.gas_con[1, w] = item.cathode.gas_con[1]
            self.gas_con[2, w] = item.cathode.gas_con[2]
            self.gas_con[3, w] = item.anode.gas_con[0]
            self.gas_con[4, w] = item.anode.gas_con[1]
            self.gas_con[5, w] = item.anode.gas_con[2]
            self.m_f[0, w] = item.cathode.mass_f[0]
            self.m_f[1, w] = item.cathode.mass_f[1]
            self.m_f[2, w] = item.cathode.mass_f[2]
            self.m_f[3, w] = item.anode.mass_f[0]
            self.m_f[4, w] = item.anode.mass_f[1]
            self.m_f[5, w] = item.anode.mass_f[2]
            self.mol_f[0, w] = item.cathode.mol_f[0]
            self.mol_f[1, w] = item.cathode.mol_f[1]
            self.mol_f[2, w] = item.cathode.mol_f[2]
            self.mol_f[3, w] = item.anode.mol_f[0]
            self.mol_f[4, w] = item.anode.mol_f[1]
            self.mol_f[5, w] = item.anode.mol_f[2]
            self.v_loss[w] = item.v_loss
            self.v_cell[w] = item.v
            self.cp[0, w] = item.cathode.cp[0]
            self.cp[1, w] = item.cathode.cp[1]
            self.cp[2, w] = item.cathode.cp[2]
            self.cp[3, w] = item.anode.cp[0]
            self.cp[4, w] = item.anode.cp[1]
            self.cp[5, w] = item.anode.cp[2]
            self.lambda_gas[0, w] = item.cathode.lambda_gas[0]
            self.lambda_gas[1, w] = item.cathode.lambda_gas[1]
            self.lambda_gas[2, w] = item.cathode.lambda_gas[2]
            self.lambda_gas[3, w] = item.anode.lambda_gas[0]
            self.lambda_gas[4, w] = item.anode.lambda_gas[1]
            self.lambda_gas[5, w] = item.anode.lambda_gas[2]
            self.visc[0, w] = item.cathode.visc[0]
            self.visc[1, w] = item.cathode.visc[1]
            self.visc[2, w] = item.cathode.visc[2]
            self.visc[3, w] = item.anode.visc[0]
            self.visc[4, w] = item.anode.visc[1]
            self.visc[5, w] = item.anode.visc[2]
            self.r_gas_cat[w] = item.cathode.r_gas
            self.r_gas_ano[w] = item.anode.r_gas
            self.cp_gas_cat[w] = item.cathode.cp_gas
            self.cp_gas_ano[w] = item.anode.cp_gas
            self.visc_gas_cat[w] = item.cathode.visc_gas
            self.visc_gas_ano[w] = item.anode.visc_gas
            self.lambda_gas_cat[w] = item.cathode.lambda_gas
            self.lambda_gas_ano[w] = item.anode.lambda_gas
            self.rho_gas_cat[w] = item.cathode.rho_gas
            self.rho_gas_cat[w] = item.anode.rho_gas
            self.u_gas_cat[w] = item.cathode.u
            self.u_gas_ano[w] = item.anode.u
            self.p_cat[w] = item.cathode.p
            self.p_ano[w] = item.anode.p
            self.ht_coef_cat[w] = item.cathode.ht_coef
            self.ht_coef_ano[w] = item.anode.ht_coef
            self.cp_fluid_cat[w] = item.cathode.cp_fluid
            self.cp_fluid_ano[w] = item.anode.cp_fluid
            self.m_flow_fluid_cat[w] = item.cathode.m_flow_fluid
            self.m_flow_fluid_ano[w] = item.anode.m_flow_fluid
            self.temp_fluid_cat[w] = item.cathode.temp_fluid
            self.temp_fluid_ano[w] = item.anode.temp_fluid
            self.stoi_cat[w] = item.cathode.stoi
            self.stoi_ano[w] = item.anode.stoi
            for q in range(5):
                self.temp_layer.append(
                    self.stack.temp_cpl_stack.temp_layer[w][q, :])
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


start = timeit.default_timer()
Simulation_runs = Simulation(sim.simulation)
Simulation_runs.update()
stop = timeit.default_timer()
print('Simulation time:', stop-start)
