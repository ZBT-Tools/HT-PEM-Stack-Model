import data.stack_dict as st_dict
import data.channel_dict as ch_dict
import data.simulation_dict as sim
import input.operating_conditions as oper_con
import system.stack as st
import numpy as np
import data.global_parameter as gpar
import system.global_functions as gfunc
import input.geometry as geo
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
        self.k_it = dict_simulation['k_it']
        self.max_it = dict_simulation['max_it']
        self.delimiter = ','
        self.csv_format = '%.9e'
        # Variables
        self.stack = None
        self.path_plot = None
        self.path_csv_data = None
        self.t_old = None
        cell_num = st_dict.stack['cell_numb']
        nodes = gpar.dict_case['nodes']
        # Arrays
        self.mdf_criteria_cat_process = []
        self.mdf_criteria_ano_process = []
        self.i_criteria_process = []
        self.t_criteria_process = []
        self.mdf_criteria_process = []
        self.t1_criteria_process = []
        self.t2_criteria_process = []
        self.t3_criteria_process = []
        self.t4_criteria_process = []
        self.t5_criteria_process = []
        self.i_criteria = []
        self.t_criteria = []
        self.tryarray = []
        self.v = []
        self.mol_flow = np.full((5, cell_num, nodes), 0.)
        self.gas_con = np.full((5, cell_num, nodes), 0.)
        self.m_f = np.full((5, cell_num, nodes), 0.)
        self.mol_f = np.full((5, cell_num, nodes), 0.)
        self.act_ov_cathode = np.full((cell_num, nodes - 1), 0.)
        self.act_ov_anode = np.full((cell_num, nodes - 1), 0.)
        self.cat_dif_los_cathode = np.full((cell_num, nodes - 1), 0.)
        self.cat_dif_los_anode = np.full((cell_num, nodes - 1), 0.)
        self.gdl_dif_los_cathode = np.full((cell_num, nodes - 1), 0.)
        self.gdl_dif_los_anode = np.full((cell_num, nodes - 1), 0.)
        self.ohm_los = np.full((cell_num, nodes - 1), 0.)
        self.v_los = np.full((cell_num, nodes - 1), 0.)
        self.v_save = np.full((cell_num, nodes - 1), 0.)
        self.cp = np.full((5, cell_num, nodes), 0.)
        self.lambda_x = np.full((5, cell_num, nodes), 0.)
        self.visc = np.full((5, cell_num, nodes), 0.)
        self.r_mix_cathode = np.full((cell_num, nodes), 0.)
        self.r_mix_anode = np.full((cell_num, nodes), 0.)
        self.cp_mix_cathode = np.full((cell_num, nodes), 0.)
        self.cp_mix_anode = np.full((cell_num, nodes), 0.)
        self.visc_mix_cathode = np.full((cell_num, nodes), 0.)
        self.visc_mix_anode = np.full((cell_num, nodes), 0.)
        self.lama_mix_cathode = np.full((cell_num, nodes), 0.)
        self.lama_mix_anode = np.full((cell_num, nodes), 0.)
        self.rho_cathode = np.full((cell_num, nodes), 0.)
        self.rho_anode = np.full((cell_num, nodes), 0.)
        self.pr_cathode = np.full((cell_num, nodes), 0.)
        self.pr_anode = np.full((cell_num, nodes), 0.)
        self.u_cathode = np.full((cell_num, nodes), 0.)
        self.u_anode = np.full((cell_num, nodes), 0.)
        self.re_cathode = np.full((cell_num, nodes), 0.)
        self.re_anode = np.full((cell_num, nodes), 0.)
        self.p_cathode = np.full((cell_num, nodes), 0.)
        self.p_anode = np.full((cell_num, nodes), 0.)
        self.ht_coef_cathode = np.full((cell_num, nodes), 0.)
        self.ht_coef_anode = np.full((cell_num, nodes), 0.)
        self.m_flow_cathode = np.full((cell_num, nodes), 0.)
        self.m_flow_anode = np.full((cell_num, nodes), 0.)
        self.cp_full_cathode = np.full((cell_num, nodes), 0.)
        self.cp_full_anode = np.full((cell_num, nodes), 0.)
        self.t_layer = []
        self.t_gas_anode = np.full((cell_num, nodes), 0.)
        self.t_gas_cathode = np.full((cell_num, nodes), 0.)
        self.t_coolant = np.full((cell_num, nodes), 0.)
        self.stoi_cathode = np.full(cell_num, 0.)
        self.stoi_anode = np.full(cell_num, 0.)

    # @do_c_profile
    def update(self):
        for q, item in enumerate(oper_con.target_current_density):
            print(q)
            gpar.dict_case['tar_cd'] = oper_con.target_current_density[q]
            self.stack = st.Stack(st_dict.stack)
            statement = True
            counter = 0
            while statement is True:
                self.save_old_value()
                self.tryarray.append(self.stack.i_ca[0, -1])
                self.stack.update()
                if self.stack.break_program is True:
                    break
                self.calc_convergence_criteria()
                if len(oper_con.target_current_density) < 1:
                    print(counter)
                counter = counter + 1
                if ((self.i_criteria < self.k_it
                    and self.t_criteria < self.k_it) and counter > 10)\
                        or counter > self.max_it:
                    statement = False
            if self.stack.break_program is False:
                self.mdf_criteria_process =\
                    (np.array(self.mdf_criteria_ano_process)
                     + np.array(self.mdf_criteria_cat_process)) * .5
                self.v.append(np.average(self.stack.v))
                self.output_plots(str(q))
                self.output_csv(str(q))
                self.tryarray = []
            else:
                oper_con.target_current_density =\
                    oper_con.target_current_density[0:-q]
                print(oper_con.target_current_density, self.v)
                break
        if len(oper_con.target_current_density) > 1:
            # comp_i = np.array([1111.11,3333.33,4444.44,5555.55,6666.66])
            # comp_v = np.array([0.675,0.582,0.5465,0.51325,0.48125])
            plt.plot(np.asarray(oper_con.target_current_density) * 1.e-4,
                     self.v, marker='.', color='k', label='Simulation')
            # plt.plot(comp_i*1e-4, comp_v, marker='^',
            #  color='r', label='Measurement')
            plt.ylabel('Voltage $[V]$', fontsize=16)
            plt.xlabel('Current Density $[A/cm²]$', fontsize=16)
            plt.tick_params(labelsize=14)
            plt.grid()
            plt.legend()
            plt.autoscale(tight=True, axis='both', enable=True)
            plt.ylim(0., 1.)
            plt.tight_layout()
            plt.savefig(os.path.join(os.path.dirname(__file__),
                                     'Polarization_curve' + '.jpg'))
            plt.close()

    def calc_convergence_criteria(self):
        self.i_criteria = np.abs(sum(((self.stack.i_ca.flatten()
                                       - self.stack.i_old.flatten())
                                      / self.stack.i_ca.flatten()) ** 2.))
        self.t_criteria =\
            np.abs(np.sum(((self.t_old
                            - self.stack.temp_cpl_stack.t_layer[0][0, 0]))
                          / self.stack.temp_cpl_stack.t_layer[0][0, 0]))
        self.t_criteria_process.append(self.t_criteria)
        self.mdf_criteria_cat_process.append(self.stack.cathode_mfd_criteria)
        self.mdf_criteria_ano_process.append(self.stack.anode_mfd_criteria)
        self.i_criteria_process.append(self.i_criteria)

    def save_old_value(self):
        self.t_old = self.stack.temp_cpl_stack.t_layer[0][0, 0]

    def plot_cell_var(self, y_var, y_label, x_label,
                      y_scale, title, xlim, x_var, y_lim):
        for l, item in enumerate(self.stack.cells):
            plt.plot(x_var, eval('self.stack.cells' +
                                 '['+str(l)+']'+'.' + y_var),
                     color=plt.cm.coolwarm(l / self.stack.cell_num),
                     marker='.')
        plt.xlabel(x_label, fontsize=16)
        plt.ylabel(y_label, fontsize=16)
        plt.yscale(y_scale)
        plt.tick_params(labelsize=14)
        plt.autoscale(tight=True, axis='both', enable=True)
        plt.xlim(xlim[0], xlim[1])
        if y_lim is not False:
            plt.ylim(y_lim[0], y_lim[1])
        plt.tight_layout()
        plt.grid()
        plt.savefig(self.path_plot + title + '.jpg')
        plt.close()

    def output_plots(self, q):
        self.path_plot = os.path.join(os.path.dirname(__file__),
                                      'output/' + 'case' + q + '/plots' + '/')
        try:
            os.makedirs(self.path_plot)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        x_node = np.linspace(0., ch_dict.cathode_channel['length'],
                             gpar.dict_case['nodes'])
        x_ele = gfunc.calc_elements_1_d(x_node)
        gfunc.output([self.mdf_criteria_process, self.i_criteria_process,
                      self.t_criteria_process], 'ERR', 'Iteration', 'log',
                     ['k', 'r', 'b'], 'Convergence', 0.,
                     len(self.t_criteria_process),
                     ['Flow Distribution', 'Current Density', 'Temperature'],
                     self.path_plot)
        self.mdf_criteria_process = []
        self.mdf_criteria_ano_process = []
        self.mdf_criteria_cat_process = []
        self.t_criteria_process = []
        self.i_criteria_process = []
        gfunc.output_x(self.stack.i_ca, x_ele, 'Current Density $[A/m²]$',
                       'Channel Location $[m]$', 'linear', 'Current Density',
                       False, [0., ch_dict.cathode_channel['length']],
                       self.path_plot)
        if self.stack.cell_num > 1:
            gfunc.output([self.stack.manifold[0].cell_stoi,
                          self.stack.manifold[1].cell_stoi],
                         'Stoichiometry', 'Cell Number', 'linear', ['k', 'r'],
                         'Stoichimetry Distribution', 0., self.stack.cell_num-1,
                         ['Cathode', 'Anode'], self.path_plot)
        self.plot_cell_var('v', 'Voltage $[V]$', 'Channel Location $[m]$',
                           'linear', 'Cell Voltage',
                           [0., ch_dict.cathode_channel['length']], x_ele,
                           [0., 1.28])
        gfunc.output_x(self.stack.temp_cpl_stack.t_cool, x_node,
                       'Coolant Temperature [K]', 'Channel Location $[m]$',
                       'linear', 'Coolant Temperature', False,
                       [0., ch_dict.cathode_channel['length']], self.path_plot)
        self.plot_cell_var('t[-1]',
                           'Anode Plate - GDE Temperature $[K]$',
                           'Channel Location $[m]$', 'linear',
                           'Anode Plate - GDE Temperature',
                           [0., ch_dict.cathode_channel['length']], x_ele,
                           False)
        self.plot_cell_var('t[-2]', 'Anode GDE - Membrane Temperature $[K]$',
                           'Channel Location $[m]$', 'linear',
                           'Anode GDE - Membrane Temperature',
                           [0., ch_dict.cathode_channel['length']], x_ele,
                           False)
        self.plot_cell_var('t[2]', 'Membrane - Cathode GDE Temperature $[K]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode GDL Temperature',
                           [0., ch_dict.cathode_channel['length']], x_ele,
                           False)
        self.plot_cell_var('cathode.t_gas', 'Cathode Channel Temperature $[K]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode_Channel_Temperature',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('t[1]', 'Cathode GDE - Plate Temperature $[K]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode GDE - Plate Temperature',
                           [0., ch_dict.cathode_channel['length']], x_ele,
                           False)
        self.plot_cell_var('anode.t_gas', 'Hydrogen Temperature $[K]$',
                           'Channel Location $[m]$', 'linear',
                           'Hydrogen Temperature',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('t[0]',
                           'Cathode Plate- Anode Plate Temperature $[K]$',
                           'Channel Location $[m]$', 'linear',
                           'Coolant Plate Temperature',
                           [0., ch_dict.cathode_channel['length']], x_ele,
                           False)
        self.plot_cell_var('cathode.mol_flow[0] * 1.e3',
                           'Oxygen Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Oxygen Molar Flow',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('cathode.mol_flow[1] * 1.e3',
                           'Water Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Water Molar Flow Cathode',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('cathode.mol_flow[2] * 1.e3',
                           'Nitrogen Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Nitrogen Molar Flow',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('anode.mol_flow[0] * 1.e3',
                           'Hydrogen Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Hydrogen Molar Flow',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('anode.mol_flow[1] * 1.e3',
                           'Water Molar Flow $[mmol/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Water Molar Flow Anode',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('cathode.mol_f[0]', 'Oxygen  Molar Fraction',
                           'Channel Location $[m]$', 'linear',
                           'Oxygen_Molar_Fraction',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('cathode.mol_f[1]', 'Gas Water  Molar Fraction',
                           'Channel Location $[m]$', 'linear',
                           'Water Molar Fraction Cathode',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('anode.mol_f[0]', 'Hydrogen Molar Fraction',
                           'Channel Location $[m]$', 'linear',
                           'Hydrogen_Molar_Fraction_Anode',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('anode.mol_f[1]', 'Gas Water  Molar Fraction',
                           'Channel Location $[m]$', 'linear',
                           'Water_Molar_Fraction_Anode',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('cathode.w * 1.e3', 'Liquid Water Flow $[mmol/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Liquid Water Flow Cathode',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('cathode.gamma * 1.e3',
                           'Water Condensation Rate $[mmol/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Water Condensation Rate Cathode',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('cathode.humidity', 'Relative Humidity',
                           'Channel Location $[m]$', 'linear',
                           'Relative Humidity Cathode',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('cathode.m_flow * 1.e6',
                           'Cathode Channel Gas Massflow $[mg/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode_Channel__Gas_Massflow',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('cathode.m_full_flow * 1.e6',
                           'Cathode Channel Massflow $[mg/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode_Channel_Full_Massflow',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('cathode.g_full * 1.e3',
                           'Cathode Capacity Flow $[mW/K]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode Capacity Flow',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('cathode.m_reac_flow * 1.e6',
                           'Oxygen Massflow $[mg/s]$', 'Channel Location $[m]$',
                           'linear', 'Oxygen_massflow',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('cathode.m_vap_water_flow * 1.e6',
                           'Vapour Massflow $[mg/s]$', 'Channel Location $[m]$',
                           'linear', 'Vapour Massflow',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('anode.m_flow * 1.e6', 'Hydrogen Massflow $[mg/s]$',
                           'Channel Location $[m]$', 'linear',
                           'Hydrogen_massflow',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('cathode.cp_full',
                           'Cathode Heat Capacity $[J/(kgK)]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode Heat Capacity',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('cathode.p', 'Cathode Channel Pressure $[Pa]$',
                           'Channel Location $[m]$', 'linear',
                           'Cathode Channel Pressure',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        self.plot_cell_var('anode.p', 'Anode Channel Pressure $[Pa]$',
                           'Channel Location $[m]$', 'linear',
                           'Anode Channel Pressure',
                           [0., ch_dict.cathode_channel['length']], x_node,
                           False)
        # Z-Axis-Temperature Plot
        x_vecz = np.array([0., geo.plate_thickness, geo.gdl_thickness,
                           geo.membrane_thickness, geo.gdl_thickness])
        x_vec_e = np.array([geo.plate_thickness, geo.plate_thickness,
                            geo.gdl_thickness, geo.membrane_thickness,
                            geo.gdl_thickness])
        x_vec_l = np.array([geo.plate_thickness, geo.plate_thickness,
                            geo.gdl_thickness, geo.membrane_thickness,
                            geo.gdl_thickness, geo.plate_thickness])
        x = []
        for l in range(self.stack.cell_num):
            if l is 0:
                x.append(x_vecz)
            elif 0 < l < self.stack.cell_num - 1:
                x.append(x_vec_e)
            else:
                x.append(x_vec_l)
        x = np.cumsum(np.block(x))
        t = self.stack.temp_cpl_stack.t_layer
        for w in range(gpar.dict_case['nodes']-1):
            t_vec = []
            for l in range(self.stack.cell_num):
                if l is not self.stack.cell_num-1:
                    t_vec.append(np.array([t[l][0, w], t[l][1, w],
                                           t[l][2, w], t[l][3, w], t[l][4, w]]))
                else:
                    t_vec.append(np.array([t[l][0, w], t[l][1, w], t[l][2, w],
                                           t[l][3, w], t[l][4, w], t[l][5, w]]))
            plt.plot(x, np.block(t_vec),
                     color=plt.cm.coolwarm((w + 1.e-20)/self.stack.cell_num))
        plt.xlim(0, x[-1])
        plt.xlabel('Stack Location $[m]$', fontsize=16)
        plt.ylabel('Temperature $[K]$', fontsize=16)
        plt.tick_params(labelsize=14)
        plt.autoscale(tight=True, axis='both', enable=True)
        plt.tight_layout()
        plt.savefig(os.path.join(self.path_plot+'Z-Cut-Temperature_'+ '.jpg'))
        plt.close()

        for q in range(self.stack.cell_num):
            print(np.average(self.stack.i_ca[q, :]))

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
            self.gas_con[0, w] = item.cathode.gas_con[0]
            self.gas_con[1, w] = item.cathode.gas_con[1]
            self.gas_con[2, w] = item.cathode.gas_con[2]
            self.gas_con[3, w] = item.anode.gas_con[0]
            self.gas_con[4, w] = item.anode.gas_con[1]
            self.m_f[0, w] = item.cathode.mass_f[0]
            self.m_f[1, w] = item.cathode.mass_f[1]
            self.m_f[2, w] = item.cathode.mass_f[2]
            self.m_f[3, w] = item.anode.mass_f[0]
            self.m_f[4, w] = item.anode.mass_f[1]
            self.mol_f[0, w] = item.cathode.mol_f[0]
            self.mol_f[1, w] = item.cathode.mol_f[1]
            self.mol_f[2, w] = item.cathode.mol_f[2]
            self.mol_f[3, w] = item.anode.mol_f[0]
            self.mol_f[4, w] = item.anode.mol_f[1]
            self.act_ov_cathode[w] = item.cathode.act_ov
            self.act_ov_anode[w] = item.anode.act_ov
            self.cat_dif_los_cathode[w] = item.cathode.cat_diff_los
            self.cat_dif_los_anode[w] = item.anode.cat_diff_los
            self.gdl_dif_los_cathode[w] = item.cathode.gdl_diff_los
            self.gdl_dif_los_anode[w] = item.anode.gdl_diff_los
            self.ohm_los[w] = item.mem_ov
            self.v_los[w] = item.v_los
            self.v_save[w] = item.v
            self.cp[0, w] = item.cathode.cp[0]
            self.cp[1, w] = item.cathode.cp[1]
            self.cp[2, w] = item.cathode.cp[2]
            self.cp[3, w] = item.anode.cp[0]
            self.cp[4, w] = item.anode.cp[1]
            self.lambda_x[0, w] = item.cathode.lambda_gas[0]
            self.lambda_x[1, w] = item.cathode.lambda_gas[1]
            self.lambda_x[2, w] = item.cathode.lambda_gas[2]
            self.lambda_x[3, w] = item.anode.lambda_gas[0]
            self.lambda_x[4, w] = item.anode.lambda_gas[1]
            self.visc[0, w] = item.cathode.visc[0]
            self.visc[1, w] = item.cathode.visc[1]
            self.visc[2, w] = item.cathode.visc[2]
            self.visc[3, w] = item.anode.visc[0]
            self.visc[4, w] = item.anode.visc[1]
            self.r_mix_cathode[w] = item.cathode.r_mix
            self.r_mix_anode[w] = item.anode.r_mix
            self.cp_mix_cathode[w] = item.cathode.cp_mix
            self.cp_mix_anode[w] = item.anode.cp_mix
            self.visc_mix_cathode[w] = item.cathode.visc_mix
            self.visc_mix_anode[w] = item.anode.visc_mix
            self.lama_mix_cathode[w] = item.cathode.lambda_mix
            self.lama_mix_anode[w] = item.anode.lambda_mix
            self.rho_cathode[w] = item.cathode.rho
            self.rho_anode[w] = item.anode.rho
            self.u_cathode[w] = item.cathode.u
            self.u_anode[w] = item.anode.u
            self.p_cathode[w] = item.cathode.p
            self.p_anode[w] = item.anode.p
            self.ht_coef_cathode[w] = item.cathode.ht_coef
            self.ht_coef_anode[w] = item.anode.ht_coef
            self.cp_full_cathode[w] = item.cathode.cp_full
            self.cp_full_anode[w] = item.anode.cp_full
            self.m_flow_cathode[w] = item.cathode.m_full_flow
            self.m_flow_anode[w] = item.anode.m_full_flow
            self.t_gas_cathode[w] = item.cathode.t_gas
            self.t_gas_anode[w] = item.anode.t_gas
            self.stoi_cathode[w] = item.cathode.stoi
            self.stoi_anode[w] = item.anode.stoi
            for q in range(5):
                self.t_layer.append(self.stack.temp_cpl_stack.t_layer[w][q, :])
        np.savetxt(self.path_csv_data + 'Temperature Layer.csv',
                   self.t_layer, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Coolant Temperature.csv',
                   self.t_coolant, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Current Density.csv', self.stack.i_ca,
                   delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Gas Temperature.csv',
                   self.t_gas_cathode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Gas Temperature.csv',
                   self.t_gas_anode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Channel Average Velocity.csv',
                   self.u_cathode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Channel Average Velocity.csv',
                   self.u_anode, delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Channel Pressure.csv',
                   self.p_cathode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Channel Pressure.csv',
                   self.p_anode, delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Oxygen Molar Flow.csv',
                   self.mol_flow[0], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Water Molar Flow.csv',
                   self.mol_flow[1], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Nitrogen Molar Flow.csv',
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
                   self.act_ov_cathode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Activation Loss.csv',
                   self.act_ov_anode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Layer Diffusion Loss.csv',
                   self.cat_dif_los_cathode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Layer Diffusion Loss.csv',
                   self.cat_dif_los_anode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode GDL Diffusion Loss.csv',
                   self.gdl_dif_los_cathode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode GDL Diffusion ´Loss.csv',
                   self.gdl_dif_los_anode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Membrane Conductivity Loss.csv',
                   self.ohm_los, delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Voltage Loss.csv', self.v_los,
                   delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cell Voltage.csv', self.v_save,
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
                   self.lambda_x[0], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Gas Water Thermal Conductivity.csv',
                   self.lambda_x[1], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Nitrogen Thermal Conductivity.csv',
                   self.lambda_x[2], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Hydrogen Thermal Conductivity.csv',
                   self.lambda_x[3], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Anode Water Gas Thermal Conductivity.csv',
                   self.lambda_x[4], delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Channel Mixture Heat Capacity.csv',
                   self.cp_mix_cathode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Anode Channel Mixture Heat Capacity.csv',
                   self.cp_mix_anode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Channel Mixture Gas Constant.csv',
                   self.r_mix_cathode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Anode Channel Mixture Gas Constant.csv',
                   self.r_mix_anode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Channel Mixture Dynamic Viscosity.csv',
                   self.visc_mix_cathode,
                   delimiter=self.delimiter, fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Anode Channel Mixture Dynamic Viscosity.csv',
                   self.visc_mix_anode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Channel Mixture Heat Conductivity.csv',
                   self.lama_mix_cathode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Anode Channel Mixture Heat Conductivity.csv',
                   self.lama_mix_anode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Channel Two Phase Heat Capacity.csv',
                   self.cp_full_cathode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Anode Channel Two Phase Heat Capacity.csv',
                   self.cp_full_anode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Cathode Channel Gas Phase Density.csv',
                   self.rho_cathode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Channel Gas Phase Density.csv',
                   self.rho_anode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Channel Two Phase Mass Flow.csv',
                   self.m_flow_cathode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Anode Channel Two Phase Mass Flow.csv',
                   self.m_flow_anode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data + 'Air Stoichiometry Distribution.csv',
                   self.stoi_cathode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Hydrogen Stoichiometry Distribution.csv',
                   self.stoi_anode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Cathode Channel Heat Convection Coefficient.csv',
                   self.ht_coef_cathode, delimiter=self.delimiter,
                   fmt=self.csv_format)
        np.savetxt(self.path_csv_data
                   + 'Anode Channel Heat Convection Coefficient .csv',
                   self.ht_coef_anode, delimiter=self.delimiter,
                   fmt=self.csv_format)


start = timeit.default_timer()
Simulation_runs = Simulation(sim.simulation)
Simulation_runs.update()
stop = timeit.default_timer()
print('Simulation time:', stop-start)
