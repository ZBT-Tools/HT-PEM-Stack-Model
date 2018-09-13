import input_parameter as i_p
import stack as st
import numpy as np
import copy
import global_parameter as gpar
import global_functions as gfunc
import cProfile
import matplotlib.pyplot as plt
import os, errno
np.set_printoptions(threshold=np.nan, linewidth=1000, precision=9, suppress=True)

def do_cprofile(func):
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

    def __init__(self, dict):
        self.k_it = dict['k_it']
        self.max_it = dict['max_it']
        self.mdf_criteria_cat_process = []
        self.mdf_criteria_ano_process = []
        self.i_criteria_process = []
        self.t1_criteria_process = []
        self.t2_criteria_process = []
        self.t3_criteria_process = []
        self.t4_criteria_process = []
        self.t5_criteria_process = []
        self.i_last_criteria_process = []
        self.tryarray = []
        self.v = []

   # @do_cprofile
    def update(self):
        for q, item in enumerate(i_p.tar_cd):
            gpar.dict_case['tar_cd'] = i_p.tar_cd[q]
            self.stack = st.Stack(i_p.stack)
            statement = True
            counter = 0
            while statement is True:
                self.save_old_value()
                self.tryarray.append(self.stack.i[0, -1])
                self.stack.update()
                if self.stack.break_programm is True:
                    break
                self.calc_convergence_criteria()
                print(counter)
                counter = counter + 1
                if ((self.i_criteria < self.k_it and counter > 100)
                    and (self.mdf_criteria_ano < self.k_it
                         and self.mfd_criteria_cat < self.k_it))\
                        or counter > self.max_it:
                    statement = False
            self.t_criteria_process = (np.array(self.t1_criteria_process)
                                       + np.array(self.t2_criteria_process)
                                       + np.array(self.t3_criteria_process)
                                       + np.array(self.t4_criteria_process)
                                       + np.array(self.t5_criteria_process)) * .2
            self.mdf_criteria_process = (np.array(self.mdf_criteria_ano_process)
                                         + np.array(self.mdf_criteria_cat_process)) * .5
            self.output(str(q))
            #self.v.append(np.average(self.stack.v))
            self.tryarray = []
        #plt.plot(i_p.tar_cd,self.v)
        #plt.ylabel('Voltage [V]')
        #plt.xlabel('Current Density [A/m²]')
        #plt.savefig(os.path.join(os.path.dirname(__file__),   'polarization curve' + '.jpg'))
        #plt.show()



    def calc_convergence_criteria(self):
        self.mfd_criteria_cat = np.abs(sum(((self.stack.q_x_cat - self.q_x_cat_old)
                                            / self.stack.q_x_cat)**2.))
        self.mdf_criteria_ano = np.abs(sum(((self.stack.q_x_ano - self.q_x_ano_old)
                                            / self.stack.q_x_ano)**2.))
        self.i_criteria = np.abs(sum(((self.stack.i.flatten()
                                       - self.stack.i_old.flatten())
                                      / self.stack.i.flatten())**2.))
        self.i_last_criteria = np.abs((self.stack.i[0, -1] - self.i_last)
                                      / self.stack.i[0, -1])
        self.t1_criteria = np.abs(sum(((self.t1_old - self.stack.cell_list[0].t1)
                                            / self.stack.cell_list[0].t1)**2.))
        self.t2_criteria = np.abs(sum(((self.t2_old - self.stack.cell_list[0].t2)
                                       / self.stack.cell_list[0].t2) ** 2.))
        self.t3_criteria = np.abs(sum(((self.t3_old - self.stack.cell_list[0].t3)
                                       / self.stack.cell_list[0].t3) ** 2.))
        self.t4_criteria = np.abs(sum(((self.t4_old - self.stack.cell_list[0].t4)
                                       / self.stack.cell_list[0].t4) ** 2.))
        self.t5_criteria = np.abs(sum(((self.t5_old - self.stack.cell_list[0].t5)
                                       / self.stack.cell_list[0].t5) ** 2.))
        self.t1_criteria_process.append(self.t1_criteria)
        self.t2_criteria_process.append(self.t2_criteria)
        self.t3_criteria_process.append(self.t3_criteria)
        self.t4_criteria_process.append(self.t4_criteria)
        self.t5_criteria_process.append(self.t5_criteria)
        self.mdf_criteria_cat_process.append(self.mfd_criteria_cat)
        self.mdf_criteria_ano_process.append(self.mdf_criteria_ano)
        self.i_criteria_process.append(self.i_criteria)
        self.i_last_criteria_process.append(self.i_last_criteria)

    def save_old_value(self):
        self.t1_old = copy.deepcopy(self.stack.cell_list[0].t1)
        self.t2_old = copy.deepcopy(self.stack.cell_list[0].t2)
        self.t3_old = copy.deepcopy(self.stack.cell_list[0].t3)
        self.t4_old = copy.deepcopy(self.stack.cell_list[0].t4)
        self.t5_old = copy.deepcopy(self.stack.cell_list[0].t5)
        self.q_x_cat_old = copy.deepcopy(self.stack.q_x_cat)
        self.q_x_ano_old = copy.deepcopy(self.stack.q_x_ano)
        self.i_last = copy.deepcopy(self.stack.i[0, -1])

    def plot_cell_var(self, y_var, y_label, x_label,
                      y_scale, color, title, q, xlim, x_var, y_lim):
        try:
            os.makedirs(os.path.join(os.path.dirname(__file__), 'Plots' + q + '/'))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        for l, item in enumerate(self.stack.cell_list):
            plt.plot(x_var, eval('self.stack.cell_list'+
                                 '['+str(l)+']'+'.' + y_var), color=plt.cm.coolwarm(l/(self.stack.cell_numb)), marker='.')

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
        plt.savefig(os.path.join(os.path.dirname(__file__),
                                 'Plots' + q + '/' + title + '.jpg'))
        plt.close()


    def output(self, q):
        x_node = np.linspace(0., i_p.channel_length, gpar.dict_case['nodes'])
        x_ele = gfunc.calc_elements(x_node)
        gfunc.output([self.mdf_criteria_process, self.i_criteria_process, self.t_criteria_process],
                     'ERR', 'Iteration', 'log', ['k', 'r', 'b'], 'Convergence',
                     q, 0., len(self.t_criteria_process),
                     ['Flow Distribution','Current Density','Temperature'])
        gfunc.output([self.tryarray], 'Current Density $[A/m²]$', 'Iteration',
                     'linear', 'k', 'Current_Density_Last_Cell', q, 0.,
                     len(self.tryarray), False)
        gfunc.output_x(self.stack.i, x_ele, 'Current Density $[A/m²]$', 'Channel Location $[m]$',
                     'linear', np.full(self.stack.cell_numb, 'k'),
                     'Current Density', q, False,[0., i_p.channel_length])
        if self.stack.cell_numb >1:
            gfunc.output([self.stack.q_x_cat / (self.stack.q_h_in_cat[-1] / self.stack.cell_numb),
                          self.stack.q_x_ano / (self.stack.q_h_in_ano[-1] / self.stack.cell_numb)],
                         'Flow Distribution', 'Cell Number', 'linear', ['k', 'r'],
                         'Flow Distribution', q , 0., self.stack.cell_numb-1,
                         ['Cathode', 'Anode'])
            gfunc.output([self.stack.q_x_cat / (self.stack.q_h_in_cat[-1]
                          / self.stack.cell_numb) * self.stack.stoi_cat,
                          self.stack.q_x_ano / (self.stack.q_h_in_ano[-1]
                          / self.stack.cell_numb) * self.stack.stoi_ano],
                         'Stoichiometry', 'Cell Number', 'linear', ['k', 'r'],
                         'Stoichimetry Distribution', q, 0., self.stack.cell_numb-1,
                         ['Cathode', 'Anode'])
            #np.savetxt('cat_flow.csv', np.flipud(self.stack.q_x_cat / (self.stack.q_h_in_cat[-1] / self.stack.cell_numb)))

        self.plot_cell_var('v', 'Voltage $[V]$', 'Channel Location $[m]$', 'linear',
                           'k', 'Cell Voltage', q, [0., i_p.channel_length], x_ele, [0., 1.28])
        self.plot_cell_var('dv', 'dV/dI $[V/A]$', 'Channel Location $[m]$', 'linear', 'k',
                           'dvdI', q, [0., i_p.channel_length], x_ele, [-1., 1.])
        self.plot_cell_var('t1', 'Cathode Temperature $[K]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Cathode Temperature', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('t4', 'Anode Temperature $[K]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Anode Temperature', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('t2', 'Cathode GDL Temperature $[K]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Cathode GDL Temperature', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.t_gas', 'Cathode Channel Temperature $[K]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Cathode_Channel_Temperature', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('t5', 'Anode GDL Temperature $[K]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Anode GDL Temperature', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('anode.t_gas', 'Hydrogen Temperature $[K]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Hydrogen Temperature', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('t1', 'Coolant Plate Temperature $[K]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Coolant Plate Temperature', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.gas_flow[0]*1.e3', 'Oxygen Molar Flow $[mmol/s]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Oxygen Molar Flow', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.gas_flow[1]*1.e3', 'Water Molar Flow $[mmol/s]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Water Molar Flow Cathode', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.gas_flow[2]*1.e3', 'Nitrogen Molar Flow $[mmol/s]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Nitrogen Molar Flow', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('anode.gas_flow[0]*1.e3', 'Hydrogen Molar Flow $[mmol/s]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Hydrogen Molar Flow', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('anode.gas_flow[1]*1.e3', 'Water Molar Flow $[mmol/s]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Water Molar Flow Anode', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.mol_f[0]', 'Oxygen  Molar Fraction', 'Channel Location $[m]$',
                           'linear', 'k', 'Oxygen_Molar_Fraction', q,[0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.mol_f[1]', 'Gas Water  Molar Fraction', 'Channel Location $[m]$',
                           'linear', 'k', 'Water Molar Fraction Cathode', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('anode.mol_f[0]', 'Hydrogen Molar Fraction', 'Channel Location $[m]$',
                           'linear', 'k', 'Hydrogen_Molar_Fraction_Anode', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('anode.mol_f[1]', 'Gas Water  Molar Fraction', 'Channel Location $[m]$',
                           'linear', 'k', 'Water_Molar_Fraction_Anode', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.w*1e3', 'Liquid Water Flow $[mmol/s]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Liquid Water Flow Cathode', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.gamma*1e3', 'Water Condensation Rate $[mmol/s]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Water Condensation Rate Cathode', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.humidity', 'Relative Humidity', 'Channel Location $[m]$',
                           'linear', 'k', 'Relative Humidity Cathode', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.m_flow*1e6', 'Cathode Channel Gas Massflow $[mg/s]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Cathode_Channel_Massflow', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.m_full_flow*1e6', 'Cathode Channel Massflow $[mg/s]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Cathode_Channel_Massflow', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.g_full*1e3', 'Cathode Capacity Flow $[mW/K]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Cathode Capacity Flow', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.m_reac_flow*1e6', 'Oxygen Massflow $[mg/s]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Oxygen_massflow', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.m_vap_water_flow*1e6', 'Vapour Massflow $[mg/s]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Vapour Massflow', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('anode.m_flow*1e6', 'Hydrogen Massflow $[mg/s]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Hydrogen_massflow', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.cp_full', 'Cathode Heat Capacity $[J/(kgK)]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Cathode Heat Capacity', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('cathode.p', 'Cathode Channel Pressure $[Pa]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Cathode Channel Pressure', q, [0., i_p.channel_length],
                           x_node, False)
        self.plot_cell_var('anode.p', 'Anode Channel Pressure $[Pa]$', 'Channel Location $[m]$',
                           'linear', 'k', 'Anode Channel Pressure', q, [0., i_p.channel_length],
                           x_node, False)

        for l, item in enumerate(self.stack.t):
            if (l > 0 and self.stack.cool_ch_bc is False) or self.stack.cool_ch_bc is True:
                plt.plot(x_node, self.stack.t[l], label=l, marker='.',
                         color=plt.cm.coolwarm((l)/(self.stack.cell_numb)))
        #plt.legend()
        plt.grid()
        plt.ylabel(r'Coolant Temperature $[K]$')
        plt.xlabel('Channel Location $[m]$')
        plt.autoscale(tight=True, axis='both', enable=True)
        plt.xlim(0, self.stack.cell_list[0].cathode.channel.length)
        plt.savefig(os.path.join(os.path.dirname(__file__),
                                 'Plots' + q + '/' + 'Coolant' + '.jpg'))
        plt.close()

        x_vecz = np.array([0., i_p.plate_thick, i_p.gde_thick,
                           i_p.mem_thick, i_p.gde_thick])
        x_vec_e = np.array([i_p.plate_thick, i_p.plate_thick,
                            i_p.gde_thick, i_p.mem_thick, i_p.gde_thick])
        x = []
        for l, item in enumerate(self.stack.cell_list):
            if q is 0:
                    x.append(x_vecz)
            else:
                    x.append(x_vec_e)
        x = np.cumsum(np.block(x))
        for w in range(gpar.dict_case['nodes']):
            t_vec = []
            for l in self.stack.cell_list:
                t_vec.append(np.array([l.t3[w], l.t2[w],
                                       l.t1[w], l.t4[w], l.t5[w]]))
            plt.plot(x, np.block(t_vec), marker='o', color='k')
            plt.xlim(0, x[-1])
            plt.xlabel('Stack Location $[m]$')
            plt.ylabel('Temperature $[K]$')
            plt.autoscale(tight=True, axis='both', enable=True)
            plt.tight_layout()
            plt.savefig(os.path.join(os.path.dirname(__file__),
                                     'Plots' + q + '/' + 'Z-Cut-Temperature_'
                                     + str(w) + '.jpg'))
            plt.close()



Simulation_runs = Simulation(i_p.simulation)
Simulation_runs.update()