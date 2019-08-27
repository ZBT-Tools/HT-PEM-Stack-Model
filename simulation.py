import input.operating_conditions as op_con
import system.stack as stack
import numpy as np
import data.global_parameters as g_par
import cProfile
import timeit
import data.input_dicts as input_dicts
import output
import os
import sys

np.set_printoptions(threshold=sys.maxsize, linewidth=10000,
                    precision=9, suppress=True)
np.seterr(all='raise')


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

        n_cells = input_dicts.dict_stack['cell_number']
        # number of stack cells
        n_nodes = g_par.dict_case['nodes']
        # node points of the x-grid

        """General variables"""
        self.temp_old = None
        # defined temperature of the last iteration
        self.mfd_cat_criteria = []
        # array of the cathodic mdf criteria over the iterations
        self.mfd_ano_criteria = []
        # array of the anodic mdf criteria over the iterations
        self.i_ca_criteria_process = []
        # array of the current density criteria over the iterations
        self.temp_criteria_process = []
        # array of the temperature criteria over the iterations
        self.i_ca_criteria = None
        # convergence criteria of the current density
        self.temp_criteria = None
        # convergence criteria of the temperature
        #self.mol_flow = np.full((6, n_cells, n_nodes), 0.)
        # molar flow of the species in the channels
        # 0: oxygen, 1: cathode water, 2: cathode nitrogen,
        # 3: hydrogen, 4: anode water, 5: anode nitrogen
        #self.gas_con = np.full((6, n_cells, n_nodes), 0.)
        # molar concentration of the species in the gas mixture
        #self.m_f = np.full((6, n_cells, n_nodes), 0.)
        # mass fraction of the species in the gas mixture
        #self.mol_f = np.full((6, n_cells, n_nodes), 0.)
        # molar fraction of the species in the gas mixture
        #self.act_loss_cat = np.full((n_cells, n_nodes - 1), 0.)
        # cathodic activation voltage loss
        #self.act_loss_ano = np.full((n_cells, n_nodes - 1), 0.)
        # anodic activation voltage loss
        #self.cl_diff_loss_cat = np.full((n_cells, n_nodes - 1), 0.)
        # cathodic catalyst layer diffusion voltage loss
        #self.cl_diff_loss_ano = np.full((n_cells, n_nodes - 1), 0.)
        # anodic catalyst layer diffusion voltage loss
        #self.gdl_diff_loss_cat = np.full((n_cells, n_nodes - 1), 0.)
        # cathodic gas diffusion layer diffusion voltage loss
        #self.gdl_diff_loss_ano = np.full((n_cells, n_nodes - 1), 0.)
        # anodic gas diffusion layer diffusion voltage loss
        #self.mem_loss = np.full((n_cells, n_nodes - 1), 0.)
        # membrane voltage loss
        #self.v_loss = np.full((n_cells, n_nodes - 1), 0.)
        # voltage loss over the stack
        #self.v_cell = np.full((n_cells, n_nodes - 1), 0.)
        # cell voltage
        #self.cp = np.full((6, n_cells, n_nodes), 0.)
        # heat capacity of the species in the gas phase
        #self.lambda_gas = np.full((6, n_cells, n_nodes), 0.)
        # heat conductivity of the species in the gas phase
        #self.visc = np.full((6, n_cells, n_nodes), 0.)
        # viscosity of the species in the gas phase
        #self.r_gas_cat = np.full((n_cells, n_nodes), 0.)
        # gas constant of the gas phase in the cathode channels
        #self.r_gas_ano = np.full((n_cells, n_nodes), 0.)
        # gas constant of the gas phase in the anode channels
        #self.cp_gas_cat = np.full((n_cells, n_nodes), 0.)
        # heat capacitiy of the gas phase in the cathode channels
        #self.cp_gas_ano = np.full((n_cells, n_nodes), 0.)
        # heat capacity of the gas phase in the anode channels
        #self.visc_gas_cat = np.full((n_cells, n_nodes), 0.)
        # viscosity of the gas phase in the cathode channels
        #self.visc_gas_ano = np.full((n_cells, n_nodes), 0.)
        # viscosity of the gas phase in the anode channels
        #self.lambda_gas_cat = np.full((n_cells, n_nodes), 0.)
        # heat conductivity of the gas phase in the cathode channels
        #self.lambda_gas_ano = np.full((n_cells, n_nodes), 0.)
        # heat conductivity of the gas phase in the anode channels
        #self.rho_gas_cat = np.full((n_cells, n_nodes), 0.)
        # density of the gas in the cathode channels
        #self.rho_gas_cat = np.full((n_cells, n_nodes), 0.)
        # density of the gas in the anode channels
        #self.pr_gas_cat = np.full((n_cells, n_nodes), 0.)
        # prandtl number of the gas in the cathode channels
        # self.pr_gas_ano = np.full((n_cells, n_nodes), 0.)
        # # prandtl number of the gas in the anode channels
        # self.u_gas_cat = np.full((n_cells, n_nodes), 0.)
        # # velocity of the fluid in the cathode channels
        # self.u_gas_ano = np.full((n_cells, n_nodes), 0.)
        # # velocity of the fluid in the anode channels
        # self.re_gas_cat = np.full((n_cells, n_nodes), 0.)
        # # reynolds number of the fluid in the cathode channels
        # self.re_gas_ano = np.full((n_cells, n_nodes), 0.)
        # # reynolds number of the fluid in the anode channels
        # self.p_cat = np.full((n_cells, n_nodes), 0.)
        # # pressure in the cathode channels
        # self.p_ano = np.full((n_cells, n_nodes), 0.)
        # # pressure in the anode channels
        # self.ht_coef_cat = np.full((n_cells, n_nodes), 0.)
        # # heat convection coefficient in the cathode channels
        # self.ht_coef_ano = np.full((n_cells, n_nodes), 0.)
        # # heat convection coefficient in the anode channels
        # self.m_flow_fluid_cat = np.full((n_cells, n_nodes), 0.)
        # # mass flow of the fluid in the cathode channels
        # self.m_flow_fluid_ano = np.full((n_cells, n_nodes), 0.)
        # # mass flow of the fluid in the anode channels
        # self.cp_fluid_cat = np.full((n_cells, n_nodes), 0.)
        # # heat capacity of the cathode fluid
        # self.cp_fluid_ano = np.full((n_cells, n_nodes), 0.)
        # # heat capacity of the anode fluid
        # self.temp_layer = []
        # # temperature layer of the stack
        # self.temp_fluid_ano = np.full((n_cells, n_nodes), 0.)
        # # temperature of the fluid in the anode channels
        # self.temp_fluid_cat = np.full((n_cells, n_nodes), 0.)
        # # temperature of the fluid in the cathode channels
        # self.temp_cool = np.full((n_cells, n_nodes), 0.)
        # # temperature of the coolant
        # self.stoi_cat = np.full(n_cells, 0.)
        # # inlet stoichiometry of the cathode channels
        # self.stoi_ano = np.full(n_cells, 0.)
        # # inlet stoichiometry of the anode channels
        # self.act_loss_ui_ano = []
        # # average activation voltage loss of the anode
        # self.act_loss_ui_cat = []
        # # average activation voltage loss of the cathode
        # self.cl_diff_loss_ui_ano = []
        # # average anode catalyst layer diffusion voltage losses
        # self.cl_diff_loss_ui_cat = []
        # # average cathode catalyst layer diffusion voltage losses
        # self.gdl_diff_loss_ui_ano = []
        # # average anode gdl diffusion voltage losses
        # self.gdl_diff_loss_ui_cat = []
        # # average cathode gdl diffusion voltage losses
        # self.mem_loss_ui = []
        # # average membrane voltage losses

        # load input dictionaries
        stack_dict = input_dicts.dict_stack
        cell_dict = input_dicts.dict_cell
        anode_dict = input_dicts.dict_anode
        cathode_dict = input_dicts.dict_cathode
        ano_channel_dict = input_dicts.dict_anode_channel
        cat_channel_dict = input_dicts.dict_cathode_channel
        ano_manifold_dict = input_dicts.dict_mfold_ano
        cat_manifold_dict = input_dicts.dict_mfold_cat
        electrical_dict = input_dicts.dict_electrical_coupling
        temperature_dict = input_dicts.dict_temp_sys
        output_dict = input_dicts.dict_output

        # initialize stack object
        self.stack = stack.Stack(stack_dict, cell_dict, anode_dict,
                                 cathode_dict, ano_channel_dict,
                                 cat_channel_dict, ano_manifold_dict,
                                 cat_manifold_dict, electrical_dict,
                                 temperature_dict)
        self.output = output.Output(output_dict)

    # @do_c_profile
    def update(self):
        """
        This function coordinates the program sequence
        """
        target_current_density = op_con.target_current_density
        if not isinstance(target_current_density, (list, tuple, np.ndarray)):
            target_current_density = [target_current_density]

        cell_voltages = []
        for i, tar_cd in enumerate(target_current_density):
            g_par.dict_case['tar_cd'] = tar_cd
            counter = 0
            while True:
                self.save_old_value()
                self.stack.update()
                if self.stack.break_program:
                    break
                self.calc_convergence_criteria()
                if len(target_current_density) < 1:
                    print(counter)
                counter += 1
                if ((self.i_ca_criteria < self.it_crit
                     and self.temp_criteria < self.it_crit) and counter > 10)\
                        or counter > self.max_it:
                    break
            if not self.stack.break_program:
                voltage_loss = self.save_voltages(self.stack)
                cell_voltages.append(np.average(self.stack.v_cell))

                mfd_criteria = \
                    (np.array(self.mfd_ano_criteria)
                     + np.array(self.mfd_cat_criteria)) * .5
                case_name = 'Case'+str(i)
                self.output.save(case_name, self.stack)
                path = os.path.join(self.output.output_dir, case_name, 'plots')
                self.output.plot([mfd_criteria,
                                  self.i_ca_criteria_process,
                                  self.temp_criteria_process],
                                 'ERR', 'Iteration', 'log', ['k', 'r', 'b'],
                                 'Convergence', 0.,
                                 len(self.temp_criteria_process),
                                 ['Flow Distribution', 'Current Density',
                                  'Temperature'], path)
            else:
                target_current_density = target_current_density[0:-i]
                break
            self.mfd_ano_criteria = []
            self.mfd_cat_criteria = []
            self.temp_criteria_process = []
            self.i_ca_criteria_process = []
        if len(target_current_density) > 1:
            self.output.plot_polarization_curve(voltage_loss, cell_voltages,
                                                target_current_density)

    @staticmethod
    def save_voltages(fc_stack):
        """
        Saves the average voltage losses of the stack
        """
        voltage_loss = \
            {'activation': {'anode': {}, 'cathode': {}},
             'diffusion': {'CL': {'anode': {}, 'cathode': {}},
                           'GDL': {'anode': {}, 'cathode': {}}},
             'membrane': {}}

        voltage_loss['activation']['anode']['cells'] = \
            np.asarray([cell.cathode.act_loss for cell in fc_stack.cells])
        voltage_loss['activation']['cathode']['cells'] = \
            np.asarray([cell.anode.act_loss for cell in fc_stack.cells])
        voltage_loss['diffusion']['CL']['anode']['cells'] = \
            np.asarray([cell.anode.cl_diff_loss for cell in fc_stack.cells])
        voltage_loss['diffusion']['CL']['cathode']['cells'] = \
            np.asarray([cell.cathode.cl_diff_loss for cell in fc_stack.cells])
        voltage_loss['diffusion']['GDL']['anode']['cells'] = \
            np.asarray([cell.anode.gdl_diff_loss for cell in fc_stack.cells])
        voltage_loss['diffusion']['GDL']['cathode']['cells'] = \
            np.asarray([cell.cathode.gdl_diff_loss for cell in fc_stack.cells])
        voltage_loss['membrane']['cells'] = \
            np.asarray([cell.mem_loss for cell in fc_stack.cells])

        voltage_loss['activation']['anode']['average'] = \
            np.average(voltage_loss['activation']['anode']['cells'])
        voltage_loss['activation']['cathode']['average'] = \
            np.average(voltage_loss['activation']['cathode']['cells'])
        voltage_loss['diffusion']['CL']['anode']['average'] = \
            np.average(voltage_loss['diffusion']['CL']['anode']['cells'])
        voltage_loss['diffusion']['CL']['cathode']['average'] = \
            np.average(voltage_loss['diffusion']['CL']['cathode']['cells'])
        voltage_loss['diffusion']['GDL']['anode']['average'] = \
            np.average(voltage_loss['diffusion']['GDL']['anode']['cells'])
        voltage_loss['diffusion']['GDL']['cathode']['average'] = \
            np.average(voltage_loss['diffusion']['GDL']['cathode']['cells'])
        voltage_loss['membrane']['average'] = \
            np.average(voltage_loss['membrane']['cells'])

        return voltage_loss

    def calc_convergence_criteria(self):
        """
        Calculates the convergence criteria according to (Koh, 2003)
        """
        self.i_ca_criteria = np.abs(np.sum(((self.stack.i_cd.flatten()
                                             - self.stack.i_cd_old.flatten())
                                            / self.stack.i_cd.flatten()) ** 2.))
        self.temp_criteria =\
            np.abs(np.sum(((self.temp_old
                            - self.stack.temp_sys.temp_layer[0][0, 0]))
                          / self.stack.temp_sys.temp_layer[0][0, 0]))

        self.temp_criteria_process.append(self.temp_criteria)
        self.mfd_cat_criteria.append(self.stack.manifold[0].criteria)
        self.mfd_ano_criteria.append(self.stack.manifold[1].criteria)
        self.i_ca_criteria_process.append(self.i_ca_criteria)

    def save_old_value(self):
        """
        Saves an defined temperature value of the current iteration
        as the old temperature value for the next iteration.
        """
        self.temp_old = self.stack.temp_sys.temp_layer[0][0, 0]


start = timeit.default_timer()
simulation = Simulation(input_dicts.simulation_dict)
simulation.update()
stop = timeit.default_timer()
print('Simulation time:', stop-start)
