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
        self.max_it = dict_simulation['maximum_iteration']
        # maximal number of iterations before force termination#
        self.min_it = dict_simulation['minimum_iteration']

        n_cells = input_dicts.dict_stack['cell_number']
        # number of stack cells
        n_nodes = g_par.dict_case['nodes']
        # node points of the x-grid

        self.timing = {'start': 0.0,
                       'initialization': 0.0,
                       'simulation': 0.0,
                       'output': 0.0}

        """General variables"""

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

        # initialize stack object
        self.stack = stack.Stack()

        output_dict = input_dicts.dict_output
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
            # g_par.dict_case['tar_cd'] = tar_cd
            counter = 0
            while True:
                if counter == 0:
                    self.stack.update(tar_cd)
                else:
                    self.stack.update()
                if self.stack.break_program:
                    break
                self.calc_convergence_criteria()
                if len(target_current_density) < 1:
                    print(counter)
                counter += 1
                if ((self.i_ca_criteria < self.it_crit
                     and self.temp_criteria < self.it_crit)
                        and counter > self.min_it) or counter > self.max_it:
                    break
            if not self.stack.break_program:
                voltage_loss = self.save_voltages(self.stack)
                cell_voltages.append(np.average([cell.v for cell in
                                                 self.stack.cells]))

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
            np.asarray([cell.membrane.v_loss for cell in fc_stack.cells])

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
        i_cd_vec = self.stack.i_cd.flatten()
        self.i_ca_criteria = \
            np.abs(np.sum(((i_cd_vec - self.stack.i_cd_old.flatten())
                           / i_cd_vec) ** 2.0))
        # self.temp_criteria =\
        #     np.abs(np.sum(((self.temp_old
        #                     - self.stack.temp_sys.temp_layer[0][0, 0]))
        #                   / self.stack.temp_sys.temp_layer[0][0, 0]))
        self.temp_criteria =\
            np.abs(np.sum(((self.stack.temp_old
                            - self.stack.temp_sys.temp_layer_vec)
                          / self.stack.temp_sys.temp_layer_vec) ** 2.0))

        self.temp_criteria_process.append(self.temp_criteria)
        # self.mfd_cat_criteria.append(self.stack.manifolds[0].criteria)
        # self.mfd_ano_criteria.append(self.stack.manifolds[1].criteria)
        self.i_ca_criteria_process.append(self.i_ca_criteria)


start_time = timeit.default_timer()
simulation = Simulation(input_dicts.simulation_dict)
simulation.timing['start'] = start_time
simulation.timing['initialization'] = timeit.default_timer()
print('Initialization time:',
      simulation.timing['initialization'] - simulation.timing['start'])
simulation.timing['start'] = start_time
simulation.update()
stop_time = timeit.default_timer()
print('Total time:', stop_time - start_time)
