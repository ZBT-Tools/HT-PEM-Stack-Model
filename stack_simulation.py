import settings.operating_conditions as op_con
import system.stack as stack
import numpy as np
import data.global_parameters as g_par
import cProfile
import timeit
import data.input_dicts as input_dicts
import system.output_object as out_obj
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

        self.timing = {'start': 0.0,
                       'initialization': 0.0,
                       'simulation': 0.0,
                       'output': 0.0}

        """General variables"""
        # initialize stack object
        self.current_control = g_par.dict_case['current_control']
        self.stack = stack.Stack(current_control=self.current_control)

        # initialize output object
        output_dict = input_dicts.dict_output
        self.output = output.Output(output_dict)

    @do_c_profile
    def update(self):
        """
        This function coordinates the program sequence
        """
        if self.current_control:
            target_value = op_con.current_density
        else:
            target_value = op_con.average_cell_voltage * self.stack.n_cells
        if not isinstance(target_value, (list, tuple, np.ndarray)):
            target_value = [target_value]
        cell_voltages = []
        current_errors = []
        temp_errors = []
        for i, tar_value in enumerate(target_value):
            # g_par.dict_case['tar_cd'] = tar_cd
            simulation_start_time = timeit.default_timer()
            counter = 0
            while True:
                if counter == 0:
                    if self.current_control:
                        self.stack.update(current_density=tar_value)
                    else:
                        self.stack.update(voltage=tar_value)
                else:
                    self.stack.update()
                if self.stack.break_program:
                    break
                current_error, temp_error = self.calc_convergence_criteria()
                current_errors.append(current_error)
                temp_errors.append(temp_error)
                if len(target_value) < 1:
                    print(counter)
                counter += 1
                if ((current_error < self.it_crit and temp_error < self.it_crit)
                        and counter > self.min_it) or counter > self.max_it:
                    break
            simulation_stop_time = timeit.default_timer()
            simulation_time = simulation_stop_time - simulation_start_time
            self.timing['simulation'] += simulation_time
            output_start_time = timeit.default_timer()

            if not self.stack.break_program:
                voltage_loss = self.get_voltage_losses(self.stack)
                cell_voltages.append(np.average([cell.v for cell in
                                                 self.stack.cells]))

                case_name = 'Case'+str(i)
                self.output.save(case_name, self.stack)
                path = os.path.join(self.output.output_dir, case_name,
                                    'plots', 'Convergence.png')
                self.output.create_figure(path, list(range(counter)),
                                          [current_errors, temp_errors],
                                          xlabels='Iteration', ylabels='Error',
                                          yscale='log',
                                          legend=['Current Density',
                                                  'Temperature'])
                # self.output.plot([current_errors, temp_errors],
                #                  'ERR', 'Iteration', 'log', ['k', 'r', 'b'],
                #                  'Convergence', 0., len(current_errors),
                #                  ['Current Density', 'Temperature'], path)
            else:
                target_value = target_value[0:-i]
                break
            output_stop_time = timeit.default_timer()
            self.timing['output'] += output_stop_time - output_start_time
        output_start_time = timeit.default_timer()
        if len(target_value) > 1:
            self.output.plot_polarization_curve(voltage_loss, cell_voltages,
                                                target_value)
        output_stop_time = timeit.default_timer()
        self.timing['output'] += output_stop_time - output_start_time

    @staticmethod
    def get_voltage_losses(fc_stack):
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
        current_error = \
            np.abs(np.sum(((i_cd_vec - self.stack.i_cd_old.flatten())
                           / i_cd_vec) ** 2.0))
        # self.temp_criteria =\
        #     np.abs(np.sum(((self.temp_old
        #                     - self.stack.temp_sys.temp_layer[0][0, 0]))
        #                   / self.stack.temp_sys.temp_layer[0][0, 0]))
        temp_error =\
            np.abs(np.sum(((self.stack.temp_old
                            - self.stack.temp_sys.temp_layer_vec)
                          / self.stack.temp_sys.temp_layer_vec) ** 2.0))
        return current_error, temp_error


start_time = timeit.default_timer()
simulation = Simulation(input_dicts.simulation_dict)
simulation.timing['start'] = start_time
simulation.timing['initialization'] = timeit.default_timer()
# simulation.timing['start'] = start_time
simulation.update()
print('Initialization time: ',
      simulation.timing['initialization'] - simulation.timing['start'])
print('Simulation time: ', simulation.timing['simulation'])
print('Output time: ', simulation.timing['output'])
stop_time = timeit.default_timer()
print('Total time:', stop_time - start_time)
print('Stack Voltage [V]: ', simulation.stack.v_stack)
print('Average Cell Voltage [V]: ',
      simulation.stack.v_stack/simulation.stack.n_cells)
print('Minimum Cell Voltage [V]: ', np.min(simulation.stack.v))
print('Average cell current density [A/m²]: ',
      [np.average(cell.i_cd, weights=cell.active_area_dx)
       for cell in simulation.stack.cells])
print('Cathode stoichiometry [-]: ', [cell.cathode.inlet_stoi
                                      for cell in simulation.stack.cells])
print('Anode stoichiometry [-]: ', [cell.anode.inlet_stoi
                                    for cell in simulation.stack.cells])
# print(out_obj.OutputObject.make_name_list())
# for obj in out_obj.OutputObject.getinstances():
#     print(obj.name_list)