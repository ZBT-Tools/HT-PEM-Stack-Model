import data.stack_dict as st_dict
import data.simulation_dict as sim
import input.operating_conditions as oper_con
import system.stack as st
import numpy as np
import data.global_parameters as g_par
import system.global_functions as g_func
import matplotlib.pyplot as plt
import os
import errno
import timeit
import data.output_variables as out_var

np.set_printoptions(threshold=np.nan, linewidth=100, precision=9, suppress=True)

class Simulation:

    def __init__(self, dict_simulation):
        # Handover
        # iteration criteria
        self.it_crit = dict_simulation['iteration_criteria']
        # maximal number of iterations before force termination
        self.max_it = dict_simulation['maximal_iteration']
        # convergence criteria of the current density
        convergence_cd = []
        # convergence criteria of the temperature system
        convergence_temp = []
        # convergence criteria of the flow distribution
        convergence_mfd = []
        # convergence criteria
        self.convergence_criteria = [convergence_cd,
                                     convergence_temp,
                                     convergence_mfd]

    def update(self):
        """
        This function coordinates the program sequence
        """

        for q, item in enumerate(oper_con.target_current_density):
            g_par.dict_case['tar_cd'] = oper_con.target_current_density[q]
            self.stack = st.Stack(st_dict.dict_stack)
            statement = True
            counter = 0
            while statement is True:
                self.save_old_value()
                self.stack.update()
                if self.stack.break_program is True:
                    break
                self.calc_convergence_criteria()
                if len(oper_con.target_current_density) < 1:
                    print(counter)
                counter = counter + 1
                if (self.convergence_criteria[0][-1] < self.it_crit
                    and self.convergence_criteria[1][-1] < self.it_crit)\
                        or counter > self.max_it:
                    statement = False
                    self.output(q)

    def calc_convergence_criteria(self):
        """
        Calculates the convergence criteria according to (Koh, 2003)
        """

        self.convergence_criteria[0].\
            append(g_func.calc_convergence_criteria(
                   self.stack.i_ca[0][-1],
                   self.stack.i_ca_old[0][-1]))
        self.convergence_criteria[1].\
            append(g_func.calc_convergence_criteria(
                   self.stack.temp_cpl_stack.temp_layer[0][0][-1],
                   self.temp_old))
        self.convergence_criteria[2].\
            append((np.average(self.stack.cathode_mfd_criteria) +
                   np.average(self.stack.anode_mfd_criteria)) * .5)
        print(self.convergence_criteria[2][-1])

    def save_old_value(self):
        """
        Saves the average layer temperature value of the current iteration
        as the old temperature value for the next iteration.
        """
        self.temp_old = self.stack.temp_cpl_stack.temp_layer[0][0][-1]

    def save_csv_data(self, path):
        """Creates csv data from input values"""
        # Save Data from cell level
        for i in range(len(out_var.cell_var)):
            file_name = str(eval(
                "out_var." + "cell_var" + "[" + str(i) + "]" + "[1]")) + ".csv"
            file_name_space = str(eval("out_var." + "cell_var"
                              + "[" + str(i) + "]" + "[0]"))
            temp = []
            for j in range(len(self.stack.cells)):
                cell_pos = "self.stack.cells" + "[" + str(j) + "]" + "."
                temp.append(eval(cell_pos + file_name_space))
            np.savetxt(path + file_name, temp)
        # Save Data from stack level
        for i in range(len(out_var.stack_var)):
            file_name = str(eval(
                "out_var." + "stack_var" + "[" + str(i) + "]" + "[1]")) + ".csv"
            file_path = "self.stack." + str(eval("out_var." + "stack_var"
                              + "[" + str(i) + "]" + "[0]"))
            np.savetxt(path + file_name, eval(file_path))

    def save_plt_data(self, path):
        x_axis = [np.linspace(0,
                             self.stack.cells[0].cathode.channel.length,
                             g_par.dict_case['nodes']),
                  np.linspace(0, self.stack.cells[0].cathode.channel.length,
                              g_par.dict_case['elements'])]
        x = []
        for i in range(len(out_var.cell_var)):
            print(i)
            file_name = str(eval(
                "out_var." + "cell_var" + "[" + str(i) + "]" + "[1]"))
            unit = str(eval(
                "out_var." + "cell_var" + "[" + str(i) + "]" + "[2]"))
            file_name_space = str(eval("out_var." + "cell_var"
                              + "[" + str(i) + "]" + "[0]"))
            if len(eval("self.stack.cells[0]." + file_name_space)) \
                    is g_par.dict_case['elements']:
                x = x_axis[1]
            else:
                x = x_axis[0]

            for j in range(len(self.stack.cells)):
                cell_pos = "self.stack.cells" + "[" + str(j) + "]" + "."
                plt.plot(x, eval(cell_pos + file_name_space),
                         color=plt.cm.coolwarm(j / (len(self.stack.cells))),
                         marker='.')

            plt.xlabel("x-Axis", fontsize=14)
            plt.ylabel(unit, fontsize=14)
            plt.title(file_name, fontsize=16)
            plt.tick_params(labelsize=14)
            plt.autoscale(tight=True, axis='both', enable=True)
            plt.tight_layout()
            plt.grid()
            plt.savefig(os.path.join(path + file_name + '.jpg'))
            plt.close()

    def output(self, q):
        path_data = [os.path.join(os.path.dirname(__file__),
                     'output/' + 'case' + str(q) + '/csv_data' + '/'),
                     os.path.join(os.path.dirname(__file__),
                     'output/' + 'case' + str(q) + '/plt_data' + '/')]
        for i, item in enumerate(path_data):
            try:
                os.makedirs(item)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
        self.save_csv_data(path_data[0])
        self.save_plt_data(path_data[1])







start = timeit.default_timer()
Simulation_runs = Simulation(sim.simulation)
Simulation_runs.update()
stop = timeit.default_timer()
print('Simulation time:', stop - start)
