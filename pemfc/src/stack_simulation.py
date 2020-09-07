import numpy as np
import timeit
from pemfc import simulation
import sys
import os

np.set_printoptions(threshold=sys.maxsize, linewidth=10000,
                    precision=9, suppress=True)
np.seterr(all='raise')


def main():
    start_time = timeit.default_timer()
    sim = simulation.Simulation()
    sim.timing['start'] = start_time
    sim.timing['initialization'] = timeit.default_timer()
    # simulation.timing['start'] = start_time

    ### return value only for testing executable
    avg_icd = sim.run()

    summary_file_path = os.path.join(sim.output.output_dir, 'pemfc_output.txt')
    with open(summary_file_path, 'w') as file:
        file.write('Initialization time: {}\n'.format(
              (sim.timing['initialization'] - sim.timing['start'])))
        file.write('Simulation time: {}\n'.format(sim.timing['simulation']))
        file.write('Output time: {}\n'.format(sim.timing['output']))
        stop_time = timeit.default_timer()
        file.write('Total time:{}\n'.format(stop_time - start_time))
        file.write('Stack Voltage [V]: {}\n'.format(sim.stack.v_stack))
        file.write('Average Cell Voltage [V]: {}\n'.format(
              sim.stack.v_stack/sim.stack.n_cells))
        file.write('Minimum Cell Voltage [V]: {}\n'.format(np.min(sim.stack.v)))
        file.write('Maximum Cell Voltage [V]: {}\n'.format(np.max(sim.stack.v)))
        average_current_density = \
            np.average([np.average(cell.i_cd, weights=cell.active_area_dx)
                        for cell in sim.stack.cells])
        file.write('Stack current density [A/m²]: {}\n'.format(
            average_current_density))
        average_current = \
            average_current_density * sim.stack.cells[0].active_area
        file.write('Stack power density [W/m²]: {}\n'.format(
              sim.stack.v_stack * average_current_density))
        file.write('Stack power [W]: {}\n'.format(sim.stack.v_stack *
                   average_current))

    return avg_icd

