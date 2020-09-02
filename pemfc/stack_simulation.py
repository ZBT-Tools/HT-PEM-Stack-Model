import numpy as np
import timeit
from pemfc import simulation
import sys

np.set_printoptions(threshold=sys.maxsize, linewidth=10000,
                    precision=9, suppress=True)
np.seterr(all='raise')


def main():
    start_time = timeit.default_timer()
    sim = simulation.Simulation()
    sim.timing['start'] = start_time
    sim.timing['initialization'] = timeit.default_timer()
    # simulation.timing['start'] = start_time
    sim.run()
    print('Initialization time: ',
          sim.timing['initialization'] - sim.timing['start'])
    print('Simulation time: ', sim.timing['simulation'])
    print('Output time: ', sim.timing['output'])
    stop_time = timeit.default_timer()
    print('Total time:', stop_time - start_time)
    print('Stack Voltage [V]: ', sim.stack.v_stack)
    print('Average Cell Voltage [V]: ',
          sim.stack.v_stack/sim.stack.n_cells)
    print('Minimum Cell Voltage [V]: ', np.min(sim.stack.v))
    print('Maximum Cell Voltage [V]: ', np.max(sim.stack.v))
    average_current_density = \
        np.average([np.average(cell.i_cd, weights=cell.active_area_dx)
                    for cell in sim.stack.cells])
    print('Stack current density [A/m²]: ', average_current_density)
    average_current = \
        average_current_density * sim.stack.cells[0].active_area
    print('Stack power [W]: ', sim.stack.v_stack * average_current)

