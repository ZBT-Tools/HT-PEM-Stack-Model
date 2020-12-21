import numpy as np
import timeit
from pemfc import simulation
import sys

np.set_printoptions(threshold=sys.maxsize, linewidth=10000,
                    precision=9, suppress=True)
np.seterr(all='raise')


start_time = timeit.default_timer()
simulation = simulation.Simulation()
simulation.timing['start'] = start_time
simulation.timing['initialization'] = timeit.default_timer()
# simulation.timing['start'] = start_time
simulation.run()
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
print('Maximum Cell Voltage [V]: ', np.max(simulation.stack.v))
average_current_density = \
    np.average([np.average(cell.i_cd, weights=cell.active_area_dx)
                for cell in simulation.stack.cells])
print('Stack current density [A/m²]: ', average_current_density)
average_current = \
    average_current_density * simulation.stack.cells[0].active_area
print('Stack power [W]: ', simulation.stack.v_stack * average_current)