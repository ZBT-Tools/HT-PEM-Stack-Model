import numpy as np
import timeit
from pemfc import simulation
import sys
import os

np.set_printoptions(threshold=sys.maxsize, linewidth=10000,
                    precision=9, suppress=True)
np.seterr(all='raise')


def main():
    np.seterr(all='raise')
    start_time = timeit.default_timer()
    sim = simulation.Simulation()
    sim.timing['start'] = start_time
    sim.timing['initialization'] = timeit.default_timer()
    # simulation.timing['start'] = start_time
    g_data, l_data = sim.run()
    return g_data, l_data, sim


if __name__ == "__main__":
    global_data, local_data, sim = main()
    summary_file_path = os.path.join(sim.output.output_dir, 'pemfc_output.txt')
    with open(summary_file_path, 'w') as file:
        file.write('Initialization time: {}\n'.format(
              (sim.timing['initialization'] - sim.timing['start'])))
        file.write('Simulation time: {}\n'.format(sim.timing['simulation']))
        file.write('Output time: {}\n'.format(sim.timing['output']))
        stop_time = timeit.default_timer()
        file.write('Total time:{}\n'.format(stop_time - sim.timing['start']))
        for k, v in global_data.items():
            file.write('{} [{}]: {}\n'.format(k, v['value'], v['units']))

