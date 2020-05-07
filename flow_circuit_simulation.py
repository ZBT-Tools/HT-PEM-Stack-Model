import settings.operating_conditions as op_con
import system.stack as stack
import numpy as np
import data.global_parameters as g_par
import cProfile
import data.input_dicts as input_dicts
import copy
import sys
import system.channel as chl
import system.fluid as fluids
import system.flow_circuit as flow_circuit
import matplotlib.pyplot as plt
import system.interpolation as ip

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


n_chl = 40
n_subchl = 1

channel_dict = {
    'name': 'Channel',
    'length': 0.65,
    'cross_sectional_shape': 'rectangular',
    'width': 4e-3,
    'height': 1e-3,
    'p_out': 101325.0,
    'temp_in': 293.15,
    'flow_direction': 1,
    'bend_number': 0,
    'bend_friction_factor': 0.1,
    'constant_friction_factor': 0.0
    }


fluid_dict = {
    'fluid_name': 'Cathode Gas',
    'fluid_components': {'O2': 'gas', 'N2': 'gas', 'H2O': 'gas-liquid'},
    'inlet_composition': [0.21, 0.79, 0.0],
    'temp_init': 293.15,
    'press_init': 101325.0,
    'nodes': 10
}

in_manifold_dict = {
    'name': 'Inlet Manifold',
    'length': 0.27,
    'p_out': channel_dict['p_out'],
    'temp_in': 293.15,
    'flow_direction': 1,
    'width': 12.5e-3,
    'height': 7.5e-3,
    'bend_number': 0,
    'bend_friction_factor': 0.0,
    'constant_friction_factor': 0.25,
    'flow_split_factor': 0.2
    }

out_manifold_dict = copy.deepcopy(in_manifold_dict)
out_manifold_dict['name'] = 'Outlet Manifold'
out_manifold_dict['constant_friction_factor'] = 0.25
out_manifold_dict['flow_split_factor'] = 0.2

flow_circuit_dict = {
    'name': 'Flow Circuit',
    'type': 'ModifiedKoh',
    'shape': 'U'
    }

flow_model = flow_circuit.factory(flow_circuit_dict, fluid_dict,
                                  channel_dict, in_manifold_dict,
                                  out_manifold_dict, n_chl,
                                  n_subchl)


x = (ip.interpolate_1d(flow_model.manifolds[0].x)
     - flow_model.manifolds[0].dx * 0.5) \
    / (flow_model.manifolds[0].length - flow_model.manifolds[0].dx[0])
# x = flow_model.manifolds[0].x / flow_model.manifolds[0].length

flow_model.update(inlet_mass_flow=8.91E-04)
q = (flow_model.normalized_flow_distribution - 1.0) * 100.0
reynolds = flow_model.manifolds[0].reynolds[0]
plt.plot(x, q, label='Re={0:.2f}'.format(reynolds), color='k')

flow_model.update(inlet_mass_flow=0.00059425)
q = (flow_model.normalized_flow_distribution - 1.0) * 100.0
reynolds = flow_model.manifolds[0].reynolds[0]
plt.plot(x, q, label='Re={0:.2f}'.format(reynolds), color='b')

flow_model.update(inlet_mass_flow=0.000297125)
q = (flow_model.normalized_flow_distribution - 1.0) * 100.0
reynolds = flow_model.manifolds[0].reynolds[0]
plt.plot(x, q, label='Re={0:.2f}'.format(reynolds), color='r')

# print('Normalized Flow Distribution: ',
#       flow_model.normalized_flow_distribution)
# np.savetxt('output/flow_distribution.txt',
#            (flow_model.normalized_flow_distribution - 1.0) * 100.0)
plt.show()
