import input.operating_conditions as op_con
import system.stack as stack
import numpy as np
import data.global_parameters as g_par
import cProfile
import data.input_dicts as input_dicts
import copy
import sys
import system.channel as chl
import system.fluid2 as fluids
import system.flow_circuit2 as flow_circuit
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


n_chl = 50
n_subchl = 1

channel_dict = {
    'name': 'Channel',
    'channel_length': 0.4,
    'cross_sectional_shape': 'circular',
    'channel_diameter': 0.00295,
    'p_out': 101325.0,
    'temp_in': 300.0,
    'hum_in': 0.1,
    'flow_direction': 1,
    'channel_width': 0.003,
    'channel_height': 0.003,
    'bend_number': 0,
    'bend_friction_factor': 0.1,
    'additional_friction_fractor': 0.0
    }


fluid_dict = {
    'fluid_name': 'Cathode Gas',
    'temp_init': op_con.temp_air_in,
    'press_init': op_con.p_manifold_cathode_out
}

in_manifold_dict = {
    'name': 'Inlet Manifold',
    'channel_length': 0.09,
    'p_out': channel_dict['p_out'],
    'temp_in': 300.0,
    'hum_in': 0.1,
    'flow_direction': 1,
    'channel_width': 0.0075,
    'channel_height': 0.0075,
    'bend_number': 0,
    'bend_friction_factor': 0.0,
    'additional_friction_fractor': 0.0
    }

out_manifold_dict = copy.deepcopy(in_manifold_dict)
out_manifold_dict['name'] = 'Outlet Manifold'

flow_circuit_dict = {
    'name': 'Flow Circuit',
    'type': 'Wang',
    'shape': 'U'
    }

flow_model = flow_circuit.flow_circuit_factory(flow_circuit_dict, fluid_dict,
                                               channel_dict, in_manifold_dict,
                                               out_manifold_dict, n_chl,
                                               n_subchl)


x = (ip.interpolate_1d(flow_model.manifolds[0].x)
     - flow_model.manifolds[0].dx * 0.5) \
    / (flow_model.manifolds[0].length - flow_model.manifolds[0].dx[0])
# x = flow_model.manifolds[0].x / flow_model.manifolds[0].length

flow_model.update(inlet_mass_flow=0.010125*3.333)
q = flow_model.channel_vol_flow / np.average(flow_model.channel_vol_flow)
reynolds = flow_model.manifolds[0].reynolds[0]
plt.plot(x, q, label='Re='+str(reynolds))

# flow_model.update(inlet_mass_flow=0.01485)
# q = flow_model.channel_vol_flow / np.average(flow_model.channel_vol_flow)
# reynolds = flow_model.manifolds[0].reynolds[0]
# plt.plot(x, q, label='Re='+str(reynolds))

plt.legend()
plt.ylim([0.4, 2.2])
plt.xlim([0.0, 1.0])

plt.yticks([0.4, 1.0, 1.6, 2.2])
plt.show()

