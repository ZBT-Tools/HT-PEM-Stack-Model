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



n_chl = 100
n_subchl = 10

channel_dict = {
    'name': 'Channel',
    'channel_length': 0.4,
    'p_out': 101325.0,
    'temp_in': 300.0,
    'hum_in': 0.1,
    'flow_direction': 1,
    'channel_width': 0.0010,
    'channel_height': 0.0010,
    'bend_number': 40,
    'bend_friction_factor': 0.2,
    'additional_friction_fractor': 0.0
    }

fluid_dict = input_dicts.dict_cathode_fluid

in_manifold_dict = {
    'name': 'Inlet Manifold',
    'channel_length': 0.1,
    'p_out': channel_dict['p_out'],
    'temp_in': 300.0,
    'hum_in': 0.1,
    'flow_direction': 1,
    'channel_width': 0.01,
    'channel_height': 0.01,
    'bend_number': 0,
    'bend_friction_factor': 0.0,
    'additional_friction_fractor': 1.7
    }

out_manifold_dict = copy.deepcopy(in_manifold_dict)
out_manifold_dict['flow_direction'] = -1
out_manifold_dict['name'] = 'Outlet Manifold'

flow_model = flow_circuit.flow_circuit_factory(fluid_dict, channel_dict,
                                               in_manifold_dict,
                                               out_manifold_dict, n_chl,
                                               n_subchl)

# temp = np.full(n_chl + 1, 300.0)
# flow_model.manifolds[0].temp = temp
# flow_model.manifolds[1].temp = temp
#
# flow_model.manifolds[0].update()
# fluid_in = flow_model.manifolds[0].fluid
# vol_flow_in = 1e-5
# total_mole_flow_in = vol_flow_in * fluid_in.density[0] \
#                      / fluid_in.mw[0] * fluid_in.mole_fraction[:, 0]
#
# ones = np.ones(n_chl)
# dmole = np.outer(total_mole_flow_in / n_chl, ones)
# flow_model.manifolds[0].update(total_mole_flow_in, -dmole)
# print(flow_model.manifolds[0].fluid.temperature)
# print(flow_model.manifolds[0].fluid.mole_fraction)
# print(flow_model.manifolds[0].fluid.density)
# print(flow_model.manifolds[0].fluid.viscosity)
# print(flow_model.manifolds[0].fluid.specific_heat)
# print(flow_model.manifolds[0].fluid.thermal_conductivity)

# flow_model.manifolds[1].update()
# print(np.sum(flow_model.manifolds[1].fluid.mole_fraction, axis=0))
# print(flow_model.manifolds[1].fluid.liquid_mole_fraction)
# print(np.sum(flow_model.manifolds[1].fluid.gas.mole_fraction, axis=0))
# print(flow_model.manifolds[1].fluid.gas.pressure)


flow_model.update(inlet_mass_flow=2e-5)
x = ip.interpolate_1d(flow_model.manifolds[0].x)
y = flow_model.channel_vol_flow / np.average(flow_model.channel_vol_flow)
plt.plot(x, y)
plt.show()






