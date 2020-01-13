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


n_chl = 20
n_subchl = 10

channel_dict = {
    'name': 'Channel',
    'channel_length': 0.651,
    'p_out': 101325.0,
    'temp_in': 300.0,
    'hum_in': 0.1,
    'flow_direction': 1,
    'channel_width': 0.0010,
    'channel_height': 0.0010,
    'bend_number': 48,
    'bend_friction_factor': 0.1,
    'additional_friction_fractor': 0.0
    }

fluid_dict = input_dicts.dict_cathode_fluid

in_manifold_dict = {
    'name': 'Inlet Manifold',
    'channel_length': 1.0,
    'p_out': channel_dict['p_out'],
    'temp_in': 300.0,
    'hum_in': 0.1,
    'flow_direction': 1,
    'channel_width': 0.0125,
    'channel_height': 0.0075,
    'bend_number': 0,
    'bend_friction_factor': 0.0,
    'additional_friction_fractor': 0.4
    }

out_manifold_dict = copy.deepcopy(in_manifold_dict)
out_manifold_dict['name'] = 'Outlet Manifold'

flow_circuit_dict = {
    'name': 'Flow Circuit',
    'shape': 'U'
    }

flow_model = flow_circuit.flow_circuit_factory(flow_circuit_dict, fluid_dict,
                                               channel_dict, in_manifold_dict,
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


flow_model.update(inlet_mass_flow=1e-3)
x = ip.interpolate_1d(flow_model.manifolds[0].x)
flow_model.manifolds[0].zeta_other = 1.0
flow_model.manifolds[1].zeta_other = 1.0
flow_model.update()
q1 = flow_model.channel_vol_flow / np.average(flow_model.channel_vol_flow)
flow_model.manifolds[0].zeta_other = 0.4
flow_model.manifolds[1].zeta_other = 0.4
flow_model.update()
q2 = flow_model.channel_vol_flow / np.average(flow_model.channel_vol_flow)
flow_model.manifolds[0].zeta_other = 1.0
flow_model.manifolds[1].zeta_other = 1.0
flow_model.update()
q3 = flow_model.channel_vol_flow / np.average(flow_model.channel_vol_flow)

plt.plot(x, q1)
plt.plot(x, q2)
plt.plot(x, q3)
plt.show()
p_in = ip.interpolate_1d(flow_model.manifolds[0].p)
p_out = ip.interpolate_1d(flow_model.manifolds[1].p)

plt.plot(p_in/np.average(p_in))
plt.plot(p_out/np.average(p_out))
plt.show()
i = 3
y = np.hstack((p_in[:i], flow_model.channels[i].p, p_out[:i]))
plt.plot(y)
i = 10
y = np.hstack((p_in[:i], flow_model.channels[i].p, p_out[:i]))
plt.plot(y)
i = 17
y = np.hstack((p_in[:i], flow_model.channels[i].p, p_out[:i]))
plt.plot(y)
plt.show()
print(p_in)
print(np.array([channel.p for channel in flow_model.channels]).transpose())
print(p_out)