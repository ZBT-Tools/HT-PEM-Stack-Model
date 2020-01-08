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


def init_flow_model(dict_fluid, dict_channel, dict_in_manifold,
                    dict_out_manifold, n_channels, channel_multiplier=1.0):
    nx = g_par.dict_case['nodes']
    fluid = \
        [fluids.fluid_factory(nx, dict_fluid['fluid_name'],
                              species_dict=dict_fluid['fluid_components'],
                              mole_fractions=
                              dict_fluid['inlet_composition'])
         for i in range(n_channels)]
    channels = [chl.Channel(dict_channel, fluid[i]) for i in range(n_channels)]
    fluid = \
        [fluids.fluid_factory(n_channels, dict_fluid['fluid_name'],
                              species_dict=dict_fluid['fluid_components'],
                              mole_fractions=
                              dict_fluid['inlet_composition'])
         for i in range(2)]
    manifolds = [chl.Channel(dict_in_manifold, fluid[0]),
                 chl.Channel(dict_out_manifold, fluid[1])]

    flow_circuit_dict = {'name': 'Flow Circuit'}
    return flow_circuit.ParallelFlowCircuit(flow_circuit_dict, manifolds,
                                            channels, channel_multiplier=
                                            channel_multiplier)


n_chl = 10
n_subchl = 10

fluid_dict = input_dicts.dict_cathode_fluid
channel_dict = input_dicts.dict_cathode_channel

in_manifold_dict = {
    'name': 'Inlet Manifold',
    'channel_length': 0.1,
    'p_out': 100000.0,
    'temp_in': 350.0,
    'hum_in': 0.5,
    'flow_direction': 1,
    'channel_width': 0.010,
    'channel_height': 0.010,
    'bend_number': 0,
    'bend_friction_factor': 0.0,
    'additional_friction_fractor': 1.7
    }

out_manifold_dict = copy.deepcopy(in_manifold_dict)
out_manifold_dict['flow_direction'] = -1

flow_model = init_flow_model(fluid_dict, channel_dict, in_manifold_dict,
                             out_manifold_dict, n_chl, n_subchl)

print(flow_model.manifolds[0].fluid.mole_fraction)

print(flow_model.manifolds[0].fluid.mole_fraction)
print(flow_model.manifolds[1].fluid.mole_fraction)
temp = np.linspace(300, 400, n_chl)
flow_model.manifolds[0].temp = temp
flow_model.manifolds[1].temp = temp
flow_model.manifolds[0].update()
flow_model.manifolds[1].update()
print(np.sum(flow_model.manifolds[1].fluid.mole_fraction, axis=0))
print(flow_model.manifolds[1].fluid.liquid_mole_fraction)
print(np.sum(flow_model.manifolds[1].fluid.gas.mole_fraction, axis=0))
print(flow_model.manifolds[1].fluid.gas.pressure)





