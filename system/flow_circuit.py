import numpy as np
import data.global_parameters as g_par
import system.interpolation as ip
import system.global_functions as g_func
import copy as copy
import system.channel as chl
from system.output_object import OutputObject
import system.fluid2 as fluids
from abc import ABC, abstractmethod


class ParallelFlowCircuit(ABC, OutputObject):
    def __new__(cls, dict_flow_circuit, manifolds, channels,
                channel_multiplier=1.0):
        fluid = manifolds[0].fluid
        if type(fluid) is fluids.IncompressibleFluid:
            return super(ParallelFlowCircuit, cls).\
                __new__(HomogeneousFluidFlowCircuit)
        elif type(fluid) is fluids.GasMixture:
            return super(ParallelFlowCircuit, cls).\
                __new__(GasMixtureFlowCircuit)
        elif type(fluid) is fluids.TwoPhaseMixture:
            return super(ParallelFlowCircuit, cls).\
                __new__(TwoPhaseMixtureFlowCircuit)
        else:
            raise NotImplementedError

    def __init__(self, dict_flow_circuit, manifolds, channels,
                 channel_multiplier=1.0):
        super().__init__()
        self.name = dict_flow_circuit['name']
        assert isinstance(dict_flow_circuit, dict)
        assert isinstance(manifolds, (list, tuple))
        assert isinstance(channels,  (list, tuple))
        err_message = 'manifolds must be tuple or list with two objects of ' \
                      'class Channel'
        if len(manifolds) != 2:
            raise ValueError(err_message)
        elif not isinstance(manifolds[0], chl.Channel):
            raise TypeError(err_message)
        if not isinstance(channels[0], chl.Channel):
            raise TypeError(err_message)
        self.manifolds = manifolds
        self.channels = channels
        self.initialize()

        self.n_channels = len(self.channels)
        self.channel_multiplier = channel_multiplier
        self.tolerance = 1e-16
        self.max_iter = 50
        self.shape = dict_flow_circuit.get('shape', 'U')
        if self.shape not in ('U', 'Z'):
            raise ValueError('shape of flow circuit must be either U or Z')

        if self.shape == 'U':
            self.manifolds[1].flow_direction = -self.manifolds[0].flow_direction
        else:
            self.manifolds[1].flow_direction = self.manifolds[0].flow_direction

        # Distribution factor
        self.alpha = np.zeros(self.n_channels)
        self.alpha.fill(1.0)

        self.mass_flow_in = self.manifolds[0].mass_flow_total[0]
        self.vol_flow_in = 0.0
        self.channel_vol_flow = \
            self.alpha * self.mass_flow_in / self.n_channels
        self.channel_mass_flow = np.zeros(self.channel_vol_flow.shape)
        self.dp_ref = self.channels[-1].p[0] - self.channels[-1].p[-1]
        self.k_perm = np.zeros(self.n_channels)

    @abstractmethod
    def update(self, vol_flow_in=None):
        pass

    def initialize(self):
        """
        Update the flow circuit
        """
        # Calculate initial properties
        for channel in self.channels:
            channel.update()
        for manifold in self.manifolds:
            manifold.update()


class HomogeneousFluidFlowCircuit(ParallelFlowCircuit, OutputObject):
    def __init__(self, dict_flow_circuit, manifolds, channels,
                 channel_multiplier=1.0):
        super().__init__(dict_flow_circuit, manifolds, channels,
                         channel_multiplier)

    def update(self, vol_flow_in=None):
        pass


class GasMixtureFlowCircuit(ParallelFlowCircuit, OutputObject):

    def __init__(self, dict_flow_circuit, manifolds, channels,
                 channel_multiplier=1.0):
        super().__init__(dict_flow_circuit, manifolds, channels,
                         channel_multiplier)

        # self.kf = dict_flow_circuit['kf']
        self.channel_length = \
            np.asarray([channel.length for channel in channels])
        self.channel_cross_area = \
            np.asarray([channel.cross_area for channel in channels])

    def single_loop(self, inlet_mass_flow=None):
        """
        Update the flow circuit
        """
        if inlet_mass_flow is not None:
            self.mass_flow_in = inlet_mass_flow

        # Outlet header update
        channel_mass_flow_out = \
            np.array([channel.mass_flow_total[-1] for channel in self.channels])
        channel_mass_flow_out *= self.channel_multiplier

        mass_fraction = np.array([channel.fluid.mass_fraction[:, -1]
                                  for channel in self.channels]).transpose()
        mass_source = channel_mass_flow_out * mass_fraction
        #mass_source = self.channel_mass_flow * mass_fraction

        self.manifolds[1].update(mass_flow_in=0.0, mass_source=mass_source)
        print(self.manifolds[1].p - self.manifolds[1].p_out)

        # Channel update
        for i, channel in enumerate(self.channels):
            channel.p_out = ip.interpolate_1d(self.manifolds[1].p)[i]
            channel.update(mass_flow_in=
                           self.channel_mass_flow[i]/self.channel_multiplier)

        # Inlet header update
        self.manifolds[0].p_out = self.channels[-1].p[0]
        mass_fraction = np.array([channel.fluid.mass_fraction[:, 0]
                                  for channel in self.channels]).transpose()
        mass_source = -self.channel_mass_flow * mass_fraction
        self.manifolds[0].update(mass_flow_in=self.mass_flow_in,
                                 mass_source=mass_source)

        self.vol_flow_in = \
            self.mass_flow_in / self.manifolds[0].fluid.density[0]

        dp_channel = np.array([channel.p[0] - channel.p[-1]
                               for channel in self.channels])
        vol_flow_channel = np.array([np.average(channel.vol_flow)
                                     for channel in self.channels])
        visc_channel = np.array([np.average(channel.fluid.viscosity)
                                 for channel in self.channels])

        self.k_perm[:] = vol_flow_channel / dp_channel * visc_channel
        self.dp_ref = dp_channel[-1]

        p_in = ip.interpolate_1d(self.manifolds[0].p)
        p_out = ip.interpolate_1d(self.manifolds[1].p)

        print(p_in - self.manifolds[1].p_out)
        print(p_out - self.manifolds[1].p_out)

        self.alpha[:] = (p_in - p_out) / self.dp_ref
        self.dp_ref = self.vol_flow_in / np.sum(self.alpha) \
            * visc_channel[-1] / self.k_perm[-1] / self.channel_multiplier

        p_in += -self.manifolds[0].p_out + self.manifolds[1].p[-1] + self.dp_ref
        self.alpha[:] = (p_in - p_out) / self.dp_ref

        self.channel_vol_flow[:] = (p_in - p_out) * self.k_perm \
            * self.channel_multiplier / visc_channel
        density = np.array([channel.fluid.density[0] for channel in
                            self.channels])
        self.channel_mass_flow[:] = self.channel_vol_flow * density

    def update(self, inlet_mass_flow=None):
        """
        Update the flow circuit
        """
        if inlet_mass_flow is not None:
            self.mass_flow_in = inlet_mass_flow
            self.channel_mass_flow[:] = self.mass_flow_in / self.n_channels

        channel_vol_flow_old = np.zeros(self.channel_vol_flow.shape)
        channel_vol_flow_old[:] = 1e3
        for i in range(self.max_iter):
            self.single_loop(self.mass_flow_in)
            error = \
                np.sum(np.divide(self.channel_vol_flow - channel_vol_flow_old,
                                 self.channel_vol_flow,
                                 where=channel_vol_flow_old != 0.0) ** 2.0)
            #print(channel_vol_flow_old)
            #print(self.channel_vol_flow)
            channel_vol_flow_old[:] = self.channel_vol_flow
            #print(error)
            if error < self.tolerance:
                break
            if i == (self.max_iter - 1):
                print('Maximum number of iterations in update() of {} '
                      'reached'.format(self))


class TwoPhaseMixtureFlowCircuit(GasMixtureFlowCircuit, OutputObject):
    def __init__(self, dict_flow_circuit, manifolds, channels,
                 channel_multiplier=1.0):
        super().__init__(dict_flow_circuit, manifolds, channels,
                         channel_multiplier)

    def update(self, inlet_mass_flow=None):
        super().update(inlet_mass_flow)


def flow_circuit_factory(dict_fluid, dict_channel, dict_in_manifold,
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
        [fluids.fluid_factory(n_channels + 1, dict_fluid['fluid_name'],
                              species_dict=dict_fluid['fluid_components'],
                              mole_fractions=
                              dict_fluid['inlet_composition'])
         for i in range(2)]
    manifolds = [chl.Channel(dict_in_manifold, fluid[0]),
                 chl.Channel(dict_out_manifold, fluid[1])]

    flow_circuit_dict = {'name': 'Flow Circuit',
                         'shape': 'Z'}
    return ParallelFlowCircuit(flow_circuit_dict, manifolds, channels,
                               channel_multiplier=channel_multiplier)
