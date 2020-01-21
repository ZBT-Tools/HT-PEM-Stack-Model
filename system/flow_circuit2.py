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
        circuit_type = dict_flow_circuit.get('type', 'Koh')
        if circuit_type == 'Koh':
            return super(ParallelFlowCircuit, cls).\
                __new__(KohFlowCircuit)
        elif circuit_type == 'Wang':
            return super(ParallelFlowCircuit, cls).\
                __new__(WangFlowCircuit)
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
        self.manifolds[0].flow_direction = 1
        self.shape = dict_flow_circuit.get('shape', 'U')
        if self.shape not in ('U', 'Z'):
            raise ValueError('shape of flow circuit must be either U or Z')
        if self.shape == 'U':
            self.manifolds[1].flow_direction = -1
        else:
            self.manifolds[1].flow_direction = 1

        self.initialize()
        self.n_channels = len(self.channels)
        self.channel_multiplier = channel_multiplier
        self.tolerance = 0.0
        self.max_iter = 10

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
        self.l_by_a = np.array([channel.length / channel.cross_area
                                for channel in self.channels])
        self.channel_length = \
            np.asarray([channel.length for channel in channels])
        self.channel_cross_area = \
            np.asarray([channel.cross_area for channel in channels])

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


class KohFlowCircuit(ParallelFlowCircuit):

    def __init__(self, dict_flow_circuit, manifolds, channels,
                 channel_multiplier=1.0):
        super().__init__(dict_flow_circuit, manifolds, channels,
                         channel_multiplier)

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
        # mass_source = self.channel_mass_flow * mass_fraction
        self.manifolds[1].update(mass_flow_in=0.0, mass_source=mass_source)

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
        velocity = np.array([np.average(channel.velocity)
                             for channel in self.channels])

        p_in = ip.interpolate_1d(self.manifolds[0].p)
        p_out = ip.interpolate_1d(self.manifolds[1].p)
        self.k_perm[:] = vol_flow_channel / dp_channel * visc_channel \
            * self.l_by_a
        self.dp_ref = dp_channel[-1]
        self.alpha[:] = (p_in - p_out) / self.dp_ref
        self.dp_ref = self.vol_flow_in / np.sum(self.alpha) * self.l_by_a \
            * visc_channel[-1] / self.k_perm[-1] / self.channel_multiplier

        p_in += -self.manifolds[0].p_out + self.manifolds[1].p[-1] + self.dp_ref
        self.alpha[:] = (p_in - p_out) / self.dp_ref

        self.channel_vol_flow[:] = (p_in - p_out) * self.k_perm / self.l_by_a \
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


class WangFlowCircuit(ParallelFlowCircuit):
    def __init__(self, dict_flow_circuit, manifolds, channels,
                 channel_multiplier=1.0):
        super().__init__(dict_flow_circuit, manifolds, channels,
                         channel_multiplier)

        self.zeta = np.zeros(self.n_channels)
        self.xsi = 1.0
        self.H = self.manifolds[0].cross_area / self.manifolds[1].cross_area
        self.M = np.sum(np.array([np.average(channel.cross_area)
                                  for channel in self.channels])) \
            / self.manifolds[0].cross_area

    def update(self, inlet_mass_flow=None):
        if inlet_mass_flow is not None:
            self.mass_flow_in = inlet_mass_flow
            self.vol_flow_in = self.mass_flow_in / self.manifolds[0].density[0]
        mfd_in = self.manifolds[0]
        mfd_out = self.manifolds[1]

        k_in_0 = 0.7
        k_out_0 = 1.0
        b_in = 0.0
        b_out = 0.0
        W_0 = self.vol_flow_in / mfd_in.cross_area

        self.zeta[:] = np.array([channel.zeta_bends * channel.n_bends
                                 for channel in self.channels]) \
            + np.array([channel.friction_factor * channel.length / channel.d_h
                        for channel in self.channels])

        self.zeta[:] += 1.0 + mfd_in.zeta_other + mfd_out.zeta_other

        for i, channel in enumerate(self.channels):
            k_in = k_in_0 + b_in * np.log(mfd_in.velocity[i] / W_0)
            k_out = k_out_0 + b_out * np.log(mfd_out.velocity[i] / W_0)
            Q_i = 2.0 / (3.0 * self.zeta[i]) * (k_in - k_out * self.H ** 2.0) \
                * self.M ** 2.0


def flow_circuit_factory(dict_circuit, dict_fluid, dict_channel,
                         dict_in_manifold, dict_out_manifold, n_channels,
                         channel_multiplier=1.0):
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

    return ParallelFlowCircuit(dict_circuit, manifolds, channels,
                               channel_multiplier=channel_multiplier)
