import numpy as np
import data.global_parameters as g_par
import system.interpolation as ip
import system.global_functions as g_func
import copy as copy
import system.channel as chl
from system.output_object import OutputObject
import system.fluid2 as fluids
import system.species as species
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

        if hasattr(self.manifolds[0].fluid, 'mass_fraction'):
            self.multi_component = True
        else:
            self.multi_component = False

        self.n_channels = len(self.channels)
        self.channel_multiplier = channel_multiplier
        self.tolerance = 1e-6
        self.max_iter = 5

        self.mass_flow_in = self.manifolds[0].mass_flow_total[0]
        self.vol_flow_in = 0.0
        self.channel_mass_flow = \
            np.ones(self.n_channels) * self.mass_flow_in / self.n_channels
        self.channel_vol_flow = np.zeros(self.channel_mass_flow.shape)

        self.channel_length = \
            np.asarray([channel.length for channel in channels])
        self.channel_cross_area = \
            np.asarray([channel.cross_area for channel in channels])
        self.initialize = True
        self.update_channels()

    def update(self, inlet_mass_flow=None):
        """
        Update the flow circuit
        """
        if inlet_mass_flow is not None:
            self.mass_flow_in = inlet_mass_flow
            self.channel_mass_flow[:] = self.mass_flow_in / self.n_channels
            self.vol_flow_in = \
                inlet_mass_flow / self.manifolds[0].fluid.density[0]
            self.initialize = True

        channel_vol_flow_old = np.zeros(self.channel_vol_flow.shape)
        channel_vol_flow_old[:] = 1e3
        for i in range(self.max_iter):
            self.single_loop()
            error = \
                np.sum(
                    np.divide(self.channel_vol_flow - channel_vol_flow_old,
                              self.channel_vol_flow,
                              where=channel_vol_flow_old != 0.0) ** 2.0)
            # print(channel_vol_flow_old)
            # print(self.channel_vol_flow)
            channel_vol_flow_old[:] = self.channel_vol_flow
            # print(error)
            if error < self.tolerance:
                break
            if i == (self.max_iter - 1):
                print('Maximum number of iterations n = {} with error = {} in '
                      'update() of {} '
                      'reached'.format(self.max_iter, error, self))

    @abstractmethod
    def single_loop(self, inlet_mass_flow=None):
        pass

    def update_channels(self):
        if self.initialize:
            channel_mass_flow_in = np.zeros(self.n_channels) \
                + self.mass_flow_in / self.n_channels
            channel_mass_flow_out = channel_mass_flow_in
        else:
            channel_mass_flow_in = self.channel_mass_flow
            channel_mass_flow_out = \
                np.array([channel.mass_flow_total[-1]
                          for channel in self.channels])
            channel_mass_flow_out *= self.channel_multiplier

        if self.multi_component:
            mass_fraction = np.array([channel.fluid.mass_fraction[:, -1]
                                      for channel in self.channels]).transpose()
        else:
            mass_fraction = 1.0
        mass_source = channel_mass_flow_out * mass_fraction
        # mass_source = self.channel_mass_flow * mass_fraction
        self.manifolds[1].update(mass_flow_in=0.0, mass_source=mass_source)

        # Channel update
        for i, channel in enumerate(self.channels):
            channel.p_out = ip.interpolate_1d(self.manifolds[1].p)[i]
            channel.update(mass_flow_in=
                           channel_mass_flow_in[i]/self.channel_multiplier)
        # Inlet header update
        self.manifolds[0].p_out = self.channels[-1].p[0]
        if self.multi_component:
            mass_fraction = np.array([channel.fluid.mass_fraction[:, 0]
                                      for channel in self.channels]).transpose()
        else:
            mass_fraction = 1.0
        mass_source = -self.channel_mass_flow * mass_fraction
        self.manifolds[0].update(mass_flow_in=self.mass_flow_in,
                                 mass_source=mass_source)
        self.vol_flow_in = \
            self.mass_flow_in / self.manifolds[0].fluid.density[0]


class KohFlowCircuit(ParallelFlowCircuit):

    def __init__(self, dict_flow_circuit, manifolds, channels,
                 channel_multiplier=1.0):
        super().__init__(dict_flow_circuit, manifolds, channels,
                         channel_multiplier)
        # Distribution factor
        self.alpha = np.ones(self.n_channels)
        self.dp_ref = self.channels[-1].p[0] - self.channels[-1].p[-1]
        self.k_perm = np.zeros(self.n_channels)
        self.l_by_a = np.array([channel.length / channel.cross_area
                                for channel in self.channels])

    def update_channels(self):
        super().update_channels()
        self.initialize = False

    def single_loop(self, inlet_mass_flow=None):
        """
        Update the flow circuit
        """
        if inlet_mass_flow is not None:
            self.mass_flow_in = inlet_mass_flow
        self.update_channels()
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


class WangFlowCircuit(ParallelFlowCircuit):
    def __init__(self, dict_flow_circuit, manifolds, channels,
                 channel_multiplier=1.0):
        super().__init__(dict_flow_circuit, manifolds, channels,
                         channel_multiplier)

        # self.zeta = np.zeros(self.n_channels)
        self.xsi = 1.0
        self.H = self.manifolds[0].cross_area / self.manifolds[1].cross_area
        F_c = np.array([np.average(channel.cross_area)
                        for channel in self.channels])
        self.M = g_func.add_source(F_c, F_c[1:], direction=-1) \
            / self.manifolds[0].cross_area
        # print('self.M: ', self.M)
        # self.M = np.sum(F_c) / self.manifolds[0].cross_area
        self.E = self.manifolds[0].length / self.manifolds[0].d_h
        self.D_star = self.manifolds[0].d_h / self.manifolds[1].d_h
        self.sqr_M = self.M ** 2.0
        self.sqr_H = self.H ** 2.0

        print('M = ', self.M)
        print('E = ', self.E)

    def update_channels(self):
        super().update_channels()
        if self.initialize:
            self.f_in = np.copy(self.manifolds[0].friction_factor)
            self.f_out = np.copy(self.manifolds[1].friction_factor)
        # if self.initialize:
        self.zeta = np.array([channel.zeta_bends * channel.n_bends
                                 for channel in self.channels]) \
            + np.array([np.sum(channel.friction_factor * channel.dx /
                               channel.d_h)
                        for channel in self.channels])
        self.zeta[:] += 1.0 + self.manifolds[0].zeta_other \
            + self.manifolds[1].zeta_other
        # self.zeta[:] = 10.0
        self.initialize = False

    def single_loop(self, inlet_mass_flow=None):
        if inlet_mass_flow is not None:
            self.mass_flow_in = inlet_mass_flow
            self.vol_flow_in = self.mass_flow_in \
                / self.manifolds[0].fluid.density[0]
        self.update_channels()

        mfd_in = self.manifolds[0]
        mfd_out = self.manifolds[1]

        k_in_0 = 0.6
        k_out_0 = 1.0
        b_in = 0.01
        b_out = 0.01
        W_0 = self.vol_flow_in / mfd_in.cross_area
        print('W_0: ', W_0)
        print('Re_0:', W_0 * mfd_in.fluid.density[0] * mfd_in.d_h /
              mfd_in.fluid.viscosity[0])
        print('mfd_in.velocity[:-1]: ', mfd_in.velocity[:-1])
        # mfd_in.velocity[0] = W_0

        print('zeta = ', self.zeta)

        f_in = mfd_in.friction_factor
        f_out = mfd_out.friction_factor
        # f_in = self.f_in
        # f_out = self.f_out
        print('f_in: ', f_in)
        print('f_out: ', f_out)
        # f_in[:] = 0.038
        # f_out[:] = 0.038
        k_in = k_in_0 + b_in * np.log(mfd_in.velocity[:-1] / W_0)
        k_out = k_out_0 + b_out * np.log(mfd_out.velocity[:-1] / W_0)

        Q = 2.0 / (3.0 * self.zeta) * (k_in - k_out * self.sqr_H) \
            * self.sqr_M
        R = - 0.25 * self.E * self.xsi / self.zeta \
            * (f_in + f_out * self.D_star * self.sqr_H) * self.sqr_M
        avg_R = np.average(R)
        avg_Q = np.average(Q)

        cube_Q = np.power(Q, 3.0)
        condition = np.square(R) + cube_Q
        avg_condition = np.square(avg_R) + np.power(avg_Q, 3.0)
        condition_0 = np.square(R[0]) + np.power(Q[0], 3.0)
        x = mfd_in.x / mfd_in.length
        one_third = 1.0 / 3.0
        print('avg_condition: ', avg_condition)
        print('condition: ', condition)
        w = 1.0
        for i in range(self.n_channels):
            # print('w_i: ', w)
            # k_in_i = k_in_0 + b_in * np.log(w)
            # k_out_i = k_out_0 + b_out * np.log(w * self.H)
            # Q_i = 2.0 / (3.0 * self.zeta[i]) * (
            #             k_in_i - k_out_i * self.sqr_H) * self.sqr_M
            # R_i = - 0.25 * self.E * self.xsi / self.zeta[i] \
            #     * (f_in[i] + f_out[i] * self.D_star * self.sqr_H) * self.sqr_M
            # cube_Q_i = np.power(Q_i, 3.0)
            # square_R_i = np.square(R_i)
            # condition_i = square_R_i + cube_Q_i
            # print('cube_Q_i: ', cube_Q_i)
            # print('square_R_i: ', square_R_i)
            condition_i = condition[i]
            R_i = R[i]
            Q_i = Q[i]
            cube_Q_i = cube_Q[i]
            # print('condition: ', condition_i)

            if condition_i < 0.0:
                theta = np.arccos(R_i/np.sqrt(-cube_Q_i))
                sqrt_Q = np.sqrt(-Q_i)
                r_1 = 2.0 * sqrt_Q * np.cos(theta * one_third)
                r_2 = 2.0 * sqrt_Q * np.cos((theta + 2.0*np.pi) * one_third)
                w = (np.exp(r_1 + r_2 * x[i+1]) - np.exp(r_2 + r_1 * x[i+1])) \
                    / (np.exp(r_1) - np.exp(r_2))
                print('i :', i, ', condition < 0,  w: ', w)
            elif condition_i == 0.0:
                r = - 0.5 * np.power(R_i, one_third)
                w = (1.0 - x[i+1]) * np.exp(r*x[i+1])
                print('i :', i, ', condition == 0,  w: ', w)
            else:
                sqrt_condition = np.sqrt(condition_i)
                term_1 = np.cbrt(R_i + sqrt_condition)
                term_2 = np.cbrt(R_i - sqrt_condition)
                B = term_1 + term_2
                J = term_1 - term_2
                sqrt3_J_by_2 = np.sqrt(3.0) * J * 0.5
                w = np.exp(-B * x[i+1] * 0.5) \
                    * np.sin(sqrt3_J_by_2 * (1.0 - x[i+1])) \
                    / np.sin(sqrt3_J_by_2)
                print('i :', i, ', condition > 0,  w: ', w)
            W = w * W_0
            mfd_in.velocity[i+1] = W
            mfd_out.velocity[i+1] = W * self.H \
                * mfd_in.fluid.density[i+1] / mfd_out.fluid.density[i+1]

        # print('condition: ', condition)
        mass_flow_in = \
            mfd_in.velocity * mfd_in.fluid.density * mfd_in.cross_area
        self.channel_mass_flow[:] = mass_flow_in[:-1] - mass_flow_in[1:]
        self.channel_vol_flow[:] = \
            self.channel_mass_flow / ip.interpolate_1d(mfd_in.fluid.density)
        # print('distribution: ', self.channel_vol_flow/(np.sum(
        #     self.channel_vol_flow)/self.n_channels))


def flow_circuit_factory(dict_circuit, dict_fluid, dict_channel,
                         dict_in_manifold, dict_out_manifold, n_channels,
                         channel_multiplier=1.0):
    nx = g_par.dict_case['nodes']
    fluid_name = dict_fluid['fluid_name']
    species_dict = dict_fluid.get('fluid_components', None)
    mole_fractions = dict_fluid.get('inlet_composition', None)
    temperature = dict_in_manifold['temp_in']
    pressure = dict_out_manifold['p_out']
    if species_dict is None:
        liquid_props = species.ConstantProperties('Test', density=997.0,
                                                  viscosity=9e-4,
                                                  thermal_conductivity=0.21,
                                                  specific_heat=4100.0)
    else:
        liquid_props = None
    fluid = \
        [fluids.fluid_factory(nx, fluid_name, liquid_props=liquid_props,
                              species_dict=species_dict,
                              mole_fractions=mole_fractions,
                              temperature=temperature, pressure=pressure)
         for i in range(n_channels)]
    channels = [chl.Channel(dict_channel, fluid[i]) for i in range(n_channels)]
    fluid = \
        [fluids.fluid_factory(n_channels + 1, fluid_name,
                              liquid_props=liquid_props,
                              species_dict=species_dict,
                              mole_fractions=mole_fractions,
                              temperature=temperature, pressure=pressure)
         for i in range(2)]
    manifolds = [chl.Channel(dict_in_manifold, fluid[0]),
                 chl.Channel(dict_out_manifold, fluid[1])]

    return ParallelFlowCircuit(dict_circuit, manifolds, channels,
                               channel_multiplier=channel_multiplier)
