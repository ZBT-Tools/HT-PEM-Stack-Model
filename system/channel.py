import numpy as np
import data.global_parameters as g_par
import system.fluid as fluid


class Channel:

    def __init__(self, dict_ch):
        self.length = dict_ch['channel_length']
        # channel length
        n_ele = g_par.dict_case['elements']
        n_nodes = n_ele + 1
        self.x = np.linspace(0.0, self.length, n_nodes)
        self.dx = np.diff(self.x)
        # element length
        self.p_out = dict_ch['p_out']
        # inlet pressure
        self.temp_in = dict_ch['temp_in']
        # inlet temperature
        self.humidity_in = dict_ch['hum_in']
        # inlet humidity
        self.flow_dir = dict_ch['flow_dir']
        # flow direction
        self.width = dict_ch['channel_width']
        # channel width
        self.height = dict_ch['channel_height']
        # channel height
        self.n_bends = dict_ch['bend_number']
        # number of channel bends
        self.bend_fri_fac = dict_ch['bend_fri_fac']
        # bend friction factor
        self.rib_width = dict_ch['rib_width']
        # rack width
        self.active_area = self.width * self.length
        # planar area of the channel
        self.active_area_dx = self.width * self.dx
        # planar area of an element of the channel
        self.cross_area = self.width * self.height
        # channel cross area
        self.circum = 2. * (self.width + self.height)
        # channel circumference
        self.d_h = 4. * self.cross_area / self.circum
        # channel hydraulic diameter
        self.velocity = np.zeros(n_nodes)
        # flow velocity
        self.fluid = \
            fluid.Fluid(n_ele, {'O2': 'gas', 'N2': 'gas', 'H2O': 'gas-liquid'},
                        mole_fractions_init=[0.205, 0.785, 0.01],
                        liquid_props=
                        {'H2O': fluid.LiquidProperties(1000.0, 1e-3,
                                                       4000.0, 0.2)})

    # def calc_flow_velocity(self):
    #     """
    #     Calculates the gas phase velocity.
    #     The gas phase velocity is taken to be the liquid water velocity as well.
    #
    #         Access to:
    #         -self.q_gas
    #         -self.temp_fluid
    #         -self.p
    #         -self.channel.cross_area
    #         -g_par.dict_uni['R']
    #
    #         Manipulate:
    #         -self.u
    #     """
    #     self.velocity = self.q_gas * g_par.dict_uni['R'] * self.temp_fluid \
    #                     / (self.p * self.cross_area)
