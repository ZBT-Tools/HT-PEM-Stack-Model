import numpy as np
import data.global_parameters as g_par
import system.global_functions as g_func
import copy as copy
import system.channel as chl
from system.output_object import OutputObject


class ParallelFlowCircuit(OutputObject):

    def __init__(self, dict_flow_circuit, manifolds, channels,
                 channel_multiplier=1.0):
        super().__init__()
        self.name = dict_flow_circuit['name']
        assert isinstance(dict_flow_circuit, dict)
        assert isinstance(manifolds, (list, tuple))
        assert isinstance(channels,  (list, tuple))
        if len(manifolds) != 2:
            raise ValueError('manifolds must be tuple or list with two objects '
                             'of class Channel')
        elif not isinstance(manifolds[0], chl.Channel):
            raise TypeError('manifolds must be tuple or list with two objects '
                            'of class Channel')
        if not isinstance(channels[0], chl.Channel):
            raise TypeError('channels must be tuple or list with objects '
                            'of class Channel')
        self.manifolds = manifolds
        self.channels = channels
        self.n_channels = len(self.channels)
        self.channel_multiplier = channel_multiplier

        # Distribution factor
        self.alpha = np.zeros(len(self.channels))
        self.alpha[:] = 1.0

        self.vol_flow_in = self.manifolds[0].vol_flow[0]
        self.channel_vol_flow = \
            self.alpha * self.vol_flow_in / self.n_channels
        # self.kf = dict_flow_circuit['kf']
        self.channel_length = \
            np.asarray([channel.length for channel in channels])
        self.channel_cross_area = \
            np.asarray([channel.cross_area for channel in channels])
        # self.head_p = np.full((2, self.n_channels), dict_flow_circuit['p_out'])
        # self.head_stoi = 1.5
        self.cell_mass_flow = None
        self.cell_mol_flow = None
        self.cell_temp = None
        self.cell_cp = None
        self.cell_visc = None
        self.cell_p = None
        self.cell_R_avg = None


        # Initialize scalar variables
        # self.cross_area = self.head_height * self.head_width
        # self.circumference = 2. * (self.head_height + self.head_width)
        # self.hydraulic_diameter = 4. * self.cross_area / self.circumference
        self.cell_ref_p_drop = 0.
        self.ref_perm = 0.
        self.cell_ref_p_drop_cor = 0.
        self.p_cor_fac = 0.
        # Initialize arrays
        self.fwd_mat = np.tril(np.full((self.n_channels, self.n_channels), 1.))
        self.fwd_mat_ele = self.fwd_mat[:-1, :-1]
        self.head_mol_flow = np.full((2, self.n_channels), 0.)
        self.head_f_mass_flow = np.full((2, self.n_channels), 0.)
        self.head_g_mass_flow = np.full((2, self.n_channels), 0.)
        self.head_temp = np.full((2, self.n_channels), 0.)
        self.head_u = np.full((2, self.n_channels), 0.)
        self.head_cp = np.full((2, self.n_channels), 0.)
        self.head_r = np.full((2, self.n_channels), 0.)
        self.head_density = np.full((2, self.n_channels), 0.)
        self.head_Re = np.full((2, self.n_channels), 0.)
        self.head_fan_fri = np.full((2, self.n_channels), 0.)
        self.p_dist_fac = np.full(self.n_channels, 0.)
        self.cell_stoi = np.full(self.n_channels, 1.5)
        self.cell_mol_flow_old = np.full(self.n_channels, 0.)
        self.criteria = 0.

    def update(self, vol_flow_in=None):
        """
        Update the flow circuit
        """
        if vol_flow_in is not None:
            self.vol_flow_in = vol_flow_in
        self.channel_vol_flow = \
            self.alpha * self.vol_flow_in / self.n_channels
        dmol_inlet_header = self.chann

    def calc_header_temperature(self):
        """
        This function mixes up the given cell outlet temperatures
        to the total header outlet temperatures over the z-axis.
        The given cell inlet temperatures
        are used as the header inlet temperatures.
        """
        self.head_temp[0] = self.cell_temp[0]
        self.head_temp[1] = np.matmul(self.fwd_mat,
                                      self.cell_f_mass_flow[1] * self.cell_cp[1] *
                                      self.cell_temp[1])\
                            / (self.head_cp[1] * self.head_f_mass_flow[1])

    def calc_ref_permeability(self):
        """"
        Calculation of the permeability of the reference cell
        and a pressure drop correction factor.
        """
        self.ref_perm = np.average(self.cell_visc[:, 0]) \
                        * self.channel_length[0] * np.average(self.cell_mol_flow[:, 0]) \
                        / (self.channel_cross_area[0] * self.cell_ref_p_drop) / self.channel_multiplier
        self.cell_ref_p_drop_cor = \
            np.average(self.cell_mol_flow[:, 0]) \
            / self.channel_multiplier * np.average(self.cell_visc[:, 0]) \
            * self.channel_length[0] / (self.channel_cross_area[0] * self.ref_perm)
        self.p_cor_fac = self.cell_ref_p_drop / self.cell_ref_p_drop_cor

    def calc_pressure_distribution_factor(self):
        """
        Calculation of the pressure distribution factor.
        """
        self.p_dist_fac = \
            (self.head_p[0] - self.head_p[1]) / self.cell_ref_p_drop

    def calc_new_ref_p_drop(self):
        """
        Calculation of the updated cell_ref_p_drop.
        """
        self.cell_ref_p_drop = self.head_mol_flow[0, -1] \
            * np.average(self.cell_visc[:, 0]) * self.channel_length[0] \
            / (self.ref_perm * np.sum(self.p_dist_fac)
               * self.channel_cross_area[0] * self.channel_multiplier)

    def calc_new_cell_flows(self):
        """
        Calculation of the new inlet cell molar flows.
        """
        self.cell_mol_flow_old = copy.deepcopy(self.cell_mol_flow[0])
        self.cell_mol_flow[0] = (self.head_p[0] - self.head_p[1]) \
            * self.ref_perm * self.channel_cross_area[0] \
            / (np.average(self.cell_visc) * self.channel_length[0] *
               self.p_cor_fac) * self.channel_multiplier

     def calc_criteria(self):
        """
        Calculation of the convergence of the flow distribution.
        """
        self.criteria = \
            np.sum(((self.cell_mol_flow[0] - self.cell_mol_flow_old)
                    / self.cell_mol_flow_old)**2)
