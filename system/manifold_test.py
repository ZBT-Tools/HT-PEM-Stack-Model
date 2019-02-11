import numpy as np
import data.global_parameters as g_par
import system.global_functions as g_func
import copy as copy
import input.geometry as geo
import input.physical_properties as phy_prop


class Manifold:

    def __init__(self, dict_manifold_const):
        self.cell_num = dict_manifold_const['cell_num']
        self.channel_num = dict_manifold_const['channel_numb']
        self.head_width = dict_manifold_const['header_width']
        self.head_height = dict_manifold_const['header_height']
        self.kf = dict_manifold_const['kf']
        self.cell_height = dict_manifold_const['cell_height']
        self.cell_ch_length = dict_manifold_const['cell_channel_length']
        self.cell_ch_ca = dict_manifold_const['cell_channel_cross_area']
        self.head_p = np.full((2, self.cell_num), dict_manifold_const['p_out'])
        self.cell_ch_d = 0.5 * (geo.channel_width * geo.channel_width)\
                         / (geo.channel_width + geo.channel_width)
        self.zeta = phy_prop.bend_pressure_loss_coefficient
        self.head_stoi = 1.5
        self.cell_mass_flow = None
        self.cell_mol_flow = None
        self.cell_temp = None
        self.cell_cp = None
        self.cell_visc = None
        self.cell_p = None
        self.cell_R_avg = None
        # Initialize scalar variables
        self.cross_area = self.head_height * self.head_width
        self.circumference = 2. * (self.head_height + self.head_width)
        self.hydraulic_diameter = 4. * self.cross_area / self.circumference
        self.cell_ref_p_drop = 0.
        self.ref_perm = 0.
        self.cell_ref_p_drop_cor = 0.
        self.p_cor_fac = 0.
        # Initialize arrays
        self.fwd_mat = np.tril(np.full((self.cell_num, self.cell_num), 1.))
        self.fwd_mat_ele = self.fwd_mat[:-1, :-1]
        self.head_mol_flow = np.full((2, self.cell_num), 0.)
        self.head_f_mass_flow = np.full((2, self.cell_num), 0.)
        self.head_g_mass_flow = np.full((2, self.cell_num), 0.)
        self.head_temp = np.full((2, self.cell_num), 0.)
        self.head_u = np.full((2, self.cell_num), 0.)
        self.head_cp = np.full((2, self.cell_num), 0.)
        self.head_r = np.full((2, self.cell_num), 0.)
        self.head_density = np.full((2, self.cell_num), 0.)
        self.head_Re = np.full((2, self.cell_num), 0.)
        self.head_fan_fri = np.full((2, self.cell_num), 0.)
        self.p_dist_fac = np.full(self.cell_num, 0.)
        self.cell_stoi = np.full(self.cell_num, 1.5)
        self.cell_mol_flow_old = np.full(self.cell_num, 0.)
        self.criteria = 0.

    def update_values(self, dict_manifold_dyn):
        self.cell_f_mass_flow = dict_manifold_dyn['f_mass_flow']
        self.cell_g_mass_flow = dict_manifold_dyn['g_mass_flow']
        self.cell_mol_flow = dict_manifold_dyn['mol_flow']
        self.cell_temp = dict_manifold_dyn['cell_temp']
        self.cell_cp = dict_manifold_dyn['cell_cp']
        self.cell_visc = dict_manifold_dyn['cell_visc']
        self.cell_p = dict_manifold_dyn['cell_p']
        self.cell_R_avg = dict_manifold_dyn['cell_r']

    def update(self):
        """
        This function coordinates the program sequence.
        """

        self.calc_header_fluid_mass_flows()
        self.calc_header_mol_flows()
        self.calc_header_heat_capacity()
        self.calc_header_temperature()
        self.calc_header_velocity()
        self.calc_header_gas_constant()
        self.calc_header_gas_density()
        self.calc_header_reynolds_numb()
        self.calc_header_fanning_friction_factor()
        self.calc_header_p_out()
        self.calc_ref_p_drop()
        self.calc_header_p_in()
        self.calc_pressure_distribution_factor()
        self.calc_new_ref_p_drop()
        self.calc_header_p_in()
        self.calc_pressure_distribution_factor()
        self.calc_new_cell_flows()
        self.calc_new_cell_stoi()
        self.calc_criteria()

    def calc_header_fluid_mass_flows(self):
        """
        This function sums up the given cell mass inlet and outlet flows
        to the total header inlet and outlet flows over the z-axis.

            Access to:
            - self.cell_mass_flow, 2-D-array, [cell inlet/outlet][cell number]

            Manipulate:
            - self.head_mass_flow, 2-D-array, [head inlet/outlet][cell number]
        """

        for q in range(2):
            self.head_f_mass_flow[q] = np.matmul(self.fwd_mat,
                                                 self.cell_f_mass_flow[q])
            self.head_g_mass_flow[q] = np.matmul(self.fwd_mat,
                                                 self.cell_g_mass_flow[q])


    def calc_header_mol_flows(self):
        """
        This function sums up the given cell molar inlet and outlet flows
        to the total header inlet and outlet flows over the z-axis.

            Access to:
            - self.cell_mol_flow, 2-D-array, [cell inlet/outlet][cell number]

            Manipulate:
            - self.head_mol_flow, 2-D-array, [head inlet/outlet][cell number]
        """

        for q in range(2):
            self.head_mol_flow[q] = np.matmul(self.fwd_mat,
                                              self.cell_mol_flow[q])

    def calc_header_heat_capacity(self):
        """
        This function mixes up the given cell outlet heat capacity
        to the total header outlet heat capacity over the z-axis.
        The given cell inlet heat capacities
        are used as the header inlet heat capacities.

            Access to:
            - self.cell_cell_cp, 2-D-array, [cell inlet/outlet][cell number]
            - self.cell_mass_flow, 2-D-array, [cell inlet/outlet][cell number]
            - self.head_mass_flow, 2-D-array, [cell inlet/outlet][cell number]

            Manipulate:
            - self.head_cp, 2-D-array, [head inlet/outlet][cell number]
        """

        self.head_cp[0] = self.cell_cp[0]
        self.head_cp[1] = np.matmul(self.fwd_mat, self.cell_f_mass_flow[1]
                                    * self.cell_cp[1]) / self.head_f_mass_flow[1]

    def calc_header_temperature(self):
        """
        This function mixes up the given cell outlet temperatures
        to the total header outlet temperatures over the z-axis.
        The given cell inlet temperatures
        are used as the header inlet temperatures.

            Access to:
            - self.cell_cell_cp, 2-D-array, [cell inlet/outlet][cell number]
            - self.cell_mass_flow, 2-D-array, [cell inlet/outlet][cell number]
            - self.head_mass_flow, 2-D-array, [cell inlet/outlet][cell number]
            - self.cell_t 2-D-array, [cell inlet/outlet][cell number]

            Manipulate:
            - self.head_t, 2-D-array, [head inlet/outlet][cell number]
        """

        self.head_temp[0] = self.cell_temp[0]
        self.head_temp[1] = np.matmul(self.fwd_mat,
                                      self.cell_f_mass_flow[1] * self.cell_cp[1] *
                                      self.cell_temp[1])\
                            / (self.head_cp[1] * self.head_f_mass_flow[1])

    def calc_header_velocity(self):
        """
        This function calculates the header inlet and outlet
        velocities over the z-axis, based on the ideal gas law

            Access to:
            - g_par.dict_uni['R], the universal gas constant
            - self.cross_area, cross area of the header channel
            - self.head_mol_flow, 2-D-array, [cell inlet/outlet][cell number]
            - self.head_t 2-D-array, [cell inlet/outlet][cell number]
            - self.head_p 2-D-array, [cell inlet/outlet][cell number]

            Manipulate:
            - self.head_u, 2-D-array, [head inlet/outlet][cell number]
        """

        for q in range(2):
            self.head_u[q] = self.head_mol_flow[q]\
                             * g_par.dict_uni['R'] * self.head_temp[q]\
                             / (self.cross_area * self.head_p[q])

    def calc_header_gas_constant(self):
        """
        This function mixes up the given cell outlet gas constant
        to the total header outlet gas constant over the z-axis.
        The given cell inlet gas constants are used as
        the header inlet gas constants.

            Access to:
            - self.cell_cell_r, 2-D-array, [cell inlet/outlet][cell number]
            - self.cell_mass_flow, 2-D-array, [cell inlet/outlet][cell number]
            - self.head_mass_flow, 2-D-array, [cell inlet/outlet][cell number]

            Manipulate:
            - self.head_r, 2-D-array, [head inlet/outlet][cell number]
        """

        self.head_r[0] = self.cell_R_avg[0]
        self.head_r[1] = np.matmul(self.fwd_mat,
                                   self.cell_g_mass_flow[1] * self.cell_R_avg[1])\
            / self.head_g_mass_flow[1]

    def calc_header_gas_density(self):
        """
        This function calculates the header gas densities over the z axis,
        based on the ideal gas law.

            Access to:
            - self.head_p, 2-D-array, [cell inlet/outlet][cell number]
            - self.head_r, 2-D-array, [cell inlet/outlet][cell number]
            - self.head_t, 2-D-array, [cell inlet/outlet][cell number]

            Manipulate:
            - self.head_density, 2-D-array, [head inlet/outlet][cell number]
        """

        for q in range(2):
            self.head_density[q] = g_func.calc_rho(self.head_p[q],
                                                   self.head_r[q],
                                                   self.head_temp[q])

    def calc_header_reynolds_numb(self):
        """
        This function calculates the header gas reynolds number.

            Access to:
            - self.head_density, 2-D-array, [cell inlet/outlet][cell number]
            - self.head_u, 2-D-array, [cell inlet/outlet][cell number]
            - self.cell_visc, 2-D-array, [cell inlet/outlet][cell number]
            - self.hydraulic_diameter, hydraulic diameter of the header channel

            Manipulate:
            - self.head_Re, 2-D-array, [head inlet/outlet][cell number]
        """
        # it should be considered to use self.head_visc
        #  here for the molar fraction and molar mass is needed
        for q in range(2):
            self.head_Re[q] = g_func.calc_Re(self.head_density[q],
                                             self.head_u[q],
                                             self.hydraulic_diameter,
                                             self.cell_visc[q])

    def calc_header_fanning_friction_factor(self):
        """
        This function calculates the header fanning friction factor.

            Access to:
            - self.head_Re, 2-D-array, [head inlet/outlet][cell number]


            Manipulate:
            - self.head_fan_fri, 2-D-array, [head inlet/outlet][cell number]
        """

        for q in range(2):
            self.head_fan_fri[q] = g_func.calc_fan_fri_fac(self.head_Re[q])

    def calc_header_p_out(self):
        """
        This function calculates the the header outlet pressure.

            Access to:
            - self.head_density, 2-D-array, [head inlet/outlet][cell number]
            - self.head_u, 2-D-array, [head inlet/outlet][cell number]
            - self.head_fan_fri, 2-D-array, [head inlet/outlet][cell number]
            - self.kf, pressure loss coefficient
            - self.cell_height, combined height of all cell_layer,
              2-D-array, [head inlet/outlet][cell number]
            - self.hydraulic_diameter, hydraulic diameter of the header channel

            Manipulate:
            self.head_p[1], 2-D-array, [head_outlet][cell_number]
        """
        drop = g_func.calc_head_p_drop(self.head_density[1][::-1][1:],
                                       self.head_u[1][::-1][:-1],
                                       self.head_u[1][::-1][1:],
                                       self.head_fan_fri[1][::-1][1:],
                                       self.kf,
                                       self.cell_height[::-1][1:],
                                       self.hydraulic_diameter)
        self.head_p[1][::-1][1:] = np.matmul(self.fwd_mat_ele, drop)\
            + self.head_p[1, -1]

    def calc_ref_p_drop(self):
        """
        This function calculates the pressure drop of the zeroth cell
        of the stack, the reference cell.

            Access to:
            - self.cell_p, 2-D-array, [cell inlet/outlet][cell number]

            Manipulate:
            self.cell_ref_p_drop, reference pressure drop in Pa
        """

        self.cell_ref_p_drop = self.cell_p[0, 0] - self.cell_p[1, 0]

    def calc_header_p_in(self):
        """
        This function calculates the pressure of the inlet header
        from the reference cell to the inlet of the inlet header.

           Access to:
           -self.cell_ref_p_drop, reference cell pressure drop
           -self.head_p, 2-D-Array, [head inlet/outlet][cell number]
           -self.head_density, 2-D-array, [head inlet/outlet][cell number]
           -self.head_u, 2-D-Array, [head inlet/outlet][cell number]
           -self.head_fan_fri, 2-D-array, [head inlet/outlet][cell number]
           -self.kf, pressure loss coefficient
           -self.cell_height, 1-a-Array, [cell number]
           self.hydraulic_diameter, hydraulic diameter of the header channel

           Manipulate:
           -self.head_p, 2-D-array, [head inlet/outlet][cell number]
        """

        self.head_p[0, 0] = self.cell_ref_p_drop + self.head_p[1, 0]
        drop = g_func.calc_head_p_drop(self.head_density[0, :-1],
                                       self.head_u[0, :-1],
                                       self.head_u[0, 1:],
                                       self.head_fan_fri[0, :-1],
                                       self.kf,
                                       self.cell_height[:-1],
                                       self.hydraulic_diameter)
        self.head_p[0, 1:] = \
            self.head_p[0, 0] + np.matmul(self.fwd_mat_ele, drop)

    def calc_pressure_distribution_factor(self):
        """
        This function calculates the pressure distributions factor.

            Access to:
            -self.head_p, 2-D-array, [head inlet/outlet][cell number]
            -self.cell_ref_p_drop, reference cell pressure drop

            Manipulate:
            - self.p_dist_fac, 1-D-Array [cell number]
        """

        self.p_dist_fac = \
            (self.head_p[0] - self.head_p[1]) / self.cell_ref_p_drop

    def calc_new_ref_p_drop(self):
        """
        This function calculates the updated cell_ref_p_drop.

            Access to:
            -self.head_mol_flow,
             2-D-array, [head inlet/outlet][cell number]
            -self.cell_visc, 2-D-array, [cell inlet/outlet][cell number]
            -self.cell_ch_length, 1-D-array [cell number]
            -self.p_cor_fac, correction factor
            -self.ref_perm, permeability of the reference cell
            -self.p_dist_fac, pressure distribution factor
            -self.cell_ch_ca, header channel cross area

            Manipulate:
            -self.cell_ref_p_drop, 1-D-array, [cell number]
        """

        self.cell_ref_p_drop = self.cell_ref_p_drop / np.sum(self.p_dist_fac) * self.cell_num

    def calc_new_cell_flows(self):
        """
        This function calculates the new inlet cell molar flows.

            Access to:
            - self.cell_mol_flow, 2-D-array, [cell inlet/outlet][cell number]
            - self.head_p,  2-D-array, [head inlet/outlet][cell number]
            -self.ref_perm, permeability of the reference cell
            -self.cell_ch_ca, header channel cross area
            -self.cell_visc, 2-D-array, [cell inlet/outlet][cell number]
            -self.cell_ch_length, 1-D-array [cell number]
            -self.p_cor_fac, correction factor

            Manipulate:
            - self.cell_mol_flow, 2-D-array, [cell inlet/outlet][cell number]
            - self.cell_mol_flow_old, 1-D-array, [cell number]
        """

        self.cell_mol_flow_old = copy.deepcopy(self.cell_mol_flow[0])
        p = (self.head_p[0] + self.head_p[1]) * .5
        pdif = self.head_p[0] - self.head_p[1]
        T = (self.cell_temp[0] + self.cell_temp[1]) * .5
        density = (self.head_density[1] + self.head_density[0]) * .5
        visc = (self.cell_visc[1]+self.cell_visc[0]) * .5
        self.zet = self.zeta * 48
        R = 8.3144598
        a = 32. * visc * self.cell_ch_length[0] / (self.cell_ch_d ** 2. * self.zeta * density)
        b = 2. * pdif / (self.zeta * density)
        u = -a + np.sqrt(a ** 2. + b)
        self.cell_mol_flow[0] = self.cell_ch_ca * p / T / R * u * self.channel_num
        print('old', self.cell_mol_flow_old )
        print('new', self.cell_mol_flow)

    def calc_new_cell_stoi(self):
        """
        This function calculates the new cell inlet stoichiometries.

            Access to:
            -self.head_stoi, stoichiometry at the header inlet
            -self.cell_num, number of cells
            -self.cell_mol_flow, 2-D-array, [cell inlet/outlet][cell number]
            -self.head_mol_flow, 2-D-array, [head inlet/outlet][cell number]

            Manipulate:
            -self.cell_stoi, 1-D-array, [cell number]
        """

        self.cell_stoi = self.cell_mol_flow[0] * self.cell_num * self.head_stoi\
            / self.head_mol_flow[0, -1]

    def calc_criteria(self):
        """
        This function calculates the convergence of the flow distribution.

            Access to:
            -self.cell_mol_flow, 2-D-array, [cell inlet/outlet][cell number]
            -self.cell_mol_flow, 1-D-array, [cell number]

            Manipulate:
            -self.criteria, convergence criteria
         """

        self.criteria = np.sum(((self.cell_mol_flow[0] - self.cell_mol_flow_old)
                                / self.cell_mol_flow_old)**2)
