import numpy as np
import data.global_parameters as g_par
import system.global_functions as g_func
import copy as copy


class FlowCircuit:

    def __init__(self, dict_flow_circuit, cells):
        self.name = dict_flow_circuit['name']
        self.n_cell = dict_flow_circuit['cell_number']
        self.n_ch = dict_flow_circuit['channel_number']
        self.head_width = dict_flow_circuit['header_width']
        self.head_height = dict_flow_circuit['header_height']
        self.kf = dict_flow_circuit['kf']
        self.cell_height = dict_flow_circuit['cell_height']
        self.cell_ch_length = dict_flow_circuit['cell_channel_length']
        self.cell_ch_ca = dict_flow_circuit['cell_channel_cross_area']
        self.head_p = np.full((2, self.n_cell), dict_flow_circuit['p_out'])
        self.head_stoi = 1.5
        self.cells = cells
        self.cell_mass_flow = None
        self.cell_mol_flow = None
        self.cell_temp = None
        self.cell_cp = None
        self.cell_visc = None
        self.cell_p = None
        self.cell_R_avg = None

        # Get name of half cells for the flow circuit (Anode or Cathode)
        half_cell_name = self.name.split()[0]
        half_cell_names = [self.cells.half_cells[0].name,
                           self.cells.half_cells[1].name]
        try:
            self.hc_id = half_cell_names.index(half_cell_name)
        except ValueError:
            raise ValueError('First part of manifold name must correspond to '
                             'any of the half cells names, '
                             'typically Anode or Cathode')

        # Initialize scalar variables
        self.cross_area = self.head_height * self.head_width
        self.circumference = 2. * (self.head_height + self.head_width)
        self.hydraulic_diameter = 4. * self.cross_area / self.circumference
        self.cell_ref_p_drop = 0.
        self.ref_perm = 0.
        self.cell_ref_p_drop_cor = 0.
        self.p_cor_fac = 0.
        # Initialize arrays
        self.fwd_mat = np.tril(np.full((self.n_cell, self.n_cell), 1.))
        self.fwd_mat_ele = self.fwd_mat[:-1, :-1]
        self.head_mol_flow = np.full((2, self.n_cell), 0.)
        self.head_f_mass_flow = np.full((2, self.n_cell), 0.)
        self.head_g_mass_flow = np.full((2, self.n_cell), 0.)
        self.head_temp = np.full((2, self.n_cell), 0.)
        self.head_u = np.full((2, self.n_cell), 0.)
        self.head_cp = np.full((2, self.n_cell), 0.)
        self.head_r = np.full((2, self.n_cell), 0.)
        self.head_density = np.full((2, self.n_cell), 0.)
        self.head_Re = np.full((2, self.n_cell), 0.)
        self.head_fan_fri = np.full((2, self.n_cell), 0.)
        self.p_dist_fac = np.full(self.n_cell, 0.)
        self.cell_stoi = np.full(self.n_cell, 1.5)
        self.cell_mol_flow_old = np.full(self.n_cell, 0.)
        self.criteria = 0.

    def update_values(self, mol_flow, temp, cp, visc, p, r, mass_flow,
                      gas_mass_flow):
        self.cell_f_mass_flow = mass_flow
        self.cell_g_mass_flow = gas_mass_flow
        self.cell_mol_flow = mol_flow
        self.cell_temp = temp
        self.cell_cp = cp
        self.cell_visc = visc
        self.cell_p = p
        self.cell_R_avg = r

    def update(self):
        """
        Coordination of the program sequence.
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
        self.calc_ref_permeability()
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
        This function sums up the cell mass inlet and outlet flows
        to the total header inlet and outlet flows over the z-axis.
        """
        for q in range(2):
            self.head_f_mass_flow[q] = np.matmul(self.fwd_mat,
                                                 self.cell_f_mass_flow[q])
            self.head_g_mass_flow[q] = np.matmul(self.fwd_mat,
                                                 self.cell_g_mass_flow[q])

    def calc_header_mol_flows(self):
        """
        Sum of the cell molar inlet and outlet flows
        to the total header inlet and outlet flows over the z-axis.
        """
        for i in range(2):
            self.head_mol_flow[i] = np.matmul(self.fwd_mat,
                                              self.cell_mol_flow[i])

    def calc_header_heat_capacity(self):
        """
        This function mixes up the given cell outlet heat capacity
        to the total header outlet heat capacity over the z-axis.
        The given cell inlet heat capacities
        are used as the header inlet heat capacities.
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
        """
        self.head_temp[0] = self.cell_temp[0]
        self.head_temp[1] = np.matmul(self.fwd_mat,
                                      self.cell_f_mass_flow[1] * self.cell_cp[1] *
                                      self.cell_temp[1])\
                            / (self.head_cp[1] * self.head_f_mass_flow[1])

    def calc_header_velocity(self):
        """
        Calculation of the header inlet and outlet
        velocities over the z-axis, based on the ideal gas law
        """
        for q in range(2):
            self.head_u[q] = self.head_mol_flow[q] \
                             * g_par.constants['R'] * self.head_temp[q] \
                             / (self.cross_area * self.head_p[q])

    def calc_header_gas_constant(self):
        """
        Mixes up the given cell outlet gas constant
        to the total header outlet gas constant over the z-axis.
        The given cell inlet gas constants are used as
        the header inlet gas constants.
        """
        self.head_r[0] = self.cell_R_avg[0]
        self.head_r[1] = np.matmul(self.fwd_mat,
                                   self.cell_g_mass_flow[1]
                                   * self.cell_R_avg[1]) \
            / self.head_g_mass_flow[1]

    def calc_header_gas_density(self):
        """
        Calculation of the header gas densities over the z axis,
        based on the ideal gas law.
        """
        for q in range(2):
            self.head_density[q] = g_func.calc_rho(self.head_p[q],
                                                   self.head_r[q],
                                                   self.head_temp[q])

    def calc_header_reynolds_numb(self):
        """
        Calculation of the header gas reynolds number.
        """
        # it should be considered to use self.head_visc
        #  here for the molar fraction and molar mass is needed
        for i in range(2):
            self.head_Re[i] = \
                g_func.calc_reynolds_number(self.head_density[i],
                                            self.head_u[i],
                                            self.hydraulic_diameter,
                                            self.cell_visc[i])

    def calc_header_fanning_friction_factor(self):
        """
        Calculation of the header fanning friction factor.
        """

        for i in range(2):
            self.head_fan_fri[i] = g_func.calc_friction_factor(self.head_Re[i])

    def calc_header_p_out(self):
        """
        Calculation of the header outlet pressure.
        """
        drop = \
            g_func.calc_pressure_drop(self.head_u[1][::-1],
                                      self.head_density[1][::-1][1:],
                                      self.head_fan_fri[1][::-1][1:],
                                      self.kf,
                                      self.cell_height[::-1][1:],
                                      self.hydraulic_diameter)
        self.head_p[1][::-1][1:] = np.matmul(self.fwd_mat_ele, drop)\
            + self.head_p[1, -1]

    def calc_ref_p_drop(self):
        """
        Calculation of the the pressure drop of the 0th cell
        of the stack, the reference cell.
        """
        self.cell_ref_p_drop = self.cell_p[0, 0] - self.cell_p[1, 0]

    def calc_ref_permeability(self):
        """"
        Calculation of the permeability of the reference cell
        and a pressure drop correction factor.
        """
        self.ref_perm = np.average(self.cell_visc[:, 0]) \
            * self.cell_ch_length[0] * np.average(self.cell_mol_flow[:, 0]) \
            / (self.cell_ch_ca[0] * self.cell_ref_p_drop) / self.n_ch
        self.cell_ref_p_drop_cor = \
            np.average(self.cell_mol_flow[:, 0]) \
            / self.n_ch * np.average(self.cell_visc[:, 0]) \
            * self.cell_ch_length[0] / (self.cell_ch_ca[0] * self.ref_perm)
        self.p_cor_fac = self.cell_ref_p_drop / self.cell_ref_p_drop_cor

    def calc_header_p_in(self):
        """
        Calculation of the pressure of the inlet header
        from the reference cell to the inlet of the inlet header.
        """
        self.head_p[0, 0] = self.cell_ref_p_drop + self.head_p[1, 0]
        drop = \
            g_func.calc_pressure_drop(self.head_u[0],
                                      self.head_density[0, :-1],
                                      self.head_fan_fri[0, :-1],
                                      self.kf,
                                      self.cell_height[:-1],
                                      self.hydraulic_diameter)
        self.head_p[0, 1:] = \
            self.head_p[0, 0] + np.matmul(self.fwd_mat_ele, drop)

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
            * np.average(self.cell_visc[:, 0]) * self.cell_ch_length[0] \
            / (self.ref_perm * np.sum(self.p_dist_fac)
               * self.cell_ch_ca[0] * self.n_ch)

    def calc_new_cell_flows(self):
        """
        Calculation of the new inlet cell molar flows.
        """
        self.cell_mol_flow_old = copy.deepcopy(self.cell_mol_flow[0])
        self.cell_mol_flow[0] = (self.head_p[0] - self.head_p[1]) \
            * self.ref_perm * self.cell_ch_ca[0] / (np.average(self.cell_visc)
                                                    * self.cell_ch_length[0]
                                                    * self.p_cor_fac)\
            * self.n_ch

    def calc_new_cell_stoi(self):
        """
        Calculation of the new cell inlet stoichiometries.
        """
        self.cell_stoi = self.cell_mol_flow[0] * self.n_cell * self.head_stoi \
                         / self.head_mol_flow[0, -1]

    def calc_criteria(self):
        """
        Calculation of the convergence of the flow distribution.
        """
        self.criteria = \
            np.sum(((self.cell_mol_flow[0] - self.cell_mol_flow_old)
                    / self.cell_mol_flow_old)**2)
