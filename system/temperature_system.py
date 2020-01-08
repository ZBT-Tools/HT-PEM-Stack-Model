import numpy as np
from scipy import linalg as sp_la
from scipy import sparse
from scipy.sparse.linalg import spsolve
import data.global_parameters as g_par
import system.global_functions as g_func
import data.water_properties as w_prop
import system.interpolation as ip
import system.matrix_functions as mtx
# from numba import jit

np.set_printoptions(linewidth=10000, threshold=None, precision=2)


class TemperatureSystem:

    def __init__(self, temp_dict, cells):
        self.dict = temp_dict
        self.cells = cells
        # Handover
        self.n_cells = temp_dict['cell_number']
        # cell number
        self.n_nodes = temp_dict['nodes']
        # node number
        self.n_ele = self.n_nodes - 1
        # element number
        n_layer = 5
        # layer number
        ch_length = temp_dict['channel_length']
        # channel length
        ch_width = temp_dict['channel_width']
        # channel width
        self.cool_ch_bc = temp_dict['cool_ch_bc']
        # coolant geometry condition
        self.n_cool = self.n_cells
        if self.cool_ch_bc:
            self.n_cool += 1
        # number of coolant channels
        self.temp_gas_in = temp_dict['temp_gas_in']
        # gas inlet temperature
        temp_cool_in = temp_dict['cool_temp_in']
        # coolant inlet temperature
        temp_layer_init = temp_dict['temp_layer_init']
        # initial temperature
        cp_cool = temp_dict['cool_cp']
        # coolant heat capacity
        m_flow_cool = temp_dict['cool_m_flow']
        # mass flow of the coolant
        rho_cool = temp_dict['cool_density']
        # density of the coolant
        visc_cool = temp_dict['cool_visc']
        # viscosity of the coolant
        height_cool = temp_dict['channel_height']
        # height of the coolant channel
        width_cool = temp_dict['channel_width']
        # width of the coolant channel
        n_cool_cell = temp_dict['cool_ch_numb']
        # number of coolant channels
        # end plate heat power
        self.lambda_cool = temp_dict['cool_lambda']
        # thermal conductivity from the channel to the coolant
        self.temp_amb = temp_dict['temp_amb']
        # ambient temperature
        self.e_tn = g_par.dict_case['v_tn']
        # thermodynamic neutral cell potential
        self.e_0 = g_par.dict_case['e_0']
        # open circuit potential

        """General values"""
        self.mtx_const = None
        # over the iterations constant conductance matrix
        self.mtx = None
        # over the iterations dynamic conductance matrix
        #self.k_gas_ch = np.full((2, self.n_cells, self.n_ele), 0.)
        # conductance of the species flow to the species channels
        # 0: cathode channels, 1: anode channels

        self.temp_layer_vec = \
            np.hstack([cell.heat_rhs for cell in self.cells])
        # unsorted result layer temperature vector
        self.rhs = np.zeros(np.shape(self.temp_layer_vec))
        # right side of the matrix system: mat T = rhs,
        # contains the power sources and explicit coupled terms
        self.temp_fluid = np.zeros((2, self.n_cells, self.n_nodes))
        self.temp_fluid[:] = self.temp_gas_in[0]
        self.temp_fluid[1, :] = self.temp_gas_in[1]
        self.temp_fluid_ele = np.zeros((2, self.n_cells, self.n_ele))
        self.temp_fluid_ele[:] = self.temp_gas_in[0]
        self.temp_fluid_ele[1, :] = self.temp_gas_in[1]
        # temperature of the fluid 0: cathode fluids, 1: anode fluids
        self.heat_fluid = np.zeros(np.shape(self.temp_fluid_ele))
        # heat transferred to the fluids
        #self.g_fluid = np.full((2, self.n_cells, self.n_ele), 0.)
        # heat capacity flow of the fluid 0: cathode fluids, 1: anode fluids
        #self.cond_rate = np.full((2, self.n_cells, self.n_nodes), 0.)
        # condensation rate of the water in the channels
        # 0: cathode fluids, 1: anode fluids
        #self.i = np.full((self.n_cells, self.n_ele), 0.)
        # electrical current in z-direction
        #self.v_loss = np.full((2, self.n_cells, self.n_ele), 0.)
        # voltage loss at the cathode side:0 and at the anode side:1
        #self.omega = np.full((self.n_cells, self.n_ele), 0.)
        # electrical resistance of the membrane
        self.g_cool = cp_cool * m_flow_cool * n_cool_cell
        # coolant heat capacity flow

        """Calculating the coolant to channel thermal conductance"""
        pr_ch = visc_cool * cp_cool / self.lambda_cool
        # prandtl number
        d_h_cool = 2. * width_cool * height_cool \
            / (width_cool + height_cool)
        # hydraulic diameter of the coolant channel
        u_ch = m_flow_cool / (width_cool * height_cool * rho_cool)
        # velocity of the coolant flow
        re_ch = rho_cool * u_ch * d_h_cool / visc_cool
        # reynolds number in the coolant channel

        nu_1 = 3.66
        nu_2 = 1.66 * np.sqrt(re_ch * pr_ch * d_h_cool / ch_length)
        nu_3 = (2. / (1. + 22. * pr_ch)) ** (1. / 6.) \
            * np.sqrt(re_ch * pr_ch * d_h_cool / ch_length)
        nu_lam = (nu_1 ** 3. + 0.7 ** 3. + (nu_2 - 0.7) ** 3.
                  + nu_3 ** 3.) ** (1. / 3.)
        # laminar nusselt number
        zeta = (1.8 * np.log(re_ch) - 1.5) ** -2.
        nu_turb = zeta / 8. * re_ch * pr_ch \
            / (1. + 12.7 * np.sqrt(zeta / 8.)
                * (pr_ch ** 2. / 3.) - 1.) \
            * (1. + (d_h_cool / ch_length) ** 2. / 3.)
        if re_ch <= 2300.:
            nu_ch = nu_lam
        elif 2300. < re_ch < 1.e4:
            gamma = (re_ch - 2300.) / 7700.
            nu_ch = (1. - gamma) * nu_lam + gamma * nu_turb
        else:
            nu_ch = nu_turb
        # turbulent nusselt number

        conv_coeff_ch = nu_ch * self.lambda_cool / d_h_cool
        # convection coefficient between the coolant and the channel wall
        conv_area = d_h_cool * np.pi * ch_length / self.n_ele
        # convection area of the channel wall
        self.k_cool = conv_coeff_ch * conv_area * n_cool_cell
        # thermal conductance between the element channel area and the coolant

        """Building up the result temperature list and arrays"""
        # List of references to cell layer temperatures
        self.temp_layer = [cell.temp_layer for cell in self.cells]

        # Coolant temperature array (cell, element
        self.temp_cool = np.zeros((self.n_cool, self.n_nodes))
        self.temp_cool[:] = temp_cool_in
        self.temp_cool_ele = np.zeros((self.n_cool, self.n_ele))
        self.temp_cool_ele[:] = temp_cool_in

        """Building up the base conductance matrix mat"""
        self.dyn_vec = np.zeros(np.shape(self.temp_layer_vec))
        # vector for dynamically changing heat conductance values

        # Add constant source and sink coefficients to heat conductance matrix
        # Heat transfer to ambient
        alpha_amb = temp_dict['alpha_amb']

        for cell in self.cells:
            cell.k_amb = cell.calc_ambient_conductance(alpha_amb)
            if cell.last_cell:
                k_amb_vector = cell.k_amb.transpose().flatten()
            else:
                k_amb_vector = cell.k_amb[:-1].transpose().flatten()

            cell.add_implicit_layer_source(cell.heat_mtx_const, -k_amb_vector)
            cell.add_explicit_layer_source(cell.heat_rhs_const,
                                           k_amb_vector * self.temp_amb)
            # Heat transfer to coolant channels
            if cell.first_cell:
                if self.cool_ch_bc:
                    cell.add_implicit_layer_source(cell.heat_mtx_const,
                                                   -self.k_cool, layer_id=0)
            else:
                cell.add_implicit_layer_source(cell.heat_mtx_const,
                                               -self.k_cool, layer_id=0)
            # mtx_2.append(cell.heat_mtx.copy())
        if self.cool_ch_bc:
            cell = self.cells[-1]
            cell.add_implicit_layer_source(cell.heat_mtx_const, -self.k_cool,
                                           layer_id=-1)
        self.rhs_const = \
            np.hstack([cell.heat_rhs_const for cell in self.cells])

        self.index_list, self.layer_index_list = \
            mtx.create_index_lists(self.cells)

        self.mtx_const = self.connect_cells()
        self.sparse_solve = True
        if self.sparse_solve:
            self.mtx_const = sparse.csr_matrix(self.mtx_const)

    def connect_cells(self):
        matrix = sp_la.block_diag(*[cell.heat_mtx_const for cell in self.cells])
        cell_ids = np.asarray([list(range(self.n_cells-1)),
                               list(range(1, self.n_cells))]).transpose()
        layer_ids = np.asarray([(-1, 0) for i in range(self.n_cells-1)])
        conductance = \
            np.asarray([self.cells[i].thermal_conductance_z[layer_ids[i][0]]
                        for i in range(self.n_cells-1)])
        mtx.connect_cells(matrix, cell_ids, layer_ids,
                          conductance, self.index_list)
        return matrix

    def update(self):
        """
        This function coordinates the program sequence
        """
        # self.change_value_shape()
        self.update_gas_channel()
        self.update_coolant_channel()
        self.update_temp_layer()

    def update_temp_layer(self):
        """
        This function coordinates the temp_layer program sequence
        """
        self.update_matrix()
        self.update_rhs()
        self.solve_system()
        self.update_cell_layer_temperatures()

    @staticmethod
    def element_heat_transfer(wall_temp, fluid_temp, g_fluid, k_fluid):
        if np.isscalar(wall_temp):
            wall_temp = np.asarray([wall_temp])
        if (len(wall_temp) + 1) != len(fluid_temp):
            raise ValueError('fluid temperature array must be node-based '
                             'and wall temperature element based')
        # Heat transfer guess based on initial temperatures
        avg_wall_temp = np.average(wall_temp)
        avg_fluid_temp = np.average(fluid_temp)
        avg_fluid_temp = np.where(fluid_temp[0] <= avg_wall_temp,
                                  np.minimum(avg_fluid_temp, avg_wall_temp),
                                  np.maximum(avg_fluid_temp, avg_wall_temp))
        avg_k_fluid = np.average(k_fluid)
        avg_heat = avg_k_fluid * (avg_wall_temp - avg_fluid_temp)
        avg_g_fluid = np.average(g_fluid)
        delta_fluid_temp = avg_heat / np.average(g_fluid)
        fluid_temp_in = fluid_temp[0]
        fluid_temp_out = fluid_temp_in + delta_fluid_temp
        fluid_temp_out = np.where(fluid_temp_in <= wall_temp,
                                  np.minimum(fluid_temp_out, avg_wall_temp),
                                  np.maximum(fluid_temp_out, avg_wall_temp))

        avg_heat = avg_g_fluid * (fluid_temp_out - fluid_temp_in) \
            / len(wall_temp)
        return fluid_temp_out, avg_heat

    def channel_heat_transfer(self, wall_temp, fluid_temp, g_fluid,
                              k_ht_coeff, flow_direction):
        if (len(wall_temp) + 1) != len(fluid_temp):
            raise ValueError('fluid temperature array must be node-based '
                             'and wall temperature element based')
        heat_fluid = np.zeros(np.shape(wall_temp))
        if np.isscalar(g_fluid):
            g_fluid = np.full_like(wall_temp, g_fluid)
        if np.isscalar(k_ht_coeff):
            k_ht_coeff = np.full_like(wall_temp, k_ht_coeff)
        n = len(wall_temp)
        if flow_direction == 1:
            for i in range(n):
                fluid_temp[i+1], heat_fluid[i] = \
                    self.element_heat_transfer(wall_temp[i],
                                               [fluid_temp[i], fluid_temp[i+1]],
                                               g_fluid[i], k_ht_coeff[i])
        elif flow_direction == -1:
            for i in reversed(range(n)):
                fluid_temp[i], heat_fluid[i] = \
                    self.element_heat_transfer(wall_temp[i],
                                               [fluid_temp[i+1], fluid_temp[i]],
                                               g_fluid[i], k_ht_coeff[i])
        else:
            raise ValueError('flow_direction must be 1 or -1')

    def update_gas_channel(self):
        """
        Calculates the fluid temperatures in the anode and cathode channels
        """
        for i, cell in enumerate(self.cells):
            self.channel_heat_transfer(cell.temp_layer[1],
                                       cell.cathode.temp_fluid,
                                       cell.cathode.g_fluid,
                                       cell.cathode.k_ht_coeff, 1)
            # for j in range(self.n_ele):
            #     self.temp_fluid[0, i, j+1], self.heat_fluid[0, i, j] = \
            #         self.element_heat_transfer(self.temp_layer[i][1, j],
            #                                    [self.temp_fluid[0, i, j],
            #                                     self.temp_fluid[0, i, j+1]],
            #                                    cell.cathode.g_fluid[j],
            #                                    cell.cathode.k_ht_coeff[j])
            cell.cathode.temp_fluid[0] = cell.cathode.channel.temp_in
            cell.cathode.temp_fluid_ele = \
                ip.interpolate_1d(cell.cathode.temp_fluid)

            self.channel_heat_transfer(cell.temp_layer[4],
                                       cell.anode.temp_fluid,
                                       cell.anode.g_fluid,
                                       cell.anode.k_ht_coeff, -1)
            # for j in reversed(range(self.n_ele)):
            #     self.temp_fluid[1, i, j], self.heat_fluid[1, i, j] = \
            #         self.element_heat_transfer(self.temp_layer[i][4, j],
            #                                    [self.temp_fluid[1, i, j+1],
            #                                     self.temp_fluid[1, i, j]],
            #                                    cell.anode.g_fluid[j],
            #                                    cell.anode.k_ht_coeff[j])
            cell.anode.temp_fluid[0] = cell.anode.channel.temp_in
            cell.anode.temp_fluid_ele = \
                ip.interpolate_1d(cell.anode.temp_fluid)

        # self.temp_fluid[0] = np.array([cell.cathode.temp_fluid
        #                                for cell in self.cells])
        # self.temp_fluid[1] = np.array([cell.anode.temp_fluid
        #                                for cell in self.cells])
        # self.temp_fluid[0, :, 0] = self.temp_gas_in[0]
        # self.temp_fluid_ele[0] = \
        #     ip.interpolate_along_axis(self.temp_fluid[0], axis=1)
        # self.temp_fluid[1, :, -1] = self.temp_gas_in[1]
        # self.temp_fluid_ele[1] = \
        #     ip.interpolate_along_axis(self.temp_fluid[1], axis=1)

    def update_coolant_channel(self):
        """
            Calculates the coolant channel temperatures.
        """
        for i, cell in enumerate(self.cells):
            self.channel_heat_transfer(cell.temp_layer[0],
                                       self.temp_cool[i],
                                       self.g_cool,
                                       self.k_cool, 1)
            # dtemp = self.k_cool / self.g_cool \
            #     * (self.temp_layer[i][0, :] - self.temp_cool_ele[i])
            # g_func.add_source(self.temp_cool[i], dtemp, 1)
            self.temp_cool_ele[i] = ip.interpolate_1d(self.temp_cool[i])
        if self.cool_ch_bc:
            cell = self.cells[-1]
            self.channel_heat_transfer(cell.temp_layer[-1],
                                       self.temp_cool[-1],
                                       self.g_cool,
                                       self.k_cool, 1)
            # dtemp = self.k_cool / self.g_cool \
            #     * (self.temp_layer[-1][-1, :] - self.temp_cool_ele[-1])
            # g_func.add_source(self.temp_cool[-1], dtemp, 1)
            self.temp_cool_ele[-1] = ip.interpolate_1d(self.temp_cool[-1])

    # @jit(nopython=True)
    def update_rhs(self):
        """
        Creates a vector with the right hand side entries,
        add explicit heat sources here.
        Sources from outside the system
        to the system must be defined negative.
        """
        for i, cell in enumerate(self.cells):
            cell.heat_rhs_dyn.fill(0.0)
            # Cathode bpp-gde source
            h_vap = w_prop.water.calc_h_vap(cell.cathode.temp_fluid[:-1])
            source = \
                cell.cathode.k_ht_coeff * cell.cathode.temp_fluid_ele \
                + h_vap * cell.cathode.cond_rate[:-1]
            cell.add_explicit_layer_source(cell.heat_rhs_dyn, source, 1)

            # Cathode gde-mem source
            current = cell.i_cd * cell.active_area_dx
            ohmic_heat_membrane = 0.5 * cell.omega * np.square(current)
            source = \
                (self.e_tn - self.e_0 + cell.cathode.v_loss) * current \
                + ohmic_heat_membrane
            cell.add_explicit_layer_source(cell.heat_rhs_dyn, source, 2)
            # Anode gde-mem source
            source = \
                cell.anode.v_loss * current \
                + ohmic_heat_membrane
            cell.add_explicit_layer_source(cell.heat_rhs_dyn, source, 3)
            # Anode bpp-gde source
            h_vap = w_prop.water.calc_h_vap(cell.anode.temp_fluid[:-1])
            source = h_vap * cell.anode.cond_rate[:-1] \
                + cell.anode.k_ht_coeff * cell.anode.temp_fluid_ele
            cell.add_explicit_layer_source(cell.heat_rhs_dyn, source, 4)

            # Cooling channels
            if cell.first_cell:
                if self.cool_ch_bc:
                    source = self.k_cool * self.temp_cool_ele[i]
                    cell.add_explicit_layer_source(cell.heat_rhs_dyn, source, 0)
            else:
                source = self.k_cool * self.temp_cool_ele[i]
                cell.add_explicit_layer_source(cell.heat_rhs_dyn, source, 0)
        if self.cool_ch_bc:
            source = self.k_cool * self.temp_cool_ele[-1]
            cell = self.cells[-1]
            cell.add_explicit_layer_source(cell.heat_rhs_dyn, source, -1)
        rhs_dyn = np.hstack([cell.heat_rhs_dyn for cell in self.cells])
        self.rhs = self.rhs_const + rhs_dyn

    # @jit(nopython=True)
    def update_matrix(self):
        """
        Updates the thermal conductance matrix
        """
        source_vectors = []
        for i, cell in enumerate(self.cells):
            cell.heat_mtx_dyn.fill(0.0)
            source_vectors.append(np.zeros(np.shape(cell.heat_rhs_dyn)))
            source = -cell.cathode.k_ht_coeff
            matrix, source_vec_1 = \
                cell.add_implicit_layer_source(cell.heat_mtx_dyn, source, 1)
            source = -cell.anode.k_ht_coeff
            matrix, source_vec_2 = \
                cell.add_implicit_layer_source(cell.heat_mtx_dyn, source, 4)
            source_vectors[i] = source_vec_1 + source_vec_2
        dyn_vec = np.hstack(source_vectors)
        if self.sparse_solve:
            self.mtx = \
                self.mtx_const + sparse.diags([dyn_vec], [0], format='csr')
        else:
            self.mtx = self.mtx_const + np.diag(dyn_vec)

    def solve_system(self):
        """
        Solves the layer temperatures.
        """
        if self.sparse_solve:
            self.temp_layer_vec[:] = spsolve(self.mtx, self.rhs)
        else:
            self.temp_layer_vec[:] = np.linalg.tensorsolve(self.mtx, self.rhs)

    def update_cell_layer_temperatures(self):
        """
        From 1D temperature vector to 2D cell temperature arrays
        """
        for i, cell in enumerate(self.cells):
            for j in range(cell.n_layer):
                index_vector = cell.index_array[j]
                cell.temp_layer[j] = self.temp_layer_vec[index_vector]

