import numpy as np
from scipy import linalg as sp_la
from scipy import sparse
from scipy.sparse.linalg import spsolve
import data.global_parameters as g_par
import system.global_functions as g_func
import data.water_properties as w_prop
import system.interpolation as ip
from numba import jit

np.set_printoptions(linewidth=10000, threshold=None, precision=2)


class TemperatureSystem:

    def __init__(self, temp_dict):
        self.dict = temp_dict
        # Handover
        self.n_cells = temp_dict['cell_numb']
        # cell number
        self.n_nodes = temp_dict['nodes']
        # node number
        self.n_ele = self.n_nodes - 1
        # element number
        self.n_layer = 5
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
        # heat conductance from the channel to the coolant
        self.k_layer = temp_dict['k_layer']
        # heat conductance array through and along the control volume
        self.k_alpha_env = temp_dict['k_alpha_env']
        # heat conductance from control volume to the environment
        self.temp_env = g_par.dict_case['temp_env']
        # environment temperature
        self.v_tn = g_par.dict_case['v_tn']
        # thermodynamic neutral cell voltage

        """General values"""
        self.mat_const = None
        # over the iterations constant conductance matrix
        self.mat_dyn = None
        # over the iterations dynamic conductance matrix
        self.k_gas_ch = np.full((2, self.n_cells, self.n_ele), 0.)
        # conductance of the species flow to the species channels
        # 0: cathode channels, 1: anode channels
        self.temp_layer_vec = \
            np.full(self.n_ele * (self.n_layer * self.n_cells + 1), 0.)
        # unsorted result layer temperature vector
        self.rhs = np.zeros_like(self.temp_layer_vec)
        # right side of the matrix system: mat T = rhs,
        # contains the power sources and explicit coupled terms
        self.temp_fluid = np.full((2, self.n_cells, self.n_nodes),
                                  self.temp_gas_in[0])
        self.temp_fluid_ele = np.full((2, self.n_cells, self.n_ele), 0.)
        # temperature of the fluid 0: cathode fluids, 1: anode fluids
        self.g_fluid = np.full((2, self.n_cells, self.n_ele), 0.)
        # heat capacity flow of the fluid 0: cathode fluids, 1: anode fluids
        self.cond_rate = np.full((2, self.n_cells, self.n_nodes), 0.)
        # condensation rate of the water in the channels
        # 0: cathode fluids, 1: anode fluids
        self.i = np.full((self.n_cells, self.n_ele), 0.)
        # electrical current in z-direction
        self.v_loss = np.full((2, self.n_cells, self.n_ele), 0.)
        # voltage loss at the cathode side:0 and at the anode side:1
        self.omega = np.full((self.n_cells, self.n_ele), 0.)
        # electrical resistance of the membrane
        self.pos_cat_ch = None
        # coordinates of the cathode channel heat conductance
        self.pos_ano_ch = None
        # coordinates of the anode channel heat conductance
        self.g_cool = cp_cool * m_flow_cool * n_cool_cell
        # coolant heat capacity flow
        self.k_cool = None
        # heat conductance between the coolant and the channel wall

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
        #print('Coolant channel convection coefficient:',
           #   conv_coeff_ch, 'W/(m²K)')
        conv_area = d_h_cool * np.pi * ch_length / self.n_ele
        # convection area of the channel wall
        self.k_cool = conv_coeff_ch * conv_area * n_cool_cell
        # thermal conductance between the element channel area and the coolant

        """Building up the result temperature list and arrays"""
        temp_layer = np.full((self.n_layer, self.n_ele), temp_layer_init)
        temp_layer_n = np.full((self.n_layer + 1, self.n_ele), temp_layer_init)
        self.temp_layer = []
        for q in range(self.n_cells - 1):
            self.temp_layer.append(temp_layer)
        self.temp_layer.append(temp_layer_n)
        # layer temperature list cell, layer, element
        #temp_cool_out = temp_cool_in + op_con.tar

        self.temp_cool = np.full((self.n_cool, self.n_nodes), temp_cool_in)
        self.temp_cool_ele = np.full((self.n_cool, self.n_ele), temp_cool_in)
        # coolant temperature array cell, element

        """Building up the base conductance matrix mat"""
        mat_base = np.full((self.n_layer, self.n_layer), 0.)
        mat_base[0, 0] = - self.k_layer[0, 2, 0]
        mat_base[0, 1] = self.k_layer[0, 2, 0]
        mat_base[1, 0] = self.k_layer[0, 2, 0]
        mat_base[1, 1] = - self.k_layer[0, 2, 0] - self.k_layer[0, 1, 0]
        mat_base[1, 2] = + self.k_layer[0, 1, 0]
        mat_base[2, 1] = + self.k_layer[0, 1, 0]
        mat_base[2, 2] = - self.k_layer[0, 1, 0] - self.k_layer[0, 0, 0]
        mat_base[2, 3] = + self.k_layer[0, 0, 0]
        mat_base[3, 2] = + self.k_layer[0, 0, 0]
        mat_base[3, 3] = - self.k_layer[0, 0, 0] - self.k_layer[0, 1, 0]
        mat_base[3, 4] = + self.k_layer[0, 1, 0]
        mat_base[4, 3] = + self.k_layer[0, 1, 0]
        mat_base[4, 4] = - self.k_layer[0, 1, 0]
        # heat conductance matrix in z-direction for the cells 0-(n-1)
        mat_n = np.full((self.n_layer + 1, self.n_layer + 1), 0.)
        mat_n[0:5, 0:5] = mat_base
        mat_n[4, 4] -= self.k_layer[0, 2, 0]
        mat_n[4, 5] = self.k_layer[0, 2, 0]
        mat_n[5, 4] = self.k_layer[0, 2, 0]
        mat_n[5, 5] = -self.k_layer[0, 2, 0]
        # heat conductance matrix in z-direction for the cell n

        list_mat = []
        for i in range(self.n_cells - 1):
            for j in range(self.n_ele):
                list_mat.append(mat_base)
        for i in range(self.n_ele):
            list_mat.append(mat_n)
        # list of all the heat conductance matrix in z-direction
        # for all cells and all elements
        # for i in range(len(list_mat)):
        #     print(list_mat[i])
        self.mat_const = sp_la.block_diag(*list_mat)
        # uncoupled heat conductance matrix in z-direction
        self.dyn_vec = np.zeros_like(self.temp_layer_vec)
        # vector for dynamically changing heat conductance values

        """Setting the coolant channel heat conductance"""
        cool_pos_n_up = \
            np.arange(self.n_ele * (self.n_cells - 1) * self.n_layer,
                      self.n_ele * (self.n_layer * self.n_cells + 1),
                      self.n_layer + 1)
        # upper cool ch pos for the n cell

        if self.cool_ch_bc:
            cool_pos_base = \
                np.arange(0, self.n_ele * (self.n_cells - 1) * self.n_layer,
                          self.n_layer)
            # cool ch pos for the 0-(n-1) cell
            cool_pos_n_down = \
                np.arange(self.n_ele * (self.n_cells - 1) * self.n_layer +
                          self.n_layer,
                          self.n_ele * (self.n_layer * self.n_cells + 1),
                          self.n_layer + 1)
            # lower cool ch pos for the n cell
            for q, item in enumerate(cool_pos_n_up):
                self.mat_const[item, item] -= self.k_cool

        else:
            cool_pos_base = \
                np.arange(self.n_ele,
                          self.n_ele * (self.n_cells - 1) * self.n_layer,
                          self.n_layer)
            # cool ch pos for the 1-(n-1) cell
        for q, item in enumerate(cool_pos_base):
            self.mat_const[item, item] -= self.k_cool
        for q, item in enumerate(cool_pos_n_down):
            self.mat_const[item, item] -= self.k_cool

        """Setting the x-axis heat conductance"""
        x_con_base = np.array([self.k_layer[1, 2, 0],
                               self.k_layer[1, 1, 0],
                               self.k_layer[1, 0, 0],
                               self.k_layer[1, 0, 0],
                               self.k_layer[1, 1, 0]])
        # heat conductance vec for one element of the 0-(n-1) cell
        x_con_n = np.hstack((x_con_base, 0.5 * self.k_layer[1, 2, 0]))
        # heat conductance vec for one element of the n cell
        x_con_zero = np.hstack((0.5 * self.k_layer[1, 2, 0], x_con_base[1:]))
        # heat conductance vec for one element of the n cell
        x_con_side_base = np.hstack((
            np.hstack((np.tile(x_con_zero, self.n_ele - 1),
                       np.zeros(self.n_layer))),
            np.tile(np.hstack((np.tile(x_con_base, self.n_ele - 1),
                               np.zeros(self.n_layer))), self.n_cells - 2),
            np.zeros((self.n_layer + 1) * (self.n_ele - 1) + 1)))
        # heat conductance vec for the right and left side
        # of the matrix main diagonal for the cells 0-(n-1)
        x_con_side_n = \
            np.hstack((np.zeros((self.n_cells - 1) * self.n_ele *
                                 self.n_layer),
                       np.tile(x_con_n, self.n_ele - 1)))
        # heat conductance vec for the right and left side
        # of the matrix main diagonal for the cell n
        x_con_mid = \
            np.hstack((np.hstack((x_con_zero,
                                  *np.tile(2. * x_con_zero, self.n_ele - 2),
                                  x_con_zero)),
                       np.tile(np.hstack((x_con_base,
                                          *np.tile(2. * x_con_base,
                                                   self.n_ele - 2),
                                          x_con_base)),
                               self.n_cells - 2),
                       np.hstack((x_con_n, *np.tile(2. * x_con_n,
                                                    self.n_ele - 2),
                                  x_con_n))))
        # heat conductance vec for
        # the diagonal of the main heat conductance matrix for the cells 0-n
        self.mat_const = self.mat_const \
            - np.diag(x_con_mid) \
            + np.diag(x_con_side_base, self.n_layer) \
            + np.diag(x_con_side_base, -self.n_layer) \
            + np.diag(x_con_side_n, self.n_layer + 1) \
            + np.diag(x_con_side_n, -(self.n_layer + 1))

        """Setting the cell connecting heat conductance, z-direction"""
        pos_r, pos_c = [], []
        for ct in range(self.n_ele * (self.n_cells - 1)):
            pos_r.append(4 + self.n_layer * ct)
            if ct <= self.n_ele * (self.n_cells - 2):
                pos_c.append(self.n_layer * self.n_ele + self.n_layer * ct)
            else:
                pos_c.append(pos_c[-1] + self.n_layer + 1)
        for ct, item in enumerate(pos_c):
            self.mat_const[pos_c[ct], pos_r[ct]] += self.k_layer[0, 2, 0]
            self.mat_const[pos_r[ct], pos_c[ct]] += self.k_layer[0, 2, 0]
        # heat conductance outside the main diagonal

        pos_base_out = \
            np.arange(4, self.n_ele * (self.n_cells - 1) * self.n_layer,
                      self.n_layer)
        # coordinates of the heat conductance
        # of the last layer of the elements for the cells 0-(n-1)
        pos_base_in = \
            np.arange(self.n_ele * self.n_layer,
                      self.n_ele * (self.n_cells - 1) * self.n_layer,
                      self.n_layer)
        # coordinates of the heat conductance
        # of the first layer of the cells 1-(n-1)
        pos_n_in = np.arange((self.n_cells - 1) * self.n_ele * 5,
                             self.n_ele * (self.n_layer * self.n_cells + 1),
                             self.n_layer + 1)
        # coordinates of the heat conductance
        # of first layer of the elements for the last cell
        pos_in = np.hstack((pos_base_in, pos_n_in))
        # coordinates inside the main diagonal
        for i, item in enumerate(pos_in):
            self.mat_const[item, item] -= self.k_layer[0, 2, 0]
            self.mat_const[pos_base_out[i],
                           pos_base_out[i]] -= self.k_layer[0, 2, 0]

        """Adding the environment heat conductance"""
        env_con_base = np.array([-self.k_alpha_env[0, 2, 0],
                                 -self.k_alpha_env[0, 1, 0],
                                 -self.k_alpha_env[0, 0, 0],
                                 -self.k_alpha_env[0, 0, 0],
                                 -self.k_alpha_env[0, 1, 0]])
        # environment heat conductance for the cells 1-(n-1)
        env_con_n = np.hstack((env_con_base,
                               - .5 * self.k_alpha_env[0, 2, 0]))
        # environment heat conductance for the cell n
        env_con_zero = np.hstack((-.5 * self.k_alpha_env[0, 2, 0],
                                  env_con_base[1:]))
        # environment heat conductance for the zeroth cell
        env_con_vec = \
            np.hstack((np.tile(env_con_zero, self.n_ele),
                       np.tile(np.tile(env_con_base, self.n_cells - 2),
                               self.n_ele),
                       np.tile(env_con_n, self.n_ele)))
        # vector of the main diagonal of the heat conductance matrix
        self.mat_const = self.mat_const + np.diag(env_con_vec)
        # self.mat_dyn = np.copy(self.mat_const)

        self.mat_const_sp = sparse.csr_matrix(self.mat_const)
        # self.mat_dyn_sp = sparse.lil_matrix(self.mat_const_sp)

        """Calculating the coordinates of the gas channel heat conductance"""
        pos_cat_ch_base = \
            np.arange(1, (self.n_cells - 1) * self.n_ele * self.n_layer,
                      self.n_layer)
        # coordinates of the cathode channels
        # heat conductance for the 0-(n-1) cells
        pos_ano_ch_base = \
            np.arange(4, (self.n_cells - 1) * self.n_ele * self.n_layer,
                      self.n_layer)
        # coordinates of the anode channels
        # heat conductance for the 0-(n-1) cells
        pos_cat_ch_n = \
            np.arange((self.n_cells - 1) * self.n_ele * self.n_layer + 1,
                      self.n_ele * (self.n_layer * self.n_cells + 1),
                      self.n_layer + 1)
        # coordinates of the cathode channels
        # heat conductance for the n cell
        pos_ano_ch_n = \
            np.arange((self.n_cells - 1) * self.n_ele * self.n_layer \
                      + (self.n_layer - 1),
                      self.n_ele * (self.n_layer * self.n_cells + 1),
                      self.n_layer + 1)
        # coordinates of the anode channels
        # heat conductance for the n cell
        self.pos_cat_ch = np.hstack((pos_cat_ch_base, pos_cat_ch_n))
        self.pos_ano_ch = np.hstack((pos_ano_ch_base, pos_ano_ch_n))

    def update_values(self, k_alpha_ch, gamma, omega, v_loss, g_gas, i):
        """
        Updates the dynamic parameters

            Access to:
            -dict_temp_sys_dyn

            Manipulate:
            -self.g_fluid
            -self.k_gas_ch
            -self.cond_rate
            -self.i
            -self.v_loss
            -self.omega
        """
        self.g_fluid[0] = ip.interpolate_along_axis(g_gas[0], axis=1)
        self.g_fluid[1] = ip.interpolate_along_axis(g_gas[1], axis=1)
        self.k_gas_ch = k_alpha_ch
        self.cond_rate = gamma
        self.i = i
        self.v_loss = v_loss
        self.omega = omega

    def update(self):
        """
        This function coordinates the program sequence
        """
        # self.change_value_shape()
        self.update_gas_channel_lin()
        self.update_coolant_channel_lin()
        self.update_temp_layer()

    def update_temp_layer(self):
        """
        This function coordinates the temp_layer program sequence
        """
        self.update_matrix()
        self.update_rhs()
        self.solve_system()
        self.sort_results()

    def update_gas_channel_lin(self):
        """
        Calculates the fluid temperatures in the anode and cathode channels

            Access to:
            -self.k_gas_ch
            -self.temp_layer
            -self.temp_fluid
            -self.g_fluid

            Manipulate:
            -self.temp_fluid
        """

        self.temp_fluid[0].fill(self.temp_gas_in[0])
        self.temp_fluid[1].fill(self.temp_gas_in[1])
        self.temp_fluid_ele[0].fill(self.temp_gas_in[0])
        self.temp_fluid_ele[1].fill(self.temp_gas_in[1])

        for i in range(self.n_cells):

            dtemp = self.k_gas_ch[0, i] / self.g_fluid[0, i] \
                * (self.temp_layer[i][1, :] - self.temp_fluid_ele[0, i])
            g_func.add_source(self.temp_fluid[0, i], dtemp, 1)
            temp_fluid_ele = ip.interpolate_1d(self.temp_fluid[0, i])
            self.temp_fluid_ele[0, i] = \
                np.minimum(temp_fluid_ele, self.temp_layer[i][0])

            dtemp = self.k_gas_ch[1, i] / self.g_fluid[1, i] \
                * (self.temp_layer[i][4, :] - self.temp_fluid_ele[1, i])
            g_func.add_source(self.temp_fluid[1, i], dtemp, -1)
            temp_fluid_ele = ip.interpolate_1d(self.temp_fluid[1, i])
            self.temp_fluid_ele[1, i] = \
                np.minimum(temp_fluid_ele, self.temp_layer[i][4])

            self.temp_fluid[0] = \
                ip.interpolate_along_axis(self.temp_fluid_ele[0], axis=1,
                                          add_edge_points=True)
            self.temp_fluid[0, :, 0] = self.temp_gas_in[0]
            self.temp_fluid[1] = \
                ip.interpolate_along_axis(self.temp_fluid_ele[1], axis=1,
                                          add_edge_points=True)
            self.temp_fluid[1, :, -1] = self.temp_gas_in[1]

    def update_coolant_channel_lin(self):
        """
                Calculates the coolant channel temperatures.

                    Access to:
                    -self.temp_layer
                    -self.g_cool
                    -self.k_cool

                    Manipulate:
                    -self.temp_cool
                    -self.temp_cool_ele
                """
        for i in range(self.n_cells):
            dtemp = self.k_cool / self.g_cool \
                * (self.temp_layer[i][0, :] - self.temp_cool_ele[i])
            g_func.add_source(self.temp_cool[i], dtemp, 1)
            self.temp_cool_ele[i] = ip.interpolate_1d(self.temp_cool[i])
        if self.cool_ch_bc:
            dtemp = self.k_cool / self.g_cool \
                * (self.temp_layer[-1][-1, :] - self.temp_cool_ele[-1])
            g_func.add_source(self.temp_cool[-1], dtemp, 1)
            self.temp_cool_ele[-1] = ip.interpolate_1d(self.temp_cool[-1])

    # @jit(nopython=True)
    def update_rhs(self):
        """
        Creates a vector with the right hand side entries,
        add explicit heat sources here.
        Sources from outside the system
        to the system must to be defined negative.

            Access to:
            -self.cond_rate
            -self.cell_numb
            -self.elements
            -self.temp_env
            -self.k_alpha_env
            -self.temp_fluid_ele
            -self.temp_cool_ele
            -self.k_gas_ch
            -self.v_loss
            -self.v_tn
            -g_par.dict_case['e_0]
            -self.omega
            -self.i
            -self.k_cool
            -self.cool_ch_bc

            Manipulate:
            -self.rhs
        """
        heat_pow = self.dict['heat_pow']

        self.rhs.fill(0.)
        rhs = self.rhs
        temp_env = self.temp_env
        k_alpha_env = self.k_alpha_env

        ct = 0
        # print(self.cond_rate)
        for i in range(self.n_cells):
            for j in range(self.n_ele):
                if i is 0:
                    rhs[ct] = -.5 * temp_env * k_alpha_env[0, 2, i]
                else:
                    rhs[ct] = - temp_env * k_alpha_env[0, 2, i]

                rhs[ct + 1] = -temp_env * k_alpha_env[0, 1, i] \
                    - self.temp_fluid_ele[0, i, j] * self.k_gas_ch[0, i, j] \
                    - w_prop.water.calc_h_vap(self.temp_fluid[0, i, j]) \
                    * self.cond_rate[0, i, j]
                rhs[ct + 2] = \
                    - temp_env * k_alpha_env[0, 0, i] \
                    - (self.v_tn - g_par.dict_case['e_0'] + self.v_loss[0, i, j]
                       + .5 * self.omega[i, j] * self.i[i, j]) * self.i[i, j]
                rhs[ct + 3] = \
                    - temp_env * k_alpha_env[0, 0, i] \
                    - (self.v_loss[1, i, j]
                       + self.omega[i, j] * self.i[i, j] * .5) * self.i[i, j]
                rhs[ct + 4] = - temp_env * k_alpha_env[0, 1, i] \
                    - self.temp_fluid_ele[1, i, j] * self.k_gas_ch[1, i, j] \
                    - w_prop.water.calc_h_vap(self.temp_fluid[1, i, j]) \
                    * self.cond_rate[1, i, j]
                if i is 0:
                    rhs[ct] -= heat_pow
                    if self.cool_ch_bc:
                        rhs[ct] -= self.k_cool * self.temp_cool_ele[0, j]
                    cr = 5
                elif 0 < i < self.n_cells - 1:
                    rhs[ct] -= self.k_cool * self.temp_cool_ele[i, j]
                    cr = 5
                else:
                    rhs[ct] -= self.k_cool * self.temp_cool_ele[i, j]
                    rhs[ct + 5] -= heat_pow \
                        - .5 * self.k_alpha_env[0, 2, 0] * temp_env
                    if self.cool_ch_bc:
                        rhs[ct + 5] -= self.k_cool * self.temp_cool_ele[-1, j]
                    cr = 6
                ct += cr
        # print(self.rhs)

    # @jit(nopython=True)
    def update_matrix(self):
        """
        Updates the thermal conductance matrix

            Access to:
            -self.elements
            -self.cell_numb
            -self.k_gas_ch
            -self.mat_const

            Manipulate:
            -self.mat_dyn
        """
        ct = 0
        for q in range(self.n_cells):
            if q is not self.n_cells - 1:
                cr = 5
            else:
                cr = 6
            for w in range(self.n_ele):
                self.dyn_vec[ct + 1] = -self.k_gas_ch[0, q, w]
                self.dyn_vec[ct + 4] = -self.k_gas_ch[1, q, w]
                ct += cr
        self.mat_dyn = \
           self.mat_const_sp + sparse.diags([self.dyn_vec], [0], format='csr')
        # self.mat_dyn = self.mat_const + np.diag(self.dyn_vec)

    def solve_system(self):
        """
        Solves the layer temperatures.

            Access to:
            -self.mat_dyn
            -self.rhs

            Manipulate:
            -self.temp_layer_vec
        """
        # self.temp_layer_vec = np.linalg.tensorsolve(self.mat_dyn, self.rhs)
        self.temp_layer_vec = spsolve(self.mat_dyn, self.rhs)

    def sort_results(self):
        """
        Sorts the temperatures in the 1-d-array self.temp_layer_vec
        to the 3-d-list self.temp_layer

            Access to:
            -self.cell_numb
            -self.elements
            -self.temp_layer_vec

            Manipulate:
            -self.temp_layer
        """
        ct = 0
        for q in range(self.n_cells):
            if q is not self.n_cells - 1:
                cr = 5
            else:
                cr = 6
            for w in range(self.n_ele):
                self.temp_layer[q][:, w] = self.temp_layer_vec[ct: ct + cr]
                ct += cr
        # print(self.temp_layer)