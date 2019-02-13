import numpy as np
import scipy.linalg as sp_l
import data.global_parameters as g_par
import system.global_functions as g_func
import data.water_properties as w_prop

np.set_printoptions(linewidth=10000, threshold=None, precision=2)


class TemperatureSystem:

    def __init__(self, temp_sys_const_dict):
        # Handover
        self.n_cells = temp_sys_const_dict['cell_numb']
        # cell number
        self.nodes = temp_sys_const_dict['nodes']
        # node number
        self.n_ele = self.nodes - 1
        # element number
        self.ch_length = temp_sys_const_dict['channel_length']
        # channel length
        self.ch_width = temp_sys_const_dict['channel_width']
        # channel width
        self.cool_ch_bc = temp_sys_const_dict['cool_ch_bc']
        # coolant geometry condition
        self.temp_gas_in = temp_sys_const_dict['temp_gas_in']
        # gas inlet temperature
        self.temp_cool_in = temp_sys_const_dict['cool_temp_in']
        # coolant inlet temperature
        self.temp_layer_init = temp_sys_const_dict['temp_layer_init']
        # initial temperature
        self.cp_cool = temp_sys_const_dict['cool_cp']
        # coolant heat capacity
        self.m_flow_cool = temp_sys_const_dict['cool_m_flow']
        # mass flow of the coolant
        self.rho_cool = temp_sys_const_dict['cool_density']
        # density of the coolant
        self.visc_cool = temp_sys_const_dict['cool_visc']
        # density of the coolant
        self.height_cool = temp_sys_const_dict['channel_height']
        # height of the coolant channel
        self.width_cool = temp_sys_const_dict['channel_width']
        # width of the coolant channel
        self.cool_numb = temp_sys_const_dict['cool_ch_numb']
        # number of coolant channels
        self.ch_numb = temp_sys_const_dict['gas_ch_numb']
        # number of gas channels
        self.heat_pow = temp_sys_const_dict['heat_pow']
        # end plate heat power
        self.lambda_cool = temp_sys_const_dict['cool_lambda']
        # heat conductance from the channel to the coolant
        self.k_layer = temp_sys_const_dict['k_layer']
        # heat conductance array through and along the control volume
        self.k_alpha_env = temp_sys_const_dict['k_alpha_env']
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
            np.full(self.n_ele * (5 * (self.n_cells - 1) + 6), 0.)
        # unsorted result layer temperature vector
        self.rhs = np.full(self.n_ele * (5 * (self.n_cells - 1) + 6), 0.)
        # right side of the matrix system: mat T = rhs,
        # contains the power sources and explicit coupled terms
        self.temp_fluid = np.full((2, self.n_cells, self.nodes),
                                  self.temp_gas_in[0])
        self.temp_fluid_ele = np.full((2, self.n_cells, self.n_ele), 0.)
        # temperature of the fluid 0: cathode fluids, 1: anode fluids
        self.g_fluid = np.full((2, self.n_cells, self.n_ele), 0.)
        # heat capacity flow of the fluid 0: cathode fluids, 1: anode fluids
        self.cond_rate = np.full((2, self.n_cells, self.nodes), 0.)
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
        self.g_cool = self.cp_cool * self.m_flow_cool * self.cool_numb
        # coolant heat capacity flow
        self.k_cool = None
        # heat conductance between the coolant and the channel wall

        """Calculating the coolant to channel thermal conductance"""
        pr_ch = self.visc_cool * self.cp_cool / self.lambda_cool
        # prandtl number
        d_h_cool = 2. * self.width_cool * self.height_cool \
            / (self.width_cool + self.height_cool)
        # hydraulic diameter of the coolant channel
        u_ch = self.m_flow_cool / (self.width_cool
                                   * self.height_cool * self.rho_cool)
        # velocity of the coolant flow
        re_ch = self.rho_cool * u_ch * d_h_cool / self.visc_cool
        # reynolds number in the coolant channel

        nu_1 = 3.66
        nu_2 = 1.66 * np.sqrt(re_ch * pr_ch * d_h_cool / self.ch_length)
        nu_3 = (2. / (1. + 22. * pr_ch)) ** (1. / 6.) \
            * np.sqrt(re_ch * pr_ch * d_h_cool / self.ch_length)
        nu_lam = (nu_1 ** 3. + 0.7 ** 3. + (nu_2 - 0.7) ** 3.
                  + nu_3 ** 3.) ** (1. / 3.)
        # laminar nusselt number
        zeta = (1.8 * np.log(re_ch) - 1.5) ** -2.
        nu_turb = zeta / 8. * re_ch * pr_ch \
            / (1. + 12.7 * np.sqrt(zeta / 8.)
                * (pr_ch ** 2. / 3.) - 1.) \
            * (1. + (d_h_cool / self.ch_length) ** 2. / 3.)
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
        conv_area = d_h_cool * np.pi * self.ch_length / self.n_ele
        # convection area of the channel wall
        self.k_cool = conv_coeff_ch * conv_area * self.cool_numb
        # thermal conductance between the element channel area and the coolant

        """Building up the result temperature list and arrays"""
        temp_layer = np.full((5, self.n_ele), self.temp_layer_init)
        temp_layer_n = np.full((6, self.n_ele), self.temp_layer_init)
        self.temp_layer = []
        for q in range(self.n_cells - 1):
            self.temp_layer.append(temp_layer)
        self.temp_layer.append(temp_layer_n)
        # layer temperature list cell, layer, element
        #temp_cool_out = self.temp_cool_in + op_con.tar
        if self.cool_ch_bc is True:
            self.temp_cool = np.full((self.n_cells + 1, self.nodes),
                                     self.temp_cool_in)
            self.temp_cool_ele = np.full((self.n_cells + 1,
                                          self.n_ele), 0.)
        else:
            self.temp_cool = np.full((self.n_cells, self.nodes),
                                     self.temp_cool_in)
            self.temp_cool_ele = np.full((self.n_cells, self.n_ele), 0.)
        # coolant temperature array cell, element

        """Building up the base conductance matrix mat"""
        mat_base = np.full((5, 5), 0.)
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
        mat_n = np.full((6, 6), 0.)
        mat_n[0:5, 0:5] = mat_base
        mat_n[4, 4] -= self.k_layer[0, 2, 0]
        mat_n[4, 5] = self.k_layer[0, 2, 0]
        mat_n[5, 4] = self.k_layer[0, 2, 0]
        mat_n[5, 5] = -self.k_layer[0, 2, 0]
        # heat conductance matrix in z-direction for the cell n

        list_mat = []
        for i in range(self.n_cells):
            for j in range(self.n_ele):
                if i is self.n_cells - 1:
                    list_mat.append(mat_n)
                else:
                    list_mat.append(mat_base)
        # list of all the heat conductance matrix in z-direction
        # for all cells and all elements
        self.mat_const = sp_l.block_diag(*list_mat)
        # uncoupled heat conductance matrix in z-direction

        """Setting the coolant channel heat conductance"""
        cool_pos_n_up = \
            np.arange(self.n_ele * (self.n_cells - 1) * 5,
                      self.n_ele * (5 * (self.n_cells - 1) + 6), 6)
        # upper cool ch pos for the n cell

        if self.cool_ch_bc is True:
            cool_pos_base = np.arange(0,
                                      self.n_ele * (self.n_cells - 1) * 5,
                                      5)
            # cool ch pos for the 0-(n-1) cell
            cool_pos_n_down = \
                np.arange(self.n_ele * (self.n_cells - 1) * 5 + 5,
                          self.n_ele * (5 * (self.n_cells - 1) + 6), 6)
            # lower cool ch pos for the n cell
            for q, item in enumerate(cool_pos_n_up):
                self.mat_const[item, item] -= self.k_cool

        else:
            cool_pos_base = \
                np.arange(self.n_ele,
                          self.n_ele * (self.n_cells - 1) * 5, 5)
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
            np.hstack((np.tile(x_con_zero, self.n_ele - 1), np.zeros(5))),
            np.tile(np.hstack((np.tile(x_con_base, self.n_ele - 1),
                               np.zeros(5))), self.n_cells - 2),
            np.zeros(6 * (self.n_ele - 1) + 1)))
        # heat conductance vec for the right and left side
        # of the matrix main diagonal for the cells 0-(n-1)
        x_con_side_n = np.hstack((np.zeros((self.n_cells - 1)
                                           * self.n_ele * 5),
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
            + np.diag(x_con_side_base, 5) \
            + np.diag(x_con_side_base, -5) \
            + np.diag(x_con_side_n, 6) \
            + np.diag(x_con_side_n, -6)

        """Setting the cell connecting heat conductance, z-direction"""
        pos_r, pos_c = [], []
        for ct in range(self.n_ele * (self.n_cells - 1)):
            pos_r.append(4 + 5 * ct)
            if ct <= self.n_ele * (self.n_cells - 2):
                pos_c.append(5 * self.n_ele + 5 * ct)
            else:
                pos_c.append(pos_c[-1] + 6)
        for ct, item in enumerate(pos_c):
            self.mat_const[pos_c[ct], pos_r[ct]] += self.k_layer[0, 2, 0]
            self.mat_const[pos_r[ct], pos_c[ct]] += self.k_layer[0, 2, 0]
        # heat conductance outside the main diagonal

        pos_base_out = np.arange(4,
                                 self.n_ele * (self.n_cells - 1) * 5, 5)
        # coordinates of the heat conductance
        # of the last layer of the elements for the cells 0-(n-1)
        pos_base_in = np.arange(self.n_ele * 5,
                                self.n_ele * (self.n_cells - 1) * 5, 5)
        # coordinates of the heat conductance
        # of the first layer of the cells 1-(n-1)
        pos_n_in = np.arange((self.n_cells - 1) * self.n_ele * 5,
                             self.n_ele * (5 * (self.n_cells - 1) + 6), 6)
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
        self.mat_dyn = self.mat_const

        """Calculating the coordinates of the gas channel heat conductance"""
        pos_cat_ch_base = np.arange(1,
                                    (self.n_cells - 1) * self.n_ele * 5, 5)
        # coordinates of the cathode channels
        # heat conductance for the 0-(n-1) cells
        pos_ano_ch_base = np.arange(4,
                                    (self.n_cells - 1) * self.n_ele * 5, 5)
        # coordinates of the anode channels
        # heat conductance for the 0-(n-1) cells
        pos_cat_ch_n = np.arange((self.n_cells - 1) * self.n_ele * 5 + 1,
                                 self.n_ele * (5 * (self.n_cells - 1) + 6),
                                 6)
        # coordinates of the cathode channels
        # heat conductance for the n cell
        pos_ano_ch_n = np.arange((self.n_cells - 1) * self.n_ele * 5 + 4,
                                 self.n_ele * (5 * (self.n_cells - 1) + 6),
                                 6)
        # coordinates of the anode channels
        # heat conductance for the n cell
        self.pos_cat_ch = np.hstack((pos_cat_ch_base, pos_cat_ch_n))
        self.pos_ano_ch = np.hstack((pos_ano_ch_base, pos_ano_ch_n))

    def update_values(self, dict_temp_sys_dyn):
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

        self.g_fluid = \
            np.array([g_func.calc_elements_2d(dict_temp_sys_dyn['g_gas'][0]),
                      g_func.calc_elements_2d(dict_temp_sys_dyn['g_gas'][1])])\
            * self.ch_numb
        self.k_gas_ch = dict_temp_sys_dyn['k_gas_ch']
        self.cond_rate = dict_temp_sys_dyn['cond_rate'] * self.ch_numb
        self.i = dict_temp_sys_dyn['i']
        self.v_loss = dict_temp_sys_dyn['v_loss']
        self.omega = dict_temp_sys_dyn['omega']

    def change_value_shape(self):
        """
        Changes the array shape

            Access to:
            -self.v_loss
            -self.k_gas_ch

            Manipulate:
            -self.v_loss
            -self.k_gas_ch
        """

        self.v_loss = np.array(self.v_loss)
        self.k_gas_ch = np.array([g_func.calc_elements_2d(self.k_gas_ch[0]),
                                  g_func.calc_elements_2d(self.k_gas_ch[1])])\
            * self.ch_numb

    def update(self):
        """
        This function coordinates the program sequence
        """

        self.change_value_shape()
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

        for q in range(self.n_cells):
            for w in range(1, self.nodes):
                self.temp_fluid[0, q, w] =\
                    g_func.calc_fluid_temp_out(self.temp_fluid[0, q, w - 1],
                                               self.temp_layer[q][1, w - 1],
                                               self.g_fluid[0, q, w - 1],
                                               self.k_gas_ch[0, q, w - 1])
            temp_fluid_ele = g_func.calc_elements_1_d(self.temp_fluid[0, q])
            self.temp_fluid_ele[0, q] = np.minimum(temp_fluid_ele,
                                                   self.temp_layer[q][0])
            for w in range(self.n_ele - 1, -1, -1):
                self.temp_fluid[1, q, w] =\
                    g_func.calc_fluid_temp_out(self.temp_fluid[1, q, w + 1],
                                               self.temp_layer[q][4, w],
                                               self.g_fluid[1, q, w],
                                               self.k_gas_ch[1, q, w])
                temp_fluid_ele = g_func.calc_elements_1_d(self.temp_fluid[1, q])
                self.temp_fluid_ele[1, q] = np.minimum(temp_fluid_ele,
                                                       self.temp_layer[q][4])
            self.temp_fluid[0] = g_func.calc_nodes_2_d(self.temp_fluid_ele[0])
            self.temp_fluid[0, :, 0] = self.temp_gas_in[0]
            self.temp_fluid[1] = g_func.calc_nodes_2_d(self.temp_fluid_ele[1])
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

        for q in range(self.n_cells):
            for w in range(1, self.nodes):
                self.temp_cool[q, w] =\
                    g_func.calc_fluid_temp_out(self.temp_cool[q, w - 1],
                                               self.temp_layer[q][0, w - 1],
                                               self.g_cool, self.k_cool)
                self.temp_cool_ele[q] = g_func.calc_elements_1_d(self.temp_cool[q])
        if self.cool_ch_bc is True:
            for w in range(1, self.nodes):
                self.temp_cool[-1, w] =\
                    g_func.calc_fluid_temp_out(self.temp_cool[-1, w - 1],
                                               self.temp_layer[-1][-1, w - 1],
                                               self.g_cool, self.k_cool)
                self.temp_cool_ele[-1] = g_func.calc_elements_1_d(self.temp_cool[-1])

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
        self.rhs = np.full(self.n_ele * (5 * (self.n_cells - 1) + 6), 0.)
        rhs = self.rhs
        temp_env = self.temp_env
        k_alpha_env = self.k_alpha_env

        ct = 0
        for q in range(self.n_cells):
            for w in range(self.n_ele):
                if q is 0:
                    rhs[ct] = -.5 * temp_env * k_alpha_env[0, 2, q]
                else:
                    rhs[ct] = - temp_env * k_alpha_env[0, 2, q]

                rhs[ct + 1] = -temp_env * k_alpha_env[0, 1, q] \
                    - self.temp_fluid_ele[0, q, w] * self.k_gas_ch[0, q, w] \
                    - w_prop.water.calc_h_vap(self.temp_fluid[0, q, w]) \
                    * self.cond_rate[0, q, w]
                rhs[ct + 2] = \
                    - temp_env * k_alpha_env[0, 0, q] \
                    - (self.v_tn - g_par.dict_case['e_0'] + self.v_loss[0, q, w]
                       + .5 * self.omega[q, w] * self.i[q, w]) * self.i[q, w]
                rhs[ct + 3] = \
                    - temp_env * k_alpha_env[0, 0, q] \
                    - (self.v_loss[1, q, w]
                       + self.omega[q, w] * self.i[q, w] * .5) * self.i[q, w]
                rhs[ct + 4] = - temp_env * k_alpha_env[0, 1, q] \
                    - self.temp_fluid_ele[1, q, w] * self.k_gas_ch[1, q, w] \
                    - w_prop.water.calc_h_vap(self.temp_fluid[1, q, w]) \
                    * self.cond_rate[1, q, w]
                if q is 0:
                    rhs[ct] -= self.heat_pow
                    if self.cool_ch_bc is True:
                        rhs[ct] -= self.k_cool * self.temp_cool_ele[0, w]
                    cr = 5
                elif 0 < q < self.n_cells - 1:
                    rhs[ct] -= self.k_cool * self.temp_cool_ele[q, w]
                    cr = 5
                else:
                    rhs[ct] -= self.k_cool * self.temp_cool_ele[q, w]
                    rhs[ct + 5] -= self.heat_pow \
                        - .5 * self.k_alpha_env[0, 2, 0] * temp_env
                    if self.cool_ch_bc is True:
                        rhs[ct + 5] -= self.k_cool * self.temp_cool_ele[-1, w]
                    cr = 6
                ct += cr

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
        dyn_vec = np.full(self.n_ele * (5 * (self.n_cells - 1) + 6), 0.)
        ct = 0
        for q in range(self.n_cells):
            if q is not self.n_cells - 1:
                cr = 5
            else:
                cr = 6
            for w in range(self.n_ele):
                dyn_vec[ct + 1] = -self.k_gas_ch[0, q, w]
                dyn_vec[ct + 4] = -self.k_gas_ch[1, q, w]
                ct += cr
        self.mat_dyn = self.mat_const + np.diag(dyn_vec)

    def solve_system(self):
        """
        Solves the layer temperatures.

            Access to:
            -self.mat_dyn
            -self.rhs

            Manipulate:
            -self.temp_layer_vec
        """

        self.temp_layer_vec = np.linalg.tensorsolve(self.mat_dyn, self.rhs)

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