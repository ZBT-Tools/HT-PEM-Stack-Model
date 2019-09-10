import numpy as np
import system.global_functions as g_func
import data.global_parameters as g_par
import system.half_cell as h_c
import system.interpolation as ip
import system.matrix_functions as mtx


class Cell:

    def __init__(self, number, cell_dict, anode_dict, cathode_dict,
                 ano_channel_dict, cat_channel_dict):
        self.cell_dict = cell_dict
        # Handover
        self.name = 'Cell ' + str(number)
        self.n_layer = 5
        if cell_dict['last_cell']:
            self.n_layer += 1
        n_nodes = g_par.dict_case['nodes']
        # number of nodes along the channel
        self.n_ele = n_nodes - 1
        n_ele = self.n_ele

        self.width = self.cell_dict['width']
        self.length = self.cell_dict['length']

        self.first_cell = cell_dict['first_cell']
        self.last_cell = cell_dict['last_cell']

        self.anode = h_c.HalfCell(anode_dict, cell_dict, ano_channel_dict)
        # anode - object of the class HalfCell
        self.cathode = h_c.HalfCell(cathode_dict, cell_dict, cat_channel_dict)
        self.half_cells = [self.cathode, self.anode]

        self.dx = self.cathode.channel.dx
        # cathode - object of the class HalfCell
        self.th_mem = cell_dict['th_mem']
        # thickness membrane
        self.lambda_bpp = [cell_dict['lambda_z_bpp'], cell_dict['lambda_x_bpp']]
        # heat conductivity of the bipolar plate
        self.lambda_gde = [cell_dict['lambda_z_gde'], cell_dict['lambda_x_gde']]
        # heat conductivity of the gas diffusion layer
        self.lambda_mem = [cell_dict['lambda_z_mem'], cell_dict['lambda_x_mem']]
        # heat conductivity of the membrane

        # reordering thermal conductivities
        self.lambda_thermal = \
            np.asarray([[cell_dict['lambda_z_bpp'],
                         cell_dict['lambda_z_gde'],
                         cell_dict['lambda_z_mem'],
                         cell_dict['lambda_z_gde'],
                         cell_dict['lambda_z_bpp']],
                        [cell_dict['lambda_x_bpp'],
                         cell_dict['lambda_x_gde'],
                         cell_dict['lambda_x_mem'],
                         cell_dict['lambda_x_gde'],
                         cell_dict['lambda_x_bpp']]])

        self.mem_base_r = cell_dict['mem_base_r']
        # basic electrical resistance of the membrane
        self.mem_acl_r = cell_dict['mem_acl_r']
        # thermal related electrical resistance gain of the membrane
        self.calc_mem_loss = cell_dict['calc_mem_loss']

        """membrane resistance parameter (Goßling)"""
        self.fac_res_fit = 0.5913
        self.fac_res_basic = 0.03
        self.res_25 = 101249.82 * self.th_mem \
            + 36.24 * self.th_mem\
            + 2805.83 * self.th_mem + 0.021
        self.res_65 = 3842453.95 * self.th_mem\
            - 0.2775 * self.th_mem\
            + 2.181 * self.th_mem + 0.029
        self.fac_m = (np.log10(self.res_65) - np.log10(self.res_25)) / (
                    (1000. / (65. + 273.15)) - (1000. / 25. + 273.15))
        self.fac_n = np.log10(self.res_65) - self.fac_m * 1000. / (65. + 273.15)

        """heat conductivity along and through the cell layers"""
        self.width_channels = \
            self.cathode.channel.width * self.cathode.n_chl \
            + self.cathode.channel.rib_width * (self.cathode.n_chl + 1)
        self.active_area_dx = self.width_channels * self.dx
        self.k_bpp_z = \
            self.lambda_bpp[0] * self.active_area_dx / self.cathode.th_bpp
        # heat conductivity through the bipolar plate
        self.k_gde_z = \
            self.lambda_gde[0] * self.active_area_dx / self.cathode.th_gde
        # heat conductivity through the gas diffusion electrode
        self.k_mem_z = \
            self.lambda_mem[0] * self.active_area_dx / self.th_mem
        # heat conductivity through the membrane
        self.k_bpp_x = self.width_channels * self.lambda_bpp[1] \
            * self.cathode.th_bpp / self.cathode.channel.dx
        # heat conductivity along the bipolar plate
        self.k_gp = (self.width_channels
                     * (self.lambda_bpp[1] * self.cathode.th_bpp
                        + self.lambda_gde[1] * self.cathode.th_gde)) \
            / (2. * self.cathode.channel.dx)
        # heat conductivity alon the bipolar plate and gas diffusion electrode
        self.k_gm = (self.width_channels
                     * (self.lambda_mem[1] * self.th_mem
                        + self.lambda_gde[1] * self.cathode.th_gde)) \
            / (2. * self.cathode.channel.dx)
        # heat conductivity along the gas diffusion electrode and membrane
        self.th_layer = \
            np.asarray([self.cathode.th_bpp,
                        self.cathode.th_gde,
                        self.th_mem,
                        self.anode.th_gde,
                        self.anode.th_bpp])

        # if cell_dict['last_cell']:
        #     k_layer_z = \
        #         np.asarray([self.k_bpp_z, self.k_gde_z, self.k_mem_z,
        #                     self. k_gde_z, self.k_bpp_z])
        #     k_layer_x = \
        #         np.asarray([self.k_bpp_x, self.k_gp, self.k_gm,
        #                     self. k_gm, self.k_gp, self.k_bpp_x])
        # else:
        #     k_layer_z = \
        #         np.asarray([self.k_bpp_z, self.k_gde_z, self.k_mem_z,
        #                     self. k_gde_z])
        #     k_layer_x = \
        #         np.asarray([self.k_bpp_x, self.k_gp, self.k_gm,
        #                     self. k_gm, self.k_gp])

        self.k_layer_z = \
            np.outer(self.lambda_thermal[0] / self.th_layer,
                     self.active_area_dx)
        self.k_layer_x = self.lambda_thermal[1] * self.th_layer
        self.k_layer_x = (self.k_layer_x + np.roll(self.k_layer_x, 1)) * 0.5
        self.k_layer_x = np.hstack((self.k_layer_x, self.k_layer_x[0]))
        self.k_layer_x = np.outer(self.k_layer_x, self.width_channels / self.dx)

        if not self.last_cell:
            self.k_layer_z = np.delete(self.k_layer_z, -1, 0)
            self.k_layer_x = np.delete(self.k_layer_x, -1, 0)

        if self.first_cell:
            self.k_layer_x[0] *= 0.5
        if self.last_cell:
            self.k_layer_x[-1] *= 0.5

        # print('k_layer')
        # print(np.sum(np.abs(k_layer_z - k_layer_z_1)))
        # print(np.sum(np.abs(k_layer_x - k_layer_x_1)))
        # print(k_layer_x)
        # print(k_layer_x_1)
        self.heat_cond_mtx = \
            mtx.build_cell_conductance_matrix(self.k_layer_x, self.k_layer_z,
                                              n_ele)
        self.heat_mtx = self.heat_cond_mtx.copy()
        self.heat_rhs = np.full(self.n_layer * n_ele, 0.0)

        # Create array for each thermal layer with indices according to
        # corresponding position in center diagonal of conductance matrix and
        # right hand side vector
        index_list = []
        for i in range(self.n_layer):
            index_list.append([(j * self.n_layer) + i for j in range(n_ele)])
        self.index_array = np.asarray(index_list)

        # Set constant thermal boundary conditions
        if self.first_cell:
            heat_bc_0 = cell_dict['heat_pow']
            np.put(self.heat_rhs, self.index_array[0], -heat_bc_0)
        if self.last_cell:
            heat_bc_0 = cell_dict['heat_pow']
            np.put(self.heat_rhs, self.index_array[-1], heat_bc_0)

        """boolean alarms"""
        self.v_alarm = False
        # True if :voltage loss > cell voltage
        self.break_program = False
        # True if the program aborts because of some critical impact

        """general parameter"""
        self.height = self.th_mem \
            + self.cathode.th_bpp \
            + self.cathode.th_gde \
            + self.anode.th_bpp \
            + self.anode.th_gde

        # height of the cell
        #self.n_nodes = g_par.dict_case['nodes']

        self.w_cross_flow = np.zeros(n_ele)
        # water cross flux through the membrane
        self.omega_ca = np.zeros(n_ele)
        # area specific membrane resistance
        self.v_loss = np.full(n_ele, 0.)
        # voltage loss
        self.temp_layer = np.full((self.n_layer, n_ele), cell_dict['temp_init'])
        # layer temperature
        #self.temp_cool_in = cell_dict['temp_cool_in']
        # coolant inlet temperature
        self.temp_cool = np.full(n_ele, cell_dict['temp_cool_in'])
        self.temp_names = ['Cathode BPP-BPP',
                           'Cathode BPP-GDE',
                           'Cathode GDE-MEM',
                           'Anode MEM-GDE',
                           'Anode GDE-BPP',
                           'Anode BPP-BPP']
        # interface names according to temperature array
        self.temp_mem = np.zeros(n_ele)
        # membrane temperature
        self.i_cd = np.full(n_ele, 1.)
        # current density
        self.omega = np.full(n_ele, 0.)
        # membrane resistance
        self.mem_loss = np.full(n_ele, 0.)
        # voltage loss at the membrane
        self.v = np.full(n_ele, 0.)
        # cell voltage
        self.resistance = np.full(n_ele, 0.)
        # cell resistance
        self.print_data = [
            {
                'Current Density': {'value': self.i_cd, 'units': 'A/m²'},
                'Coolant Temperature': {'value': self.temp_cool, 'units': 'K'}
            },
            {
                'Temperature': {self.temp_names[i]:
                                {'value': self.temp_layer[i], 'units': 'K'}
                                for i in range(len(self.temp_layer))}
            }]

    def calc_ambient_conductance(self, alpha_amb):
        """
        :param alpha_amb: heat transfer coefficient for free or forced
        convection to the ambient
        :return: discretized conductance for each layer of cell based on
        external surface
        """
        # geometry model
        ext_surface_factor = (self.width + self.length) \
            / (self.cathode.channel.length * self.width_channels)
        k_amb = np.full((self.n_layer, self.n_ele), 0.)
        # convection conductance to the environment
        th_layer_amb = (self.th_layer + np.roll(self.th_layer, 1)) * 0.5
        if self.cell_dict['last_cell']:
            th_layer_amb = np.hstack((th_layer_amb, th_layer_amb[0]))
        # k_amb[0, :] = self.cathode.th_bpp
        # k_amb[1, :] = (self.cathode.th_bpp + self.cathode.th_gde)
        # k_amb[1] = 0.5 * alpha_amb * self.dx\
        #     * (self.cathode.th_bpp + self.cathode.th_gde) / ext_surface_factor
        # k_amb[0] = 0.5 * alpha_amb * self.dx \
        #     * (self.cathode.th_bpp + self.th_mem) / ext_surface_factor
        # k_amb[2] = alpha_amb * self.dx * self.cathode.th_bpp / ext_surface_factor
        k_amb = np.outer(th_layer_amb, self.dx) * alpha_amb / ext_surface_factor
        if not self.last_cell:
            k_amb = np.delete(k_amb, -1, 0)
        if self.first_cell:
            k_amb[0] *= 0.5
        if self.last_cell:
            k_amb[-1] *= 0.5
        return k_amb

    def add_explicit_layer_source(self, source_term, layer_id=None):
        if not layer_id:
            if np.isscalar(source_term):
                source_vector = np.full_like(self.heat_rhs, -source_term)
            else:
                source_vector = np.asarray(-source_term)
        else:
            source_vector = np.zeros_like(self.heat_rhs)
            np.put(source_vector, self.index_array[layer_id], -source_term)
        self.heat_rhs += source_vector
        return source_vector

    def add_implicit_layer_source(self, coefficient, layer_id=None):
        if layer_id is None:
            if np.isscalar(coefficient):
                source_vector = np.full_like(self.heat_rhs, coefficient)
            else:
                source_vector = np.asarray(coefficient)
        else:
            source_vector = np.zeros_like(self.heat_rhs)
            np.put(source_vector, self.index_array[layer_id], coefficient)
        self.heat_mtx += np.diag(source_vector)

    def update(self):
        """
        This function coordinates the program sequence
        """
        is_ht_pem = self.cell_dict['is_ht_pem']
        self.temp_mem[:] = .5 * (self.temp_layer[2] + self.temp_layer[3])
        if not is_ht_pem:
            self.cathode.is_ht_pem = False
            self.anode.is_ht_pem = False
            self.cathode.w_cross_flow[:] = self.w_cross_flow
            self.anode.w_cross_flow[:] = self.w_cross_flow
        self.cathode.i_cd[:] = self.i_cd
        self.anode.i_cd[:] = self.i_cd
        # self.cathode.set_layer_temperature([self.temp[2], self.temp[3],
        #                                     self.temp[4]])
        # self.anode.set_layer_temperature([self.temp[0], self.temp[1]])
        self.cathode.update()
        self.anode.update()
        if self.anode.break_program or self.cathode.break_program:
            self.break_program = True
        else:
            if not is_ht_pem:
                self.calc_cross_water_flux()
                self.calc_mem_resistivity_springer()
                # self.calc_mem_resistivity_gossling()
            else:
                self.calc_mem_resistivity_kvesic()
            self.calc_membrane_loss()
            self.calc_voltage()
            self.calc_resistance()

    def calc_cross_water_flux(self):
        """
        Calculates the water cross flux through the membrane
        according to (Springer, 1991).
        """
        vap_coeff = g_par.dict_case['vap_m_temp_coeff']

        dw = 2.1e-7 * np.exp(-2436. / self.temp_mem)
        humidity = np.asarray([self.cathode.humidity, self.anode.humidity])
        humidity_ele = \
            np.array([ip.interpolate_1d(humidity[0]),
                      ip.interpolate_1d(humidity[1])])
        water_content = 0.043 + 17.81 * humidity_ele \
            - 39.85 * humidity_ele ** 2. + 36. * humidity_ele ** 3.
        zeta_plus = water_content[0] + water_content[1] \
            + self.i_cd / (2. * vap_coeff * g_par.dict_case['mol_con_m']
                           * g_par.dict_uni['F'])
        zeta_negative = \
            (water_content[0] - water_content[1]
             + 5. * self.i_cd / (2. * vap_coeff * g_par.dict_case['mol_con_m']
             * g_par.dict_uni['F'])) \
            / (1. + dw * zeta_plus / (self.th_mem * vap_coeff))
        m_c = 0.5 * (zeta_plus + zeta_negative)
        m_a = 0.5 * (zeta_plus - zeta_negative)
        self.w_cross_flow = \
            self.i_cd / g_par.dict_uni['F'] + g_par.dict_case['mol_con_m'] \
            * dw * (m_a ** 2. - m_c ** 2.) / (2. * self.th_mem)

    def calc_mem_resistivity_kvesic(self):
        """
        Calculates the membrane resistivity and resistance
        according to (Kvesic, 2013).
        """
        self.omega_ca[:] = (self.mem_base_r
                         - self.mem_acl_r * self.temp_mem) * 1.e-4
        self.omega[:] = self.omega_ca / self.active_area_dx

    def calc_mem_resistivity_gossling(self):
        """
        Calculates the membrane resitace for NT-PEMFC according to Goßling
        """
        res_t = np.exp(self.fac_m * 1.e3 / self.temp_mem + self.fac_n)
        r_avg = (self.cathode.humidity + self.anode.humidity) * 0.5
        lambda_x = np.where(r_avg > 0,
                            0.3 + 6. * r_avg * (1. - np.tanh(r_avg - 0.5))
                            + 3.9 * np.sqrt(r_avg)
                            * (1. + np.tanh((r_avg - 0.89) / 0.23)),
                            -1. / (r_avg - (3. + 1. / 3.)))
        res_lambda = -0.007442 + 0.006053 * lambda_x \
            + 0.0004702 * np.exp(1.144 * (lambda_x - 8.))
        res = res_lambda * res_t / 0.01415
        rp = self.th_mem / res
        self.omega_ca[:] = 1.e-4 * (self.fac_res_basic
                                 + rp * self.fac_res_fit)
        self.omega[:] = self.omega_ca / self.active_area_dx

    def calc_mem_resistivity_springer(self):
        """
        Calculates the membrane resistivity
        for NT-PEMFC according to (Springer, 1991).
        """
        humidity = (self.cathode.humidity + self.anode.humidity) * 0.5
        humidity_ele = ip.interpolate_1d(humidity)
        lambda_springer = \
            np.where(humidity_ele < 1.0,
                     0.043 + 17.81 * humidity_ele
                     - 39.85 * humidity_ele ** 2. + 36. * humidity_ele ** 3.,
                     14.0 + 1.4 * (humidity_ele - 1.0))
        lambda_springer[lambda_springer < 1.0] = 1.0
        mem_cond = (0.005139 * lambda_springer - 0.00326) \
            * np.exp(1268 * (0.0033 - 1. / self.temp_mem))
        self.omega_ca[:] = self.th_mem / mem_cond * 1.e-4

    def calc_membrane_loss(self):
        """
        Calculates the voltage loss at the membrane.
        """
        if not self.calc_mem_loss:
            self.mem_loss[:] = 0.
        else:
            self.mem_loss[:] = self.omega_ca * self.i_cd

    def calc_voltage(self):
        """
        Calculates the cell voltage loss. If the cell voltage loss is larger
        than the open circuit cell voltage, the cell voltage is set to zero.
        """
        self.v_loss[:] = \
            self.mem_loss + self.cathode.v_loss + self.anode.v_loss
        if np.any(self.v_loss) >= g_par.dict_case['e_0']:
            self.v_alarm = True
        self.v_loss[:] = np.minimum(self.v_loss, g_par.dict_case['e_0'])
        self.v[:] = g_par.dict_case['e_0'] - self.v_loss

    def calc_resistance(self):
        """
            Calculates the electrical resistance of the element in z-direction
        """
        self.resistance[:] = self.v_loss / self.i_cd + 2. \
            * g_par.dict_case['bpp_resistivity'] * self.cathode.th_bpp
