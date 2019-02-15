import numpy as np
import system.global_functions as g_func
import data.global_parameters as g_par
import system.half_cell as h_c
import data.half_cell_dict as hc_dict


class Cell:

    def __init__(self, dict_cell):
        self.dict_cell = dict_cell
        # Handover
        self.anode = h_c.HalfCell(hc_dict.dict_anode, dict_cell)
        # anode - object of the class HalfCell
        self.cathode = h_c.HalfCell(hc_dict.dict_cathode, dict_cell)
        # cathode - object of the class HalfCell
        self.th_mem = dict_cell['th_mem']
        # thickness membrane
        self.lambda_bpp = [dict_cell['lambda_z_bpp'], dict_cell['lambda_x_bpp']]
        # heat conductivity of the bipolar plate
        self.lambda_gde = [dict_cell['lambda_z_gde'], dict_cell['lambda_x_gde']]
        # heat conductivity of the gas diffusion layer
        self.lambda_mem = [dict_cell['lambda_z_mem'], dict_cell['lambda_x_mem']]
        # heat conductivity of the membrane
        self.temp_cool_in = dict_cell['temp_cool_in']
        # coolant inlet temperature
        self.mem_base_r = dict_cell['mem_base_r']
        # basic electrical resistance of the membrane
        self.mem_acl_r = dict_cell['mem_acl_r']
        # thermal related electrical resistance gain of the membrane
        self.calc_mem_loss = dict_cell['calc_mem_loss']

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
        self.width_channels = self.cathode.channel.width\
            * self.cathode.channel_numb\
            + self.cathode.channel.rack_width\
            * (self.cathode.channel_numb + 1)
        self.active_area_dx = self.width_channels * self.cathode.channel.dx
        self.k_bpp_z = self.lambda_bpp[0] * self.active_area_dx \
            / self.cathode.th_bpp
        # heat conductivity through the bipolar plate
        self.k_gde_z = self.lambda_gde[0]\
            * self.active_area_dx \
            / self.cathode.th_gde
        # heat conductivity through the gas diffusion electrode
        self.k_mem_z = self.lambda_mem[0]\
            * self.active_area_dx \
            / self.th_mem
        # heat conductivity through the membrane
        self.k_bpp_x = self.width_channels * self.lambda_bpp[1] \
            * self.cathode.th_bpp / self.cathode.channel.dx
        # heat conductivity along the bipolar plate
        self.k_gp = (self.width_channels
                     * (self.lambda_bpp[1] * self.cathode.th_bpp
                        + self.lambda_gde[1] * self.cathode.th_gde))\
            / (2. * self.cathode.channel.dx)
        # heat conductivity alon the bipolar plate and gas diffusion electrode
        self.k_gm = (self.width_channels
                     * (self.lambda_mem[1] * self.th_mem
                        + self.lambda_gde[1] * self.cathode.th_gde))\
            / (2. * self.cathode.channel.dx)
        # heat conductivity alon the gas diffusion electrode and membrane

        """boolean alarms"""
        self.v_alarm = False
        # True if :voltage loss > cell voltage
        self.break_program = False
        # True if the program aborts because of some critical impact

        """general parameter"""
        self.height = self.th_mem\
            + 2. * self.cathode.th_bpp\
            + 2. * self.cathode.th_gde
        # height of the cell
        n_nodes = g_par.dict_case['nodes']
        # number of nodes along the channel
        n_ele = n_nodes - 1
        self.w_cross_flow = np.zeros(n_ele)
        # water cross flux through the membrane
        self.omega_ca = np.zeros(n_nodes)
        # area specific membrane resistance
        self.v_loss = np.full(n_ele, 0.)
        # voltage loss
        self.temp = np.full((5, n_ele), dict_cell['temp_init'])
        # layer temperature
        self.temp_mem = np.zeros(n_nodes)
        # membrane temperature
        self.i_cd = np.full(n_ele, 0.)
        # current density
        self.omega = np.full(n_ele, 0.)
        # membrane resistance
        self.mem_loss = np.full(n_ele, 0.)#
        # voltage loss at the membrane
        self.v = np.full(n_ele, 0.)
        # cell voltage
        self.resistance = np.full(n_ele, 0.)
        # cell resistance

    def update(self):
        """
        This function coordinates the program sequence
        """
        is_ht_pem = self.dict_cell['is_ht_pem']
        self.temp_mem = .5 * (self.temp[2] + self.temp[3])
        if not is_ht_pem:
            self.cathode.is_ht_pem = False
            self.anode.is_ht_pem = False
            self.cathode.w_cross_flow = self.w_cross_flow
            self.anode.w_cross_flow = self.w_cross_flow
        self.cathode.i_ca = self.i_cd
        self.anode.i_ca = self.i_cd
        self.cathode.set_layer_temperature([self.temp[2],
                                            self.temp[3],
                                            self.temp[4]])
        self.anode.set_layer_temperature([self.temp[0], self.temp[1]])
        self.cathode.update()
        self.anode.update()
        if self.anode.break_program or self.cathode.break_program:
            self.break_program = True
        else:
            if not is_ht_pem:
                self.calc_cross_water_flux()
                self.calc_mem_resistivity_springer()
            else:
                self.calc_mem_resistivity_kvesic()
            self.calc_membrane_loss()
            self.calc_voltage()
            self.calc_resistance()

    def calc_cross_water_flux(self):
        """
        Calculates the water cross flux through the membrane
        according to (Springer, 1991).

            Access to:
            -g_par.dict_case['vap_m_t_coef']
            -self.cathode.free_w
            -self.anode.free_w
            -g_par.dict_case['mol_con_m']
            -g_par.dict_uni['F']
            -self.i_ca
            -self.temp_mem
            -self.th_mem

            Manipulate:
            -self.w_cross_flow
        """
        vap_coeff = g_par.dict_case['vap_m_temp_coeff']
        humidity = [self.cathode.humidity, self.anode.humidity]
        humidity_ele = np.array([g_func.calc_elements_1_d(humidity[0]),
                                 g_func.calc_elements_1_d(humidity[1])])
        a = 0.043 + 17.81 * humidity_ele
        b = -39.85 * humidity_ele ** 2. + 36. * humidity_ele ** 3.
        free_w_content = a + b
        zeta_plus = free_w_content[0] + free_w_content[1] \
                    + self.i_cd / (2. * vap_coeff
                                   * g_par.dict_case['mol_con_m']
                                   * g_par.dict_uni['F'])
        zeta_negative =\
            (free_w_content[0]
             - free_w_content[1]
             + 5. * self.i_cd / (2. * vap_coeff * g_par.dict_case['mol_con_m']
                                 * g_par.dict_uni['F'])) \
            / (1. + g_func.dw(self.temp_mem) * zeta_plus
               / (self.th_mem * vap_coeff))
        m_c = 0.5 * (zeta_plus + zeta_negative)
        m_a = 0.5 * (zeta_plus - zeta_negative)
        self.w_cross_flow = \
            self.i_cd / g_par.dict_uni['F'] + g_par.dict_case['mol_con_m'] \
            * g_func.dw(self.temp_mem) * (m_a ** 2. - m_c ** 2.) \
            / (2. * self.th_mem)

    def calc_mem_resistivity_kvesic(self):
        """
        Calculates the membrane resistivity and resistance
        according to (Kvesic, 2013).

            Access to:
            -self.mem_base_r
            -self.mem_acl_r
            -self.temp_mem
            -self.omega_ca
            -self.cathode.cathode.active_area_dx_ch

            Manipulate:
            -self.omega_ca
            -self.omega
        """
        self.omega_ca = (self.mem_base_r
                         - self.mem_acl_r * self.temp_mem) * 1.e-4
        self.omega = self.omega_ca / self.active_area_dx

    def calc_mem_resistivity_gossling(self):
        """
        Calculates the membrane resitace for NT-PEMFC according to Goßling
        """
        lambda_x = np.full(g_par.dict_case['nodes'], 0.)
        res_t = np.exp(self.fac_m * 1.e3 / self.temp_mem + self.fac_n)
        r_avg = (self.cathode.humidity + self.anode.humidity) * 0.5
        for q in range(g_par.dict_case['nodes']):
            if r_avg[q] > 0:
                lambda_x[q] = 0.3 + 6. * r_avg[q] * (
                            1. - np.tanh(r_avg[q] - 0.5)) + 3.9 * np.sqrt(
                    r_avg[q]) * (1. + np.tanh((r_avg[q] - 0.89) / 0.23))
            else:
                lambda_x[q] = -1. / (r_avg[q] - (3. + 1. / 3.))
        a = -0.007442
        b = 0.006053
        c = 0.0004702
        d = 1.144
        e = 8.
        res_lambda = a + b * lambda_x + c * np.exp(d * (lambda_x - e))
        res = res_lambda * res_t / 0.01415
        rp = self.th_mem / res
        self.omega_ca = 1.e-4 * (self.fac_res_basic
                                 + rp * self.fac_res_fit)
        self.omega = self.omega_ca / self.active_area_dx

    def calc_mem_resistivity_springer(self):
        """
        Calculates the membrane resistivity
        for NT-PEMFC according to (Springer, 1991).

            Access to:
            -self.cathode.humidity
            -self.anode.humidity
            -self.temp_mem
            -self.th_mem
            -self.cathode.cathode.active_area_dx_ch

            Manipulate:
            -self.omega_ca
            -self.omega
        """
        humidity = (self.cathode.humidity + self.anode.humidity) * 0.5
        humidity_ele = g_func.calc_elements_1_d(humidity)
        a = 0.043 + 17.81 * humidity_ele
        b = -39.85 * humidity_ele ** 2. + 36. * humidity_ele ** 3.
        free_water_content = a + b
        mem_el_con = 0.005139 * free_water_content - 0.00326
        mem_el_con_temp =\
            np.exp(1268 * (0.0033 - 1. / self.temp_mem)) * mem_el_con
        self.omega_ca = self.th_mem / mem_el_con_temp * 1.e-4

    def calc_membrane_loss(self):
        """
        Calculates the voltage loss at the membrane.

            Access to:
            -self.omega_ca
            -self.i_ca

            Manipulate:
            -self.mem_los
        """
        #self.mem_loss = self.omega_ca * self.i_cd
        if self.calc_mem_loss is False:
            self.mem_loss = 0.
        else:
            self.mem_loss = self.omega_ca * self.i_cd

    def calc_voltage(self):
        """
            Calculates the cell voltage loss.
            Calculates the cell voltage
            If the cell voltage loss is larger than the
            open circuit cell voltage, the cell voltage is set to zero.

            Access to:
            -self.mem_loss
            -self.cathode.v_loss
            -self.anode.v_loss
            -g_par.dict_case['e_0']

            Manipulate:
            -self.v_loss
            -self.v
            -self.v_alarm
        """
        self.v_loss = self.mem_loss + self.cathode.v_loss + self.anode.v_loss
        if any(self.v_loss) >= g_par.dict_case['e_0']:
            self.v_alarm = True
        self.v_loss = np.minimum(self.v_loss, g_par.dict_case['e_0'])
        self.v = g_par.dict_case['e_0'] - self.v_loss

    def calc_resistance(self):
        """
            Calculates the electrical resistance of the element in z-direction

            Access to:
            -self.v_loss
            -self.i_ca
            -self.cathode.th_bpp
            -g_par.dict_case['plate_resistivity']

            Manipulate:
            -self.resistance
        """
        self.resistance = self.v_loss / self.i_cd + 2. \
                          * g_par.dict_case['bpp_resistivity'] \
                          * self.cathode.th_bpp
