import numpy as np
import system.global_functions as g_func
import data.global_parameter as g_par
import system.half_cell as h_c
import data.half_cell_dict as hc_dict


class Cell:

    def __init__(self, dict_cell):
        # Handover
        self.anode = h_c.HalfCell(hc_dict.anode)
        self.cathode = h_c.HalfCell(hc_dict.cathode)
        self.th_mem = dict_cell['th_mem']
        self.lambda_p = [dict_cell['plate_h_con'], dict_cell['plate_hi_con']]
        self.lambda_g = [dict_cell['gde_h_con'], dict_cell['gde_hi_con']]
        self.lambda_m = [dict_cell['mem_h_con'], dict_cell['mem_hi_con']]
        self.t_cool_in = dict_cell['t_cool_in']
        self.mem_bas_r = dict_cell['mem_bas_r']
        self.mem_acl_r = dict_cell['mem_acl_r']
        # Scalar variables
        self.v_alarm = False
        self.fac_res_fit = 0.5913
        self.fac_res_basic = 0.03
        self.ko = 6.2
        self.ho = -52300.
        self.break_program = False
        self.res_25 = 101249.82 * self.th_mem \
                      + 36.24 * self.th_mem\
                      + 2805.83 * self.th_mem + 0.021
        self.res_65 = 3842453.95 * self.th_mem\
                      - 0.2775 * self.th_mem\
                      + 2.181 * self.th_mem + 0.029
        self.fac_m = (np.log10(self.res_65) - np.log10(self.res_25)) / (
                    (1000. / (65. + 273.15)) - (1000. / 25. + 273.15))
        self.fac_n = np.log10(self.res_65) - self.fac_m * 1000. / (65. + 273.15)
        self.k_p = self.lambda_p[0] * self.cathode.channel.plane_dx\
            / self.cathode.th_plate
        self.k_g = self.lambda_g[0] * self.cathode.channel.plane_dx\
            / self.cathode.th_gde
        self.k_m = self.lambda_m[0] * self.cathode.channel.plane_dx\
            / self.th_mem
        self.k_pp = self.cathode.channel.width * self.lambda_p[1]\
            * self.cathode.th_plate / self.cathode.channel.d_x
        self.k_gp = (self.cathode.channel.width
                     * (self.lambda_p[1] * self.cathode.th_plate
                        + self.lambda_g[1] * self.cathode.th_gde))\
            / (2. * self.cathode.channel.d_x)
        self.k_gm = (self.cathode.channel.width
                     * (self.lambda_m[1] * self.th_mem
                        + self.lambda_g[1] * self.cathode.th_gde))\
            / (2. * self.cathode.channel.d_x)
        self.height = self.th_mem\
            + 2. * self.cathode.th_plate\
            + 2. * self.cathode.th_gde
        nodes = g_par.dict_case['nodes']
        self.t_mem = np.full(nodes, 0.)
        self.j = np.full(nodes, 0.)
        self.m_a = np.full(nodes, 0.)
        self.m_c = np.full(nodes, 0.)
        self.m_0 = np.full(nodes, 0.)
        self.m_1 = np.full(nodes, 0.)
        self.m_pos_0 = np.full(nodes, 0.)
        self.m_pos_1 = np.full(nodes, 0.)
        self.omega_ca = np.full(nodes, 0.)
        self.con_ov = np.full(nodes, 0.)
        self.v_los = np.full(nodes - 1, 0.)
        self.t = np.full((5, nodes-1), dict_cell['t_init'])
        self.t_membrane = np.full(nodes, 0.)
        self.i_ca = np.full((nodes - 1), 0.)
        self.i_nodes = np.full(nodes, 0.)
        self.zeta_plus = np.full(nodes, 0.)
        self.zeta_negative = np.full(nodes, 0.)
        self.ke = np.full(nodes, 0.)
        self.omega_ca = np.full((nodes - 1), 0.)
        self.omega = np.full(nodes-1, 0.)
        self.mem_ov = np.full((nodes - 1), 0.)
        self.omega_gling = np.full(nodes, 0.)
        self.v = np.full((nodes - 1), 0.)
        self.resistance = np.full((nodes - 1), 0.)

    def update(self):
        self.t_membrane = .5 * (self.t[2] + self.t[1])
        if g_par.dict_case['pem_type'] is False:
            self.cathode.set_pem_type(False)
            self.anode.set_pem_type(False)
            self.cathode.set_j(self.j)
            self.anode.set_j(self.j)
        self.cathode.set_i(self.i_ca)
        self.anode.set_i(self.i_ca)
        self.cathode.set_t([self.t[2], self.t[3], self.t[4]])
        self.anode.set_t([self.t[0], self.t[1]])
        self.cathode.update()
        self.anode.update()
        if self.anode.break_program is True\
                or self.cathode.break_program is True:
            self.break_program = True
        else:
            self.calc_mem_resistivity()
            self.calc_membrane_losses()
            if g_par.dict_case['pem_type'] is False:
                self.calc_mem_block_1()
                self.calc_cross_over_water_flux()
                self.calc_mem_block_2()
                self.calc_con_ov()
            self.calc_voltage()
            self.calc_resistance()

    def set_i(self, i_ca):
        self.i_ca = i_ca
        self.i_nodes = g_func.calc_nodes_1_d(self.i_ca)

    def calc_mem_block_1(self):
        vap_coef = g_par.dict_case['vap_m_t_coef']
        self.zeta_plus = self.cathode.free_water\
                         + self.anode.free_water\
                         + self.i_nodes / (2.
                                           * vap_coef
                                           * g_par.dict_case['mol_con_m']
                                           * g_par.dict_uni['F'])
        self.zeta_negative = (self.cathode.free_water
                              - self.anode.free_water
                              + 5. * self.i_nodes
                              / (2.
                                 * vap_coef
                                 * g_par.dict_case['mol_con_m']
                                 * g_par.dict_uni['F'])) \
                                 / (1. + g_func.dw(self.t_membrane)
                                * self.zeta_plus / (self.th_mem * vap_coef))
        self.m_c = 0.5 * (self.zeta_plus + self.zeta_negative)
        self.m_a = 0.5 * (self.zeta_plus - self.zeta_negative)

    def calc_cross_over_water_flux(self):
        self.j = self.i_nodes / g_par.dict_uni['F'] + g_par.dict_case[
            'mol_con_m'] * g_func.dw(self.t_membrane) * (
                             self.m_a ** 2. - self.m_c ** 2.) / (
                             2. * self.th_mem)

    def calc_mem_block_2(self):
        self.ke = self.ko * np.exp(-self.ho / g_par.dict_uni['R'] * (
                    1. / self.t_membrane - 1. / g_par.dict_case['t_u']))
        a = self.m_c ** 2.
        b = self.m_a ** 2. - self.m_c ** 2.
        self.m_0 = np.sqrt(a - b)  # z = 0
        self.m_1 = np.sqrt(a)  # z = 1
        self.m_pos_0 = -self.ke * self.m_0 * 0.5 + np.sqrt(
            (self.ke * self.m_0 * 0.5) ** 2. + self.ke * self.m_0)
        self.m_pos_1 = -self.ke * self.m_1 * 0.5 + np.sqrt(
            (self.ke * self.m_1 * 0.5) ** 2. + self.ke * self.m_1)

    def calc_con_ov(self):
        self.con_ov = g_par.dict_uni['R'] * self.t_membrane * np.log(
            self.m_pos_1 / self.m_pos_0) / g_par.dict_uni['F']

    def calc_mem_resistivity(self):
        self.omega_ca = (self.mem_bas_r
                         - self.mem_acl_r * self.t_membrane) * 1.e-4
        self.omega = self.omega_ca / self.cathode.channel.plane_dx

    def calc_membrane_losses(self):
        self.mem_ov = self.omega_ca * self.i_ca

    def calc_mem_resistivity_gling(self):
        lambda_x = np.full(g_par.dict_case['nodes'], 0.)
        res_t = np.exp(self.fac_m * 1.e3 / self.t_membrane + self.fac_n)
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
        res_lama = a + b * lambda_x + c * np.exp(d * (lambda_x - e))
        res = res_lama * res_t / 0.01415
        rp = self.th_mem / res
        self.omega_gling = 1.e4 * (self.fac_res_basic
                                   + rp * self.fac_res_fit)

    def calc_voltage(self):
        self.v_los = self.mem_ov + self.cathode.v_los + self.anode.v_los
        if any(self.v_los) >= g_par.dict_case['e_0']:
            self.v_alarm = True
        self.v_los = np.minimum(self.v_los, g_par.dict_case['e_0'])
        self.v = g_par.dict_case['e_0'] - self.v_los

    def calc_resistance(self):
        self.resistance = self.v_los / self.i_ca + 2.\
            * g_par.dict_case['plate_resistivity'] * self.cathode.th_plate
