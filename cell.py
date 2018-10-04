import numpy as np
import global_functions as g_func
import global_parameter as g_par
import half_cell as h_c
import input_parameter as i_p


class Cell:

    def __init__(self, dict_cell):
        self.anode = h_c.HalfCell(i_p.anode)
        self.cathode = h_c.HalfCell(i_p.cathode)
        self.thick_mem = dict_cell['mem_thick']
        self.k_p = dict_cell['plate_h_con']
        self.k_g = dict_cell['gde_h_con']
        self.k_m = dict_cell['mem_h_con']
        self.ki_p = dict_cell['plate_hi_con']
        self.ki_g = dict_cell['gde_hi_con']
        self.ki_m = dict_cell['mem_hi_con']
        self.t_cool_in = dict_cell['t_cool_in']
        self.width = dict_cell['width']
        self.length = dict_cell['length']
        self.mem_bas_r = dict_cell['mem_bas_r']
        self.mem_acl_r = dict_cell['mem_acl_r']
        self.init_param()
        self.init_arrays()
        self.init_func()

    def init_func(self):
        self.t0 = (self.anode.t1 + self.cathode.t1) / 2.
        self.extent = 2. * (self.length + self.width)
        self.plane = self.length * self.width
        self.res_25 = 1.6658 ** -9. * 1.e7 * self.thick_mem \
                      + 4.3093 ** -7 * 1.e6 * self.thick_mem \
                      + 1.2894 ** -5 * 1.e4 * self.thick_mem + 0.021
        self.res_65 = 1.1270 ** -8 * 1.e7 * self.thick_mem \
                      - 8.6430 ** -7 * 1.e6 * self.thick_mem \
                      + 5.3983 ** -5 * 1.e4 * self.thick_mem + 0.029
        self.fac_m = (np.log10(self.res_65) - np.log10(self.res_25)) \
                     / ((1000. / (65. + 273.15)) - (1000. / 25. + 273.15))
        self.fac_n = np.log10(self.res_65) - self.fac_m * 1000. / (65. + 273.15)
        self.r_p = self.cathode.thickness_plate\
                   / (self.k_p * self.cathode.channel.plane_dx)
        self.r_g = self.cathode.thickness_gde\
                   / (self.k_g * self.cathode.channel.plane_dx)
        self.r_m = self.thick_mem / (self.k_m * self.cathode.channel.plane_dx)
        self.r_pp = self.cathode.channel.d_x / (self.cathode.channel.width
                                                * self.ki_p
                                                * self.cathode.thickness_plate)
        self.r_gp = 2. * self.cathode.channel.d_x \
                    / (self.cathode.channel.width
                       * (self.ki_p * self.cathode.thickness_plate
                          + self.ki_g * self.cathode.thickness_gde))
        self.r_gm = 2. * self.cathode.channel.d_x \
                    / (self.cathode.channel.width
                       * (self.ki_m * self.thick_mem
                          + self.ki_g * self.cathode.thickness_gde))
        self.height = self.thick_mem + 2. * self.cathode.thickness_plate \
                      + 2. * self.cathode.thickness_gde

    def init_param(self):
        self.fac_res_fit = 0.5913
        self.fac_res_basic = 0.03
        self.t1 = self.cathode.t1
        self.t2 = self.cathode.t2
        self.t2e = g_func.calc_elements(self.cathode.t2)
        self.t3 = self.cathode.t3
        self.t4 = self.anode.t1
        self.t5 = self.anode.t2
        self.t5e = g_func.calc_elements(self.anode.t2)
        self.ko = 6.2
        self.ho = -52300.
        self.break_program = False

    def init_arrays(self):
        self.j = np.full(g_par.dict_case['nodes'], 0.)
        self.m_a = np.full(g_par.dict_case['nodes'], 0.)
        self.m_c = np.full(g_par.dict_case['nodes'], 0.)
        self.m0 = np.full(g_par.dict_case['nodes'], 0.)
        self.m1 = np.full(g_par.dict_case['nodes'], 0.)
        self.m_pos0 = np.full(g_par.dict_case['nodes'], 0.)
        self.m_pos1 = np.full(g_par.dict_case['nodes'], 0.)
        self.omega = np.full(g_par.dict_case['nodes'], 0.)
        self.psi = np.full(g_par.dict_case['nodes'], 0.)
        self.t = np.full(g_par.dict_case['nodes'], self.t_cool_in)
        self.v_los = np.full(g_par.dict_case['nodes']-1, 0.)

    def update(self):
        self.t0 = .5 * (self.anode.t1 + self.cathode.t1)
        self.calc_temperature_elements()
        if g_par.dict_case['pem_type'] is False:
            self.cathode.set_pem_type(False)
            self.cathode.set_i(self.i)
            self.cathode.set_j(self.j)
            self.cathode.set_t([self.t1, self.t2, self.t3])
            self.cathode.update()
            self.anode.set_pem_type(False)
            self.anode.set_i(self.i)
            self.anode.set_j(self.j)
            self.anode.set_t([self.t4, self.t5])
            self.anode.update()
            if self.anode.break_program is True\
                    or self.cathode.break_program is True:
                self.break_program = True
            else:
                self.calc_mem_block_1()
                self.calc_cross_over_water_flux()
                self.calc_mem_block_2()
                self.calc_con_overpotential()
                self.calc_mem_resistivity()
                self.calc_voltage()
        else:
            self.cathode.set_pem_type(True)
            self.cathode.set_i(self.i)
            self.cathode.set_t([self.t1, self.t2, self.t3])
            self.cathode.update()
            self.anode.set_pem_type(True)
            self.anode.set_i(self.i)
            self.anode.set_t([self.t4, self.t5])
            self.anode.update()
            if self.anode.break_program is True\
                    or self.cathode.break_program is True:
                self.break_program = True
            else:
                self.calc_mem_resistivity()
                self.psi = np.full(g_par.dict_case['nodes'], 0.)
                self.calc_voltage()
                self.calc_resistance()

    def set_i(self, i):
        self.i = i
        self.i_n = g_func.calc_nodes_1d(self.i)

    def calc_temperature_elements(self):
        self.t2e = g_func.calc_elements(self.t2)
        self.t3e = g_func.calc_elements(self.t3)
        self.t5e = g_func.calc_elements(self.t5)

    def calc_mem_block_1(self):
        self.zeta_plus = self.cathode.free_water + self.anode.free_water + \
                         self.i_n / (2. * g_par.dict_case['vap_m_t_coef']
                                   * g_par.dict_case['mol_con_m'] *
                                   g_par.dict_uni['f'])
        self.zeta_negative = (self.cathode.free_water - self.anode.free_water
                              + 5. * self.i_n / (
                                      2. * g_par.dict_case['vap_m_t_coef']
                                      * g_par.dict_case['mol_con_m']
                                      * g_par.dict_uni['f']))\
                             / (1. + g_func.dw(self.t0) * self.zeta_plus
                                / (self.thick_mem
                                   * g_par.dict_case['vap_m_t_coef']))
        self.m_c = 0.5 * (self.zeta_plus + self.zeta_negative)
        self.m_a = 0.5 * (self.zeta_plus - self.zeta_negative)

    def calc_cross_over_water_flux(self):
        self.j = self.i_n / g_par.dict_uni['f'] \
                 + g_par.dict_case['mol_con_m'] * g_func.dw(self.t0) \
                 * (self.m_a ** 2. - self.m_c ** 2.) / (2. * self.thick_mem)

    def calc_mem_block_2(self):
        self.ke = self.ko * np.exp(-self.ho / g_par.dict_uni['r']
                                   * (1. / self.t0
                                      - 1. / g_par.dict_case['t_u']))
        a = self.m_c ** 2.
        b = self.m_a ** 2. - self.m_c ** 2.
        self.m0 = np.sqrt(a - b)  # z = 0
        self.m1 = np.sqrt(a)  # z = 1
        self.m_pos0 = -self.ke * self.m0 * 0.5 \
                      + np.sqrt((self.ke * self.m0 * 0.5) ** 2.
                                + self.ke * self.m0)
        self.m_pos1 = -self.ke * self.m1 * 0.5 \
                      + np.sqrt((self.ke * self.m1 * 0.5) ** 2.
                                + self.ke * self.m1)

    def calc_con_overpotential(self):  # nodewise
        self.psi = g_par.dict_uni['r'] * self.t0 / g_par.dict_uni['f'] \
                   * np.log(self.m_pos1 / self.m_pos0)

    def calc_mem_resistivity(self):  # nodewise
        self.omega = (self.mem_bas_r - self.mem_acl_r * self.t0) * 1.e-4
        self.omega_a = self.omega / self.cathode.channel.plane_dx
        self.omega_ele = g_func.calc_elements(self.omega)

    def calc_mem_resistivity_goessling(self):
        lama = np.full(g_par.dict_case['nodes'], 0.)
        res_t = np.exp(self.fac_m * 1000. / self.t0 + self.fac_n)
        r_avg = (self.cathode.humidity + self.anode.humidity) * 0.5
        for q in range(g_par.dict_case['nodes']):
            if r_avg[q] > 0:
                lama[q] = 0.3 + 6. * r_avg[q] * (1. - np.tanh(r_avg[q] - 0.5))\
                          + 3.9 * np.sqrt(r_avg[q])\
                          * (1. + np.tanh((r_avg[q] - 0.89) / 0.23))
            else:
                lama[q] = -1. / (r_avg[q] - (3. + 1. / 3.))
        a = -0.007442
        b = 0.006053
        c = 0.0004702
        d = 1.144
        e = 8.
        res_lama = a + b * lama + c * np.exp(d * (lama - e))
        res = res_lama * res_t / (0.005738 * 15. - 0.07192)
        rp = self.thick_mem / res
        self.omega_goessling = 1.e4 * (self.fac_res_basic
                               + rp * self.fac_res_fit)
        # correction not implemented

    def calc_voltage(self):
        self.v_los = + self.omega_ele * self.i\
                     + self.cathode.act_ov\
                     + self.cathode.gde_dif_los\
                     + self.cathode.cat_dif_los

        self.v_los = np.minimum(self.v_los, g_par.dict_case['e_0'])
        self.v = g_par.dict_case['e_0'] - self.v_los
        self.v_th = g_func.calc_nodes_1d(self.v)

    def calc_resistance(self):
        self.resistance = self.v_los / self.i +\
                          2. * g_par.dict_case['plate_resistivity']\
                          * self.cathode.thickness_plate


