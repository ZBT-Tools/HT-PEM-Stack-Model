import numpy as np
import copy
import matrix_database as m_d
import saturation_pressure_vapour as p_sat
import input_parameter as i_p
import global_functions as g_func
import global_parameter as g_par
import cell as cl


class Stack:

    def __init__(self, dict):
        self.cell_numb = dict['cell_numb']
        self.heat_pow = dict['heat_power']
        self.heigth = dict['heigth']
        self.width = dict['width']
        self.kf = dict['dis_dis_fac']
        self.stoi_cat = dict['stoi_cat']
        self.stoi_ano = dict['stoi_ano']
        self.alpha_f_conv = g_par.dict_case['conv_coef']
        self.cool_ch_bc = dict['cool_ch_bc']
        self.h_col = dict['h_col']
        self.m_col = dict['m_flow_col']
        self.cp_col = dict['cp_col']
        self.alpha_cool = dict['alpha_cool']
        self.cell_list =[]
        for i in range(self.cell_numb):
            x = cl.Cell(i_p.cell)
            self.cell_list.append(x)
        self.init_func()
        self.init_arrays()
        self.init_param()

    def init_func(self):
        self.g_cool = self.m_col * self.cp_col
        if self.cool_ch_bc is False:
            self.t = np.full((self.cell_numb, g_par.dict_case['nodes']),
                                 i_p.t_cool_in)
            self.a_cool = np.full(self.cell_numb, self.alpha_cool)
            self.a_cool[0] = 1.e-20
        else:
            self.a_cool = np.full(self.cell_numb+1, self.alpha_cool)
            self.t = np.full((self.cell_numb+1, g_par.dict_case['nodes']),
                                 i_p.t_cool_in)
        self.extent = 2. * (self.width + self.heigth)
        self.cross_area = self.width * self.heigth
        self.h_d = 4. * self.cross_area / self.extent
        self.q_h_in_cat = np.full(self.cell_numb, 0.)
        self.q_h_in_ano = np.full(self.cell_numb, 0.)
        self.q_h_in_cat[-1] = self.stoi_cat \
                             * self.cell_numb \
                             * g_par.dict_case['tar_cd'] \
                             * np.average(self.cell_list[0].cathode.channel.width) \
                             * np.average(self.cell_list[0].cathode.channel.length) \
                             / (4. * g_par.dict_uni['f']) \
                             * (1. + self.cell_list[-1].cathode.n2o2ratio) \
                             * (1. + self.cell_list[-1].cathode.channel.phi
                                * p_sat.water.calc_psat(self.cell_list[-1].cathode.channel.t_in)
                                / self.cell_list[-1].cathode.channel.p_in
                                   - self.cell_list[-1].cathode.channel.phi
                                   * p_sat.water.calc_psat(self.cell_list[-1].cathode.channel.t_in))
        self.q_h_in_ano[-1] = self.stoi_ano \
                             * self.cell_numb \
                             * g_par.dict_case['tar_cd'] \
                             * self.cell_list[-1].anode.channel.width \
                             * self.cell_list[-1].anode.channel.length \
                             / (2. * g_par.dict_uni['f']) \
                             * (1. + self.cell_list[-1].anode.channel.phi
                                * p_sat.water.calc_psat(self.cell_list[-1].anode.channel.t_in)
                                / (self.cell_list[-1].anode.channel.p_in
                                   - self.cell_list[-1].anode.channel.phi
                                   * p_sat.water.calc_psat(self.cell_list[-1].anode.channel.t_in)))
        self.set_stoi(np.full(self.cell_numb, self.stoi_cat),
                      np.full(self.cell_numb, self.stoi_ano))

    def init_param(self):
        self.d_col = 4. * (self.cell_list[0].cathode.channel.width * self.h_col)\
                     / (2. * (self.h_col + self.cell_list[0].cathode.channel.width))
        self.r_alpha_col = 1./(self.a_cool * np.pi * self.cell_list[0].cathode.channel.d_x * self.d_col)
        self.r_alpha_col = 1. / (1. / self.r_alpha_col)
        self.break_programm = False

    def init_arrays(self):
        x = np.full(g_par.dict_case['nodes']-1, g_par.dict_case['tar_cd'])
        #x = np.linspace(g_par.dict_case['tar_cd']*1.1,g_par.dict_case['tar_cd']*0.9, g_par.dict_case['nodes']-1)
        y = []
        for q in range(self.cell_numb): y.append(x)
        self.i = np.array(y)
        self.q_h_out_cat = np.full(self.cell_numb, 0.)
        self.q_h_out_ano = np.full(self.cell_numb, 0.)
        self.t_h_out_cat = np.full(self.cell_numb, 0.)
        self.t_h_out_ano = np.full(self.cell_numb, 0.)
        self.r_mix_h_in_cat = np.full(self.cell_numb, 0.)
        self.r_mix_h_out_cat = np.full(self.cell_numb, 0.)
        self.r_mix_h_in_ano = np.full(self.cell_numb, 0.)
        self.r_mix_h_out_ano = np.full(self.cell_numb, 0.)
        self.p_h_in_cat = np.full(self.cell_numb, i_p.p_cat_in)
        self.p_h_out_cat = np.full(self.cell_numb, i_p.p_cat_in)
        self.p_h_out_ano = np.full(self.cell_numb, i_p.p_ano_in)
        self.p_h_in_ano = np.full(self.cell_numb, i_p.p_ano_in)
        self.q_x_cat = np.full(self.cell_numb,
                               self.q_h_in_cat[-1] / self.cell_numb)
        self.q_x_ano = np.full(self.cell_numb,
                               self.q_h_in_ano[-1] / self.cell_numb)
        self.b = m_d.b(g_par.dict_case['nodes']-1,
                       self.cell_numb,
                       self.cell_list[0].cathode.channel.d_x)
        if self.cell_numb >=3:
            self.c = m_d.c(g_par.dict_case['nodes'] - 1,
                           self.cell_numb)
        self.stack_r()
        self.stack_resistance()
        self.zero = np.full(self.cell_numb, 0.)
        self.r_zero = np.full(self.cell_numb, 1.e50)
        self.r_alpha_gp = np.full(self.cell_numb, 0.)
        self.r_alpha_pp = np.full(self.cell_numb, 0.)
        self.r_alpha_gm = np.full(self.cell_numb, 0.)
        self.r_alpha_gegm = np.full(self.cell_numb, 1.e50)
        self.r_alpha_gegp = np.full(self.cell_numb, 1.e50)
        self.r_alpha_gepp = np.full(self.cell_numb, 1.e50)
        self.cool_m1 = m_d.col_mat_m1(g_par.dict_case['nodes'])
        fac = 1.
        for q, item in enumerate(self.cell_list):
            self.r_alpha_gp[q] = 2. / (self.alpha_f_conv
                                    * self.cell_list[q].cathode.channel.d_x
                                    * (self.cell_list[q].cathode.thickness_plate
                                       + self.cell_list[q].cathode.thickness_gde)) * fac
            self.r_alpha_gm[q] = 2. / (self.alpha_f_conv
                                    * self.cell_list[q].cathode.channel.d_x
                                    * (self.cell_list[q].cathode.thickness_plate
                                       + self.cell_list[q].thickness_mea)) * fac
            self.r_alpha_pp[q] = 1. / (self.alpha_f_conv
                                    * self.cell_list[q].cathode.channel.d_x
                                    * self.cell_list[q].cathode.thickness_plate) * fac
            self.r_alpha_gegm[q] = 2. / (self.alpha_f_conv
                                         * self.cell_list[q].cathode.channel.width
                                         * (self.cell_list[q].cathode.thickness_gde
                                            + self.cell_list[q].thickness_mea))
            self.r_alpha_gegp[q] = 2. / (self.alpha_f_conv
                                         * self.cell_list[q].cathode.channel.width
                                         * (self.cell_list[q].cathode.thickness_gde
                                            + self.cell_list[q].cathode.thickness_plate))
            self.r_alpha_gepp[q] = 1. / (self.alpha_f_conv
                                         * self.cell_list[q].cathode.channel.width
                                         * self.cell_list[q].cathode.thickness_plate)
        self.r_alpha_gegp = np.full(self.cell_numb, 1.e50)
        self.r_alpha_gegm = np.full(self.cell_numb, 1.e50)
        self.r_alpha_gepp = np.full(self.cell_numb, 1.e50)

    def update(self):
        for j in range(self.cell_numb):
            self.cell_list[j].set_i(self.i[j, :])
            self.cell_list[j].update()
            if self.cell_list[j].break_programm is True:
                self.break_programm = True
                break
        if self.break_programm is False:
            self.update_temperatur_coupling()
            if self.cell_numb > 1:
                self.update_flows()
            self.i_old = copy.deepcopy(self.i)
            self.update_electrical_coupling()

    def update_flows(self):
        self.sum_header_flows()
        self.calc_header_temp()
        self.calc_header_velocity()
        self.calc_header_r_mix()
        self.calc_header_roh()
        self.calc_header_reynolds_numb()
        self.calc_header_fanning_friction_factor()
        self.calc_header_p_out()
        self.calc_ref_p_drop()
        self.calc_ref_perm()
        self.calc_header_p_in()
        self.calc_flow_distribution_factor()
        self.calc_new_ref_p_drop()
        self.calc_header_p_in()
        self.calc_flow_distribution_factor()
        self.calc_new_cell_flows()
        self.calc_new_cell_stoi()
        self.set_stoi(self.new_stoi_cat, self.new_stoi_ano)
        self.set_p()

    def update_electrical_coupling(self):
        self.stack_v()
        self.stack_v_los()
        self.stack_dv()
        #self.calc_i()
        self.calc_i_v_interface()

    def update_temperatur_coupling(self):
        self.stack_alpha()
        self.calc_coolant_channel_t()
        self.calc_gas_channel_t()
        self.calc_layer_t()

    def set_i(self, i):
        self.i = i

    def stack_resistance(self):
        d_p = []
        for q, item in enumerate(self.cell_list):
               d_p.append(self.cell_list[q].cathode.thickness_plate)
        self.resistance = g_par.dict_case['plate_resistivity'] * np.average(d_p)

    def stack_r(self):
        r_p, r_g, r_m = [], [], []
        r_pp, r_gp, r_gm = [], [], []
        for i, item in enumerate(self.cell_list):
            r_p = np.hstack((r_p, self.cell_list[i].r_p))
            r_g = np.hstack((r_g, self.cell_list[i].r_g))
            r_m = np.hstack((r_m, self.cell_list[i].r_m))
            r_pp = np.hstack((r_pp, self.cell_list[i].r_pp))
            r_gp = np.hstack((r_gp, self.cell_list[i].r_gp))
            r_gm = np.hstack((r_gm, self.cell_list[i].r_gm))
        self.r_m = r_m
        self.r_g = r_g
        self.r_p = r_p
        self.r_pp = r_pp
        self.r_gp = r_gp
        self.r_gm = r_gm

    def stack_alpha(self):
        r_alpha_cat, r_alpha_ano = [], []
        for i, item in enumerate(self.cell_list):
            r_alpha_cat = np.hstack((r_alpha_cat, self.cell_list[i].cathode.r_ht_coef_a))
            r_alpha_ano = np.hstack((r_alpha_ano, self.cell_list[i].anode.r_ht_coef_a))
        self.r_alpha_cat = r_alpha_cat
        #self.r_alpha_cat = np.full(self.cell_numb, 1.e50)
        self.r_alpha_ano = r_alpha_ano
        #self.r_alpha_ano = np.full(self.cell_numb, 1.e50)

    def stack_v(self):
        var = []
        for i, item in enumerate(self.cell_list):
            var = np.hstack((var, self.cell_list[i].v))
        self.v = var
        # running

    def stack_v_los(self):
        var = []
        for i, item in enumerate(self.cell_list):
            var = np.hstack((var, self.cell_list[i].v_los))
        self.v_los = var

    def stack_c_r(self):
        var = []
        for i, item in enumerate(self.cell_list):
            var = np.hstack((var, self.cell_list[i].resistance))
        self.stack_cell_r = var

    def stack_dv(self):
        var = []
        for i, item in enumerate(self.cell_list):
            var = np.hstack((var, self.cell_list[i].dv))
        self.dv = var

    def set_stoi(self, stoi_cat, stoi_ano):
        for q in range(self.cell_numb):
            self.cell_list[q].cathode.stoi = stoi_cat[q]
            self.cell_list[q].anode.stoi = stoi_ano[q]

    def calc_stoi(self, q, val):
        if val is 4:
            b = 1 + self.cell_list[-1].cathode.n2o2ratio
            c = 1 + self.cell_list[-1].cathode.channel.phi \
                * p_sat.water.calc_psat(self.cell_list[-1].cathode.channel.t_in) \
                / (self.cell_list[-1].cathode.channel.p_in
                   - self.cell_list[-1].cathode.channel.phi
                   * p_sat.water.calc_psat(self.cell_list[-1].cathode.channel.t_in))
            a = g_par.dict_case['tar_cd']\
                * self.cell_list[-1].cathode.channel.width\
                * self.cell_list[-1].cathode.channel.heigth\
                / (val * g_par.dict_uni['f'])
        else:
            c = 1. + self.cell_list[-1].anode.channel.phi \
                * p_sat.water.calc_psat(self.cell_list[-1].anode.channel.t_in) \
                / (self.cell_list[-1].anode.channel.p_in
                   - self.cell_list[-1].anode.channel.phi
                   * p_sat.water.calc_psat(self.cell_list[-1].anode.channel.t_in))
            a = g_par.dict_case['tar_cd']\
                * self.cell_list[-1].anode.channel.width\
                * self.cell_list[-1].anode.channel.heigth\
                / (val * g_par.dict_uni['f'])
        return q /(a * b * c)

    def set_p(self):
        for q in range(self.cell_numb):
            self.cell_list[q].cathode.channel.p_in = self.p_h_in_cat[q]
            self.cell_list[q].anode.channel.p_in = self.p_h_in_ano[q]

    def sum_header_flows(self): # ok change here for parallel flow
        self.q_h_in_cat[0] = self.cell_list[0].cathode.q_sum[0]
        self.q_h_out_cat[0] = self.cell_list[0].cathode.q_sum[-1]
        self.q_h_in_ano[0] = self.cell_list[0].anode.q_sum[-1]
        self.q_h_out_ano[0] = self.cell_list[0].anode.q_sum[0]
        for q in range(1, self.cell_numb):
            self.q_h_in_cat[q] = self.q_h_in_cat[q - 1] \
                                 + self.cell_list[q].cathode.q_sum[0]
            self.q_h_out_cat[q] = self.q_h_out_cat[q - 1] \
                                  + self.cell_list[q].cathode.q_sum[-1]
            self.q_h_in_ano[q] = self.q_h_in_ano[q - 1] \
                                 + self.cell_list[q].anode.q_sum[-1]
            self.q_h_out_ano[q] = self.q_h_out_ano[q - 1] \
                                  + self.cell_list[q].anode.q_sum[0]

    def calc_header_temp(self):
        self.t_h_out_cat[0] = self.cell_list[0].cathode.t_gas[-1]
        self.t_h_out_ano[0] = self.cell_list[0].anode.t_gas[0]
        for q in range(1, self.cell_numb):
            self.t_h_out_cat[q] = (self.t_h_out_cat[q - 1] * self.q_h_out_cat[q - 1]
                                   + self.cell_list[q].cathode.q_sum[-1]
                                   * self.cell_list[q].cathode.t2[-1]) \
                                  / self.q_h_out_cat[q]
            self.t_h_out_ano[q] = (self.t_h_out_ano[q - 1] * self.q_h_out_ano[q - 1]
                                   + self.cell_list[q].anode.q_sum[0]
                                   * self.cell_list[q].anode.t2[0]) \
                                  / self.q_h_out_ano[q]

    def calc_header_velocity(self):
        self.v_h_in_cat = self.q_h_in_cat \
                          * g_par.dict_uni['r'] \
                          * self.cell_list[-1].cathode.channel.t_in \
                          / (self.p_h_in_cat * self.cross_area)
        self.v_h_in_ano = self.q_h_in_ano \
                          * g_par.dict_uni['r'] \
                          * self.cell_list[-1].anode.channel.t_in \
                          / (self.p_h_in_ano * self.cross_area)
        self.v_h_out_cat = self.q_h_out_cat \
                           * g_par.dict_uni['r'] \
                           * self.t_h_out_cat \
                           / (self.p_h_out_cat * self.cross_area)
        self.v_h_out_ano = self.q_h_out_ano \
                           * g_par.dict_uni['r'] \
                           * self.t_h_out_ano \
                           / (self.p_h_out_ano * self.cross_area)

    def calc_header_r_mix(self):
        self.r_mix_h_in_cat = np.full(self.cell_numb, self.cell_list[0].cathode.r_mix[0])
        self.r_mix_h_in_ano = np.full(self.cell_numb, self.cell_list[0].anode.r_mix[-1])
        self.r_mix_h_out_cat[0] = self.cell_list[0].cathode.r_mix[-1]
        self.r_mix_h_out_ano[0] = self.cell_list[0].anode.r_mix[0]
        for q in range(1, self.cell_numb):
            self.r_mix_h_out_cat[q] = (self.r_mix_h_out_cat[q - 1] * self.q_h_out_cat[q - 1]
                                       + self.cell_list[q].cathode.q_sum[-1]
                                       * self.cell_list[q].cathode.r_mix[-1]) \
                                      / self.q_h_out_cat[q]
            self.r_mix_h_out_ano[q] = (self.r_mix_h_out_ano[q - 1] * self.q_h_out_ano[q - 1]
                                       + self.cell_list[q].anode.q_sum[0]
                                       * self.cell_list[q].anode.r_mix[0]) \
                                      / self.q_h_out_ano[q]

    def calc_header_roh(self):
        self.roh_h_in_cat = g_func.calc_rho(self.p_h_in_cat,
                                            self.r_mix_h_in_cat,
                                            self.cell_list[-1].cathode.channel.t_in)
        self.roh_h_in_ano = g_func.calc_rho(self.p_h_in_ano,
                                            self.r_mix_h_in_ano,
                                            self.cell_list[-1].anode.channel.t_in)
        self.roh_h_out_cat = g_func.calc_rho(self.p_h_out_cat,
                                             self.r_mix_h_out_cat,
                                             self.t_h_out_cat)
        self.roh_h_out_ano = g_func.calc_rho(self.p_h_out_ano,
                                             self.r_mix_h_out_ano,
                                             self.t_h_out_ano)

    def calc_header_reynolds_numb(self):
        for q in range(self.cell_numb):
            self.re_h_in_cat = g_func.calc_re(self.roh_h_in_cat,
                                              self.v_h_in_cat,
                                              self.h_d,
                                              self.cell_list[q].cathode.visc_mix[0])
            self.re_h_in_ano = g_func.calc_re(self.roh_h_in_ano,
                                              self.v_h_in_ano,
                                              self.h_d,
                                              self.cell_list[q].anode.visc_mix[-1])
            self.re_h_out_cat = g_func.calc_re(self.roh_h_out_cat,
                                               self.v_h_out_cat,
                                               self.h_d,
                                               self.cell_list[q].cathode.visc_mix[-1])
            self.re_h_out_ano = g_func.calc_re(self.roh_h_out_ano,
                                               self.v_h_out_ano,
                                               self.h_d,
                                               self.cell_list[q].anode.visc_mix[0])

    def calc_header_fanning_friction_factor(self):
        self.f_h_in_cat = g_func.calc_fanning_friction_factor(self.re_h_in_cat)
        self.f_h_in_ano = g_func.calc_fanning_friction_factor(self.re_h_in_ano)
        self.f_h_out_cat = g_func.calc_fanning_friction_factor(self.re_h_out_cat)
        self.f_h_out_ano = g_func.calc_fanning_friction_factor(self.re_h_out_ano)

    def calc_header_p_out(self):  # change for parallel flow
        self.p_h_out_cat[self.cell_numb - 1] = g_par.dict_case['header_p_in_cat']
        self.p_h_out_ano[self.cell_numb - 1] = g_par.dict_case['header_p_in_ano']
        for q in range(self.cell_numb - 1, 0, -1):
            self.p_h_out_cat[q - 1] = self.p_h_out_cat[q] \
                                      + g_func.calc_header_pressure_drop(self.roh_h_out_cat[q-1],
                                                                         self.v_h_out_cat[q],
                                                                         self.v_h_out_cat[q-1],
                                                                         self.f_h_out_cat[q-1],
                                                                         self.kf,
                                                                         self.cell_list[q].heigth,
                                                                         self.h_d)
            self.p_h_out_ano[q - 1] = self.p_h_out_ano[q] \
                                      + g_func.calc_header_pressure_drop(self.roh_h_out_ano[q-1],
                                                                         self.v_h_out_ano[q],
                                                                         self.v_h_out_ano[q-1],
                                                                         self.f_h_out_ano[q-1],
                                                                         self.kf,
                                                                         self.cell_list[q].heigth,
                                                                         self.h_d)

    def calc_ref_perm(self):
        self.perm_ano = np.average(self.cell_list[0].anode.visc_mix)\
                        / self.cell_list[0].anode.channel.cross_area\
                        * self.cell_list[0].anode.channel.length\
                        * np.average(self.cell_list[0].anode.q_sum)\
                        / self.ref_p_drop_ano
        self.perm_cat = np.average(self.cell_list[0].cathode.visc_mix)\
                        / self.cell_list[0].cathode.channel.cross_area\
                        * self.cell_list[0].cathode.channel.length\
                        * np.average(self.cell_list[0].cathode.q_sum)\
                        / self.ref_p_drop_cat
        self.ref_p_new_drop_cat = np.average(self.cell_list[0].cathode.q_sum)\
                         / self.cell_list[0].cathode.channel.cross_area\
                         * np.average(self.cell_list[0].cathode.visc_mix)\
                         * self.cell_list[0].cathode.channel.length\
                         / self.perm_cat
        self.p_adjusting_factor_cat = self.ref_p_drop_cat / self.ref_p_new_drop_cat
        self.ref_p_new_drop_ano = np.average(self.cell_list[0].anode.q_sum)\
                         / self.cell_list[0].anode.channel.cross_area\
                         * np.average(self.cell_list[0].anode.visc_mix)\
                         * self.cell_list[0].anode.channel.length\
                         / self.perm_ano
        self.p_adjusting_factor_ano = self.ref_p_drop_ano / self.ref_p_new_drop_ano

    def calc_ref_p_drop(self):
        self.ref_p_drop_cat = self.cell_list[0].cathode.p[0]\
                              - self.cell_list[0].cathode.p[-1]
        self.ref_p_drop_ano = self.cell_list[0].anode.p[-1]\
                              - self.cell_list[0].anode.p[0]

    def calc_new_ref_p_drop(self):
        self.ref_p_drop_cat = self.q_h_in_cat[self.cell_numb-1]\
                              / np.sum(self.alpha_cat)\
                              * np.average(self.cell_list[0].cathode.visc_mix)\
                              * self.cell_list[0].cathode.channel.length\
                              / self.cell_list[0].cathode.channel.cross_area\
                              / self.perm_cat\
                              * self.p_adjusting_factor_cat
        self.ref_p_drop_ano = self.q_h_in_ano[self.cell_numb-1]\
                              / np.sum(self.alpha_ano)\
                              * np.average(self.cell_list[0].anode.visc_mix)\
                              * self.cell_list[0].anode.channel.length\
                              / self.cell_list[0].anode.channel.cross_area\
                              / self.perm_ano\
                              * self.p_adjusting_factor_ano

    def calc_header_p_in(self):
        self.p_h_in_cat[0] = self.ref_p_drop_cat + self.p_h_out_cat[0]
        self.p_h_in_ano[0] = self.ref_p_drop_ano + self.p_h_out_ano[0]
        for q in range(1, self.cell_numb):
            self.p_h_in_cat[q] = self.p_h_in_cat[q - 1] \
                                 + g_func.calc_header_pressure_drop(self.roh_h_in_cat[q],
                                                                    self.v_h_in_cat[q - 1],
                                                                    self.v_h_in_cat[q],
                                                                    self.f_h_in_cat[q],
                                                                    self.kf,
                                                                    self.cell_list[q].heigth,
                                                                    self.h_d)
            self.p_h_in_ano[q] = self.p_h_in_ano[q - 1] \
                                 + g_func.calc_header_pressure_drop(self.roh_h_in_ano[q],
                                                                    self.v_h_in_ano[q - 1],
                                                                    self.v_h_in_ano[q],
                                                                    self.f_h_in_ano[q],
                                                                    self.kf,
                                                                    self.cell_list[q].heigth,
                                                                    self.h_d)

    def calc_flow_distribution_factor(self):
        self.alpha_cat = (self.p_h_in_cat - self.p_h_out_cat)\
                         / self.ref_p_drop_cat
        self.alpha_ano = (self.p_h_in_ano - self.p_h_out_ano)\
                         / self.ref_p_drop_ano

    def calc_new_cell_flows(self):
        self.q_x_cat = (self.p_h_in_cat - self.p_h_out_cat)\
                        * self.perm_cat\
                        * self.cell_list[0].cathode.channel.cross_area\
                        / np.average(self.cell_list[0].cathode.visc_mix)\
                        / self.cell_list[0].cathode.channel.length\
                        / self.p_adjusting_factor_cat
        self.q_x_ano = (self.p_h_in_ano - self.p_h_out_ano) \
                       * self.perm_ano \
                       * self.cell_list[0].anode.channel.cross_area \
                       / np.average(self.cell_list[0].anode.visc_mix) \
                       / self.cell_list[0].anode.channel.length \
                       / self.p_adjusting_factor_ano

    def calc_new_cell_stoi(self):
        self.new_stoi_cat = self.q_x_cat /\
                            (self.q_h_in_cat[-1] / self.cell_numb)\
                            * self.stoi_cat
        self.new_stoi_ano = self.q_x_ano\
                            / (self.q_h_in_ano[-1] / self.cell_numb)\
                            * self.stoi_ano
        min_cat = np.amin(self.new_stoi_cat)
        min_ano = np.amin(self.new_stoi_ano)
        if min_cat <= 1.15:
            print('Warning, raise the cathode stoichiometry')
            self.break_programm = True
        if min_ano <= 1.15:
            print('Warning, raise the anode stoichiometry')
            self.break_programm = True
        #print(min_cat, min_ano, min_min, 'min_min', g_par.dict_case['tar_cd'])
        #if min_min <= 1.1:
          #  stoi_fac = i_p.stoi_min / min_min
         #   g_par.dict_case['tar_cd'] = g_par.dict_case['tar_cd'] / stoi_fac
           # self.i = self.i / stoi_fac
           # print(self.new_stoi_cat, 'stoi_before_change')
           # self.new_stoi_cat = self.new_stoi_cat * stoi_fac
           # print(self.new_stoi_cat, 'stoi_after_change')
           # self.new_stoi_ano = self.new_stoi_ano * stoi_fac
           # self.stoi_ano = self.stoi_ano * stoi_fac
           # self.stoi_cat = self.new_stoi_cat * stoi_fac

    def calc_i(self):
        self.stack_c_r()
        r_cell_ele = (.5 * (self.stack_cell_r[:-(g_par.dict_case['nodes']-1)]
                            + self.stack_cell_r[(g_par.dict_case['nodes']-1):]))
        r_mat = m_d.electrical_mat(g_par.dict_case['nodes']-1,
                                   self.cell_numb,
                                   self.stack_cell_r,
                                   r_cell_ele,
                                   g_par.dict_case['plate_resistivity']
                                   * self.cell_list[0].cathode.channel.d_x * 0.5)
        self.v_plate = np.sum(self.v_los)/(g_par.dict_case['nodes']-1)
        self.neuman_vec = np.hstack((-2.*self.v_plate / self.stack_cell_r[:g_par.dict_case['nodes']-1],
                                     np.full((self.cell_numb-1)*(g_par.dict_case['nodes']-1), 0.)))
        self.v_neu = np.array(np.linalg.tensorsolve(r_mat, self.neuman_vec))
        self.v_neu = np.hstack((np.full((g_par.dict_case['nodes']-1), self.v_plate),
                                self.v_neu[0],
                                np.full((g_par.dict_case['nodes']-1),0.)))
        self.v_dif = self.v_neu[:-(g_par.dict_case['nodes']-1)] - self.v_neu[(g_par.dict_case['nodes']-1):]
        self.r_vec_1d = np.hstack((self.stack_cell_r[:g_par.dict_case['nodes']-1]*.5,
                                   r_cell_ele,self.stack_cell_r[-g_par.dict_case['nodes']+1:]*.5))
        self.i_vec = self.v_dif /self.r_vec_1d
        self.i_ele_vec = .5 * (self.i_vec[:-(g_par.dict_case['nodes']-1)]
                          + self.i_vec[(g_par.dict_case['nodes']-1):])
        self.i = g_func.toarray(self.i_ele_vec, self.cell_numb, g_par.dict_case['nodes']-1)
        self.i = self.i / np.average(self.i) * g_par.dict_case['tar_cd']

    def calc_i_v_interface(self):
        self.stack_c_r()
        r_mat = m_d.electrical_mat_interface_v(g_par.dict_case['nodes']-1,
                                   self.cell_numb,
                                   self.stack_cell_r,
                                   g_par.dict_case['plate_resistivity']
                                   * self.cell_list[0].cathode.channel.d_x * 0.5)
        self.v_plate = np.sum(self.v_los) / (g_par.dict_case['nodes']-1)
        self.neuman_vec = np.hstack((-self.v_plate / self.stack_cell_r[:g_par.dict_case['nodes']-1],
                                     np.full((self.cell_numb-2)*(g_par.dict_case['nodes']-1), 0.)))
        self.v_neu = np.array(np.linalg.tensorsolve(r_mat, self.neuman_vec))
        self.v_neu = np.hstack((np.full((g_par.dict_case['nodes']-1), self.v_plate),
                                self.v_neu[0],
                                np.full((g_par.dict_case['nodes']-1), 0.)))
        self.v_dif = self.v_neu[:-(g_par.dict_case['nodes'] - 1)]\
                     - self.v_neu[(g_par.dict_case['nodes'] - 1):]
        self.i_vec = self.v_dif / self.stack_cell_r
        self.i = g_func.toarray(self.i_vec, self.cell_numb, g_par.dict_case['nodes']-1)
        self.i = self.i / np.average(self.i) * g_par.dict_case['tar_cd']

    def calc_gas_channel_t(self):
        for q, item in enumerate(self.cell_list):
            t_gas = np.full(g_par.dict_case['nodes'], 0.)
            t_gas[0] = self.cell_list[q].cathode.channel.t_in
            for w in range(1, g_par.dict_case['nodes']):
                self.cell_list[q].cathode.t_gas[w] = (self.cell_list[q].cathode.t_gas[w-1]
                                                      - self.cell_list[q].t2e[w-1])\
                                                     * np.exp(-1./(self.r_alpha_cat[q]
                                                                   * self.cell_list[q].cathode.g_full_e[w-1]))\
                                                     + self.cell_list[q].t2e[w-1]\
                                                     + self.cell_list[q].anode.cpe[w - 1]\
                                                     * self.cell_list[q].anode.m_reac_flow_delta[w - 1]\
                                                     * self.cell_list[q].anode.t_gas_e[w - 1]
            for w in range(g_par.dict_case['nodes']-2, -1, -1):
                self.cell_list[q].anode.t_gas[w] = (self.cell_list[q].anode.t_gas[w+1]
                                                    - self.cell_list[q].t5e[w])\
                                                   * np.exp(-1./(self.r_alpha_ano[q]
                                                                 * self.cell_list[q].anode.g_full_e[w]))\
                                                   + self.cell_list[q].t5e[w]\
                                                   - self.cell_list[q].anode.cpe[w]\
                                                   * self.cell_list[q].anode.m_reac_flow_delta[w]\
                                                   * self.cell_list[q].anode.t_gas_e[w]

    def calc_coolant_channel_t(self):
        for q, item in enumerate(self.cell_list):
            var1 = 1./(self.r_alpha_col[q]* self.g_cool)
            for w in range(1, g_par.dict_case['nodes']):
                self.t[q, w] = (self.t[q, w-1] - self.cell_list[q].t3e[w-1])\
                              * np.exp(-1. * var1) + self.cell_list[q].t3e[w-1]
        if self.cool_ch_bc is True:
            var1 = 0.5 / (self.r_alpha_col[self.cell_numb-1] * self.g_cool)
            for w in range(1,g_par.dict_case['nodes']):
                self.t[self.cell_numb, w] = (self.t[self.cell_numb, w-1]
                                             - self.cell_list[q].t5e[w-1])\
                                            * np.exp(-1. * var1)\
                                            + self.cell_list[q].t5e[w-1]

    def calc_layer_t(self):
        self.I = g_func.iepolate_nodes(self.i) * self.cell_list[0].cathode.channel.plane_dx
        temp_mat = m_d.t_mat_no_bc_col(g_par.dict_case['nodes'],
                                         self.cell_numb,
                                         self.r_g,
                                         self.r_m,
                                         self.r_p,
                                         self.r_gp,
                                         self.r_gm,
                                         self.r_pp,
                                         self.r_alpha_col,
                                         self.r_alpha_cat,
                                         self.r_alpha_ano,
                                         self.r_alpha_gm,
                                         self.r_alpha_gp,
                                         self.r_alpha_pp,
                                         self.r_alpha_gegm,
                                         self.r_alpha_gegp,
                                         self.r_alpha_gepp,
                                         self.cool_ch_bc)
        r_side = []
        if self.cell_numb > 1:
            for q, item in enumerate(self.cell_list):
                for w in range(g_par.dict_case['nodes']):
                    if w is 0 or w is g_par.dict_case['nodes']-1: # bc nodes
                        if q is 0:
                            r_side.append(+ .5 /self.r_alpha_ano[q] * self.cell_list[q].anode.t_gas[w]
                                          + .5 /self.r_alpha_gp[q] * g_par.dict_case['tu']
                                          + 1./self.r_alpha_gegp[q] * g_par.dict_case['tu'])
                            r_side.append(self.cell_list[q].omega_a[w] * self.I[q, w] ** 2 * 0.25
                                          +.5 /self.r_alpha_gm[q] * g_par.dict_case['tu']
                                          + 1./self.r_alpha_gegm[q] * g_par.dict_case['tu'])
                            r_side.append((g_par.dict_case['vtn'] - self.cell_list[q].v_th[w] - self.cell_list[q].omega_a[w]
                                           * self.I[q,w] * 0.5) * self.I[q,w] * 0.5
                                          +.5/self.r_alpha_gm[q] * g_par.dict_case['tu']
                                          + 1./self.r_alpha_gegm[q] * g_par.dict_case['tu'])
                            r_side.append(+ .5/self.r_alpha_gp[q] * g_par.dict_case['tu']
                                          + 1./self.r_alpha_gegp[q] * g_par.dict_case['tu']
                                          + 1./self.r_alpha_cat[q] * self.cell_list[q].cathode.t_gas[w] * 0.5
                                          + 0.5 * self.cell_list[q].cathode.gamma[w] * g_par.dict_uni['h_vap'])
                            if self.cool_ch_bc is True:
                                r_side.append(-self.heat_pow * .5
                                              - .5 / self.r_alpha_pp[q] * g_par.dict_case['tu']
                                              - 1. / self.r_alpha_gepp[q] * g_par.dict_case['tu']
                                              -.5 / self.r_alpha_col[q] * self.t[q, w])
                            else:
                                r_side.append(-self.heat_pow * .5
                                              - .5/self.r_alpha_pp[q] * g_par.dict_case['tu']
                                              - 1./self.r_alpha_gepp[q] * g_par.dict_case['tu'])
                        elif q is self.cell_numb-1:
                            if self.cool_ch_bc is True:
                                r_side.append(+ self.heat_pow * .5
                                              + .5 / self.r_alpha_ano[q] * self.cell_list[q].anode.t_gas[w]
                                              + .5 / self.r_alpha_gp[q] * g_par.dict_case['tu']
                                              + 1. / self.r_alpha_gegp[q] * g_par.dict_case['tu']
                                              + .5 / self.r_alpha_col[q] * self.t[q+1, w])
                            else:
                                r_side.append(+ self.heat_pow * .5
                                              + .5 / self.r_alpha_ano[q] * self.cell_list[q].anode.t_gas[w]
                                              + .5 / self.r_alpha_gp[q] * g_par.dict_case['tu']
                                              + 1. / self.r_alpha_gegp[q] * g_par.dict_case['tu'])
                            r_side.append(self.cell_list[q].omega_a[w] * self.I[q, w] ** 2 * 0.25
                                          + .5/self.r_alpha_gm[q] * g_par.dict_case['tu']
                                          + 1./self.r_alpha_gegm[q] * g_par.dict_case['tu'])
                            r_side.append((g_par.dict_case['vtn'] - self.cell_list[q].v_th[w] - self.cell_list[q].omega_a[w] *
                                           self.I[q, w] * 0.5) * self.I[q, w] * 0.5
                                          + .5/self.r_alpha_gm[q] * g_par.dict_case['tu']
                                          + 1./self.r_alpha_gegm[q] * g_par.dict_case['tu'])
                            r_side.append(+ .5/self.r_alpha_gp[q] * g_par.dict_case['tu']
                                          + 1./self.r_alpha_gegp[q] * g_par.dict_case['tu']
                                          + 1./self.r_alpha_cat[q] * self.cell_list[q].cathode.t_gas[w] * 0.5
                                          + 0.5 * self.cell_list[q].cathode.gamma[w] * g_par.dict_uni['h_vap'])
                            r_side.append(-.5 / self.r_alpha_col[q] * self.t[q, w]
                                          - .5/self.r_alpha_pp[q] * g_par.dict_case['tu']
                                          - 1./self.r_alpha_gepp[q] * g_par.dict_case['tu'])
                        else:
                            r_side.append(+ .5/self.r_alpha_ano[q] * self.cell_list[q].anode.t_gas[w]
                                          + .5/self.r_alpha_gp[q] * g_par.dict_case['tu']
                                          + 1./self.r_alpha_gegp[q] * g_par.dict_case['tu'])
                            r_side.append(self.cell_list[q].omega_a[w] * self.I[q, w] ** 2 * 0.25
                                          + .5/self.r_alpha_gm[q] * g_par.dict_case['tu']
                                          + 1./self.r_alpha_gegm[q] * g_par.dict_case['tu'])
                            r_side.append((g_par.dict_case['vtn'] - self.cell_list[q].v_th[w] - self.cell_list[q].omega_a[w] *
                                           self.I[q, w] * 0.5) * self.I[q, w] * 0.5
                                          +.5/self.r_alpha_gm[q] * g_par.dict_case['tu']
                                          + 1./self.r_alpha_gegm[q] * g_par.dict_case['tu'])
                            r_side.append(+ 1./self.r_alpha_cat[q] * self.cell_list[q].cathode.t_gas[w] * 0.5
                                          +.5 /self.r_alpha_gp[q] * g_par.dict_case['tu']
                                          + 1./self.r_alpha_gegp[q] * g_par.dict_case['tu']
                                          + 0.5 * self.cell_list[q].cathode.gamma[w] * g_par.dict_uni['h_vap'])
                            r_side.append(-.5 / self.r_alpha_col[q] * self.t[q, w]
                                          - .5/self.r_alpha_pp[q] * g_par.dict_case['tu']
                                          - 1. / self.r_alpha_gepp[q] * g_par.dict_case['tu'])
                    else:
                        if q is 0:
                            r_side.append(+ 1./self.r_alpha_ano[q] * self.cell_list[q].anode.t_gas[w]
                                          + 1./self.r_alpha_gp[q] * g_par.dict_case['t_u'])
                            r_side.append(self.cell_list[q].omega_a[w] * self.I[q, w] ** 2 * 0.5
                                          + 1./self.r_alpha_gm[q] * g_par.dict_case['tu'])
                            r_side.append((g_par.dict_case['vtn'] - self.cell_list[q].v_th[w] - self.cell_list[q].omega_a[w] *
                                           self.I[q, w] * 0.5) * self.I[q, w]
                                          + 1./self.r_alpha_gm[q] * g_par.dict_case['tu'])
                            r_side.append(+ 1./self.r_alpha_gp[q] * g_par.dict_case['tu'] + 1./self.r_alpha_cat[q]
                                          * self.cell_list[q].cathode.t_gas[w]
                                          + self.cell_list[q].cathode.gamma[w] * g_par.dict_uni['h_vap'])
                            if self.cool_ch_bc is True:
                                r_side.append(-self.heat_pow - 1. / self.r_alpha_pp[q] * g_par.dict_case['tu']
                                              - self.t[q,w] * 1./self.r_alpha_col[q])
                            else:
                                r_side.append(-self.heat_pow - 1./self.r_alpha_pp[q] * g_par.dict_case['tu'])
                        elif q is self.cell_numb-1:
                            if self.cool_ch_bc is True:
                                r_side.append(+ self.heat_pow
                                              + 1./self.r_alpha_ano[q] * self.cell_list[q].anode.t_gas[w]
                                              + 1./self.r_alpha_gp[q] * g_par.dict_case['t_u'] + self.t[q+1,w] * 1./self.r_alpha_col[q])
                            else:
                                r_side.append(+ self.heat_pow
                                              + 1. / self.r_alpha_ano[q] * self.cell_list[q].anode.t_gas[w]
                                              + 1. / self.r_alpha_gp[q] * g_par.dict_case['t_u'])
                            r_side.append(self.cell_list[q].omega_a[w] * self.I[q, w] ** 2 * 0.5
                                          + 1./self.r_alpha_gm[q] * g_par.dict_case['tu'])
                            r_side.append((g_par.dict_case['vtn'] - self.cell_list[q].v_th[w] - self.cell_list[q].omega_a[w] *
                                           self.I[q, w] * 0.5) * self.I[q, w]
                                          + 1./self.r_alpha_gm[q] * g_par.dict_case['tu'])
                            r_side.append(+ 1./self.r_alpha_gp[q] * g_par.dict_case['tu']
                                          + 1./self.r_alpha_cat[q] * self.cell_list[q].cathode.t_gas[w]
                                          + self.cell_list[q].cathode.gamma[w] * g_par.dict_uni['h_vap'])
                            r_side.append(-1. / self.r_alpha_col[q] * self.t[q, w]
                                          - 1./self.r_alpha_pp[q] * g_par.dict_case['tu'])
                        else:
                            r_side.append(+ 1./self.r_alpha_ano[q] * self.cell_list[q].anode.t_gas[w]
                                          + 1./self.r_alpha_gp[q] * g_par.dict_case['t_u'])
                            r_side.append(self.cell_list[q].omega_a[w] * self.I[q, w]**2 * 0.5
                                          + 1./self.r_alpha_gm[q] * g_par.dict_case['tu'])
                            r_side.append((g_par.dict_case['vtn'] - self.cell_list[q].v_th[w] - self.cell_list[q].omega_a[w] *
                                           self.I[q, w] * 0.5) * self.I[q, w]
                                          + 1./self.r_alpha_gm[q] * g_par.dict_case['tu'])
                            r_side.append(+ 1. / self.r_alpha_gp[q] * g_par.dict_case['tu']
                                          + 1. / self.r_alpha_cat[q] * self.cell_list[q].cathode.t_gas[w]
                                          + self.cell_list[q].cathode.gamma[w] * g_par.dict_uni['h_vap'])
                            r_side.append(-1. / self.r_alpha_col[q] * self.t[q, w]
                                          - 1./self.r_alpha_pp[q] * g_par.dict_case['tu'])
        else:
            for q, item in enumerate(self.cell_list):
                for w in range(g_par.dict_case['nodes']):
                    #print(w)
                    if w is 0 or w is g_par.dict_case['nodes']-1:
                        r_side.append(self.heat_pow * 0.5
                                      + .5 / self.r_alpha_ano[q] * self.cell_list[q].anode.t_gas[w]
                                      + .5 / self.r_alpha_gp[q] * g_par.dict_case['tu']
                                      + 1. / self.r_alpha_gegp[q] * g_par.dict_case['tu']
                                      + .5 / self.r_alpha_col[q] * self.t[q+1, w])
                        #print(r_side[-1],'ano_gdl 0.5')
                        r_side.append(self.cell_list[q].omega_a[w] * self.I[q, w] ** 2 * 0.25
                                      + .5 / self.r_alpha_gm[q] * g_par.dict_case['tu']
                                      + 1. / self.r_alpha_gegm[q] * g_par.dict_case['tu'])
                        #print(r_side[-1], 'ano 0.5')
                        r_side.append((g_par.dict_case['vtn'] - self.cell_list[q].v_th[w] - self.cell_list[q].omega_a[w]
                                       * self.I[q, w] * 0.5) * self.I[q, w] * 0.5
                                      + .5 / self.r_alpha_gm[q] * g_par.dict_case['tu']
                                      + 1. / self.r_alpha_gegm[q] * g_par.dict_case['tu'])
                        #print(r_side[-1], 'cat 0.5')
                        r_side.append(+ .5 / self.r_alpha_gp[q] * g_par.dict_case['tu']
                                      + 1. / self.r_alpha_gegp[q] * g_par.dict_case['tu']
                                      + 1. / self.r_alpha_cat[q] * self.cell_list[q].cathode.t_gas[w] * 0.5
                                      + 0.5 * self.cell_list[q].cathode.gamma[w] * g_par.dict_uni['h_vap'])
                        #print(r_side[-1], 'cat_gdl 0.5')
                        r_side.append(-self.heat_pow * .5 - .5 / self.r_alpha_pp[q] * g_par.dict_case['tu']
                                      - 1. / self.r_alpha_gepp[q] * g_par.dict_case['tu']
                                      - .5 / self.r_alpha_col[q] * self.t[q, w])
                        #print(r_side[-1], 'cat_plate 0.5')
                    else:
                        r_side.append(self.heat_pow +
                                      1. / self.r_alpha_ano[q] * self.cell_list[q].anode.t_gas[w]
                                      + 1. / self.r_alpha_gp[q] * g_par.dict_case['tu']
                                      + 1. / self.r_alpha_gegp[q] * g_par.dict_case['tu']
                                      + 1. / self.r_alpha_col[q] * self.t[q+1, w])
                        #print(r_side[-1], 'ano_gdl')
                        r_side.append(self.cell_list[q].omega_a[w] * self.I[q, w] ** 2. * 0.5
                                      + 1. / self.r_alpha_gm[q] * g_par.dict_case['tu']
                                      + 1. / self.r_alpha_gegm[q] * g_par.dict_case['tu'])
                       # print(r_side[-1], 'ano')
                        r_side.append((g_par.dict_case['vtn'] - self.cell_list[q].v_th[w] - self.cell_list[q].omega_a[w]
                                       * self.I[q, w]*.5) * self.I[q, w]
                                      + 1. / self.r_alpha_gm[q] * g_par.dict_case['tu']
                                      + 1. / self.r_alpha_gegm[q] * g_par.dict_case['tu'])
                       # print(r_side[-1], 'cat')
                        r_side.append(+ 1. / self.r_alpha_gp[q] * g_par.dict_case['tu']
                                      + 1. / self.r_alpha_gegp[q] * g_par.dict_case['tu']
                                      + 1. / self.r_alpha_cat[q] * self.cell_list[q].cathode.t_gas[w] * 1.
                                      + 1. * self.cell_list[q].cathode.gamma[w] * g_par.dict_uni['h_vap'])
                       # print(r_side[-1], 'cat_gdl')
                        r_side.append(-self.heat_pow * 1. - 1. / self.r_alpha_pp[q] * g_par.dict_case['tu']
                                      - 1. / self.r_alpha_gepp[q] * g_par.dict_case['tu']
                                      - 1. / self.r_alpha_col[q] * self.t[q, w])
                       # print(r_side[-1], 'cat_plate')


        #t_vec = np.linalg.tensorsolve(np.asarray(temp_mat), r_side)
        #t_vec = np.linalg.lstsq(np.asarray(temp_mat), r_side, rcond=None)
        t_vec = np.linalg.solve(np.asarray(temp_mat), r_side)
        counter = 0
        for q, item in enumerate(self.cell_list):
            for w in range(g_par.dict_case['nodes']):
                self.cell_list[q].t5[w] = t_vec[counter]
                self.cell_list[q].t4[w] = t_vec[counter + 1]
                self.cell_list[q].t1[w] = t_vec[counter + 2]
                self.cell_list[q].t2[w] = t_vec[counter + 3]
                self.cell_list[q].t3[w] = t_vec[counter + 4]
                counter = counter + 5