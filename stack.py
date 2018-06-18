import numpy as np
import copy
import matrix_database as m_d
import saturation_pressure_vapour as p_sat
import input_parameter as i_p
import global_functions as g_func
import matplotlib.pyplot as plt


class Stack:

    def __init__(self, cell, cell_numb, end_plate_t, alpha, resistance,
                 adiabat, heigth, width, kf, stoi_cat, stoi_ano):
        self.cell = cell
        self.cell_numb = cell_numb
        self.t_w = end_plate_t
        self.alpha = alpha
        self.resistance = resistance
        self.adiabat = adiabat
        self.heigth = heigth
        self.width = width
        self.kf = kf
        self.stoi_cat = stoi_cat
        self.stoi_ano = stoi_ano
        self.init_param()
        self.init_func()
        self.init_arrays()

    def init_func(self):
        for i in range(self.cell_numb):
            self.cell_list.append(copy.deepcopy(self.cell))
        if self.cell.coolant_channel is True:
            self.t = np.full((self.cell_numb,
                             self.cell.cathode.channel.nodes),
                             self.cell.T_cool_in)
            self.g = 529.  # coolant flow
            self.a = 4000.  # heat_coef
        else:
            self.t = 0.
            self.g = 0.
            self.a = 0.
        self.extent = 2. * (self.width + self.heigth)
        self.cross_area = self.width * self.heigth
        self.h_d = 4. * self.cross_area / self.extent
        self.cell_heigth = self.cell.thickness_mea \
                           + 2 * self.cell.thickness_plate \
                           + 2 * self.cell.thickness_gde
        self.q_h_in_cat = np.full(self.cell_numb, 0.)
        self.q_h_in_ano = np.full(self.cell_numb, 0.)
        self.q_h_in_cat[-1] = self.stoi_cat \
                             * self.cell_numb \
                             * self.cell.cathode.tar_cd \
                             * self.cell.cathode.channel.width \
                             * self.cell.cathode.channel.length \
                             / (4. * self.cell.cathode.f_const) \
                             * (1. + 0.79 / 0.21) \
                             * (1. + self.cell.cathode.channel.phi
                                * p_sat.water.calc_psat(self.cell.cathode.channel.t_in)
                                / (self.cell.cathode.channel.p_in
                                   - self.cell.cathode.channel.phi
                                   * p_sat.water.calc_psat(self.cell.cathode.channel.t_in)))
        self.q_h_in_ano[-1] = self.stoi_ano \
                             * self.cell_numb \
                             * self.cell.cathode.tar_cd \
                             * self.cell.cathode.channel.width \
                             * self.cell.cathode.channel.length \
                             / (2. * self.cell.cathode.f_const) \
                             * (1. + self.cell.anode.channel.phi
                                * p_sat.water.calc_psat(self.cell.anode.channel.t_in)
                                / (self.cell.anode.channel.p_in
                                   - self.cell.anode.channel.phi
                                   * p_sat.water.calc_psat(self.cell.anode.channel.t_in)))
        self.set_stoi(np.full(self.cell_numb, self.stoi_cat),
                      np.full(self.cell_numb, self.stoi_ano))
        # print(self.q_h_in_cat[-1])

    def init_param(self):
        self.alpha_1 = 2. * self.alpha * self.cell.extend \
                       / (self.cell.thickness_gde + self.cell.thickness_plate)
        self.alpha_2 = 2. * self.alpha * self.cell.extend \
                       / (self.cell.thickness_mea + self.cell.thickness_gde)
        self.alpha_3 = self.alpha * self.cell.extend / self.cell.thickness_plate
        self.h_vap = 45400.
        self.cell_list = []

    def init_arrays(self):
        x = np.linspace(self.cell.cathode.tar_cd * 1.5,
                        self.cell.cathode.tar_cd * 0.5,
                        self.cell.cathode.channel.nodes)
        y = []
        for q in range(self.cell_numb): y.append(x)
        self.i = np.array(y)
        self.temp_mat_inv = m_d.temperature_matrix_conv(self.cell.cathode.channel.nodes,
                                                        self.cell_numb, self.cell.mu_p,
                                                        self.cell.mu_g, self.cell.mu_m,
                                                        self.a, self.alpha_1,
                                                        self.alpha_2, self.alpha_3)
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

    def update(self):
        for j in range(self.cell_numb):
            self.cell_list[j].set_i(self.i[j, :])
            self.cell_list[j].update()
        if self.cell.coolant_channel is True:
            self.calc_coolant_t_fit()
            # self.calc_coolant_t()
        self.calc_layer_t()
        self.stack_v()
        self.stack_dv()
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

    def set_i(self, i):
        self.i = i

    def stack_v(self):
        var = []
        for i, item in enumerate(self.cell_list):
            var = np.hstack((var, self.cell_list[i].v))
        self.v = var
        # running

    def stack_dv(self):
        var = []
        for i, item in enumerate(self.cell_list):
            var = np.hstack((var, self.cell_list[i].dv))
        self.dv = var
        # running

    def calc_coolant_t(self):
        for q, item in enumerate(self.cell_list):
            if q is 0:
                for w in range(self.cell.cathode.channel.nodes):
                    if w is 0:
                        self.cell_list[q].t[w] = self.cell_list[q].t_cool_in
                    else:
                        if self.endplate_fix_heatflow is True:
                            var2 = self.cell_list[q].t2[
                                       w] - 1 * self.endplate_heat_flow.density
                        else:
                            var2 = self.cell_list[q].t2[w] - 1 * self.cell_list[q].t3[w]
                        var1 = self.cell.mu_p / self.g * self.cell.cathode.channel.d_x
                        self.cell_list[q].t[w] = var1 * var2 + self.cell_list[q].t[w - 1]
                        self.cell_list[q].t[w] = self.cell.cathode.channel.t_in \
                                                 + w / self.cell.cathode.channel.nodes
            else:
                for w in range(self.cell.cathode.channel.nodes):
                    if w is 0:
                        self.cell_list[q].t[w] = self.cell_list[q].t_cool_in
                    else:
                        var1 = self.cell.mu_p / self.g * self.cell.cathode.channel.d_x
                        var2 = self.cell_list[q].t2[w] + self.cell_list[q - 1].t5[w] \
                               - 2 * self.cell_list[q].t3[w]
                        self.cell_list[q].t[w] = var1 * var2 + self.cell_list[q].t[w - 1]

    def calc_coolant_t_fit(self):
        for q, item in enumerate(self.cell_list):
            self.cell_list[q].t = np.linspace(self.cell.t_cool_in,
                                              self.cell.cathode.channel.t_in,
                                              self.cell.cathode.channel.nodes)

    def calc_layer_t(self):
        res_vec = []
        res = np.full(5, 0.)
        for q, item in enumerate(self.cell_list):
            for w in range(self.cell.cathode.channel.nodes):
                res[0] = self.h_vap * self.cell_list[q].cathode.gamma[w] \
                         - self.alpha_1 * self.cell.t_u
                var1 = self.cell.vtn - self.cell_list[q].v[w]
                var2 = self.cell_list[q].omega[w] * self.cell_list[q].i[w] / 2.
                res[1] = (var1 - var2) * self.cell_list[q].i[w] \
                         - self.alpha_2 * self.cell.t_u
                res[2] = 0.5 * self.cell_list[q].omega[w] * self.i[q][w] ** 2. \
                         - self.alpha_2 * self.cell.t_u

                if self.adiabat is True:
                    if q is 0:
                        res[3] = self.h_vap * self.cell_list[q].anode.gamma[w] \
                                 - self.alpha_1 * self.cell.t_u \
                                 + self.cell.mu_p * self.cell_list[q + 1].t3[w]
                        res[4] = -self.cell.mu_p * self.cell_list[q].t3[w] \
                                 - self.alpha_3 * self.cell.t_u \
                                 - self.a * self.cell_list[q].t[w]
                    elif q is self.cell_numb - 1:
                        res[3] = self.h_vap * self.cell_list[q].anode.gamma[w] \
                                 - self.alpha_1 * self.cell.t_u \
                                 + self.cell.mu_p * self.cell_list[q].t5[w]
                        res[4] = - self.cell.mu_p * self.cell_list[q - 1].t5[w] \
                                 - self.alpha_3 * self.cell.t_u \
                                 - self.cell_list[q].t[w] * self.a
                    else:
                        res[3] = self.h_vap * self.cell_list[q].anode.gamma[w] \
                                 - self.alpha_1 * self.cell.t_u \
                                 + self.cell.mu_p * self.cell_list[q + 1].t3[w]
                        res[4] = - self.cell.mu_p * self.cell_list[q - 1].t5[w] \
                                 - self.alpha_3 * self.cell.t_u \
                                 - self.cell_list[q].t[w] * self.a
                else:  # fixed Bctemp
                    if q is 0:
                        res[3] = self.h_vap * self.cell_list[q].anode.gamma[w] \
                                 - self.alpha_1 * self.cell.t_u \
                                 + self.cell.mu_p * self.cell_list[q + 1].t3[w]
                        res[4] = -self.cell.mu_p * self.t_w \
                                 - self.alpha_3 * self.cell.t_u \
                                 - self.a * self.cell_list[q].t[w]
                    elif q is self.cell_numb - 1:
                        res[3] = self.h_vap * self.cell_list[q].anode.gamma[w] \
                                 - self.alpha_1 * self.cell.t_u \
                                 + self.cell.mu_p * self.t_w
                        res[4] = - self.cell.mu_p * self.cell_list[q - 1].t5[w] \
                                 - self.alpha_3 * self.cell.t_u \
                                 - self.cell_list[q].t[w] * self.a
                    else:
                        res[3] = self.h_vap * self.cell_list[q].anode.gamma[w] \
                                 - self.alpha_1 * self.cell.t_u \
                                 + self.cell.mu_p * self.cell_list[q + 1].t3[w]
                        res[4] = - self.cell.mu_p * self.cell_list[q - 1].t5[w] \
                                 - self.alpha_3 * self.cell.t_u \
                                 - self.cell_list[q].t[w] * self.a
                res_vec = np.hstack((res_vec, res))

        t_vec = np.matmul(self.temp_mat_inv, res_vec)
        counter = 0
        for q, item in enumerate(self.cell_list):
            for w in range(self.cell.cathode.channel.nodes):
                self.cell_list[q].t1[w] = t_vec[counter]
                self.cell_list[q].t2[w] = t_vec[counter + 1]
                self.cell_list[q].t3[w] = t_vec[counter + 2]
                self.cell_list[q].t4[w] = t_vec[counter + 3]
                self.cell_list[q].t5[w] = t_vec[counter + 4]
            counter = counter + 5

    def calc_layer_t_off(self):
        for q, item in enumerate(self.cell_list):
            self.cell_list[q].t1 = np.full(self.cell.cathode.channel.nodes, 300.)
            self.cell_list[q].t2 = np.full(self.cell.cathode.channel.nodes, 300.)
            self.cell_list[q].t3 = np.full(self.cell.cathode.channel.nodes, 300.)
            self.cell_list[q].t4 = np.full(self.cell.cathode.channel.nodes, 300.)
            self.cell_list[q].t5 = np.full(self.cell.cathode.channel.nodes, 300.)

    def set_stoi(self, stoi_cat, stoi_ano):
        for q in range(self.cell_numb):
            self.cell_list[q].cathode.stoi = stoi_cat[q]
            self.cell_list[q].anode.stoi = stoi_ano[q]

    def calc_stoi(self, q, val):
        if val is 4:
            b = 1 + 0.79 / 0.21
            c = 1 + i_p.phi_cat \
                * p_sat.water.calc_psat(i_p.t_cat_in) \
                / (i_p.p_cat_in - i_p.phi_cat * p_sat.water.calc_psat(i_p.t_cat_in))
        else:
            c = 1 + i_p.phi_ano \
                * p_sat.water.calc_psat(i_p.t_ano_in) \
                / (i_p.p_ano_in - i_p.phi_ano * p_sat.water.calc_psat(i_p.t_ano_in))
        a = i_p.tar_cd * i_p.channel_width * i_p.channel_heigth\
            / (val * self.cell.cathode.f_const)
        return q/(a * b * c)

    def set_p(self):
        for q in range (self.cell_numb):
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
        # print(self.q_h_in_cat)
        # print(self.q_h_out_cat)
        # print(self.q_h_in_ano)
        # print(self.q_h_out_ano)

    def calc_header_temp(self):
        self.t_h_out_cat[0] = self.cell_list[0].cathode.t2[-1]
        self.t_h_out_ano[0] = self.cell_list[0].anode.t2[0]
        for q in range(1, self.cell_numb):
            self.t_h_out_cat[q] = (self.t_h_out_cat[q - 1] * self.q_h_out_cat[q - 1]
                                   + self.cell_list[q].cathode.q_sum[-1]
                                   * self.cell_list[q].cathode.t2[-1]) \
                                  / self.q_h_out_cat[q]
            self.t_h_out_ano[q] = (self.t_h_out_ano[q - 1] * self.q_h_out_ano[q - 1]
                                   + self.cell_list[q].anode.q_sum[0]
                                   * self.cell_list[q].anode.t2[0]) \
                                  / self.q_h_out_ano[q]
        # print(self.t_h_out_cat)
        # print(self.t_h_out_ano)

    def calc_header_velocity(self):
        self.v_h_in_cat = self.q_h_in_cat \
                          * self.cell.cathode.r \
                          * self.cell.cathode.channel.t_in \
                          / (self.p_h_in_cat * self.cross_area)
        self.v_h_in_ano = self.q_h_in_ano \
                          * self.cell.cathode.r \
                          * self.cell.anode.channel.t_in \
                          / (self.p_h_in_ano * self.cross_area)
        self.v_h_out_cat = self.q_h_out_cat \
                           * self.cell.cathode.r \
                           * self.t_h_out_cat \
                           / (self.p_h_out_cat * self.cross_area)
        self.v_h_out_ano = self.q_h_out_ano \
                           * self.cell.anode.r \
                           * self.t_h_out_ano \
                           / (self.p_h_out_ano * self.cross_area)
        # print(self.v_h_in_cat)
        # print(self.v_h_out_cat)
        # print(self.v_h_in_ano)
        # print(self.v_h_out_ano)

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
        # print(self.r_mix_h_in_cat)
        # print(self.r_mix_h_out_cat)
        # print(self.r_mix_h_in_ano)
        # print(self.r_mix_h_out_ano)

    def calc_header_roh(self):
        self.roh_h_in_cat = g_func.calc_roh(self.p_h_in_cat,
                                            self.r_mix_h_in_cat,
                                            self.cell.cathode.channel.t_in)
        self.roh_h_in_ano = g_func.calc_roh(self.p_h_in_ano,
                                            self.r_mix_h_in_ano,
                                            self.cell.anode.channel.t_in)
        self.roh_h_out_cat = g_func.calc_roh(self.p_h_out_cat,
                                             self.r_mix_h_out_cat,
                                             self.t_h_out_cat)
        self.roh_h_out_ano = g_func.calc_roh(self.p_h_out_ano,
                                             self.r_mix_h_out_ano,
                                             self.t_h_out_ano)
        # print(self.roh_h_in_cat)
        # print(self.roh_h_out_cat)
        # print(self.roh_h_in_ano)
        # print(self.roh_h_out_ano)

    def calc_header_reynolds_numb(self):
        self.re_h_in_cat = g_func.calc_re(self.roh_h_in_cat,
                                          self.v_h_in_cat,
                                          self.h_d,
                                          self.cell.cathode.visc)
        self.re_h_in_ano = g_func.calc_re(self.roh_h_in_ano,
                                          self.v_h_in_ano,
                                          self.h_d,
                                          self.cell.anode.visc)
        self.re_h_out_cat = g_func.calc_re(self.roh_h_out_cat,
                                           self.v_h_out_cat,
                                           self.h_d,
                                           self.cell.cathode.visc)
        self.re_h_out_ano = g_func.calc_re(self.roh_h_out_ano,
                                           self.v_h_out_ano,
                                           self.h_d,
                                           self.cell.anode.visc)
        # print(self.re_h_in_cat)
        # print(self.re_h_out_cat)
        # print(self.re_h_in_ano)
        # print(self.re_h_out_ano)

    def calc_header_fanning_friction_factor(self):
        self.f_h_in_cat = g_func.calc_fanning_friction_factor(self.re_h_in_cat)
        self.f_h_in_ano = g_func.calc_fanning_friction_factor(self.re_h_in_ano)
        self.f_h_out_cat = g_func.calc_fanning_friction_factor(self.re_h_out_cat)
        self.f_h_out_ano = g_func.calc_fanning_friction_factor(self.re_h_out_ano)
        # print(self.f_h_in_cat)
        # print(self.f_h_out_cat)
        # print(self.f_h_in_ano)
        # print(self.f_h_out_ano)

    def calc_header_p_out(self):  # change for parallel flow
        self.p_h_out_cat[self.cell_numb - 1] = self.cell.cathode.channel.p_in
        self.p_h_out_ano[self.cell_numb - 1] = self.cell.anode.channel.p_in
        for q in range(self.cell_numb - 1, 0, -1):
            self.p_h_out_cat[q - 1] = self.p_h_out_cat[q] \
                                      + g_func.calc_header_pressure_drop(self.roh_h_out_cat[q-1],
                                                                         self.v_h_out_cat[q],
                                                                         self.v_h_out_cat[q-1],
                                                                         self.f_h_out_cat[q-1],
                                                                         self.kf,
                                                                         self.cell_heigth,
                                                                         self.h_d)
            self.p_h_out_ano[q - 1] = self.p_h_out_ano[q] \
                                      + g_func.calc_header_pressure_drop(self.roh_h_out_ano[q-1],
                                                                         self.v_h_out_ano[q],
                                                                         self.v_h_out_ano[q-1],
                                                                         self.f_h_out_ano[q-1],
                                                                         self.kf,
                                                                         self.cell_heigth,
                                                                         self.h_d)
        # print(self.p_h_out_cat)
        # print(self.p_h_out_ano)
        # plt.plot(self.p_h_out_cat)
        # plt.plot(self.p_h_out_ano)
        # plt.show()

    def calc_ref_perm(self):
        self.perm_ano = self.cell_list[0].anode.visc\
                        / self.cell_list[0].anode.channel.cross_area\
                        * self.cell_list[0].anode.channel.length\
                        * np.median(self.cell_list[0].anode.q_sum)\
                        / self.ref_p_drop_ano
        self.perm_cat = self.cell_list[0].cathode.visc\
                        / self.cell_list[0].cathode.channel.cross_area\
                        * self.cell_list[0].cathode.channel.length\
                        * np.median(self.cell_list[0].cathode.q_sum)\
                        / self.ref_p_drop_cat
        self.ref_p_new_drop_cat = self.cell_list[0].cathode.q_sum[0]\
                         / self.cell.cathode.channel.cross_area\
                         * self.cell.cathode.visc\
                         * self.cell.cathode.channel.length\
                         / self.perm_cat
        self.p_adjusting_factor_cat = self.ref_p_drop_cat / self.ref_p_new_drop_cat
        self.ref_p_new_drop_ano = self.cell_list[0].anode.q_sum[-1]\
                         / self.cell.anode.channel.cross_area\
                         * self.cell.anode.visc\
                         * self.cell.anode.channel.length\
                         / self.perm_ano
        self.p_adjusting_factor_ano = self.ref_p_drop_ano / self.ref_p_new_drop_ano
        # print(self.p_adjusting_factor_cat)
        # print(self.ref_p_new_drop_cat)
        # print(self.p_adjusting_factor_ano)
        # print(self.ref_p_new_drop_ano)
        # print(self.ref_p_new_drop_ano * self.p_adjusting_factor_ano)
        # print(self.perm_cat, self.perm_ano)

    def calc_ref_p_drop(self):
        self.ref_p_drop_cat = self.cell_list[0].cathode.p[0]\
                              - self.cell_list[0].cathode.p[-1]
        self.ref_p_drop_ano = self.cell_list[0].anode.p[-1]\
                              - self.cell_list[0].anode.p[0]

    def calc_new_ref_p_drop(self):
        self.ref_p_drop_cat = self.q_h_in_cat[self.cell_numb-1]\
                              / np.sum(self.alpha_cat) * self.cell.cathode.visc\
                              * self.cell.cathode.channel.length\
                              / self.cell.cathode.channel.cross_area\
                              / self.perm_cat\
                              * self.p_adjusting_factor_cat
        self.ref_p_drop_ano = self.q_h_in_ano[self.cell_numb-1]\
                              / np.sum(self.alpha_ano) * self.cell.anode.visc\
                              * self.cell.anode.channel.length\
                              / self.cell.anode.channel.cross_area\
                              / self.perm_ano\
                              * self.p_adjusting_factor_ano
        # print(self.ref_p_drop_cat)
        # print(self.ref_p_drop_ano)

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
                                                                    self.cell_heigth,
                                                                    self.h_d)
            self.p_h_in_ano[q] = self.p_h_in_ano[q - 1] \
                                 + g_func.calc_header_pressure_drop(self.roh_h_in_ano[q],
                                                                    self.v_h_in_ano[q - 1],
                                                                    self.v_h_in_ano[q],
                                                                    self.f_h_in_ano[q],
                                                                    self.kf,
                                                                    self.cell_heigth,
                                                                    self.h_d)
        # print(self.p_h_in_cat)
        # print(self.p_h_out_cat)
        # print(self.p_h_in_ano)
        # plt.plot(self.p_h_in_cat)
        # plt.plot(self.p_h_out_cat)
        # plt.show()

    def calc_flow_distribution_factor(self):
        self.alpha_cat = (self.p_h_in_cat - self.p_h_out_cat)\
                         / self.ref_p_drop_cat
        self.alpha_ano = (self.p_h_in_ano - self.p_h_out_ano)\
                         / self.ref_p_drop_ano
        # plt.plot(self.alpha_cat)
        # plt.plot(self.alpha_ano)
        # plt.xlim(self.cell_numb, 0)
        # plt.show()
        # print(self.alpha_cat)
        # print(self.alpha_ano)

    def calc_new_cell_flows(self):
        self.q_x_cat = (self.p_h_in_cat - self.p_h_out_cat)\
                        * self.perm_cat\
                        * self.cell.cathode.channel.cross_area\
                        / self.cell.cathode.visc\
                        / self.cell.cathode.channel.length\
                        / self.p_adjusting_factor_cat
        self.q_x_ano = (self.p_h_in_ano - self.p_h_out_ano) \
                       * self.perm_ano \
                       * self.cell.anode.channel.cross_area \
                       / self.cell.anode.visc \
                       / self.cell.anode.channel.length \
                       / self.p_adjusting_factor_ano
            # print(self.q_x_ano[q])
            # print(self.cell_list[q].anode.q_sum[0])

    def calc_new_cell_stoi(self):
        self.new_stoi_cat = self.q_x_cat /\
                            (self.q_h_in_cat[-1] / self.cell_numb)\
                            * self.stoi_cat
        self.new_stoi_ano = self.q_x_ano\
                            / (self.q_h_in_ano[-1] / self.cell_numb)\
                            * self.stoi_ano
        # print(np.median(self.new_stoi_cat))
        # print(np.median(self.new_stoi_ano))
