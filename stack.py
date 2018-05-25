import numpy as np
import copy
import matrix_database as m_d


class Stack:

    def __init__(self, cell, cell_numb, end_plate_t, alpha, resistance, adiabat):
        self.cell = cell
        self.cell_numb = cell_numb
        self.t_w = end_plate_t
        self.alpha = alpha
        self.alpha_1 = 2. * self.alpha * self.cell.extend \
                       / (self.cell.thickness_gde + self.cell.thickness_plate)
        self.alpha_2 = 2. * self.alpha * self.cell.extend \
                       / (self.cell.mem_thic + self.cell.thickness_gde)
        self.alpha_3 = self.alpha * self.cell.extend / self.cell.thickness_plate
        # print(self.alpha_1, self.alpha_2, self.alpha_3)
        self.cell_list = []
        for i in range(self.cell_numb):
            self.cell_list.append(copy.deepcopy(cell))
        self.resistance = resistance
        self.adiabat = adiabat
        self.i = np.full((cell_numb, cell.cathode.channel.nodes), self.cell.cathode.tar_cd)
        self.h_vap = 45400.
        if self.cell.coolant_channel is True:
            self.t = np.full((cell_numb, cell.cathode.channel.nodes), self.cell.T_cool_in)
            self.g = 529.  # coolant flow
            self.a = 4000.  # heat_coef
        else:
            self.t = 0.
            self.g = 0.
            self.a = 0.
        self.temp_mat_inv = m_d.temperature_matrix_conv(self.cell.cathode.channel.nodes,
                                                        self.cell_numb, self.cell.mu_p,
                                                        self.cell.mu_g, self.cell.mu_m,
                                                        self.a, self.alpha_1,
                                                        self.alpha_2, self.alpha_3)

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
                        self.cell_list[q].t[w] = self.cell.cathode.channel.t_in\
                                                 + w / self.cell.cathode.channel.nodes
            else:
                for w in range(self.cell.cathode.channel.nodes):
                    if w is 0:
                        self.cell_list[q].t[w] = self.cell_list[q].t_cool_in
                    else:
                        var1 = self.cell.mu_p / self.g * self.cell.cathode.channel.d_x
                        var2 = self.cell_list[q].t2[w] + self.cell_list[q - 1].t5[w]\
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
                res[0] = self.h_vap * self.cell_list[q].cathode.gamma[w]\
                         - self.alpha_1 * self.cell.t_u
                var1 = self.cell.vtn - self.cell_list[q].v[w]
                var2 = self.cell_list[q].omega[w] * self.cell_list[q].i[w] / 2.
                res[1] = (var1 - var2) * self.cell_list[q].i[w]\
                         - self.alpha_2 * self.cell.t_u
                res[2] = 0.5 * self.cell_list[q].omega[w] * self.i[q][w] ** 2.\
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
                else:
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
