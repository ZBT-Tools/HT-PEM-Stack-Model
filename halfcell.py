import matrix_database
import saturation_pressure_vapour
import numpy as np


class Halfcell:
    def __init__(self, channel, stoi, spec_numb, val, tar_cd):
        self.channel = channel
        self.stoi = stoi
        self.spec_numb = spec_numb
        self.val = val
        self.tar_cd = tar_cd
        self.f_const = 96485.33289
        self.r = 8.314459848
        self.pem_type = True
        self.j = []
        self.w = np.full(self.channel.nodes, 0.)
        self.p = np.full(self.channel.nodes, self.channel.p_in)
        self.gamma = np.full(self.channel.nodes, 0.)
        self.humidity = np.full(self.channel.nodes, 0.)
        self.free_water = np.full(self.channel.nodes, 0.)
        self.i = np.full(self.channel.nodes, self.tar_cd)
        self.node_forward = \
            matrix_database.forward_matrix_reac_flow_new(self.channel.nodes)
        self.node_backward = \
            matrix_database.backward_matrix_reac_flow_new(self.channel.nodes)
        for l in range(self.spec_numb):  # Init flow and con arrays (1-n)
            exec("self.q%d = np.full(self.channel.nodes,0.)" % (l + 1))
            # 1 reactant(hydrogen/oxygen), 2 water, 3 (if oxygen) nitrogen
            exec("self.c%d = np.full(self.channel.nodes,0.)" % (l + 1))
            # 1 reactant(hydrogen/oxygen), 2 water, 3 (if oxygen) nitrogen
            exec('self.t%d = np.full(self.channel.nodes,self.channel.t_in)'
                 % (l + 1))
            # 1 catalyst layer, 2 channel layer, 3 coolant plate layer

    def update(self):
        if self.pem_type is True:  # HT
            self.calc_reac_flow()
            self.calc_water_flow()  # stabel with j=0?
            self.calc_pressure()
            self.calc_con()
            self.calc_fluid_water()
            self.calc_cond_rates()
        else:  # NT
            self.calc_reac_flow()
            self.calc_water_flow()
            self.calc_pressure()
            self.calc_con()
            self.calc_fluid_water()
            self.calc_cond_rates()
            self.calc_rel_humidity()
            self.calc_free_water()

    def set_i(self, i):
        self.i = i

    def set_j(self, j):
        self.j = j

    def set_pem_type(self, pem_type):
        self.pem_type = pem_type

    def set_t(self, var):  # Ã¼bergabe der Layertemperaturen als Array(T1-T2/3)
        for l, item in enumerate(var):
            exec('self.t%d = var[l]' % (l + 1))

    def calc_reac_flow(self):
        if self.channel.flow_dir is True:
            var1 = self.stoi * self.tar_cd * self.channel.length / \
                   (self.val * self.f_const)
            self.q1 = var1 + np.matmul(self.node_forward, self.i) \
                      * self.channel.d_x / (self.val * self.f_const)

        else:
            var1 = self.stoi * self.tar_cd * self.channel.length \
                   / (self.val * self.f_const)
            self.q1 = var1 + np.matmul(self.node_backward, self.i) \
                      * self.channel.d_x / (self.val * self.f_const)

    def calc_water_flow(self):  # water molar flux in the channel#elementwise
        if self.pem_type is True:
            if self.channel.flow_dir is True:
                a = self.channel.d_x / (self.val * self.f_const * 0.5) \
                    * np.matmul(-self.node_forward, self.i)
                # production
                b = saturation_pressure_vapour.water.calc_psat(self.channel.t_in) \
                    * self.channel.phi  # richtig?
                c = (0.79 / 0.21 + 1.) * self.q1[0] * b / (self.channel.p_in - b)
                # Value at q2[0]
                self.q2 = a + c
                self.q3 = np.full(self.channel.nodes, self.q1[0] * 0.79 / 0.21)

            elif self.channel.flow_dir is False:
                a = saturation_pressure_vapour.water.calc_psat(self.channel.t_in) \
                    * self.channel.phi
                b = a / (self.channel.p_in - a)

                self.q2 = np.full(self.channel.nodes, self.q1[-1] * b
                                  / (self.channel.p_in - b))
        else:
            if self.channel.flow_dir is True:
                a = self.channel.d_x / (self.val * self.f_const * 0.5) \
                    * np.matmul(-self.node_forward, self.i)
                # production
                b = saturation_pressure_vapour.water.calc_psat(self.channel.t_in) \
                    * self.channel.phi
                c = (0.79 / 0.21 + 1.) * self.q1[0] * b \
                    / (self.channel.p_in - b)  # Value at q2[0]
                d = self.channel.d_x * \
                    np.matmul(-self.node_forward, self.j), self.channel.nodes
                # crossover
                self.q2 = a + c + d
                self.q3 = np.full(self.channel.nodes, self.q1[0] * 0.79 / 0.21)

            elif self.channel.flow_dir is False:
                a = saturation_pressure_vapour.water.calc_psat(self.channel.t_in) \
                    * self.channel.phi
                b = a / (self.channel.p_in - a)
                c = self.channel.d_x * np.matmul(-self.node_backward, self.j)
                self.q2 = c + b

    def calc_pressure(self):  # channel pressure #qo qc qn / qh qa#elementwise
        a = self.channel.p_in
        b = self.channel.d_x * self.channel.k
        if self.spec_numb is 3:
            c = np.matmul(-self.node_forward, self.q1 + self.q2 - self.w)\
                + self.q3
            d = self.r * self.t1 / self.p
            self.p = a - b * c * d

        else:
            c = self.q1 + self.q2 - self.w
            d = self.r * self.t1 / self.p
            self.p = a - b * np.matmul(-self.node_backward, c) * d

    def calc_con(self):  # Gas cocentrations in the channel [moles/mÂ³]#nodewise
        if self.spec_numb is 3:  # cathode
            for w in range(self.channel.nodes):
                var1 = self.p[w] / (self.r * self.t2[w])
                var2 = self.q2[w] / (self.q1[w] + self.q2[w] + self.q3[w])
                self.c2[w] = var1 * var2
                if self.c2[w] >= saturation_pressure_vapour.water.calc_psat(self.t2[w]) \
                        / (self.r * self.t2[w]):  # saturated
                    a = (self.p[w]
                         - saturation_pressure_vapour.water.calc_psat(self.t2[w]))\
                        / (self.r * self.t2[w])
                    b = self.q1[w] / (self.q1[w] + self.q3[w])
                    self.c1[w] = a * b
                    self.c2[w] = saturation_pressure_vapour.water.calc_psat(self.t2[w]) \
                                 / (self.r * self.t2[w])
                else:  # not saturated
                    var3 = self.p[w] / (self.r * self.t2[w])
                    var4 = self.q1[w] / (self.q1[w] + self.q2[w] + self.q3[w])
                    self.c1[w] = var3 * var4
        else:
            for w in range(self.channel.nodes):
                var1 = self.p[w] / (self.r * self.t2[w])
                var2 = self.q2[w] / (self.q1[w] + self.q2[w])
                self.c2[w] = var1 * var2
                if self.c2[w] >= saturation_pressure_vapour.water.calc_psat(self.t2[w]) \
                        / (self.r * self.t2[w]):  # saturated
                    a = (self.p[w]
                         - saturation_pressure_vapour.water.calc_psat(self.t2[w]))\
                        / (self.r * self.t2[w])
                    self.c1[w] = a
                    self.c2[w] = saturation_pressure_vapour.water.calc_psat(self.t2[w]) \
                                 / (self.r * self.t2[w])
                else:  # not saturated
                    var3 = self.p[w] / (self.r * self.t2[w])
                    var4 = self.q1[w] / (self.q1[w] + self.q2[w])
                    self.c1[w] = var3 * var4

    def calc_fluid_water(self):  # fluid water in the channel#nodewise
        self.w = self.q2 - self.c2 / self.c1 * self.q1

    def calc_cond_rates(self):  # condensation rates of the vapour in the channel
        if self.spec_numb is 3:
            self.gamma = np.gradient(self.w, self.channel.d_x)
        else:
            self.gamma = -np.gradient(self.w, self.channel.d_x)

    def calc_rel_humidity(self):  # relative humidity in the channel
        self.humidity = self.c2 * self.r * self.t2 \
                        / saturation_pressure_vapour.water.calc_psat(self.t2)
        # ok

    def calc_free_water(self):  # membrane free water content#nodewise
        a = 0.043 + 17.81 * self.humidity
        b = -39.85 * self.humidity ** 2 + 36. * self.humidity ** 3
        self.free_water = a + b
