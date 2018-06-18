import matrix_database
import global_functions as g_func
import saturation_pressure_vapour as p_sat
import numpy as np
import matplotlib.pyplot as plt
# from numba import jit


class Halfcell:
    def __init__(self, channel, spec_numb, val, tar_cd):
        self.channel = channel
        self.spec_numb = spec_numb
        self.val = val
        self.tar_cd = tar_cd
        self.init_param()
        self.init_arrays()
        self.init_func()

    def init_param(self):
        self.f_const = 96485.33289
        self.r = 8.314459848
        self.m_h2 = 2.
        self.m_o2 = 32.
        self.m_h2o = 18.
        self.m_n2 = 28.
        self.pem_type = True
        self.j = []
        self.re = []
        self.roh = []
        self.stoi = 1.5

    def init_arrays(self):
        self.w = np.full(self.channel.nodes, 0.)
        self.p = np.full(self.channel.nodes, self.channel.p_in)
        self.gamma = np.full(self.channel.nodes, 0.)
        self.humidity = np.full(self.channel.nodes, 0.)
        self.free_water = np.full(self.channel.nodes, 0.)
        self.i = np.full(self.channel.nodes, self.tar_cd)
        self.u = np.full(self.channel.nodes,    0.)
        self.node_forward = \
            matrix_database.forward_matrix_reac_flow_new(self.channel.nodes)
        self.node_backward = \
            matrix_database.backward_matrix_reac_flow_new(self.channel.nodes)
        for l in range(self.spec_numb):  # Init flow and con arrays (1-n)
            exec('self.t%d = np.full(self.channel.nodes,self.channel.t_in)'
                 % (l + 1))
            # 1 catalyst layer, 2 channel layer, 3 coolant plate layer
        self.perm = np.full(self.channel.nodes, 0.)
        self.gas_flow = np.full((self.spec_numb, self.channel.nodes), 0.)
        self.gas_con = np.full((self.spec_numb, self.channel.nodes), 0.)

    def init_func(self):
        if self.spec_numb is 3:
            self.visc = 17.1e-6
        else:
            self.visc = 8.4e-6
        self.dh = 4 * self.channel.cross_area/self.channel.extent
        self.r_h2 = self.r/self.m_h2*1e3
        self.r_o2 = self.r/self.m_o2*1e3
        self.r_h2o = self.r/self.m_h2o*1e3
        self.r_n2 = self.r/self.m_n2*1e3
        self.r_dry_air = 287.058
        self.r_vap = 461.523

    def update(self):
        self.calc_reac_flow()
        self.calc_water_flow()  # stabel with j=0?
        #self.calc_pressure()
        self.calc_con()
        self.calc_fluid_water()
        self.sum_flows()
        self.calc_cond_rates()
        self.calc_mass_fraction()
        self.calc_r_mix()
        self.calc_roh()
        self.calc_rel_humidity()
        self.calc_flow_velocity()
        self.calc_re()
        self.calc_pressure_scheuer_lam()
        #self.calc_channel_permeability()
        if self.pem_type is False:  # NT
            self.calc_free_water()

    def set_i(self, i):
        self.i = i

    def set_j(self, j):
        self.j = j

    def set_stoi(self, stoi):
        self.stoi = stoi

    def set_pem_type(self, pem_type):
        self.pem_type = pem_type

    def set_t(self, var):  # Ã¼bergabe der Layertemperaturen als Array(T1-T2/3)
        for l, item in enumerate(var):
            exec('self.t%d = var[l]' % (l + 1))

    def calc_reac_flow(self):
        var1 = self.stoi * self.tar_cd * self.channel.length \
               * self.channel.width / (self.val * self.f_const)
        if self.channel.flow_dir is True:
            self.gas_flow[0] = var1 + np.matmul(self.node_forward, self.i) \
                      * self.channel.d_x * self.channel.width\
                      / (self.val * self.f_const)

        else:
            self.gas_flow[0] = var1 + np.matmul(self.node_backward, self.i) \
                      * self.channel.d_x * self.channel.width\
                      / (self.val * self.f_const)

    def calc_water_flow(self):  # water molar flux in the channel#elementwise
        if self.pem_type is True:
            if self.channel.flow_dir is True:
                q_0_water = self.gas_flow[0][0] \
                            * (1. + 0.79 / 0.21) \
                            * p_sat.water.calc_psat(self.channel.t_in)\
                            * self.channel.phi\
                            / (self.channel.p_in
                               - self.channel.phi
                               * p_sat.water.calc_psat(self.channel.t_in))
                a = self.channel.d_x * self.channel.width\
                    / (self.val * self.f_const * 0.5) \
                    * np.matmul(-self.node_forward, self.i)
                # production
                self.gas_flow[1] = a + q_0_water
                self.gas_flow[2] = np.full(self.channel.nodes,
                                         self.gas_flow[0][0] * 0.79 / 0.21)

            elif self.channel.flow_dir is False:
                q_0_water = self.gas_flow[0][-1] \
                            * p_sat.water.calc_psat(self.channel.t_in)\
                            * self.channel.phi\
                            / (self.channel.p_in
                               - self.channel.phi
                               * p_sat.water.calc_psat(self.channel.t_in))
                self.gas_flow[1] = np.full(self.channel.nodes, q_0_water)
        else:
            if self.channel.flow_dir is True:
                q_0_water = self.gas_flow[0][0] \
                            * (1. + 0.79 / 0.21) \
                            * p_sat.water.calc_psat(self.channel.t_in)\
                            * self.channel.phi\
                            / (self.channel.p_in
                               - self.channel.phi
                               * p_sat.water.calc_psat(self.channel.t_in))
                a = self.channel.d_x * self.channel.width\
                    / (self.val * self.f_const * 0.5) \
                    * np.matmul(-self.node_forward, self.i)
                # production
                b = self.channel.d_x * self.channel.width\
                    * np.matmul(-self.node_forward, self.j)
                # crossover
                self.gas_flow[1] = a + b + q_0_water
                self.gas_flow[3] = np.full(self.channel.nodes, self.gas_flow[0][0] * 0.79 / 0.21)

            elif self.channel.flow_dir is False:
                q_0_water = self.gas_flow[0][-1] \
                            * p_sat.water.calc_psat(self.channel.t_in)\
                            * self.channel.phi\
                            / (self.channel.p_in
                               - self.channel.phi
                               * p_sat.water.calc_psat(self.channel.t_in))
                a = self.channel.d_x * self.channel.width\
                    * np.matmul(-self.node_backward, self.j)
                self.gas_flow[1] = a + q_0_water

    def calc_pressure(self):  # channel pressure #qo qc qn / qh qa#elementwise
        a = self.channel.p_in
        b = self.channel.d_x * self.channel.k
        if self.spec_numb is 3:
            c = np.matmul(-self.node_forward, self.gas_flow[0] + self.gas_flow[1] - self.w)\
                + self.gas_flow[2]
            d = self.r * self.t1 / (self.p * self.channel.width)
            self.p = a - b * c * d

        else:
            c = self.gas_flow[0] + self.gas_flow[1] - self.w
            d = self.r * self.t1 / (self.p * self.channel.width)
            self.p = a - b * np.matmul(-self.node_backward, c) * d

    def calc_pressure_scheuer_lam(self):
        if self.spec_numb is 3:
            mat = self.node_forward
        else:
            mat = self.node_backward
        self.p = self.channel.p_in + 32./self.dh\
                 * np.matmul(mat, self.roh * self.u**2. / self.re)\
                 * self.channel.d_x

    def calc_con(self):  # Gas cocentrations in the channel [moles/mÂ³]#nodewise
        if self.spec_numb is 3:  # cathode
            for w in range(self.channel.nodes):
                var1 = self.p[w] / (self.r * self.t2[w])
                var2 = self.gas_flow[1][w] / (self.gas_flow[0][w] + self.gas_flow[1][w] + self.gas_flow[2][w])
                self.gas_con[1][w] = var1 * var2
                if self.gas_con[1][w] >= p_sat.water.calc_psat(self.t2[w]) \
                        / (self.r * self.t2[w]):  # saturated
                    a = p_sat.water.calc_psat(self.t2[w])
                    b = self.gas_flow[0][w] + self.gas_flow[2][w]
                    c = self.gas_flow[0][w] / b
                    d = self.gas_flow[2][w] / b
                    e = self.r * self.t2[w]
                    self.gas_con[0][w] = (self.p[w] - a) / e * c
                    self.gas_con[1][w] = a / e
                    self.gas_con[2][w] = (self.p[w] - a) / e * d
                else:  # not saturated
                    var3 = self.p[w] / (self.r * self.t2[w])
                    var4 = (self.gas_flow[0][w] + self.gas_flow[1][w] + self.gas_flow[2][w])
                    var5 = var3 / var4
                    self.gas_con[0][w] = var5 * self.gas_flow[0][w]
                    self.gas_con[2][w] = var5 * self.gas_flow[2][w]
        else:
            for w in range(self.channel.nodes):
                var1 = self.p[w] / (self.r * self.t2[w])
                var2 = self.gas_flow[1][w] / (self.gas_flow[0][w] + self.gas_flow[1][w])
                self.gas_con[1][w] = var1 * var2
                if self.gas_con[1][w] >= p_sat.water.calc_psat(self.t2[w]) \
                        / (self.r * self.t2[w]):  # saturated
                    a = (self.p[w]
                         - p_sat.water.calc_psat(self.t2[w]))\
                        / (self.r * self.t2[w])
                    self.gas_con[0][w] = a
                    self.gas_con[1][w] = p_sat.water.calc_psat(self.t2[w]) \
                                 / (self.r * self.t2[w])
                else:  # not saturated
                    var3 = self.p[w] / (self.r * self.t2[w])
                    var4 = self.gas_flow[0][w] / (self.gas_flow[0][w] + self.gas_flow[1][w])
                    self.gas_con[0][w] = var3 * var4

    def calc_mass_fraction(self):
        if self.spec_numb is 2:
            a = self.gas_con[0] * self.m_h2
            b = self.gas_con[1] * self.m_h2o
            self.mf_1 = a / (a + b)
            self.mf_2 = b / (a+b)
        else:
            a = self.gas_con[0] * self.m_o2
            b = self.gas_con[1] * self.m_h2o
            c = self.gas_con[2] * self.m_n2
            self.mf_1 = a / (a + b + c)
            self.mf_2 = b / (a + b + c)
            self.mf_3 = c / (a + b + c)

    def calc_r_mix(self):
        if self.spec_numb is 2:
            self.r_mix = self.mf_1 * self.r_h2 + self.mf_2 * self.r_h2o
        else:
            self.r_mix = self.mf_1 * self.r_o2\
                         + self.mf_2 * self.r_h2o\
                         + self.mf_3 * self.r_n2

    def calc_flow_velocity(self):
        self.u = self.q_sum * self.r * self.t2 / (self.p * self.channel.cross_area)

    def calc_roh(self):
        self.roh = g_func.calc_roh(self.p, self.r_mix, self.t2)
        #print(self.roh)

    def calc_re(self):
        self.re = g_func.calc_re(self.roh, self.u, self.dh, self.visc)
        #print(self.re)

    def calc_fluid_water(self):  # fluid water in the channel#nodewise
        self.w = self.gas_flow[1] - self.gas_con[1] / self.gas_con[0] * self.gas_flow[0]

    def calc_cond_rates(self):  # condensation rates of the vapour in the channel
        if self.spec_numb is 3:
            self.gamma = np.gradient(self.w, self.channel.d_x)
        else:
            self.gamma = -np.gradient(self.w, self.channel.d_x)

    def calc_rel_humidity(self):  # relative humidity in the channel
        self.humidity = self.gas_con[1] * self.r * self.t2 \
                        / p_sat.water.calc_psat(self.t2)
        # ok

    def calc_free_water(self):  # membrane free water content#
        a = 0.043 + 17.81 * self.humidity
        b = -39.85 * self.humidity ** 2 + 36. * self.humidity ** 3
        self.free_water = a + b

    def sum_flows(self):
        if self.spec_numb is 2:
            self.q_sum = self.gas_flow[0] + self.gas_flow[1] - self.w
        else:
            self.q_sum = self.gas_flow[0] + self.gas_flow[1] + self.gas_flow[2] - self.w
        # plt.plot(self.q_sum)
        # plt.show()

    def calc_channel_permeability(self):
        if self.spec_numb is 2:
            counter = 1
            for q in range(self.channel.nodes-2, -1, -1):
                self.perm[q] = self.visc / self.channel.cross_area\
                               * counter * self.channel.d_x * self.q_sum[q]\
                               / (-self.p[q]+self.p[self.channel.nodes-1])
                counter = counter + 1
        else:
            for q in range(1, self.channel.nodes):
                self.perm[q] = self.visc / self.channel.cross_area\
                               * q * self.channel.d_x * self.q_sum[q]\
                               / (-self.p[q] + self.p[0])


