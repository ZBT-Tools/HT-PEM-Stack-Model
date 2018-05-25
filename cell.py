import numpy as np
import global_functions as g_func


class Cell:

    def __init__(self, anode, cathode, gamma, alpha, mem_thic, mu_p, mu_g,
                 mu_m, e0, c_ref, i_ref, delta, pem_type, coolant_channel,
                 t_cool_in, width, thickness_gde, thickness_plate, t_u):
        self.anode = anode
        self.cathode = cathode
        self.gamma = gamma
        self.alpha = alpha
        self.mem_thic = mem_thic
        self.mu_p = mu_p
        self.mu_g = mu_g
        self.mu_m = mu_m
        self.e0 = e0
        self.c_ref = c_ref
        self.i_ref = i_ref
        self.delta = delta
        self.pem_type = pem_type
        self.coolant_channel = coolant_channel
        self.t_cool_in = t_cool_in
        self.width = width
        self.thickness_gde = thickness_gde
        self.thickness_plate = thickness_plate
        self.t_u = t_u
        self.extend = 2. * (self.cathode.channel.length + self.width)
        self.plane = self.cathode.channel.length * self.width
        self.j = np.full(self.cathode.channel.elements + 1, 0.)
        self.ke = 0.
        self.ho = -52300.
        self.ko = 6.2
        self.m_a = np.full(self.cathode.channel.nodes, 0.)
        self.m_c = np.full(self.cathode.channel.nodes, 0.)
        self.m0 = np.full(self.cathode.channel.nodes, 0.)
        self.m1 = np.full(self.cathode.channel.nodes, 0.)
        self.m_pos0 = np.full(self.cathode.channel.nodes, 0.)
        self.m_pos1 = np.full(self.cathode.channel.nodes, 0.)
        self.omega = np.full(self.cathode.channel.nodes, 0.)
        self.psi = np.full(self.cathode.channel.nodes, 0.)
        self.t = np.full(self.cathode.channel.nodes, self.t_cool_in)
        self.t0 = (self.anode.t1 + self.cathode.t1) / 2.
        self.t1 = self.cathode.t1
        self.t2 = self.cathode.t2
        self.t3 = self.cathode.t3
        self.t4 = self.anode.t1
        self.t5 = self.anode.t2
        self.vtn = 1.28

    def update(self):
        if self.pem_type is False:
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
            self.calc_mem_block_1()
            self.calc_cross_over_water_flux()
            self.calc_mem_block_2()
            self.calc_con_overpotential()
            self.calc_mem_resistivity()
            self.calc_voltage()
            self.calc_dvdi()
            self.t0 = 0.5 * self.anode.t1 + 0.5 * self.cathode.t1
        else:
            self.cathode.set_pem_type(True)
            self.cathode.set_i(self.i)
            self.cathode.set_t([self.t1, self.t2, self.t3])
            self.cathode.update()
            self.anode.set_pem_type(True)
            self.anode.set_i(self.i)
            self.anode.set_t([self.t4, self.t5])
            self.anode.update()
            # self.calc_con_overpotential() need a different function here
            self.psi = np.full(self.cathode.channel.nodes, 0.)
            self.calc_mem_resistivity()
            self.calc_voltage()
            self.calc_dvdi()
            self.t0 = 0.5 * self.anode.t1 + 0.5 * self.cathode.t1

    def set_i(self, i):
        self.i = i

    def calc_mem_block_1(self):  # nodewise
        a = self.cathode.free_water + self.anode.free_water
        b = np.reshape(self.i
                       / (2. * self.gamma * self.alpha * self.cathode.f_const), self.cathode.channel.nodes)
        self.zeta_plus = a + b
        c = self.cathode.free_water - self.anode.free_water
        d = np.reshape(5. * self.i
                       / (2. * self.gamma * self.alpha * self.cathode.f_const), self.cathode.channel.nodes)
        e = 1. + g_func.dw(self.t0) * self.zeta_plus \
            / (self.mem_thic * self.gamma)
        self.zeta_negative = (c + d) / e
        self.m_c = 0.5 * (self.zeta_plus + self.zeta_negative)
        self.m_a = 0.5 * (self.zeta_plus - self.zeta_negative)
        # ok

    def calc_cross_over_water_flux(self):  # nodewise
        a = np.reshape(self.i / self.cathode.f_const, self.cathode.channel.nodes)
        aa = self.alpha * g_func.dw(self.t0)
        b = (self.m_a ** 2. - self.m_c ** 2.) / (2. * self.mem_thic)
        self.j = a + aa * b
        # ok

    def calc_mem_block_2(self):  # nodewise
        self.ke = self.ko * np.exp(-self.ho / self.cathode.r
                                   * (1. / self.t0 - 1. / self.t_u))
        a = self.m_c ** 2.
        b = self.m_a ** 2. - self.m_c ** 2.
        self.m0 = np.sqrt(a - b)  # z = 0
        self.m1 = np.sqrt(a)  # z = 1
        self.m_pos0 = -self.ke * self.m0 * 0.5\
                      + np.sqrt((self.ke * self.m0 * 0.5) ** 2 + self.ke * self.m0)
        self.m_pos1 = -self.ke * self.m1 * 0.5\
                      + np.sqrt((self.ke * self.m1 * 0.5) ** 2 + self.ke * self.m1)
        # ok

    def calc_con_overpotential(self):  # nodewise
        self.psi = self.cathode.r * self.t0 / self.cathode.f_const\
                   * np.log(self.m_pos1 / self.m_pos0)
        # ok

    def calc_mem_resistivity(self):  # nodewise
        self.omega = (0.4025 - 0.0007 * self.t0) * 1.e-4
        # ok

    def calc_voltage(self):  # nodewise
        a = self.e0 - self.i * self.omega
        b = self.cathode.r * self.t1 / self.cathode.f_const
        c = np.log(self.i * self.c_ref
                   / (self.i_ref * (self.cathode.c1 - self.delta * self.i)))
        self.v = a + self.psi - b * c
        # ok

    def calc_dvdi(self):  # nodewise
        a = -self.omega
        b = -self.cathode.r * self.t1 / (self.cathode.f_const)
        c = self.cathode.c1 / (self.i * (self.cathode.c1 - self.delta * self.i))
        self.dv = a + b * c
        # ok
