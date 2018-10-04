import warnings
import copy
import matrix_database as m_d
import global_functions as g_func
import saturation_pressure_vapour as p_sat
import gas_properties_fit as g_fit
import input_parameter as i_p
import numpy as np
import global_parameter as g_par
import channel as ch

warnings.filterwarnings("ignore")

class HalfCell:

    def __init__(self, dict_hc):
        if dict_hc['type'] is True:
            self.channel = ch.Channel(i_p.channel_cat)
            self.o2_con_in = i_p.o2_con_in
            self.n2o2ratio = (1. - self.o2_con_in) / self.o2_con_in

        else:
            self.channel = ch.Channel(i_p.channel_ano)
        self.spec_numb = dict_hc['spec_numb']
        self.val = dict_hc['val']
        self.type = dict_hc['type']  # Anode is False; Cathode is True
        self.thickness_gde = dict_hc['thick_gde']
        self.thickness_plate = dict_hc['plate_thick']
        self.thickness_cat = dict_hc['cat_thick']
        self.m_m = dict_hc['m_mass']
        self.t_init = dict_hc['t_init']
        self.vol_ex_cd = dict_hc['vol_ex_cd']
        self.prot_con = dict_hc['prot_con']
        self.dif_coef_cat = dict_hc['dif_coef_cat']
        self.dif_coef_gdl = dict_hc['dif_coef_gdl']
        self.tafel_slope = dict_hc['tafel_slope']
        self.init_param()
        self.init_arrays()
        self.init_func()

    def init_param(self):
        self.break_program = False
        self.pem_type = False
        self.j = []
        self.re = []
        self.i_theta = np.sqrt(2. * self.vol_ex_cd * self.prot_con
                               * self.tafel_slope)
        self.index_cat = g_par.dict_case['nodes']-1
        self.p_drop_bends = 0.

    def init_arrays(self):
        self.zero_ele = np.full(g_par.dict_case['nodes']-1, 1.e-50)
        self.zero_node = np.full(g_par.dict_case['nodes'], 1.e-50)
        self.w = np.full(g_par.dict_case['nodes'], 0.)
        self.p = np.full(g_par.dict_case['nodes'], self.channel.p_in)
        self.gamma = np.full(g_par.dict_case['nodes'], 0.)
        self.humidity = np.full(g_par.dict_case['nodes'], 0.)
        self.free_water = np.full(g_par.dict_case['nodes'], 0.)
        self.i = np.full(g_par.dict_case['nodes']-2, g_par.dict_case['tar_cd'])
        self.u = np.full(g_par.dict_case['nodes'], 0.)
        self.node_for = m_d.for_mat_reac_flow(g_par.dict_case['nodes']-1)
        self.node_back = m_d.back_mat_reac_flow(g_par.dict_case['nodes']-1)

        if self.type is True:
            self.t1 = np.full(g_par.dict_case['nodes'], self.t_init)
            self.t2 = np.full(g_par.dict_case['nodes'], self.t_init)
            self.t3 = np.full(g_par.dict_case['nodes'], self.t_init)
        else:
            self.t1 = np.full(g_par.dict_case['nodes'], self.t_init)
            self.t2 = np.full(g_par.dict_case['nodes'], self.t_init)
            # 1 catalyst layer, 2 channel layer, 3 coolant plate layer
        self.perm = np.full(g_par.dict_case['nodes'], 0.)
        self.gas_flow = np.full((self.spec_numb, g_par.dict_case['nodes']), 0.)
        self.gas_con = np.full((self.spec_numb, g_par.dict_case['nodes']), 0.)
        self.gas_con_ele = g_func.calc_elements(self.gas_con[0])
        self.mixture = [None] * g_par.dict_case['nodes']
        self.t_gas = np.full(g_par.dict_case['nodes'], self.channel.t_in)
        self.rho = np.full(g_par.dict_case['nodes'], 1.)
        self.visc_mix = np.full(g_par.dict_case['nodes'], 1.e-5)
        self.nu = np.full(g_par.dict_case['nodes'], 0.)
        self.mol_f = np.full((self.spec_numb, g_par.dict_case['nodes']), 0.)
        self.mf = np.full((self.spec_numb, g_par.dict_case['nodes']), 0.)
        self.r_mix = np.full((self.spec_numb, g_par.dict_case['nodes']), 0.)
        self.r_el = np.full(self.spec_numb, 0.)
        self.cp = np.full((self.spec_numb, g_par.dict_case['nodes']), 0.)
        self.lam = np.full((self.spec_numb, g_par.dict_case['nodes']), 0.)
        self.visc = np.full((self.spec_numb, g_par.dict_case['nodes']), 0.)

    def init_func(self):
        self.dh = 4. * self.channel.cross_area / self.channel.extent
        for q in range(len(self.m_m)):
            self.r_el[q] = g_par.dict_uni['r'] / self.m_m[q] * 1.e3

    def update(self):
        self.calc_t_gas_e()
        self.calc_reac_flow()
        self.calc_water_flow()
        self.calc_con()
        self.calc_voltage_losses_parameter()
        if self.break_program is False:
            # print(self.stoi,'stoi')
            self.calc_activation_losses()
            # print(self.act_ov,'act_ov')
            self.calc_transport_losses_catalyst_layer()
            # print(self.cat_dif_los,'cat_dif')
            self.calc_transport_losses_diffusion_layer()
            # print(self.gde_dif_los,'gde_dif')
            self.calc_fluid_water()
            self.sum_flows()
            self.calc_cond_rates()
            self.calc_mass_fraction()
            self.calc_mol_fraction()
            self.calc_gas_properties()
            self.calc_gas_mix_properties()
            self.calc_rel_humidity()
            self.calc_flow_velocity()
            self.calc_mass_flow()
            self.calc_m_mix_properties()
            self.calc_re()
            self.calc_nu()
            self.calc_heat_transfer_coef()
            self.calc_pressure_drop_bends()
            self.calc_pressure()
            if self.pem_type is False:  # NT
                self.calc_free_water()

    def set_i(self, i):
        self.i = i

    def set_j(self, j):
        self.j = j[1:]

    def set_stoi(self, stoi):
        self.stoi = stoi

    def set_pem_type(self, pem_type):
        self.pem_type = pem_type

    def set_t(self, var):
        for l, item in enumerate(var):
            exec('self.t%d = var[l]' % (l + 1))

    def calc_t_gas_e(self):
        self.t_gas_copy = copy.deepcopy(self.t_gas)
        self.t_gas_e = g_func.calc_elements(self.t_gas_copy)

    def calc_reac_flow(self):
        var1 = self.stoi * g_par.dict_case['tar_cd']\
               * self.channel.plane / (self.val * g_par.dict_uni['f'])
        if self.type is True:
            self.gas_flow[0, 0] = var1
            self.gas_flow[0, 1:] = var1 - np.matmul(self.node_for, self.i) \
                                   * self.channel.plane_dx\
                                   / (self.val * g_par.dict_uni['f'])

        else:
            self.gas_flow[0, -1] = var1
            self.gas_flow[0, :-1] = var1 - np.matmul(self.node_back, self.i) \
                                    * self.channel.plane_dx \
                                    / (self.val * g_par.dict_uni['f'])
        self.gas_flow[0] = np.maximum(self.gas_flow[0], self.zero_node)

    def calc_water_flow(self):
        b = 0
        if self.type is True:
            q_0_water = self.gas_flow[0][0] \
                        * (1. + self.n2o2ratio) \
                        * p_sat.water.calc_psat(self.channel.t_in) \
                        * self.channel.phi \
                        / (self.channel.p_in
                           - self.channel.phi
                           * p_sat.water.calc_psat(self.channel.t_in))
            a = self.channel.plane_dx / (self.val * g_par.dict_uni['f'] * 0.5) \
                * np.matmul(self.node_for, self.i)
            # production
            if self.pem_type is False:
                b = self.channel.plane_dx \
                   * np.matmul(self.node_for, self.j)
                # crossover
            self.gas_flow[1, 0] = q_0_water
            self.gas_flow[1, 1:] = a + b + q_0_water
            self.gas_flow[2] = np.full(g_par.dict_case['nodes'],
                                       self.gas_flow[0][0] * self.n2o2ratio)
        else:
            q_0_water = self.gas_flow[0][-1] \
                        * p_sat.water.calc_psat(self.channel.t_in) \
                        * self.channel.phi \
                        / (self.channel.p_in
                           - self.channel.phi
                           * p_sat.water.calc_psat(self.channel.t_in))
            if self.pem_type is False:
                b = self.channel.plane_dx\
                    * np.matmul(self.node_back, self.j)
            self.gas_flow[1, -1] = q_0_water
            self.gas_flow[1, :-1] = b + q_0_water

        self.gas_flow[1] = np.choose(self.gas_flow[0] > 1.e-50,
                                     [self.zero_node, self.gas_flow[1]])
        if self.type is True:
            for w in range(1, g_par.dict_case['nodes']):
                if self.gas_flow[1, w] < 1.e-49:
                    self.index_cat = w - 1
                    self.gas_flow[1, -self.index_cat - 1:] =\
                        self.gas_flow[1, self.index_cat]
                    break

    def calc_pressure_drop_bends(self):
        self.p_drop_bends = self.channel.bend_fri_fac * np.average(self.rho)\
                            * np.average(self.u)**2. * self.channel.bend_numb\
                            / (g_par.dict_case['nodes']-1) * .5

    def calc_pressure(self):
        self.rho_ele = g_func.calc_elements(self.rho)
        self.u_ele = g_func.calc_elements(self.u)
        self.re_ele = g_func.calc_elements(self.re)
        if self.type is True:
            mat = self.node_for
            self.p[0] = self.channel.p_in
            self.p[1:] = self.channel.p_in + 32. / self.dh\
                         * np.matmul(-mat, self.rho_ele * self.u_ele ** 2.
                                     / self.re_ele) \
                         * self.channel.d_x - self.p_drop_bends
        else:
            mat = self.node_back
            self.p[-1] = self.channel.p_in
            self.p[:-1] = self.channel.p_in + 32. / self.dh\
                          * np.matmul(-mat, self.rho_ele * self.u_ele**2.
                                      / self.re_ele)\
                          * self.channel.d_x - self.p_drop_bends

    def calc_con(self):
        for w in range(g_par.dict_case['nodes']):
            id_lw = self.p[w] / (g_par.dict_uni['r'] * self.t2[w])
            var4 = np.sum(self.gas_flow[:, w])
            var2 = self.gas_flow[1][w] / var4
            self.gas_con[1][w] = id_lw * var2
            a = p_sat.water.calc_psat(self.t2[w])
            e = g_par.dict_uni['r'] * self.t_gas[w]
            if self.gas_con[1][w] >= a / e:  # saturated
                if self.type is True:
                    b = self.gas_flow[0][w] + self.gas_flow[2][w]
                    c = self.gas_flow[0][w] / b
                    d = self.gas_flow[2][w] / b
                    self.gas_con[0][w] = (self.p[w] - a) / e * c
                    self.gas_con[2][w] = (self.p[w] - a) / e * d
                else:
                    self.gas_con[0][w] = (self.p[w] - a) / e
                self.gas_con[1][w] = a / e
            else:  # not saturated
                var5 = id_lw / var4
                self.gas_con[0][w] = var5 * self.gas_flow[0][w]
                if self.type is True:
                    self.gas_con[2][w] = var5 * self.gas_flow[2][w]
        if self.type is False:
            self.gas_con[1] = self.zero_node
        self.gas_con_ele = g_func.calc_elements(self.gas_con[0])

    def calc_mass_fraction(self):
        temp_var_1 = []
        for q, item in enumerate(self.gas_con):
            temp_var_1.append(self.gas_con[q] * self.m_m[q])
        for q, item in enumerate(self.gas_con):
            self.mf[q] = temp_var_1[q] / sum(temp_var_1)

    def calc_mol_fraction(self):
        for q in range(len(self.gas_con)):
            self.mol_f[q] = self.gas_con[q] / sum(self.gas_con)

    def calc_gas_properties(self):
        if self.type is True:
            self.cp[0] = g_fit.oxygen.calc_cp(self.t_gas)
            self.cp[2] = g_fit.nitrogen.calc_cp(self.t_gas)
            self.lam[0] = g_fit.oxygen.calc_lambda(self.t_gas, self.p)
            self.lam[2] = g_fit.nitrogen.calc_lambda(self.t_gas, self.p)
            self.visc[0] = g_fit.oxygen.calc_visc(self.t_gas)
            self.visc[2] = g_fit.nitrogen.calc_visc(self.t_gas)
            self.cp[1] = g_fit.water.calc_cp(self.t_gas)
            self.lam[1] = g_fit.water.calc_lambda(self.t_gas, self.p)
            self.visc[1] = g_fit.water.calc_visc(self.t_gas)
        else:
            self.cp[0] = g_fit.hydrogen.calc_cp(self.t_gas)
            self.lam[0] = g_fit.hydrogen.calc_lambda(self.t_gas, self.p)
            self.visc[0] = g_fit.hydrogen.calc_visc(self.t_gas)
            self.cp[1] = g_fit.water.calc_cp(self.t_gas)
            self.lam[1] = g_fit.water.calc_lambda(self.t_gas, self.p)
            self.visc[1] = g_fit.water.calc_visc(self.t_gas)
        self.cpe = g_func.calc_elements(self.cp[0])

    def calc_gas_mix_properties(self):
        temp1, temp2 = [], []
        for q in range(self.spec_numb):
            temp1.append(self.mf[q] * self.r_el[q])
            temp2.append(self.mf[q] * self.cp[q])
        self.r_mix = sum(temp1)
        self.cp_mix = sum(temp2)
        self.cp_mix_ele = g_func.calc_elements(self.cp_mix)
        self.visc_mix = g_func.calc_visc_mix(self.visc, self.mol_f, self.m_m)
        self.lambda_mix = g_func.calc_lambda_mix(self.lam, self.mol_f,
                                                 self.visc, self.m_m)
        self.rho = g_func.calc_rho(self.p, self.r_mix, self.t_gas)
        self.pr = self.visc_mix * self.cp_mix / self.lambda_mix

    def calc_flow_velocity(self):
        self.u = self.q_sum * g_par.dict_uni['r'] * self.t_gas\
                 / (self.p * self.channel.cross_area)

    def calc_mass_flow(self):
        self.m_flow = self.u * self.rho * self.channel.cross_area
        self.m_reac_flow = self.gas_flow[0] * self.m_m[0] * 1.e-3
        self.m_liq_water = self.w * self.m_m[1] * 1.e-3
        self.m_vap_water_flow = (self.gas_flow[1]-self.w) * self.m_m[1] * 1.e-3
        self.m_reac_flow_delta = abs(g_func.calc_dif(self.m_reac_flow))
        self.m_vap_water_flow_delta = abs(g_func.calc_dif(self.m_vap_water_flow))

    def calc_m_mix_properties(self):
        self.m_full_flow = self.m_flow + self.m_liq_water
        self.cp_full = (self.m_flow * self.cp_mix
                        + self.m_liq_water * g_par.dict_uni['cp_liq'])\
                        / self.m_full_flow
        self.g_full = self.m_full_flow * self.cp_full
        self.g_full_e = g_func.calc_elements(self.g_full)

    def calc_re(self):
        self.re = g_func.calc_re(self.rho, self.u, self.dh, self.visc_mix)

    def calc_nu(self):
        for w in range(0, g_par.dict_case['nodes'], g_par.dict_case['nodes']-1):
            if 0. < self.re[w] < 2300.:
                self.nu[w] = 3.66
            elif 2300. <= self.re[w]:
                lama = (self.re[w] - 2300.)/(1.e4 - 2300.)
                self.nu[w] = (1. - lama) * 3.66 \
                             + lama\
                             * g_func.calc_turbo_nu_numb(1.e4,
                                                         self.pr[w],
                                                         self.dh,
                                                         self.channel.d_x * w)
            else:
                self.nu[w] = g_func.calc_turbo_nu_numb(self.re[w],
                                                       self.pr[w],
                                                       self.dh,
                                                       self.channel.d_x * w)

    def calc_heat_transfer_coef(self):
        self.ht_coef = self.lambda_mix * self.nu / self.dh
        self.ht_coef = np.sum(self.ht_coef) / 2.
        self.r_ht_coef_a = 1. / (self.ht_coef * np.pi * self.channel.d_x * self.dh)

    def calc_fluid_water(self):  # fluid water in the channel
        self.w = self.gas_flow[1] - self.gas_con[1] / self.gas_con[0]\
                 * self.gas_flow[0]

    def calc_cond_rates(self):  # condensation rates of the vapour in the channel
        if self.type is True:
            self.gamma = g_func.calc_nodes_1d(np.ediff1d(self.w))
        else:
            self.gamma = -g_func.calc_nodes_1d(np.ediff1d(self.w))

    def calc_rel_humidity(self):  # relative humidity in the channel
        self.humidity = self.gas_con[1] * g_par.dict_uni['r'] * self.t_gas \
                        / p_sat.water.calc_psat(self.t_gas)

    def calc_free_water(self):  # membrane free water content
        a = 0.043 + 17.81 * self.humidity
        b = -39.85 * self.humidity ** 2 + 36. * self.humidity ** 3
        self.free_water = a + b

    def sum_flows(self):
        self.q_sum = sum(self.gas_flow) - self.w

    def calc_voltage_losses_parameter(self):
        self.i_star = self.prot_con\
                      * self.tafel_slope\
                      / self.thickness_cat
        self.i_lim = 4. * g_par.dict_uni['f'] * self.gas_con[0, :-1]\
                     * self.dif_coef_gdl / self.thickness_gde
        self.i_hat = self.i / self.i_star
        short_save = np.sqrt(2. * self.i_hat)
        self.beta = short_save / (1. + np.sqrt(1.12 * self.i_hat)
                                  * np.exp(short_save))\
                    + np.pi * self.i_hat/(2. + self.i_hat)
        self.var = 1. - self.i\
                   / (self.i_lim * self.gas_con_ele / self.gas_con[0, :-1])
        self.i_square = self.i ** 2.

    def calc_activation_losses_tafel(self):
        self.act_ov = g_par.dict_uni['r'] * g_func.calc_elements(self.t2)\
                      / (g_par.dict_uni['f'] * 0.54)\
                      * np.log10(self.i * self.gas_con[0, :-1]
                                 / (self.vol_ex_cd
                                    * self.thickness_cat
                                    * self.gas_con_ele))

    def calc_activation_losses(self):
        self.act_ov = self.tafel_slope\
                      * np.arcsinh((self.i/self.i_theta)**2.
                                   / (2. * (self.gas_con_ele
                                            / self.gas_con[0, :-1])
                                      * (1. - np.exp(-self.i /
                                                     (2. * self.i_star)))))

    def calc_transport_losses_catalyst_layer(self):
        self.cat_dif_los = ((self.prot_con * self.tafel_slope ** 2.)
                            / (4. * g_par.dict_uni['f'] * self.dif_coef_cat
                               * self.gas_con_ele)
                            * (self.i / self.i_star
                               - np.log10(1. + self.i_square
                                          / (self.i_star ** 2.
                                             * self.beta ** 2.))))\
                           / self.var

    def calc_transport_losses_diffusion_layer(self):
        self.gde_dif_los = -self.tafel_slope * np.log10(self.var)
        nan_list = np.isnan(self.gde_dif_los)
        bol = nan_list.any()
        if bol == True:
            print('True')
            self.gde_dif_los[np.argwhere(nan_list)[0, 0]:] = 1.e20
