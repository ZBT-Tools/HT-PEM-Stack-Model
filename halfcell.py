import matrix_database
import global_functions
import saturation_pressure_vapour as p_sat
import gas_properties_fit as g_fit
import input_parameter as i_p
import numpy as np
import global_parameter as g_par
import channel as ch
# from numba import jit
import matplotlib.pyplot as plt


class Halfcell:

    def __init__(self, dict):
        if dict['type'] is True:
            self.channel = ch.Channel(i_p.channel_cat)
            self.o2_con_in = i_p.o2_con_in
            self.n2o2ratio = (1. - self.o2_con_in) / self.o2_con_in

        else:
            self.channel = ch.Channel(i_p.channel_ano)
        self.spec_numb = dict['spec_numb']
        self.val = dict['val']
        self.type = dict['type']  # Anode is False Cathode is True
        self.gas_con_ref = dict['gas_con_ref']
        self.act_energy = dict['act_energy']
        self.sym = dict['sym_fac']
        self.thickness_gde = dict['thick_gde']
        self.thickness_plate = dict['plate_thick']
        self.thickness_cat = dict['cat_thick']
        self.m_m = dict['m_mass']
        self.t_init = dict['t_init']
        self.vol_ex_cd = dict['vol_ex_cd']
        self.proton_conductivity = dict['prot_con']
        self.dif_coef_cat = dict['dif_coef_cat']
        self.dif_coef_gdl = dict['dif_coef_gdl']
        self.tafel_slope = dict['tafel_slope']
        self.init_param()
        self.init_arrays()
        self.init_func()

    def init_param(self):
        self.pem_type = False
        self.j = []
        self.re = []
        self.i_theta = np.sqrt(2. * self.vol_ex_cd * self.proton_conductivity
                               * self.tafel_slope)
        self.index_cat = g_par.dict_case['nodes']-1

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
        self.node_forward = matrix_database.forward_matrix_reac_flow(g_par.dict_case['nodes']-1)
        self.node_backward = matrix_database.backward_matrix_reac_flow(g_par.dict_case['nodes']-1)

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
        self.gas_con_ele = global_functions.calc_elements(self.gas_con[0])
        self.mixture = [None] * g_par.dict_case['nodes']
        self.t_gas = np.full(g_par.dict_case['nodes'], self.channel.t_in)
        self.t_gas_test = np.full(g_par.dict_case['nodes'], self.channel.t_in)
        self.rho = np.full(g_par.dict_case['nodes'], 1.)
        self.visc_mix = np.full(g_par.dict_case['nodes'], 1.e-5)
        self.nu = np.full(g_par.dict_case['nodes'], 0.)
        self.mol_f = np.full((self.spec_numb, g_par.dict_case['nodes']), 0.)
        self.mf = np.full((self.spec_numb, g_par.dict_case['nodes']), 0.)
        self.r_mix = np.full((self.spec_numb, g_par.dict_case['nodes']), 0.)
        self.r_el = np.full(self.spec_numb,0.)
        self.cp = np.full((self.spec_numb, g_par.dict_case['nodes']), 0.)
        self.lam = np.full((self.spec_numb, g_par.dict_case['nodes']), 0.)
        self.visc = np.full((self.spec_numb, g_par.dict_case['nodes']), 0.)

    def init_func(self):
        self.dh = 4. * self.channel.cross_area / self.channel.extent
        for q in range(len(self.m_m)):
            self.r_el[q] = g_par.dict_uni['r'] / self.m_m[q] * 1.e3

    def update(self):
        #print(self.i)
        self.calc_reac_flow()
        self.calc_water_flow()
        #print(self.gas_flow,'g_flow')
        self.calc_con()
        #print(self.gas_con_ele,'g_con')
        self.calc_activation_losses()
        #print(self.act_ov, 'act_ov')
        #print(self.act_dov, 'act_ov_dv')
        self.calc_transport_losses_catalyst_layer()
        #print(self.cat_dif_los,'cat_dif_los')
        #print(self.cat_dif_los_dv,'cat_dif_los_dv')
        self.calc_transport_losses_diffusion_layer()
        #print(self.gde_dif_los, 'gde_dif_los')
        #print(self.gde_dif_los_dv, 'gde_dif_los_dv')
        #if self.type is True:
            #plt.plot(self.gde_dif_los_dv)
            #self.cat_dif_los[-1] = 2.*self.cat_dif_los[-2]-self.cat_dif_los[-3]
            #self.cat_dif_los_dv[-1] = 2. * self.cat_dif_los_dv[-2] - self.cat_dif_los_dv[-3]
            #self.gde_dif_los[-1] = 2. * self.gde_dif_los[-2] - self.gde_dif_los[-3]
            #self.gde_dif_los_dv[-1] = 2. * self.gde_dif_los_dv[-2] - self.gde_dif_los_dv[-3]
            #self.cat_dif_los_dv[-1] = self.cat_dif_los_dv[-2]
            #self.cat_dif_los[-1] = self.cat_dif_los[-2]
            #self.gde_dif_los_dv[-1] = 0. #self.gde_dif_los_dv[-2]
            #self.gde_dif_los[-1] = self.gde_dif_los[-2]
            #plt.plot(self.gde_dif_los, color='k')
         #   plt.plot(self.cat_dif_los, color='r')
            #plt.plot(self.gde_dif_los_dv, color='y')
            #plt.show()
           # print(self.gde_dif_los_dv, 'gde_dif_los_dv')
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
        self.calc_re()
        self.calc_nu()
        self.calc_heat_transfer_coef()
        self.calc_heat_flow()
        self.calc_pressure_zimmer_lam()
        if self.pem_type is False:  # NT
            self.calc_free_water()

    def set_i(self, i):
        self.i = i

    def set_j(self, j):
        self.j = j[1:]
        #print(self.j)

    def set_stoi(self, stoi):
        self.stoi = stoi

    def set_pem_type(self, pem_type):
        self.pem_type = pem_type

    def set_t(self, var):
        for l, item in enumerate(var):
            exec('self.t%d = var[l]' % (l + 1))

    #@jit
    def calc_reac_flow(self):
        var1 = self.stoi * g_par.dict_case['tar_cd']\
               * self.channel.plane / (self.val * g_par.dict_uni['f'])
        if self.type is True:
            self.gas_flow[0, 0] = var1
            self.gas_flow[0, 1:] = var1 - np.matmul(self.node_forward, self.i)\
                                   * self.channel.plane_dx\
                                   / (self.val * g_par.dict_uni['f'])

        else:
            self.gas_flow[0, -1] = var1
            self.gas_flow[0, :-1] = var1 - np.matmul(self.node_backward, self.i)\
                                    * self.channel.plane_dx\
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
                * np.matmul(self.node_forward, self.i)
            # production
            if self.pem_type is False:
                b = self.channel.plane_dx \
                   * np.matmul(self.node_forward, self.j)
                # crossover
            self.gas_flow[1,0] = q_0_water
            self.gas_flow[1,1:] = a + b + q_0_water
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
                    * np.matmul(self.node_backward, self.j)
            self.gas_flow[1, -1] = q_0_water
            self.gas_flow[1, :-1] = b + q_0_water

        self.gas_flow[1] = np.choose(self.gas_flow[0] > 1.e-50,
                                     [self.zero_node, self.gas_flow[1]])
        if self.type is True:
            for w in range(1, g_par.dict_case['nodes']):
                if self.gas_flow[1, w] < 1.e-49:
                    self.index_cat = w - 1
                    self.gas_flow[1, -self.index_cat - 1:] = self.gas_flow[1, self.index_cat]
                    break

    def calc_pressure_zimmer_lam(self):
        self.rho_ele = global_functions.calc_elements(self.rho)
        self.u_ele = global_functions.calc_elements(self.u)
        self.re_ele = global_functions.calc_elements(self.re)
        if self.type is True:
            mat = self.node_forward
            self.p[0] = self.channel.p_in
            self.p[1:] = self.channel.p_in + 32. / self.dh\
                         * np.matmul(-mat, self.rho_ele * self.u_ele ** 2. / self.re_ele) \
                         * self.channel.d_x
        else:
            mat = self.node_backward
            self.p[-1] = self.channel.p_in
            self.p[:-1] = self.channel.p_in + 32. / self.dh\
                     * np.matmul(-mat, self.rho_ele * self.u_ele**2. /self.re_ele)\
                     * self.channel.d_x

    def calc_con(self):  # Gas concentrations in the channel [moles/mÂ³]
        for w in range(g_par.dict_case['nodes']):
            id_lw = self.p[w] / (g_par.dict_uni['r'] * self.t_gas[w])
            var4 = np.sum(self.gas_flow[:, w])
            var2 = self.gas_flow[1][w] / var4
            self.gas_con[1][w] = id_lw * var2
            a = p_sat.water.calc_psat(self.t_gas[w])
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
        self.gas_con_ele = global_functions.calc_elements(self.gas_con[0])

    def calc_mass_fraction(self):
        temp_var_1 = []
        for q, item in enumerate(self.gas_con):
            temp_var_1.append(self.gas_con[q] * self.m_m[q])
        for q, item in enumerate(self.gas_con):
            self.mf[q] = temp_var_1[q] / sum(temp_var_1)

    def calc_mol_fraction(self):
        for q in range(len(self.gas_con)):  # Init flow and con arrays (1-n)
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

    def calc_gas_mix_properties(self):
        temp1, temp2 = [], []
        for q in range(self.spec_numb):
            temp1.append(self.mf[q] * self.r_el[q])
            temp2.append(self.mf[q] * self.cp[q])
        self.r_mix = sum(temp1)
        self.cp_mix = sum(temp2)
        self.visc_mix = global_functions.calc_visc_mix(self.visc, self.mol_f, self.m_m)
        self.lambda_mix = global_functions.calc_lambda_mix(self.lam, self.mol_f, self.visc, self.m_m)
        self.rho = global_functions.calc_rho(self.p, self.r_mix, self.t_gas)
        self.pr = self.visc_mix * self.cp_mix / self.lambda_mix

    def calc_flow_velocity(self):
        self.u = self.q_sum * g_par.dict_uni['r'] * self.t_gas\
                 / (self.p * self.channel.cross_area)

    def calc_mass_flow(self):
        self.m_flow = self.u * self.rho * self.channel.cross_area

    def calc_heat_flow(self):
        self.h_flow = self.cp_mix * self.m_flow #/ self.thickness_plate
        self.h_flow = np.average(self.h_flow)

    def calc_re(self):
        self.re = global_functions.calc_re(self.rho, self.u, self.dh, self.visc_mix)

    def calc_nu(self):
        for w in range(0, g_par.dict_case['nodes'], g_par.dict_case['nodes']-1):
            if 0. < self.re[w] < 2300.:
                self.nu[w] = 3.66
            elif 2300. <= self.re[w]:
                lama = (self.re[w] - 2300.)/(1.e4 - 2300.)
                self.nu[w] = (1. - lama) * 3.66 \
                             + lama\
                             * global_functions.calc_turbo_nu_numb(1.e4,
                                                                   self.pr[w],
                                                                   self.dh,
                                                                   self.channel.d_x * w)
            else:
                self.nu[w] = global_functions.calc_turbo_nu_numb(self.re[w],
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
            self.gamma = np.gradient(self.w, self.channel.d_x)
        else:
            self.gamma = -np.gradient(self.w, self.channel.d_x)

    def calc_rel_humidity(self):  # relative humidity in the channel
        self.humidity = self.gas_con[1] * g_par.dict_uni['r'] * self.t_gas \
                        / p_sat.water.calc_psat(self.t_gas)

    def calc_free_water(self):  # membrane free water content
        a = 0.043 + 17.81 * self.humidity
        b = -39.85 * self.humidity ** 2 + 36. * self.humidity ** 3
        self.free_water = a + b

    def sum_flows(self):
        self.q_sum = sum(self.gas_flow) - self.w

    def calc_activation_losses(self):
        self.i_star = self.proton_conductivity\
                      * self.tafel_slope\
                      / self.thickness_cat
        self.act_ov = self.tafel_slope\
                      * np.arcsinh((self.i/self.i_theta)**2.
                                   /(2. * (self.gas_con_ele / self.gas_con_ref)
                                     * (1. - np.exp(-self.i / (2. * self.i_star)))))
        #print(len(self.act_ov), self.type)
        self.act_dov = -.25 * self.tafel_slope *\
                       (4. * self.gas_con_ref * self.i / (self.gas_con_ele * self.i_theta**2
                                                          * (np.exp(-0.5 * self.i * self.i_star) - 1.))
                        + self.gas_con_ref * self.i**2.
                        * np.exp(-0.5 * self.i / self.i_star)
                        / (self.gas_con_ele * self.i_star * self.i_theta**2.
                           * (np.exp(-.5 * self.i / self.i_star) - 1.)**2.))\
                       / np.sqrt(0.25 * self.gas_con_ref**2. *
                                 self.i**4. / (self.gas_con_ele**2.
                                               * self.i_theta**4. *
                                               (np.exp(-0.5 * self.i / self.i_star)
                                                - 1.)**2.) + 1.)


    def calc_transport_losses_catalyst_layer(self):
        self.i_lim = 4. * g_par.dict_uni['f'] * self.gas_con[0]\
                     * self.dif_coef_gdl / self.thickness_gde
        self.i_lim = self.i_lim[:-1]
        i_hat = self.i / self.i_star
        beta = np.sqrt(2. * i_hat) / (1. + np.sqrt(1.12 * i_hat)
                                      * np.exp(2. * i_hat))\
               + np.pi * i_hat/(2. + i_hat)
        self.var = 1. - self.i/(self.i_lim * self.gas_con_ele / self.gas_con_ref)
        self.cat_dif_los = ((self.proton_conductivity * self.tafel_slope**2.)
                            / (4. * g_par.dict_uni['f'] * self.thickness_cat
                               * self.gas_con_ele)
                            * (self.i / self.i_star
                               - np.log10(1. + self.i**2.
                                          / (self.i_star**2. * beta**2.))))\
                           / self.var
        self.cat_dif_los_dv = -1./.4 * self.tafel_slope**2.\
                              * self.proton_conductivity\
                              * ((self.i**2. * (2. * np.sqrt(2.)
                                                * (2.11660104885167
                                                   * np.sqrt(self.i / self.i_star)
                                                   * np.exp(2. * self.i / self.i_star)
                                                   / self.i_star
                                                   + 0.529150262212918
                                                   * np.exp(2. * self.i / self.i_star)
                                                   / (self.i_star
                                                     * np.sqrt(self.i / self.i_star)))
                                                * np.sqrt(self.i / self.i_star)
                                                / (1.05830052442584
                                                   * np.sqrt(self.i / self.i_star)
                                                   * np.exp(2. * self.i / self.i_star) + 1.) ** 2.
                                                - 2. * np.pi / (self.i_star
                                                                * (self.i / self.i_star + 2.))
                                                + 2. * np.pi * self.i
                                                / (self.i_star ** 2. * (self.i / self.i_star + 2.)**2.)
                                                - np.sqrt(2.)
                                                / ((1.05830052442584 * np.sqrt(self.i / self.i_star)
                                                    * np.exp(2.*self.i / self.i_star) + 1.)
                                                   * self.i_star * np.sqrt(self.i / self.i_star)))
                                  / (self.i_star**2. * (np.pi * self.i / (self.i_star * (self.i / self.i_star + 2.))
                                                        + np.sqrt(2.) * np.sqrt(self.i / self.i_star)
                                                        / (1.05830052442584 * np.sqrt(self.i / self.i_star)
                                                           * np.exp(2. * self.i/self.i_star) + 1.))**3.)
                                  + 2. * self.i/(self.i_star**2. *(np.pi*self.i
                                                                  / (self.i_star*(self.i/self.i_star + 2.))
                                                                  + np.sqrt(2.) * np.sqrt(self.i/self.i_star)
                                                                  /(1.05830052442584 * np.sqrt(self.i/self.i_star)
                                                                    * np.exp(2.*self.i/self.i_star) + 1.))**2.))
                                 /(self.i**2. / (self.i_star**2. * (np.pi * self.i/(self.i_star*(self.i/self.i_star + 2.))
                                                                + np.sqrt(2.) * np.sqrt(self.i / self.i_star)
                                                                /(1.05830052442584 * np.sqrt(self.i / self.i_star)
                                                                  * np.exp(2.*self.i / self.i_star) + 1.))**2.) + 1.)
                                 - 1. / self.i_star)/(self.gas_con_ele * self.thickness_cat * g_par.dict_uni['f'])\
                              - self.gas_con_ref / (self.gas_con_ele * self.i_lim)

    def calc_transport_losses_diffusion_layer(self):
        self.gde_dif_los = -self.tafel_slope * np.log10(self.var)
        self.gde_dif_los_dv = - self.tafel_slope\
                              * self.gas_con_ref\
                              /(self.gas_con_ele * self.i_lim
                                * (self.gas_con_ref *
                                   self.i / (self.gas_con_ele * self.i_lim) - 1))



###
# error for x = 0 in g_func.calc_nu_turb
# dp/dx for re >2300 is not implemented
# implementation of the Butler-Volmer eq. instead of the tafel eq. still open
# calc_i_zero deactivated values for self.act_energy open
#


