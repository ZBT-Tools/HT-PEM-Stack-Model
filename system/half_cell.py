import warnings
import system.matrix_database as m_d
import system.global_functions as g_func
import data.saturation_pressure_vapour as p_sat
import data.gas_properties_fit as g_fit
import numpy as np
import data.global_parameter as g_par
import system.channel as ch
import data.channel_dict as ch_dict
import input.physical_property as phy_prop

warnings.filterwarnings("ignore")


class HalfCell:

    def __init__(self, dict_hc):
        # Handover
        # Check if the object is an anode or a cathode
        if dict_hc['type'] is True:
            self.channel = ch.Channel(ch_dict.cathode_channel)
            self.o2_con_in = phy_prop.oxygen_inlet_concentration
            self.n2o2ratio = (1. - self.o2_con_in) / self.o2_con_in
            self.spec_num = 3
            self.val_num = 4.
            self.mol_mass = np.array([32., 18., 28.]) * 1.e-3
        else:
            self.channel = ch.Channel(ch_dict.anode_channel)
            self.h2_con_in = phy_prop.hydrogen_inlet_concentration
            self.n2h2ratio = (1. - self.h2_con_in) / self.h2_con_in
            self.spec_num = 3
            self.val_num = 2.
            self.mol_mass = np.array([2., 18., 28.]) * 1.e-3
        self.type = dict_hc['type']  # Anode is False; Cathode is True
        self.th_gdl = dict_hc['th_gdl']
        self.th_plate = dict_hc['th_plate']
        self.th_cat = dict_hc['th_electrode']
        self.th_gde = self.th_gdl + self.th_cat
        self.vol_ex_cd = dict_hc['vol_ex_cd']
        self.prot_con = dict_hc['prot_con']
        self.diff_coef_cat = dict_hc['diff_coef_cat']
        self.diff_coef_gdl = dict_hc['diff_coef_gdl']
        self.tafel_slope = dict_hc['tafel_slope']
        # Scalar variables
        self.break_program = False
        self.pem_type = True
        self.stoi = 0.
        self.i_theta = np.sqrt(2. * self.vol_ex_cd * self.prot_con
                               * self.tafel_slope)
        self.dh = 4. * self.channel.cross_area / self.channel.extent
        self.index_cat = g_par.dict_case['nodes'] - 1
        self.p_drop_bends = 0.
        self.i_star = self.prot_con * self.tafel_slope / self.th_cat
        nodes = g_par.dict_case['nodes']
        # Arrays
        self.j = np.full(nodes-1, 0.)
        self.g_full = np.full(nodes, 0.)
        self.cp_full = np.full(nodes, 0.)
        self.Re = np.full(nodes, 0.)
        self.zero_ele = np.full(nodes - 1, 1.e-50)
        self.zero_node = np.full(nodes, 1.e-50)
        self.w = np.full(nodes, 0.)
        self.p = np.full(nodes, self.channel.p_in)
        self.gamma = np.full(nodes, 0.)
        self.humidity = np.full(nodes, 0.)
        self.free_water = np.full(nodes, 0.)
        self.i_ca = np.full(nodes - 1, g_par.dict_case['tar_cd'])
        self.u = np.full(nodes, 0.)
        self.node_fwd = m_d.fwd_mat_reac_flow(nodes - 1)
        self.node_bwd = m_d.bwd_mat_reac_flow(nodes - 1)
        self.m_flow = np.full(nodes, 0.)
        self.m_reac_flow = np.full(nodes, 0.)
        self.m_liq_water = np.full(nodes, 0.)
        self.m_vap_water_flow = np.full(nodes, 0.)
        self.m_reac_flow_delta = np.full(nodes, 0.)
        self.m_vap_water_flow_delta = np.full(nodes, 0.)
        self.m_full_flow = np.full(nodes, 0.)
        self.act_ov = np.full(nodes - 1, 0.)
        self.gdl_diff_los = np.full(nodes - 1, 0.)
        self.cat_diff_los = np.full(nodes - 1, 0.)
        self.v_los = np.full(nodes-1, 0.)
        self.i_lim = np.full(nodes - 1, 0.)
        self.i_hat = np.full(nodes - 1, 0.)
        self.beta = np.full(nodes - 1, 0.)
        self.var = np.full(nodes - 1, 0.)
        self.i_square = np.full(nodes - 1, 0.)
        self.q_gas = np.full(nodes, 0.)
        self.perm = np.full(nodes, 0.)
        self.mol_flow = np.full((self.spec_num, nodes), 0.)
        self.gas_con = np.full((self.spec_num, nodes), 0.)
        self.gas_con_ele = np.full((nodes-1), 0.)
        self.mixture = [None] * nodes
        self.t_gas = np.full(nodes, self.channel.t_in)
        self.rho = np.full(nodes, 1.)
        self.visc_mix = np.full(nodes, 1.e-5)
        self.Nu = np.full(nodes, 0.)
        self.mol_f = np.full((self.spec_num, nodes), 0.)
        self.mass_f = np.full((self.spec_num, nodes), 0.)
        self.r_mix = np.full((self.spec_num, nodes), 0.)
        self.r_el = np.full(self.spec_num, 0.)
        self.cp = np.full((self.spec_num, nodes), 0.)
        self.lambda_gas = np.full((self.spec_num, nodes), 0.)
        self.visc = np.full((self.spec_num, nodes), 0.)
        self.t_gas_ele = np.full((nodes-1), 0.)
        self.cp_ele = np.full((nodes-1), 0.)
        self.cp_mix = np.full(nodes, 0.)
        self.ht_coef = np.full(nodes, 0.)
        self.k_ht_coef_ca = np.full(nodes, 0.)
        self.cp_mix_ele = np.full((nodes-1), 0.)
        self.lambda_mix = np.full(nodes, 0.)
        self.Pr = np.full(nodes, 0.)
        if self.type is True:
            self.t_nodes = np.full((3, nodes), 0.)
        else:
            self.t_nodes = np.full((2, nodes-1), 0.)
            # 1 catalyst layer, 2 channel layer, 3 coolant plate layer
        for q, item in enumerate(self.mol_mass):
            self.r_el[q] = g_par.dict_uni['R'] / item

    def update(self):
        self.calc_t_gas_e()
        self.calc_reac_flow()
        self.calc_water_flow()
        self.calc_con()
        self.calc_voltage_losses_parameter()
        if self.break_program is False:
            self.calc_activation_losses()
            self.calc_transport_losses_catalyst_layer()
            self.calc_transport_losses_diffusion_layer()
            self.calc_electrode_losses()
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

    def set_i(self, i_ca):
        self.i_ca = i_ca

    def set_j(self, j):
        self.j = j[1:]

    def set_stoi(self, stoi):
        self.stoi = stoi

    def set_pem_type(self, pem_type):
        self.pem_type = pem_type

    def set_t(self, var):
        var = g_func.calc_nodes_2_d(np.array(var))
        if self.type is True:
            self.t = np.array([var[0], var[1], var[2]])
        else:
            self.t = np.array([var[0], var[1]])

    def calc_t_gas_e(self):
        self.t_gas_ele = g_func.calc_elements_1_d(self.t_gas)

    def calc_reac_flow(self):
        """Calculates the reactant molar flow [mol/s]

            Access too:
            -self.stoi
            -g_par.dict_case['tar_cd']
            -g_par.dict_case['F']
            -self.channel.plane
            -self.channel.plane_dx
            -self.i_ca
            -self.node_fwd
            -self.node_bwd
            -self.val_num
            -self.zero_node

            Manipulate:
            -self.mol_flow
        """
        f = g_par.dict_uni['F']
        var1 = self.stoi * g_par.dict_case['tar_cd'] \
            * self.channel.plane / (self.val_num * f)
        if self.type is True:
            self.mol_flow[0, 0] = var1
            self.mol_flow[0, 1:] = var1 - np.matmul(self.node_fwd, self.i_ca)\
                * self.channel.plane_dx / (self.val_num * f)

        else:
            self.mol_flow[0, -1] = var1
            self.mol_flow[0, :-1] = var1 - np.matmul(self.node_bwd, self.i_ca) \
                * self.channel.plane_dx \
                / (self.val_num * f)
        self.mol_flow[0] = np.maximum(self.mol_flow[0], self.zero_node)

    def calc_water_flow(self):
        """" Calculates the water and nitrogen molar flows [mol/s]

            Access to:
            -self.channel.t_in
            -self.channel.plane_dx
            -self.mol_flow
            -self.n2o2ratio
            -self.channel.humidity_in
            -self.channel.p_in
            -self.val_num
            -g_par.dict_uni['F']
            -self.node_fwd
            -self.i_ca
            -self.node_bwd
            -self.j
            -self.p
            -self.pem_type
            -self.zero

            Manipulate:
            -self.mol_flow
            -self.index_cat
        """
        sat_p = p_sat.water.calc_p_sat(self.channel.t_in)
        plane_dx = self.channel.plane_dx
        b = 0.
        if self.type is True:
            q_0_water = self.mol_flow[0][0] \
                        * (1. + self.n2o2ratio) \
                        * sat_p \
                        * self.channel.humidity_in \
                        / (self.channel.p_in
                           - self.channel.humidity_in
                           * sat_p)
            a = plane_dx \
                / (self.val_num * g_par.dict_uni['F'] * 0.5) \
                * np.matmul(self.node_fwd, self.i_ca)
            # production
            if self.pem_type is False:
                b = plane_dx \
                    * np.matmul(self.node_fwd, self.j)
                # crossover
            self.mol_flow[1, 0] = q_0_water
            self.mol_flow[1, 1:] = a + b + q_0_water
            self.mol_flow[2] = np.full(g_par.dict_case['nodes'],
                                       self.mol_flow[0][0] * self.n2o2ratio)
        else:
            q_0_water = self.mol_flow[0][0] \
                        * (1. + self.n2h2ratio) \
                        * sat_p \
                        * self.channel.humidity_in \
                        / (self.channel.p_in
                           - self.channel.humidity_in
                           * sat_p)
            if self.pem_type is False:
                b = plane_dx \
                    * np.matmul(self.node_bwd, self.j)
            self.mol_flow[1, -1] = q_0_water
            self.mol_flow[1, :-1] = b + q_0_water
            self.mol_flow[2] = np.full(g_par.dict_case['nodes'],
                                       self.mol_flow[0][-1] * self.n2h2ratio)

        self.mol_flow[1] = np.choose(self.mol_flow[0] > 1.e-50,
                                     [self.zero_node, self.mol_flow[1]])
        if self.type is True:
            for w in range(1, g_par.dict_case['nodes']):
                if self.mol_flow[1, w] < 1.e-49:
                    self.index_cat = w - 1
                    self.mol_flow[1, -self.index_cat - 1:] = \
                        self.mol_flow[1, self.index_cat]
                    break

    def calc_pressure_drop_bends(self):
        """Calculates the channel pressure
           drop through the channel bends for each element

            Access to:
            -self.channel.bend_fri_fac
            -self.rho
            -self.u
            -self.channel.bend_num
            -gpar.dict_case['nodes']

            Manipulate:
            -self.p_drop_bends
        """
        self.p_drop_bends = self.channel.bend_fri_fac * np.average(self.rho) \
            * np.average(self.u) ** 2. * self.channel.bend_num \
            / (g_par.dict_case['nodes'] - 1) * .5

    def calc_pressure(self):
        """ Calculates the total channel pressure for each element

            Access to:
            -self.channel.p_in
            -self.u
            -self.Re
            -self.node_fwd
            -self.dh
            -self.channel.d_x
            -self.p_drop_bends

            Manipulate:
            -self.p
        """
        p_in = self.channel.p_in
        rho_ele = g_func.calc_elements_1_d(self.rho)
        u_ele = g_func.calc_elements_1_d(self.u)
        Re_ele = g_func.calc_elements_1_d(self.Re)
        if self.type is True:
            mat = self.node_fwd
            self.p[0] = p_in
            self.p[1:] = p_in + 32. / self.dh\
                * np.matmul(-mat, rho_ele * u_ele ** 2. / Re_ele)\
                * self.channel.d_x - self.p_drop_bends
        else:
            mat = self.node_bwd
            self.p[-1] = p_in
            self.p[:-1] = p_in + 32. / self.dh \
                * np.matmul(-mat, rho_ele * u_ele ** 2. / Re_ele) \
                * self.channel.d_x - self.p_drop_bends

    def calc_con(self):
        for w in range(g_par.dict_case['nodes']):
            id_lw = self.p[w] / (g_par.dict_uni['R'] * self.t[1, w])
            var4 = np.sum(self.mol_flow[:, w])
            var2 = self.mol_flow[1][w] / var4
            self.gas_con[1][w] = id_lw * var2
            a = p_sat.water.calc_p_sat(self.t[1, w])
            e = g_par.dict_uni['R'] * self.t_gas[w]
            if self.gas_con[1][w] >= a / e:  # saturated
                b = self.mol_flow[0][w] + self.mol_flow[2][w]
                c = self.mol_flow[0][w] / b
                d = self.mol_flow[2][w] / b
                self.gas_con[0][w] = (self.p[w] - a) / e * c
                self.gas_con[2][w] = (self.p[w] - a) / e * d
                self.gas_con[1][w] = a / e
            else:  # not saturated
                var5 = id_lw / var4
                self.gas_con[0][w] = var5 * self.mol_flow[0][w]
                self.gas_con[2][w] = var5 * self.mol_flow[2][w]
        self.gas_con_ele = g_func.calc_elements_1_d(self.gas_con[0])

    def calc_mass_fraction(self):
        temp_var_1 = []
        for q, item in enumerate(self.gas_con):
            temp_var_1.append(item * self.mol_mass[q])
        for q in range(len(self.gas_con)):
            self.mass_f[q] = temp_var_1[q] / sum(temp_var_1)
        # print(self.mass_f, 'mass_f')

    def calc_mol_fraction(self):
        for q, item in enumerate(self.gas_con):
            self.mol_f[q] = item / sum(self.gas_con)
        # print(self.mol_f, 'mol_f')

    def calc_gas_properties(self):
        if self.type is True:
            self.cp[0] = g_fit.oxygen.calc_cp(self.t_gas)
            self.lambda_gas[0] = g_fit.oxygen.calc_lambda(self.t_gas, self.p)
            self.visc[0] = g_fit.oxygen.calc_visc(self.t_gas)
        else:
            self.cp[0] = g_fit.hydrogen.calc_cp(self.t_gas)
            self.lambda_gas[0] = g_fit.hydrogen.calc_lambda(self.t_gas, self.p)
            self.visc[0] = g_fit.hydrogen.calc_visc(self.t_gas)
        self.cp[1] = g_fit.water.calc_cp(self.t_gas)
        self.cp[2] = g_fit.nitrogen.calc_cp(self.t_gas)
        self.lambda_gas[1] = g_fit.water.calc_lambda(self.t_gas, self.p)
        self.lambda_gas[2] = g_fit.nitrogen.calc_lambda(self.t_gas, self.p)
        self.visc[1] = g_fit.water.calc_visc(self.t_gas)
        self.visc[2] = g_fit.nitrogen.calc_visc(self.t_gas)
        self.cp_ele = g_func.calc_elements_1_d(self.cp[0])

    def calc_gas_mix_properties(self):
        temp1, temp2 = [], []
        for q in range(self.spec_num):
            temp1.append(self.mass_f[q] * self.r_el[q])
            temp2.append(self.mass_f[q] * self.cp[q])
        self.r_mix = sum(temp1)
        self.cp_mix = sum(temp2)
        self.cp_mix_ele = g_func.calc_elements_1_d(self.cp_mix)
        self.visc_mix = g_func.calc_visc_mix(self.visc,
                                             self.mol_f,
                                             self.mol_mass)
        self.lambda_mix = g_func.calc_lambda_mix(self.lambda_gas, self.mol_f,
                                                 self.visc, self.mol_mass)
        self.rho = g_func.calc_rho(self.p, self.r_mix, self.t_gas)
        self.Pr = self.visc_mix * self.cp_mix / self.lambda_mix
        # print(self.mass_f,'mass_f')
        # print(self.mol_f, 'mol_f')
        # print(self.r_mix, 'r_mix')
        # print(self.cp_mix, 'cp_mix')
        # print(self.visc_mix, 'visc_mix')
        # print(self.lambda_mix, 'lambda_mix')
        # print(self.rho, 'rho')
        # print(self.Pr, 'Pr')

    def calc_flow_velocity(self):
        self.u = self.q_gas * g_par.dict_uni['R'] * self.t_gas \
                 / (self.p * self.channel.cross_area)
        # print(self.u, 'u')

    def calc_mass_flow(self):
        self.m_flow = self.u * self.rho * self.channel.cross_area
        self.m_reac_flow = self.mol_flow[0] * self.mol_mass[0]
        self.m_liq_water = self.w * self.mol_mass[1]
        self.m_vap_water_flow = (self.mol_flow[1] - self.w) * self.mol_mass[1]
        self.m_reac_flow_delta = abs(g_func.calc_dif(self.m_reac_flow))
        self.m_vap_water_flow_delta = abs(g_func.calc_dif(
            self.m_vap_water_flow))
        self.m_full_flow = self.m_flow + self.m_liq_water
        # print(self.m_flow)
        # print(self.m_reac_flow)
        # print(self.m_liq_water)
        # print(self.m_vap_water_flow)
        # print(self.m_reac_flow_delta)
        # print(self.m_reac_flow_delta)
        # print(self.m_full_flow)

    def calc_m_mix_properties(self):
        self.cp_full = (self.m_flow * self.cp_mix
                        + self.m_liq_water * g_par.dict_uni['cp_liq']) \
                       / self.m_full_flow
        self.g_full = self.m_full_flow * self.cp_full
        # print(self.cp_full, 'cp_full')
        # print(self.g_full, 'g_full')

    def calc_re(self):
        self.Re = g_func.calc_Re(self.rho, self.u, self.dh, self.visc_mix)
        # print(self.Re, 'Re')

    def calc_nu(self):
        self.Nu = np.full(g_par.dict_case['nodes'], 3.66)
        # print(self.Nu, 'Nu')

    def calc_heat_transfer_coef(self):
        self.ht_coef = self.lambda_mix * self.Nu / self.dh
        self.k_ht_coef_ca = self.ht_coef * np.pi * self.channel.d_x * self.dh

    def calc_fluid_water(self):  # fluid water in the channel
        self.w = self.mol_flow[1] - self.gas_con[1] / self.gas_con[0] \
                 * self.mol_flow[0]
        # print(self.w, 'w')

    def calc_cond_rates(self):  # condensation rates of vapour in the channel
        if self.type is True:
            self.gamma = g_func.calc_nodes_1_d(np.ediff1d(self.w))
        else:
            self.gamma = -g_func.calc_nodes_1_d(np.ediff1d(self.w))
        # print(self.gamma, 'gamma')

    def calc_rel_humidity(self):  # relative humidity in the channel
        self.humidity = self.gas_con[1] * g_par.dict_uni['R'] \
                        * self.t_gas / p_sat.water.calc_p_sat(self.t_gas)
        # print(self.humidity, 'humidity')

    def calc_free_water(self):  # membrane free water content
        a = 0.043 + 17.81 * self.humidity
        b = -39.85 * self.humidity ** 2 + 36. * self.humidity ** 3
        self.free_water = a + b
        # print(self.free_water, 'free_water')

    def sum_flows(self):
        self.q_gas = sum(self.mol_flow) - self.w
        # print(self.q_gas, 'q_gas')

    def calc_voltage_losses_parameter(self):
        self.i_lim = 4. * g_par.dict_uni['F'] * self.gas_con[0, :-1] \
                     * self.diff_coef_gdl / self.th_gdl
        self.i_hat = self.i_ca / self.i_star
        short_save = np.sqrt(2. * self.i_hat)
        self.beta =\
            short_save / (1. + np.sqrt(1.12 * self.i_hat) * np.exp(short_save))\
            + np.pi * self.i_hat / (2. + self.i_hat)
        self.var =\
            1. - self.i_ca\
            / (self.i_lim * self.gas_con_ele / self.gas_con[0, :-1])
        self.i_square = self.i_ca ** 2.
        # print(self.i_lim)
        # print(self.i_hat)
        # print(short_save)
        # print(self.beta)
        # print(self.var)
        # print(self.i_square)

    def calc_activation_losses(self):
        self.act_ov = self.tafel_slope \
                      * np.arcsinh((self.i_ca / self.i_theta) ** 2.
                                   / (2. * (self.gas_con_ele
                                            / self.gas_con[0, :-1])
                                      * (1. - np.exp(-self.i_ca /
                                                     (2. * self.i_star)))))

    def calc_transport_losses_catalyst_layer(self):
        self.cat_diff_los = ((self.prot_con * self.tafel_slope ** 2.)
                             / (4. * g_par.dict_uni['F'] * self.diff_coef_cat
                                * self.gas_con_ele)
                             * (self.i_ca / self.i_star
                                - np.log10(1. + self.i_square
                                           / (self.i_star ** 2.
                                              * self.beta ** 2.)))) / self.var

    def calc_transport_losses_diffusion_layer(self):
        self.gdl_diff_los = -self.tafel_slope * np.log10(self.var)
        nan_list = np.isnan(self.gdl_diff_los)
        bol = nan_list.any()
        if bol == True:
            self.gdl_diff_los[np.argwhere(nan_list)[0, 0]:] = 1.e20

    def calc_electrode_losses(self):
        self.v_los = self.act_ov + self.cat_diff_los + self.gdl_diff_los
