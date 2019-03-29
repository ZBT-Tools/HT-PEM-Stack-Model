import warnings
import system.global_functions as g_func
import data.water_properties as p_sat
import data.gas_properties as g_fit
import numpy as np
import data.global_parameters as g_par
import system.channel as ch
import data.channel_dict as ch_dict
import input.physical_properties as phy_prop


warnings.filterwarnings("ignore")


class HalfCell:

    def __init__(self, dict_hc):
        # number of nodes along the channel
        nodes = g_par.dict_case['nodes']

        self.species_list = dict_hc['species']
        # check if the object is an anode or a cathode
        # catalyst layer specific handover
        if dict_hc['cl_type'] is True:
            # object cathode gas channel
            self.channel = ch.Channel(ch_dict.dict_cathode_channel)
            # nitrogen to oxygen molar flow ratio
            self.n2o2ratio = 1. / phy_prop.oxygen_inlet_concentration - 1.
            # number of species of the gas mixture
            self.spec_num = 3
            # electrical charge number
            self.val_num = 4.
            # molar mass
            self.mol_mass = np.array([32., 18., 28.]) * 1.e-3
        else:
            # object anode gas channel
            self.channel = ch.Channel(ch_dict.dict_anode_channel)
            # volumetric hydrogen to nitrogen ratio
            self.n2h2ratio = 1. / phy_prop.hydrogen_inlet_concentration - 1.
            # number of species of the gas mixture
            self.spec_num = 3
            # electrical charge number
            self.val_num = 2.
            # molar mass
            self.mol_mass = np.array([2., 18., 28.]) * 1.e-3
        # set anode or cathode Anode is false; Cathode is True
        self.cl_type = dict_hc['cl_type']

        """geometry"""
        # number of channels in each cell
        self.channel_numb = dict_hc['channel_numb']
        # cell width [m]
        self.cell_width = dict_hc['cell_width']
        # cell length [m]
        self.cell_length = dict_hc['cell_length']
        # thickness of the gas diffusion layer [m]
        self.th_gdl = dict_hc['th_gdl']
        # thickness of the bipolar plate [m]
        self.th_bpp = dict_hc['th_bpp']
        # thickness of the catalyst layer [m]
        self.th_cl = dict_hc['th_cl']
        # thickness gas diffusion electrode (gdl + cl) [m]
        self.th_gde = self.th_gdl + self.th_cl

        """voltage loss parameter, (Kulikovsky, 2013)"""
        # exchange current density [A/m³]
        self.vol_ex_cd = dict_hc['vol_ex_cd']
        # proton conductivity of the catalyst layer [1/(Ohm * m)]
        self.prot_con_cl = dict_hc['prot_con_cl']
        # diffusion coefficient of the reactant in catalyst layer [m²/s]
        self.diff_coeff_cl = dict_hc['diff_coeff_cl']
        # diffusion coefficient of the reactant in gas diffusion layer [m²/s]
        self.diff_coeff_gdl = dict_hc['diff_coeff_gdl']
        # tafel slope of the electrode [V]
        self.tafel_slope = dict_hc['tafel_slope']
        # Kulikovsky variable
        self.i_sigma = np.sqrt(2. * self.vol_ex_cd * self.prot_con_cl
                               * self.tafel_slope)
        self.i_ca_char = self.prot_con_cl * self.tafel_slope / self.th_cl
        # voltage loss along the channel [V]
        self.v_loss = np.zeros(nodes-1)
        # activation voltage loss [V]
        self.act_loss = []
        # catalyst layer diffusion voltage loss [V]
        self.cl_diff_loss = []
        # gas diffusion layer loss [V]
        self.gdl_diff_loss = []
        # Kulikovsky variable
        self.beta = np.zeros(nodes - 1)
        # Kulikovsky variable
        self.var = np.zeros(nodes - 1)
        # current density² [A²/m^4]
        self.i_ca_square = np.zeros(nodes - 1)

        """general parameter"""
        # factor active area with ribs / active channel area
        area_fac = self.cell_length * self.cell_width\
            / (self.channel.active_area * self.channel_numb)
        # active area belonging to the channel plan area dx [m²]
        self.active_area_dx_ch = area_fac * self.channel.active_area\
            / g_par.dict_case['elements']
        # active area belonging to the channel plan area [m²]
        self.active_area_ch = area_fac * self.channel.active_area
        # boolean to hint if the cell voltage runs below zero
        self.break_program = False
        # if HT-PEMFC True; if NT-PEMFC False
        self.ht_pem = True
        # stoichiometry of the reactant at the channel inlet
        self.stoi = None
        # cross water flux through the membrane [mol/s]
        self.w_cross_flow = np.zeros(nodes - 1)
        # heat capacity flow of the species mixture including fluid water [W/K]
        self.g_fluid = np.zeros(nodes)
        # heat capacity of the species mixture including fluid water [J/(K mol)]
        self.cp_fluid = np.zeros(nodes)
        # molar liquid water flux [mol/s]
        self.liq_w_flow = np.zeros(nodes)
        # channel pressure [Pa]
        self.p = np.full(nodes, self.channel.p_out)
        # condensation rate of water [mol/s]
        self.cond_rate = np.zeros(nodes)
        # gas mixture humidity
        self.humidity = np.zeros(nodes)
        # current density [A/m²]
        self.i_ca = np.full(nodes - 1, g_par.dict_case['tar_cd'])
        # channel velocity [m/s]
        self.u = np.zeros(nodes)
        # forward matrix
        self.fwd_mat = np.tril(np.full((nodes - 1, nodes - 1), 1.))
        # backward matrix
        self.bwd_mat = np.triu(np.full((nodes - 1, nodes - 1), 1.))
        # mass flow of the gas mixture [kg/s]
        self.m_flow_gas = np.zeros(nodes)
        # reactant mass flow [kg/s]
        self.m_flow_reac = np.zeros(nodes)
        # liquid water mass flow [kg/s]
        self.m_flow_liq_w = np.zeros(nodes)
        # steam mass flow [kg/s]
        self.m_flow_vap_w = np.zeros(nodes)
        # fluid mass flow [kg/s]
        self.m_flow_fluid = np.zeros(nodes)
        # molar flux of the gas phase [mol/s]
        self.q_gas = np.zeros(nodes)
        # molar flow of each species, 0: Reactant, 1: Water 2: Nitrogen [mol/s]
        self.mol_flow = np.full((self.spec_num, nodes), 0.)
        # molar concentration of each species,
        # 0: Reactant, 1: Water 2: Nitrogen [mol/s]
        self.gas_con = np.full((self.spec_num, nodes), 0.)
        # element based molar concentration of the reactant [mol/m³]
        self.gas_con_ele = np.full((nodes-1), 0.)
        # temperature of the fluid in the channel [K]
        self.temp_fluid = np.full(nodes, self.channel.temp_in)
        # density of the gas phase [kg/m³]
        self.rho_gas = np.full(nodes, 1.)
        # dynamic viscosity of the gas phase [kg/s/m]
        self.visc_gas = np.full(nodes, 1.e-5)
        # nusselt number
        self.Nu = 3.66
        # molar fraction of all species in the gas phase
        self.mol_f = np.full((self.spec_num, nodes), 0.)
        # mass fraction of all species in the gas phase
        self.mass_f = np.full((self.spec_num, nodes), 0.)
        # gas constant of the gas phase [J/mol/K]
        self.r_gas = np.full((self.spec_num, nodes), 0.)
        # gas constant of all species [J/mol/K]
        self.r_species = np.full(self.spec_num, 0.)
        # heat capacity of all species in the gas phase [J/mol/K]
        self.cp = np.full((self.spec_num, nodes), 0.)
        # heat conductivity of the species in the gas phase [W/m/K]
        self.lambdas = np.full((self.spec_num, nodes), 0.)
        # viscosity of all species in the gas phase [kg/s/m]
        self.visc = np.full((self.spec_num, nodes), 0.)
        # heat capacity of the gas phase [J/kg/K]
        self.cp_gas = np.zeros(nodes)
        # heat conductivity between the gas phase and the channel [W/K]
        self.k_ht_coef_ca = np.zeros(nodes)
        # heat conductivity of the gas phase [W/ K^-1]
        self.lambda_gas = np.zeros(nodes)
        # prandtl number of the gas phase
        self.Pr = np.zeros(nodes)
        for q, item in enumerate(self.mol_mass):
            self.r_species[q] = g_par.dict_uni['R'] / item

    def update(self):
        """
        This function coordinates the program sequence
        """
        self.calc_mass_balance()
        self.calc_voltage_losses_parameter()
        if self.break_program is False:
            self.update_voltage_loss()
            self.calc_liquid_water_flow()
            self.sum_flows()
            self.calc_cond_rates()
            self.calc_mass_fraction()
            self.calc_mol_fraction()
            self.calc_species_properties()
            self.calc_gas_properties()
            self.calc_rel_humidity()
            self.calc_flow_velocity()
            self.calc_mass_flow()
            self.calc_fluid_flow_properties()
            self.calc_heat_transfer_coef()
            self.calc_pressure()

    def calc_mass_balance(self):
        self.calc_reac_flow()
        self.calc_water_flow()
        self.calc_con()

    def update_voltage_loss(self):
        self.calc_electrode_loss()

    def set_current_density(self, i_ca):
        """
        This function sets the current density.
        The current density can be obtained
        from the electrical coupling.

            Manipulate:
            - self.i_ca, scalar
        """

        self.i_ca = i_ca

    def set_water_cross_flux(self, j):
        """
        This function sets the water cross flux through the membran.
        The water cruss flux can be obtained from the cell level.

            Manipulate:
            - self.water_cross_flux, 1-D-array, [cell elements]
        """

        self.w_cross_flow = j

    def set_stoichiometry(self, stoi):
        """
        This function sets the reactant stoichiometry.
        The stoichiometry can be obtained from the manifold methods.

            Manipulate:
            - self.stoi, scalar
        """

        self.stoi = stoi

    def set_pem_type(self, pem_type):
        """
        This function defines the PEM-Type.

            Manipulate:
            - self.pem_type, boolean
        """

        self.ht_pem = pem_type

    def set_layer_temperature(self, var):
        """
        This function sets the layer Temperatures,
        they can be obtained from the temperature system.

            Manipulate:
            - self.temp, 2-D-array, [layer][elements]
        """

        var = g_func.calc_nodes_2_d(np.array(var))
        if self.cl_type is True:
            self.temp = np.array([var[0], var[1], var[2]])
        else:
            self.temp = np.array([var[0], var[1]])

    def calc_reac_flow(self):
        """
        Calculates the reactant molar flow [mol/s]

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
            * self.active_area_ch / (self.val_num * f)
        if self.cl_type is True:
            self.mol_flow[0, 0] = var1
            self.mol_flow[0, 1:] = var1 - np.matmul(self.fwd_mat, self.i_ca) \
                * self.active_area_dx_ch / (self.val_num * f)

        else:
            self.mol_flow[0, -1] = var1
            self.mol_flow[0, :-1] = var1 - np.matmul(self.bwd_mat, self.i_ca) \
                * self.active_area_dx_ch \
                / (self.val_num * f)
        self.mol_flow[0] = np.maximum(self.mol_flow[0],
                                      np.zeros(g_par.dict_case['elements'] + 1))

    def calc_water_flow(self):
        """"
        Calculates the water and nitrogen molar flows [mol/s]

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

        sat_p = p_sat.water.calc_p_sat(self.channel.temp_in)
        plane_dx = self.active_area_dx_ch
        b = 0.
        if self.cl_type is True:
            q_0_water = self.mol_flow[0][0] \
                        * (1. + self.n2o2ratio) \
                        * sat_p \
                        * self.channel.humidity_in \
                        / (self.channel.p_out
                           - self.channel.humidity_in
                           * sat_p)
            a = plane_dx \
                / (self.val_num * g_par.dict_uni['F'] * 0.5) \
                * np.matmul(self.fwd_mat, self.i_ca)
            # production
            if self.ht_pem is False:
                b = plane_dx \
                    * np.matmul(self.fwd_mat, self.w_cross_flow)
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
                        / (self.channel.p_out
                           - self.channel.humidity_in
                           * sat_p)
            if self.ht_pem is False:
                b = plane_dx \
                    * np.matmul(-self.bwd_mat, self.w_cross_flow)
            self.mol_flow[1, -1] = q_0_water
            self.mol_flow[1, :-1] = b + q_0_water
            self.mol_flow[2] = np.full(g_par.dict_case['nodes'],
                                       self.mol_flow[0][-1] * self.n2h2ratio)
        self.mol_flow[1] = np.maximum(self.mol_flow[1], 0.)
        self.mol_flow[1] = np.choose(self.mol_flow[0] > 1.e-50,
                                     [np.zeros(g_par.dict_case['nodes']),
                                      self.mol_flow[1]])
        if self.cl_type is True:
            for w in range(1, g_par.dict_case['nodes']):
                if self.mol_flow[1, w] < 1.e-49:
                    index_cat = w - 1
                    self.mol_flow[1, -index_cat - 1:] = \
                        self.mol_flow[1, index_cat]
                    break

    def calc_pressure_drop_bends(self):
        """
        Calculates the channel pressure drop
        through the channel bends for each element.

            Access to:
            -self.channel.bend_fri_fac
            -self.rho
            -self.u
            -self.channel.bend_num
            -gpar.dict_case['nodes']

            Manipulate:
            -self.p_drop_bends
        """

        p_drop_bends = self.channel.bend_fri_fac\
            * np.average(self.rho_gas) * np.average(self.u) ** 2.\
            * self.channel.bend_numb / (g_par.dict_case['nodes'] - 1) * .5
        return p_drop_bends

    def calc_pressure(self):
        """
        Calculates the total channel pressure for each element.

            Access to:
            -self.channel.p_in
            -self.u
            -self.Re
            -self.rho_gas
            -self.fwd_mat
            -self.bwd_mat
            -self.channel.d_h
            -self.channel.dx
            -self.p_drop_bends

            Manipulate:
            -self.p
        """
        p_drop_bends = self.calc_pressure_drop_bends()
        p_out = self.channel.p_out
        rho_ele = g_func.calc_elements_1_d(self.rho_gas)
        u_ele = g_func.calc_elements_1_d(self.u)
        Re_ele = g_func.calc_elements_1_d(self.calc_Re())
        if self.cl_type is True:
            mat = self.bwd_mat
            self.p[-1] = p_out
            self.p[:-1] = p_out + 32. / self.channel.d_h \
                * np.matmul(mat, rho_ele * u_ele ** 2. / Re_ele) \
                * self.channel.dx\
                + np.linspace(p_drop_bends * (g_par.dict_case['nodes']),0,
                              g_par.dict_case['nodes']-1)
        else:
            mat = self.fwd_mat
            self.p[0] = p_out
            self.p[1:] = p_out + 32. / self.channel.d_h \
                * np.matmul(mat, rho_ele * u_ele ** 2. / Re_ele) \
                * self.channel.dx\
                + np.linspace(0, p_drop_bends * (g_par.dict_case['nodes']),
                              g_par.dict_case['nodes']-1)

    def calc_con(self):
        """
        Calculates the gas phase molar concentrations.

            Access to:
            -self.p
            -self.mol_flow
            -self.temp_fluid
            -g.par.dict_uni['R']

            Manipulate:
            -self.gas_con
        """

        for w in range(g_par.dict_case['nodes']):
            id_lw = self.p[w] / (g_par.dict_uni['R'] * self.temp_fluid[w])
            var4 = np.sum(self.mol_flow[:, w])
            var2 = self.mol_flow[1][w] / var4
            self.gas_con[1][w] = id_lw * var2
            a = p_sat.water.calc_p_sat(self.temp_fluid[w])
            e = g_par.dict_uni['R'] * self.temp_fluid[w]
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
        """
        Calculates the gas phase mass fractions.

            Access to:
            -self.gas_con
            -self.mol_mass

            Manipulate:
            -self.mass_f
        """

        temp_var_1 = []
        for q, item in enumerate(self.gas_con):
            temp_var_1.append(item * self.mol_mass[q])
        for q in range(len(self.gas_con)):
            self.mass_f[q] = temp_var_1[q] / sum(temp_var_1)

    def calc_mol_fraction(self):
        """
        Calculates the gas phase molar fractions.

            Access to:
            -self.gas_con

            Manipulate:
            -self.mol_f
        """

        for q, item in enumerate(self.gas_con):
            self.mol_f[q] = item / sum(self.gas_con)

    def calc_species_properties(self):
        """
        Calculates the properties of the species in the gas phase

            Access to:
            -self.temp_fluid
            -self.p
            -self.cl_type
            -g_fit.object.methods

            Manipulate:
            -self.cp
            -self.lambdas
            -self.visc
            -self.cp_ele
        """

        if self.cl_type is True:
            self.cp[0] = g_fit.oxygen.calc_cp(self.temp_fluid)
            self.lambdas[0] = g_fit.oxygen.calc_lambda(self.temp_fluid, self.p)
            self.visc[0] = g_fit.oxygen.calc_visc(self.temp_fluid)
        else:
            self.cp[0] = g_fit.hydrogen.calc_cp(self.temp_fluid)
            self.lambdas[0] = g_fit.hydrogen.calc_lambda(self.temp_fluid,
                                                         self.p)
            self.visc[0] = g_fit.hydrogen.calc_visc(self.temp_fluid)
        self.cp[1] = g_fit.water.calc_cp(self.temp_fluid)
        self.cp[2] = g_fit.nitrogen.calc_cp(self.temp_fluid)
        self.lambdas[1] = g_fit.water.calc_lambda(self.temp_fluid, self.p)
        self.lambdas[2] = g_fit.nitrogen.calc_lambda(self.temp_fluid, self.p)
        self.visc[1] = g_fit.water.calc_visc(self.temp_fluid)
        self.visc[2] = g_fit.nitrogen.calc_visc(self.temp_fluid)

    def calc_gas_properties(self):
        """
        Calculates the properties of the gas phase

            Access to:
            -self.spec_num
            -self.mass_f
            -self.r_species
            -self.cp
            -self.visc
            -self.lambdas
            -self.mol_f
            -self.temp_fluid

            Manipulate:
            -self.r_gas
            -self.cp_gas
            -self.cp_gas_ele
            -self.cp_ele
            -self.visc_gas
            -self.lambda_gas
            -self.rho_gas
            -self.Pr
        """

        temp1, temp2 = [], []
        for q in range(self.spec_num):
            temp1.append(self.mass_f[q] * self.r_species[q])
            temp2.append(self.mass_f[q] * self.cp[q])
        self.r_gas = sum(temp1)
        self.cp_gas = sum(temp2)
        self.visc_gas = g_func.calc_visc_mix(self.visc,
                                             self.mol_f,
                                             self.mol_mass)
        self.lambda_gas = g_func.calc_lambda_mix(self.lambdas, self.mol_f,
                                                 self.visc, self.mol_mass)
        self.rho_gas = g_func.calc_rho(self.p, self.r_gas, self.temp_fluid)
        self.Pr = self.visc_gas * self.cp_gas / self.lambda_gas

    def calc_flow_velocity(self):
        """
        Calculates the gas phase velocity.
        The gas phase velocity is taken to be the liquid water velocity as well.

            Access to:
            -self.q_gas
            -self.temp_fluid
            -self.p
            -self.channel.cross_area
            -g_par.dict_uni['R']


            Manipulate:
            -self.u
        """

        self.u = self.q_gas * g_par.dict_uni['R'] * self.temp_fluid \
            / (self.p * self.channel.cross_area)

    def calc_mass_flow(self):
        """
        Calculates the relevant mass flows

            Access to:
            -self.u
            -self.rho_gas
            -self.channel_cross_area
            -self.mol_flow
            -self.mol_mass
            -self.l_w_flow

            Manipulate:
            -self.m_flow_gas
            -self.m_flow_reac
            -self.m_flow_liq_w
            -self.m_flow_vap_w
            -self.m_flow_reac_delta
            -self.m_flow_vap_water_delta
            -self.m_flow_fluid
        """

        self.m_flow_gas = self.u * self.rho_gas * self.channel.cross_area
        self.m_flow_reac = self.mol_flow[0] * self.mol_mass[0]
        self.m_flow_liq_w = self.liq_w_flow * self.mol_mass[1]
        self.m_flow_vap_w = (self.mol_flow[1] - self.liq_w_flow)\
                            * self.mol_mass[1]
        self.m_flow_fluid = self.m_flow_gas + self.m_flow_liq_w

    def calc_fluid_flow_properties(self):
        """
        Calculate the fluid flow properties

            Access to:
            -self.m_flow_gas
            -self.cp_gas
            -self.m_flow_liq_w
            -self.m_flow_fluid
            -g_par.dict_uni['cp_liq']
        """

        self.cp_fluid = (self.m_flow_gas * self.cp_gas + self.m_flow_liq_w
                         * g_par.dict_uni['cp_liq']) / self.m_flow_fluid
        self.g_fluid = self.m_flow_fluid * self.cp_fluid

    def calc_Re(self):
        """
        Calculates the reynolds number of the gas flow in the channel

            Access to:
            -self.rho_gas
            -self.u
            -self.visc_gas
            -self.channel.d_h
        """

        Re = g_func.calc_Re(self.rho_gas, self.u,
                            self.channel.d_h, self.visc_gas)
        return Re

    def calc_heat_transfer_coef(self):
        """
        Calculates the convection coefficient between the channel and
        the gas phase.
        Deputy for the convection coefficient between the fluid and the channel.

        Calculates the heat conductivity based on the convection coefficient
        and the element area of the channel.

            Access to:
            -self.lambda_gas
            -self.Nu
            -self.channel.d_h
            -self.channel.dx

            Manipulate:
            -self.ht_coef
            self.k_ht_coef_ca
        """

        ht_coef = self.lambda_gas * self.Nu / self.channel.d_h
        self.k_ht_coef_ca = ht_coef * np.pi * self.channel.dx * self.channel.d_h

    def calc_liquid_water_flow(self):
        """
        Calculates the liquid water flow.

            Access to:
            -self.mol_flow
            -self.gas_con

            Manipulate:
            -self.liq_w_flow
        """

        self.liq_w_flow = self.mol_flow[1]\
            - self.gas_con[1] / self.gas_con[0] * self.mol_flow[0]

    def calc_cond_rates(self):
        """
        Calculates the molar condensation rate of water in the channel.

            Access to:
            -self.cl_type
            -self.liq_w_flow

            Manipulate:
            -self.cond_rate
        """

        if self.cl_type is True:
            self.cond_rate = g_func.calc_nodes_1_d(np.ediff1d(self.liq_w_flow))
        else:
            self.cond_rate = -g_func.calc_nodes_1_d(np.ediff1d(self.liq_w_flow))

    def calc_rel_humidity(self):
        """
        Calculates the relative humidity of the fluid.

            Access to:
            -self.gas_con
            -self.temp_fluid
            -g_par.dict_uni['R']

            Manipulate:
            self.humidity
        """

        self.humidity = self.gas_con[1] * g_par.dict_uni['R'] \
            * self.temp_fluid / p_sat.water.calc_p_sat(self.temp_fluid)

    def sum_flows(self):
        """
        Calculates the molar flow of the gas phase

        Access to:
        -self.mol_flow
        -self.liq_w_flow

        Manipulate:
        -self.q_gas
        """

        self.q_gas = sum(self.mol_flow) - self.liq_w_flow

    def calc_voltage_losses_parameter(self):
        """
        Calculates multiply used supporting parameters
        to calculate the voltage loss according to (Kulikovsky, 2013).

        Access to:
        - self.gas_con
        -self.diff_coeff_gdl
        -self.th_gdl
        -self.i_ca
        -self.gas_con_ele
        -g_par.dict_uni['F']

        Manipulate:
        -self.var
        -self.i_ca_square
        """

        i_lim = 4. * g_par.dict_uni['F'] * self.gas_con[0, :-1] \
            * self.diff_coeff_gdl / self.th_gdl
        self.var = 1.\
            - self.i_ca / (i_lim * self.gas_con_ele / self.gas_con[0, :-1])
        self.i_ca_square = self.i_ca ** 2.

    def calc_activation_loss(self):
        """
        Calculates the activation voltage loss,
        according to (Kulikovsky, 2013).

        Access to:
        -self.tafel_slope
        -self.i_ca
        -self.i_sigma
        -self.gas_con
        -self.gas_con_ele
        -self.i_ca_char

        Manipulate:
        -self.act_ov
        """

        act_loss = self.tafel_slope\
            * np.arcsinh((self.i_ca / self.i_sigma) ** 2.
                         / (2. * (self.gas_con_ele / self.gas_con[0, :-1])
                            * (1. - np.exp(-self.i_ca /
                                           (2. * self.i_ca_char)))))
        return act_loss

    def calc_transport_loss_catalyst_layer(self):
        """
        Calculates the diffusion voltage loss in the catalyst layer
        according to (Kulikovsky, 2013).

        Access to:
        -self.i_ca
        -self.i_ca_char
        -self.prot_con_cl
        -self.tafel_slope
        -self.diff_coeff_cl
        -self.i_ca_char
        -self.gas_con_ele
        -self.i_ca_square
        -self.var
        -g_par.dict_uni['F']

        Manipulate:
        -self.cl_diff_loss
        """

        i_hat = self.i_ca / self.i_ca_char
        short_save = np.sqrt(2. * i_hat)
        beta = short_save / (1. + np.sqrt(1.12 * i_hat) * np.exp(short_save))\
            + np.pi * i_hat / (2. + i_hat)
        cl_diff_loss = \
            ((self.prot_con_cl * self.tafel_slope ** 2.)
             / (4. * g_par.dict_uni['F']
                * self.diff_coeff_cl * self.gas_con_ele)
                * (self.i_ca / self.i_ca_char
                   - np.log10(1. + self.i_ca_square /
                              (self.i_ca_char ** 2. * beta ** 2.)))) / self.var
        return cl_diff_loss

    def calc_transport_loss_diffusion_layer(self):
        """
        Calculates the diffusion voltage loss in the gas diffusion layer
        according to (Kulikovsky, 2013).

        Access to:
        -self.tafel_slope
        -self.var

        Manipulate:
        -self.gdl_diff_loss
        """

        gdl_diff_loss = -self.tafel_slope * np.log10(self.var)
        nan_list = np.isnan(gdl_diff_loss)
        bol = nan_list.any()
        if bol == True:
            gdl_diff_loss[np.argwhere(nan_list)[0, 0]:] = 1.e50
        return gdl_diff_loss

    def calc_electrode_loss(self):
        """
        Calculates the full volatege losses of the electrode

            Access to:
            -self.act_loss
            -self.cl_diff_loss
            -self.gdl_diff_loss

            Manipulate:
            -self.v_loss
        """
        if g_par.dict_case['calc_act_loss'] is False:
            self.act_loss = 0.
        else:
            self.act_loss = self.calc_activation_loss()
        if g_par.dict_case['calc_gdl_diff_loss'] is False:
            self.gdl_diff_loss = 0.
        else:
            self.gdl_diff_loss = self.calc_transport_loss_diffusion_layer()
        if g_par.dict_case['calc_cl_diff_loss'] is False:
            self.cl_diff_loss = 0.
        else:
            self.cl_diff_loss = self.calc_transport_loss_catalyst_layer()

        self.v_loss = self.act_loss + self.cl_diff_loss + self.gdl_diff_loss
