import warnings
import system.global_functions as g_func
import data.water_properties as w_prop
import data.gas_properties as g_fit
import numpy as np
import data.global_parameters as g_par
import system.channel as ch
import data.channel_dict as ch_dict
import input.physical_properties as phy_prop
import sys
import system.interpolation as ip


warnings.filterwarnings("ignore")


class HalfCell:

    def __init__(self, dict_hc, dict_cell):
        self.name = dict_hc['name']
        self.n_nodes = g_par.dict_case['nodes']
        n_nodes = self.n_nodes
        n_ele = n_nodes - 1
        self.n_ele = n_ele
        # discretization in elements and nodes along the x-axis (flow axis)

        self.flow_direction = dict_hc['flow_direction']
        if self.flow_direction not in (-1, 1):
            raise sys.exit('Member variable flow_direction of class HalfCell '
                           'must be either 1 or -1')

        self.id_reac = 0
        self.id_h2o = 1
        self.id_inert = 2
        # check if the object is an anode or a cathode
        # catalyst layer specific handover
        if dict_hc['is_cathode']:
            self.channel = ch.Channel(ch_dict.dict_cathode_channel)
            reac_con_in = phy_prop.oxygen_inlet_concentration
            # volumetric inlet oxygen ratio
            self.inert_reac_ratio = (1. - reac_con_in) / reac_con_in
            # volumetric nitrogen to oxygen ratio
            self.n_species = 3
            # number of  species in the gas mixture
            self.n_charge = 4.
            self.n_stoi = np.zeros(self.n_species)
            self.n_stoi[self.id_reac] = -1.0
            self.n_stoi[self.id_h2o] = 2.0
            self.n_stoi[self.id_inert] = 0.0
            # electrical charge number
            self.mol_mass = np.array([32., 18., 28.]) * 1.e-3
            # molar mass
        else:
            self.channel = ch.Channel(ch_dict.dict_anode_channel)
            reac_con_in = phy_prop.hydrogen_inlet_concentration
            # volumetric inlet hydrogen ratio
            self.inert_reac_ratio = (1. - reac_con_in) / reac_con_in
            # volumetric hydrogen to oxygen ratio
            self.n_species = 3
            # number of species in the gas mixture
            self.n_charge = 4.
            self.n_stoi = np.zeros(self.n_species)
            self.n_stoi[self.id_reac] = -2.0
            self.n_stoi[self.id_h2o] = 0.0
            self.n_stoi[self.id_inert] = 0.0
            # electrical charge number
            self.mol_mass = np.array([2., 18., 28.]) * 1.e-3
            # molar mass

        self.is_cathode = dict_hc['is_cathode']
        # anode is false; Cathode is true
        self.calc_act_loss = dict_hc['calc_act_loss']
        self.calc_cl_diff_loss = dict_hc['calc_cl_diff_loss']
        self.calc_gdl_diff_loss = dict_hc['calc_gdl_diff_loss']

        """geometry"""
        self.n_chl = dict_hc['channel_numb']
        # number of channels of each cell
        self.cell_width = dict_hc['cell_width']
        # height of the cell
        self.cell_length = dict_hc['cell_length']
        # length of the cell
        self.th_gdl = dict_hc['th_gdl']
        # thickness of the gas diffusion layer
        self.th_bpp = dict_hc['th_bpp']
        # thickness of the bipolar plate
        self.th_cl = dict_hc['th_cl']
        # thickness of the catalyst layer
        self.th_gde = self.th_gdl + self.th_cl
        # thickness gas diffusion electrode

        """voltage loss parameter, (Kulikovsky, 2013)"""
        self.vol_ex_cd = dict_hc['vol_ex_cd']
        # exchange current density
        self.prot_con_cl = dict_hc['prot_con_cl']
        # proton conductivity of the catalyst layer
        self.diff_coeff_cl = dict_hc['diff_coeff_cl']
        # diffusion coefficient of the reactant in the catalyst layer
        self.diff_coeff_gdl = dict_hc['diff_coeff_gdl']
        # diffusion coefficient of the reactant in the gas diffusion layer
        self.tafel_slope = dict_hc['tafel_slope']
        # tafel slope of the electrode
        self.i_sigma = np.sqrt(2. * self.vol_ex_cd * self.prot_con_cl
                               * self.tafel_slope)
        # could use a better name see (Kulikovsky, 2013) not sure if 2-D
        # exchange current densisty
        self.index_cat = self.n_nodes - 1
        # index of the first element with negative cell voltage
        self.i_ca_char = self.prot_con_cl * self.tafel_slope / self.th_cl
        # not sure if the name is ok, i_ca_char is the characteristic current
        # densisty, see (Kulikovsky, 2013)
        self.act_loss = np.zeros(n_ele)
        # activation voltage loss
        self.gdl_diff_loss = np.zeros(n_ele)
        # diffusion voltage loss at the gas diffusion layer
        self.cl_diff_loss = np.zeros(n_ele)
        # diffusion voltage loss at the catalyst layer
        self.v_loss = np.zeros(n_ele)
        # sum of the activation and diffusion voltage loss
        self.beta = np.zeros(n_ele)
        # dimensionless parameter
        self.var = np.zeros(n_ele)
        # term used in multiple functions
        self.i_ca_square = np.zeros(n_ele)
        # current density squared

        """general parameter"""
        area_fac = self.cell_length * self.cell_width\
            / (self.channel.active_area * self.n_chl)
        # factor active area with racks / active channel area
        # self.active_area_dx_ch = self.channel.active_area_dx
        # active area belonging to the channel plan area dx
        # self.active_area_ch = area_fac * self.channel.active_area
        # active area belonging to the channel plan area
        self.break_program = False
        # boolean to hint if the cell voltage runs below zero
        self.is_ht_pem = dict_cell['is_ht_pem']
        # if HT-PEMFC True; if NT-PEMFC False
        self.stoi = None
        # stoichiometry of the reactant at the channel inlet
        self.p_drop_bends = 0.
        # pressure drop in the channel through bends
        self.w_cross_flow = np.zeros(n_ele)
        # cross water flux through the membrane
        self.g_fluid = np.zeros(n_nodes)
        # heat capacity flow of the species mixture including fluid water
        self.cp_fluid = np.zeros(n_nodes)
        # heat capacity of the species mixture including fluid water
        # self.mol_flow_liq_w = np.zeros(n_nodes)
        # molar liquid water flux
        self.p = np.full(n_nodes, self.channel.p_out)
        # channel pressure
        self.cond_rate = np.zeros(n_nodes)
        # condensation rate of water
        self.humidity = np.zeros(n_nodes)
        # gas mixture humidity
        self.i_cd = np.full(n_ele, g_par.dict_case['tar_cd'])
        # current density
        self.u = np.zeros(n_nodes)
        # channel velocity
        self.fwd_mat = np.tril(np.full((n_ele, n_ele), 1.))
        # forward matrix
        self.bwd_mat = np.triu(np.full((n_ele, n_ele), 1.))
        # backward matrix
        self.mol_flow_total = np.zeros(n_nodes)
        self.mass_flow_total = np.zeros(n_nodes)
        self.mol_flow_gas_total = np.zeros(n_nodes)
        self.mass_flow_gas_total = np.zeros(n_nodes)
        # fluid mass flow
        self.q_gas = np.zeros(n_nodes)
        # molar flux of the gas phase
        self.mol_flow = np.full((self.n_species, n_nodes), 0.)
        self.mol_flow_liq = np.array(self.mol_flow)
        self.mol_flow_gas = np.array(self.mol_flow)
        self.mass_flow = np.array(self.mol_flow)
        self.mass_flow_gas = np.array(self.mol_flow)
        self.mass_flow_liq = np.array(self.mol_flow)

        # molar and mass flows of each species
        # (0: Reactant, 1: Water, 2: Inert Species

        self.gas_conc = np.array(self.mol_flow)
        # molar concentration of each species
        # 0: Reactant, 1: Water, 2: Inert Species
        self.reac_con_ele = np.full(n_ele, 0.)
        # element based molar concentration of the reactant
        self.temp_fluid = np.full(n_nodes, self.channel.temp_in)
        # temperature of the fluid in the channel
        self.rho_gas = np.full(n_nodes, 1.)
        # density of the gas phase
        self.visc_gas = np.full(n_nodes, 1.e-5)
        # viscosity of the gas phase
        self.mol_fraction = np.array(self.mol_flow)
        self.mol_fraction_gas = np.array(self.mol_fraction)
        # molar fraction of the species in the gas phase
        self.mass_fraction = np.array(self.mol_fraction)
        self.mass_fraction_gas = np.array(self.mass_fraction)

        # mass fraction of the species in the gas phase
        self.r_gas = np.array(self.mol_flow)
        # gas constant of the gas phase
        self.r_species = np.full(self.n_species, 0.)
        # gas constant of the species
        self.cp = np.array(self.mol_flow)
        # heat capacity of the species in the gas phase
        self.lambdas = np.array(self.mol_flow)
        # heat conductivity of the species in the gas phase
        self.visc = np.array(self.mol_flow)
        # viscosity of the species in the gas phase
        # self.temp_fluid_ele = np.zeros(n_ele)
        # element based temperature of the gas phase
        self.cp_gas = np.zeros(n_nodes)
        # heat capacity of the gas phase
        self.ht_coeff = np.zeros(n_ele)
        # convection coefficient between the gas phase and the channel
        self.k_ht_coeff_ca = np.zeros(n_ele)
        # heat conductivity between the gas phase and the channel
        self.cp_gas_ele = np.zeros(n_ele)
        # element based heat capacity
        self.lambda_gas = np.zeros(n_ele)
        # heat conductivity of the gas phase
        # self.Pr = np.zeros(n_ele)
        # prandtl number of the gas phase
        for i, item in enumerate(self.mol_mass):
            self.r_species[i] = g_par.dict_uni['R'] / item

    def update(self):
        """
        This function coordinates the program sequence
        """
        # self.calc_temp_fluid_ele()
        self.calc_mass_balance()
        if not self.break_program:
            self.calc_cond_rates()
            self.calc_species_properties()
            self.calc_gas_properties()
            self.calc_flow_velocity()
            self.calc_pressure()
            self.calc_humidity()
            self.calc_fluid_properties()
            self.calc_heat_transfer_coeff()
            self.update_voltage_loss()

    def calc_mass_balance(self):
        self.calc_reac_flow()
        self.calc_water_flow()
        self.mol_flow = np.maximum(self.mol_flow, 0.)
        self.mol_flow_total = np.sum(self.mol_flow, axis=0)
        self.mol_fraction = self.calc_fraction(self.mol_flow)
        self.mass_fraction = \
            self.molar_to_mass_fraction(self.mol_fraction, self.mol_mass)
        self.calc_mass_flow()
        self.calc_concentrations()
        self.calc_two_phase_flow()

    def update_voltage_loss(self):
        self.calc_electrode_loss()

    def set_layer_temperature(self, var):
        """
        This function sets the layer Temperatures,
        they can be obtained from the temperature system.
        """
        var = ip.interpolate_along_axis(np.array(var), axis=1,
                                        add_edge_points=True)
        if self.is_cathode:
            self.temp = np.array([var[0], var[1], var[2]])
        else:
            self.temp = np.array([var[0], var[1]])

    def calc_reac_flow(self):
        """
        Calculates the reactant molar flow [mol/s]
        """
        faraday = g_par.dict_uni['F']
        chl = self.channel
        tar_cd = g_par.dict_case['tar_cd']
        self.mol_flow[self.id_reac] = self.stoi * tar_cd * chl.active_area * \
            abs(self.n_stoi[self.id_reac]) / (self.n_charge * faraday)
        dmol = self.i_cd * chl.active_area_dx * self.n_stoi[self.id_reac] / \
            (self.n_charge * faraday)
        g_func.add_source(self.mol_flow[self.id_reac], dmol,
                          self.flow_direction)

    def calc_water_flow(self):
        """"
        Calculates the water and nitrogen molar flows [mol/s]
        """
        sat_p = w_prop.water.calc_p_sat(self.channel.temp_in)
        i_cd = self.i_cd
        area = self.channel.active_area_dx
        chl = self.channel
        q_0_water = \
            self.mol_flow[self.id_reac][0] * (1. + self.inert_reac_ratio) * \
            sat_p * chl.humidity_in / (chl.p_out - chl.humidity_in * sat_p)
        h2o_in = q_0_water
        h2o_source = np.zeros_like(i_cd)
        h2o_prod = area * self.n_stoi[self.id_h2o] * i_cd \
            / (self.n_charge * g_par.dict_uni['F'])
        h2o_source += h2o_prod
        h2o_cross = area * self.w_cross_flow * self.flow_direction
        h2o_source += h2o_cross
        self.mol_flow[self.id_h2o] = h2o_in
        g_func.add_source(self.mol_flow[self.id_h2o],
                          h2o_source, self.flow_direction)

        if self.is_cathode:
            self.mol_flow[self.id_inert] = \
                np.full(self.n_nodes,
                        self.mol_flow[self.id_reac][0] * self.inert_reac_ratio)
        else:
            self.mol_flow[self.id_inert] = \
                np.full(self.n_nodes,
                        self.mol_flow[self.id_reac][-1] * self.inert_reac_ratio)

    def calc_mass_flow(self):
        """
        Calculates the relevant mass flows
        """
        self.mass_flow = (self.mol_flow.transpose() * self.mol_mass).transpose()
        self.mass_flow_total = np.sum(self.mass_flow, axis=0)
        self.mass_flow_total = np.sum(self.mass_flow, axis=0)

    def calc_pressure(self):
        """
        Calculates the total channel pressure for each element.
        """
        chl = self.channel
        rho_ele = ip.interpolate_1d(self.rho_gas)
        u_ele = ip.interpolate_1d(self.u)
        reynolds = self.rho_gas * chl.d_h * self.u / self.visc_gas
        reynolds_ele = ip.interpolate_1d(reynolds)
        zeta_bends = chl.bend_fri_fac * chl.n_bends / self.n_ele
        friction_factor = 64.0 / reynolds_ele
        dp = (friction_factor * chl.dx / chl.d_h + zeta_bends) \
            * rho_ele * 0.5 * u_ele ** 2.0
        pressure_direction = -self.flow_direction
        self.p.fill(chl.p_out)
        g_func.add_source(self.p, dp, pressure_direction)

    @staticmethod
    def molar_to_mass_fraction(mol_fraction, mol_mass):
        """
        Calculates the mass fraction from molar fractions and molar masses.
        Molar fractions must be a (multi-dimensional) array with species
        along the first (0th) axis. Molar masses must be a one-dimensional
        array with the different species molar masses.
        """
        x_mw = mol_fraction.transpose() * mol_mass
        return x_mw.transpose() / np.sum(x_mw, axis=1)

    @staticmethod
    def calc_fraction(species_flow):
        """
        Calculates the species mixture fractions based on a multi-dimensional
        array with different species along the first (0th) axis.
        """
        return species_flow / np.sum(species_flow, axis=0)

    def calc_concentrations(self):
        """
        Calculates the gas phase molar concentrations.
        """
        id_reac = 0
        id_h2o = 1
        id_inert = 2

        gas_constant = g_par.dict_uni['R']
        total_mol_conc = self.p / (gas_constant * self.temp_fluid)
        conc = total_mol_conc * self.mol_fraction
        p_sat = w_prop.water.calc_p_sat(self.temp_fluid)
        sat_conc = p_sat / (gas_constant * self.temp_fluid)
        dry_air_mol_flow = np.copy(self.mol_flow)
        dry_air_mol_flow[id_h2o] = 0.0
        dry_air_fraction = self.calc_fraction(dry_air_mol_flow)
        print(self.name, 'dry_air_fraction: ', dry_air_fraction)
        self.gas_conc = conc
        self.gas_conc[id_reac] = \
            np.where(self.gas_conc[id_h2o] > sat_conc,
                     (total_mol_conc - sat_conc) * dry_air_fraction[id_reac],
                     self.gas_conc[id_reac])
        self.gas_conc[id_inert] = \
            np.where(self.gas_conc[id_h2o] > sat_conc,
                     (total_mol_conc - sat_conc) * dry_air_fraction[id_inert],
                     self.gas_conc[id_inert])
        self.gas_conc[id_h2o] = \
            np.where(self.gas_conc[id_h2o] > sat_conc,
                     sat_conc, self.gas_conc[id_h2o])
        self.reac_con_ele = ip.interpolate_1d(self.gas_conc[id_reac])
        print(self.name, 'gas_conc: ', self.gas_conc)

    def calc_two_phase_flow(self):
        """
        Calculates the condensed phase flow and updates mole and mass fractions
        """
        id_reac = 0
        id_h2o = 1
        self.mol_flow_liq = np.zeros_like(self.mol_flow)
        self.mol_flow_liq[id_h2o] = self.mol_flow[id_h2o] \
            - self.gas_conc[id_h2o] / \
            self.gas_conc[id_reac] * self.mol_flow[id_reac]
        self.mass_flow_liq = np.zeros_like(self.mass_flow)
        self.mass_flow_liq[id_h2o] = \
            self.mol_flow_liq[id_h2o] * self.mol_mass[id_h2o]
        self.mol_flow_gas = self.mol_flow - self.mol_flow_liq
        self.mass_flow_gas = self.mass_flow - self.mass_flow_liq
        self.mol_flow_gas_total = np.sum(self.mol_flow_gas, axis=0)
        self.mass_flow_gas_total = np.sum(self.mass_flow_gas, axis=0)
        self.q_gas = self.mol_flow_gas_total
        self.mol_fraction_gas = self.calc_fraction(self.mol_flow_gas)
        self.mass_fraction_gas = self.calc_fraction(self.mass_flow_gas)

    def calc_species_properties(self):
        """
        Calculates the properties of the species in the gas phase
        """
        id_reac = 0
        id_h2o = 1
        id_inert = 2
        if self.is_cathode:
            self.cp[id_reac] = g_fit.oxygen.calc_cp(self.temp_fluid)
            self.lambdas[id_reac] = \
                g_fit.oxygen.calc_lambda(self.temp_fluid, self.p)
            self.visc[id_reac] = g_fit.oxygen.calc_visc(self.temp_fluid)
        else:
            self.cp[id_reac] = g_fit.hydrogen.calc_cp(self.temp_fluid)
            self.lambdas[id_reac] = \
                g_fit.hydrogen.calc_lambda(self.temp_fluid, self.p)
            self.visc[id_reac] = g_fit.hydrogen.calc_visc(self.temp_fluid)
        self.cp[id_h2o] = g_fit.water.calc_cp(self.temp_fluid)
        self.cp[id_inert] = g_fit.nitrogen.calc_cp(self.temp_fluid)
        self.lambdas[id_h2o] = g_fit.water.calc_lambda(self.temp_fluid, self.p)
        self.lambdas[id_inert] = g_fit.nitrogen.calc_lambda(self.temp_fluid, self.p)
        self.visc[id_h2o] = g_fit.water.calc_visc(self.temp_fluid)
        self.visc[id_inert] = g_fit.nitrogen.calc_visc(self.temp_fluid)

    def calc_gas_properties(self):
        """
        Calculates the properties of the gas phase
        """
        self.r_gas = \
            np.sum(self.mass_fraction_gas.transpose() * self.r_species, axis=1)
        self.cp_gas = \
            np.sum(self.mass_fraction_gas * self.cp, axis=0)
        self.cp_gas_ele = ip.interpolate_1d(self.cp_gas)
        self.visc_gas = \
            g_func.calc_visc_mix(self.visc, self.mol_fraction_gas,
                                 self.mol_mass)
        self.lambda_gas = \
            g_func.calc_lambda_mix(self.lambdas, self.mol_fraction_gas,
                                   self.visc, self.mol_mass)
        self.lambda_gas_ele = ip.interpolate_1d(self.lambda_gas)
        self.rho_gas = g_func.calc_rho(self.p, self.r_gas, self.temp_fluid)
        # self.Pr = self.visc_gas * self.cp_gas / self.lambda_gas

    def calc_flow_velocity(self):
        """
        Calculates the gas phase velocity.
        The gas phase velocity is taken to be the liquid water velocity as well.
        """
        self.u = self.mass_flow_gas_total / self.rho_gas \
            * self.channel.cross_area

    def calc_fluid_properties(self):
        """
        Calculate the fluid flow properties
        """
        cp_liq = g_par.dict_case['cp_liq']
        self.cp_fluid = \
            ((self.mass_flow_total - self.mass_flow_gas_total) * cp_liq
             + self.mass_flow_gas_total * self.cp_gas) / self.mass_flow_total
        self.g_fluid = self.mass_flow_total * self.cp_fluid

    def calc_heat_transfer_coeff(self):
        """
        Calculates the convection coefficient between the channel and
        the gas phase.
        Deputy for the convection coefficient between the fluid and the channel.
        Calculates the heat conductivity based on the convection coefficient
        and the element area of the channel.
        """
        nusselt = 3.66
        chl = self.channel
        self.ht_coeff = self.lambda_gas_ele * nusselt / chl.d_h
        self.k_ht_coeff_ca = \
            self.ht_coeff * chl.dx * 2.0 * (chl.width + chl.height)

    def calc_cond_rates(self):
        """
        Calculates the molar condensation rate of water in the channel.
        """
        id_h2o = 1
        cond_rate_ele = np.ediff1d(self.mol_flow_liq[id_h2o])
        self.cond_rate = \
            self.flow_direction * ip.interpolate_1d(cond_rate_ele,
                                                    add_edge_points=True)

    def calc_humidity(self):
        """
        Calculates the relative humidity of the fluid.
        """
        # self.humidity = self.gas_con[1] * g_par.dict_uni['R'] \
        #     * self.temp_fluid / w_prop.water.calc_p_sat(self.temp_fluid)
        p_sat = w_prop.water.calc_p_sat(self.temp_fluid)
        self.humidity = self.mol_fraction_gas[1] * self.p / p_sat

    def calc_voltage_losses_parameter(self):
        """
        Calculates multiply used supporting parameters
        to calculate the voltage loss according to (Kulikovsky, 2013).
        """
        i_lim = 4. * g_par.dict_uni['F'] * self.gas_conc[0, :-1] \
            * self.diff_coeff_gdl / self.th_gdl
        # print(self.name, i_lim)
        # print(self.name, self.gas_con[0, :])
        self.var = \
            1. - self.i_cd / (i_lim * self.reac_con_ele / self.gas_conc[0, :-1])
        self.i_ca_square = np.square(self.i_cd)
        # print(self.name, self.var)

    def calc_activation_loss(self):
        """
        Calculates the activation voltage loss,
        according to (Kulikovsky, 2013).
        """
        self.act_loss = self.tafel_slope \
            * np.arcsinh((self.i_cd / self.i_sigma) ** 2.
                         / (2. * (self.reac_con_ele / self.gas_conc[0, :-1])
                            * (1. - np.exp(-self.i_cd /
                                           (2. * self.i_ca_char)))))

    def calc_transport_loss_catalyst_layer(self):
        """
        Calculates the diffusion voltage loss in the catalyst layer
        according to (Kulikovsky, 2013).
        """
        i_hat = self.i_cd / self.i_ca_char
        short_save = np.sqrt(2. * i_hat)
        beta = short_save / (1. + np.sqrt(1.12 * i_hat) * np.exp(short_save))\
            + np.pi * i_hat / (2. + i_hat)
        self.cl_diff_loss = \
            ((self.prot_con_cl * self.tafel_slope ** 2.)
             / (4. * g_par.dict_uni['F']
                * self.diff_coeff_cl * self.reac_con_ele)
                * (self.i_cd / self.i_ca_char
                   - np.log10(1. + self.i_ca_square /
                              (self.i_ca_char ** 2. * beta ** 2.)))) / self.var

    def calc_transport_loss_diffusion_layer(self):
        """
        Calculates the diffusion voltage loss in the gas diffusion layer
        according to (Kulikovsky, 2013).
        """
        self.gdl_diff_loss = -self.tafel_slope * np.log10(self.var)
        nan_list = np.isnan(self.gdl_diff_loss)
        if nan_list.any():
            self.gdl_diff_loss[np.argwhere(nan_list)[0, 0]:] = 1.e50

    def calc_electrode_loss(self):
        """
        Calculates the full voltage losses of the electrode
        """
        self.calc_voltage_losses_parameter()
        self.calc_activation_loss()
        self.calc_transport_loss_catalyst_layer()
        self.calc_transport_loss_diffusion_layer()
        if not self.calc_gdl_diff_loss:
            self.gdl_diff_loss = 0.
        if not self.calc_cl_diff_loss:
            self.cl_diff_loss = 0.
        if not self.calc_act_loss:
            self.act_loss = 0.
        self.v_loss = self.act_loss + self.cl_diff_loss + self.gdl_diff_loss
        # print(self.name, 'gdl loss: ', self.gdl_diff_loss)
        # print(self.name, 'cl loss: ', self.cl_diff_loss)
        # print(self.name, 'activation loss: ', self.act_loss)
        # print(self.name, 'electrode loss: ', self.v_loss)
