import warnings
import data.gas_properties as g_prop
import numpy as np
import data.global_parameters as g_par
import system.global_functions as g_func
import system.channel as chl
import system.fluid2 as fluids
import system.layers as layers
import sys
import system.interpolation as ip

warnings.filterwarnings("ignore")


class HalfCell:

    # Class variables constant across all instances of the class
    # (under construction)

    def __init__(self, halfcell_dict, cell_dict, channel):
        self.name = halfcell_dict['name']
        self.n_nodes = g_par.dict_case['nodes']
        n_ele = self.n_nodes - 1
        self.n_ele = n_ele
        # Discretization in elements and nodes along the x-axis (flow_circuit.py axis)

        """half cell geometry parameter"""
        self.width = cell_dict["width"]
        self.length = cell_dict["length"]

        # Get reference to channel object
        self.channel = channel
        # number of channels of each half cell
        self.n_channel = halfcell_dict['channel_number']
        area_factor = self.length * self.width \
            / (self.channel.base_area * self.n_channel)
        if area_factor < 1.0:
            raise ValueError('width and length of cell result in a cell '
                             'surface area  smaller than the area covered by '
                             'channels')

        self.rib_width = self.channel.width * (area_factor - 1.0)
        self.width_straight_channels = self.channel.width * self.n_channel \
            + self.rib_width * (self.n_channel + 1)
        self.length_straight_channels = (self.length * self.width) \
            / self.width_straight_channels
        self.active_area = area_factor * self.channel.base_area
        # self.active_area = area_factor * self.channel.base_area
        # factor active area with ribs / active channel area
        self.active_area_dx = area_factor * self.channel.base_area_dx
        # self.active_area_dx = area_factor * self.channel.base_area_dx

        self.id_fuel = 0
        self.id_h2o = 2
        self.id_inert = 1
        # self.species_names = halfcell_dict['species_names']
        # self.species = []
        # for i, name in enumerate(self.species_names):
        #     for j in range(len(g_prop.species)):
        #         if name == g_prop.species[j].name:
        #             self.species.append(g_prop.species[j])
        #             break
        # self.n_species = len(self.species_names)
        self.n_charge = halfcell_dict['charge_number']
        self.n_stoi = np.asarray(halfcell_dict['reaction_stoichiometry'])
        # self.mol_mass = np.asarray(halfcell_dict['molar_mass'])
        # check if the object is an anode or a cathode
        # catalyst layer specific handover
        self.inlet_composition = halfcell_dict['inlet_composition']
        self.inert_reac_ratio = \
            self.inlet_composition[self.id_inert] \
            / self.inlet_composition[self.id_fuel]

        self.faraday = g_par.constants['F']
        self.target_cd = g_par.dict_case['target_current_density']

        self.is_cathode = halfcell_dict['is_cathode']
        # anode is false; Cathode is true
        self.calc_act_loss = halfcell_dict['calc_act_loss']
        self.calc_cl_diff_loss = halfcell_dict['calc_cl_diff_loss']
        self.calc_gdl_diff_loss = halfcell_dict['calc_gdl_diff_loss']

        self.th_gdl = halfcell_dict['th_gdl']
        # thickness of the gas diffusion layer
        self.th_bpp = halfcell_dict['th_bpp']
        # thickness of the bipolar plate
        self.th_cl = halfcell_dict['th_cl']
        # thickness of the catalyst layer
        self.th_gde = self.th_gdl + self.th_cl
        # thickness gas diffusion electrode

        bpp_layer_dict = \
            {'thickness': halfcell_dict['th_bpp'],
             'width': self.width_straight_channels,
             'length': self.length_straight_channels,
             'electrical conductivity':
                 cell_dict['electrical conductivity bpp'],
             'thermal conductivity':
                 cell_dict['thermal conductivity bpp']}
        # 'porosity': self.channel.cross_area * self.n_channel / (
        #             self.th_bpp * self.width)}
        self.bpp = layers.SolidLayer(bpp_layer_dict, self.channel.dx)
        gde_layer_dict = \
            {'thickness': halfcell_dict['th_gdl'] + halfcell_dict['th_cl'],
             'width': self.width_straight_channels,
             'length': self.length_straight_channels,
             'electrical conductivity':
                 cell_dict['electrical conductivity gde'],
             'thermal conductivity':
                 cell_dict['thermal conductivity gde']}
        # 'porosity':
        #    (self.th_gdl * halfcell_dict['porosity gdl']
        #     + self.th_cl * halfcell_dict['porosity cl'])
        #    / (self.th_gde + self.th_cl)}
        self.gde = layers.SolidLayer(gde_layer_dict, self.channel.dx)

        self.thickness = self.bpp.thickness + self.gde.thickness

        """voltage loss parameter, (Kulikovsky, 2013)"""
        vol_ex_cd = halfcell_dict['vol_ex_cd']
        # exchange current density
        self.prot_con_cl = halfcell_dict['prot_con_cl']
        # proton conductivity of the catalyst layer
        self.diff_coeff_cl = halfcell_dict['diff_coeff_cl']
        # diffusion coefficient of the reactant in the catalyst layer
        self.diff_coeff_gdl = halfcell_dict['diff_coeff_gdl']
        # diffusion coefficient of the reactant in the gas diffusion layer
        self.tafel_slope = halfcell_dict['tafel_slope']
        # tafel slope of the electrode
        self.i_sigma = np.sqrt(2. * vol_ex_cd * self.prot_con_cl
                               * self.tafel_slope)
        # could use a better name see (Kulikovsky, 2013) not sure if 2-D
        # exchange current densisty
        self.index_cat = self.n_nodes - 1
        # index of the first element with negative cell voltage
        self.i_cd_char = self.prot_con_cl * self.tafel_slope / self.th_cl
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

        self.break_program = False
        # boolean to hint if the cell voltage runs below zero
        # if HT-PEMFC True; if NT-PEMFC False
        self.stoi = halfcell_dict['stoichiometry']
        # # stoichiometry of the reactant at the channel inlet
        # self.p_drop_bends = 0.
        # pressure drop in the channel through bends
        self.w_cross_flow = np.zeros(n_ele)
        # cross water flux through the membrane
        # self.g_fluid = np.zeros(n_nodes)
        # heat capacity flow_circuit.py of the species mixture including fluid water
        # self.cp_fluid = np.zeros(n_nodes)
        # heat capacity of the species mixture including fluid water
        # self.mol_flow_liq_w = np.zeros(n_nodes)
        # molar liquid water flux
        # self.p = np.full(n_nodes, self.channel.p_out)
        # # channel pressure
        # self.cond_rate = np.zeros(n_nodes)
        # self.cond_rate_ele = np.zeros(n_ele)
        # # condensation rate of water
        # self.humidity = np.zeros(n_nodes)
        # # gas mixture humidity
        # self.u = np.zeros(n_nodes)
        # channel velocity
        # self.fwd_mat = np.tril(np.full((n_ele, n_ele), 1.))
        # forward matrix
        # self.bwd_mat = np.triu(np.full((n_ele, n_ele), 1.))
        # backward matrix
        # self.mol_flow_total = np.zeros(n_nodes)
        # self.mass_flow_total = np.zeros(n_nodes)
        # self.mol_flow_gas_total = np.zeros(n_nodes)
        # self.mass_flow_gas_total = np.zeros(n_nodes)

        # # fluid mass flow_circuit.py
        # self.vol_flow_gas = np.zeros(n_nodes)
        # # molar flux of the gas phase
        # (0: Reactant, 1: Water, 2: Inert Species
        # self.mol_flow = np.full((self.n_species, n_nodes), 0.)
        # self.mol_flow_liq = np.array(self.mol_flow)
        # self.mol_flow_gas = np.array(self.mol_flow)
        # self.mass_flow = np.array(self.mol_flow)
        # self.mass_flow_gas = np.array(self.mol_flow)
        # self.mass_flow_liq = np.array(self.mol_flow)
        # self.gas_conc = np.array(self.mol_flow)
        # molar concentration of each species
        # self.mol_fraction = np.array(self.mol_flow)
        # self.mol_fraction_gas = np.array(self.mol_fraction)
        # # molar fraction of the species in the gas phase
        # self.mass_fraction = np.array(self.mol_fraction)
        # self.mass_fraction_gas = np.array(self.mass_fraction)
        # self.temp_fluid = np.full(n_nodes, self.channel.temp_in)
        # self.temp_fluid_ele = np.full(n_ele, self.channel.temp_in)
        # # temperature of the fluid in the channel
        # self.rho_gas = np.full(n_nodes, 1.)
        # # density of the gas phase
        # self.visc_gas = np.full(n_nodes, 1.e-5)
        # # viscosity of the gas phase
        #
        #
        # # mass fraction of the species in the gas phase
        # self.r_gas = np.full(n_nodes, 0.)
        # # gas constant of the gas phase
        # # self.r_species = np.full(self.n_species, 0.)
        # # gas constant of the species
        # self.cp = np.array(self.mol_flow)
        # # heat capacity of the species in the gas phase
        # self.lambdas = np.array(self.mol_flow)
        # # heat conductivity of the species in the gas phase
        # self.visc = np.array(self.mol_flow)
        # # viscosity of the species in the gas phase
        # # self.temp_fluid_ele = np.zeros(n_ele)
        # # element based temperature of the gas phase
        # self.cp_gas = np.zeros(n_nodes)
        # # heat capacity of the gas phase
        # self.ht_coeff = np.zeros(n_ele)
        # # convection coefficient between the gas phase and the channel
        # self.k_ht_coeff = np.zeros(n_ele)
        # # heat conductivity between the gas phase and the channel
        # self.cp_gas_ele = np.zeros(n_ele)
        # # element based heat capacity
        # self.lambda_gas = np.zeros(n_nodes)
        # self.lambda_gas_ele = np.zeros(n_ele)
        # # heat conductivity of the gas phase
        # self.Pr = np.zeros(n_ele)
        # # prandtl number of the gas phase
        # for i, item in enumerate(self.mol_mass):
        #     self.r_species[i] = g_par.constants['R'] / item

    def update(self, current_density):
        """
        This function coordinates the program sequence
        """
        # self.calc_temp_fluid_ele()
        # mole_flow_in, mole_source = self.calc_mass_balance(current_density)
        if not self.break_program:
            # self.channel.update(mole_flow_in, mole_source)
            # self.channel.mole_flow[:] = mole_flow_in
            self.channel.mass_source[:], self.channel.mole_source[:] = \
                self.calc_mass_source(current_density)
            self.update_voltage_loss(current_density)

    # def calc_mass_balance(self, current_density, stoi=None):
    #     n_species = self.channel.fluid.n_species
    #     mole_flow_in = np.zeros((n_species, self.n_nodes))
    #     mole_source = np.zeros((n_species, self.n_ele))
    #     mole_flow_in[self.id_fuel, :], mole_source[self.id_fuel, :] = \
    #         self.calc_fuel_flow(current_density, stoi)
    #     mole_flow_in[self.id_inert, :] = \
    #         mole_flow_in[self.id_fuel, self.channel.id_in] \
    #         * self.inert_reac_ratio
    #     air_flow_in = np.sum(mole_flow_in[:, self.channel.id_in])
    #     mole_flow_in[self.id_h2o, :], mole_source[self.id_h2o, :] = \
    #         self.calc_water_flow(current_density, air_flow_in)
    #     return mole_flow_in, mole_source

    def calc_mass_balance(self, current_density, stoi=None):
        mass_flow_in, mole_flow_in = self.calc_inlet_flow(stoi)
        mass_flow_in = g_func.fill_transposed(mass_flow_in,
                                              self.channel.mass_flow.shape)
        mole_flow_in = g_func.fill_transposed(mole_flow_in,
                                              self.channel.mole_flow.shape)
        mass_source, mole_source = self.calc_mass_source(current_density)
        return mass_flow_in, mole_flow_in, mass_source, mole_source

    def calc_inlet_flow(self, stoi=None):
        if stoi is None:
            stoi = self.stoi
        mole_flow_in = np.zeros(self.channel.fluid.n_species)
        mole_flow_in[self.id_fuel] = self.target_cd * self.active_area * stoi \
            * abs(self.n_stoi[self.id_fuel]) / (self.n_charge * self.faraday)

        inlet_composition = \
            self.channel.fluid.mole_fraction[:, self.channel.node_in]
        for i in range(len(mole_flow_in)):
            if i != self.id_fuel:
                mole_flow_in[i] = mole_flow_in[self.id_fuel] \
                    * inlet_composition[i] / inlet_composition[self.id_fuel]
        mass_flow_in = mole_flow_in / self.channel.fluid.species.mw
        return mass_flow_in, mole_flow_in

    def calc_mass_source(self, current_density):
        mole_source = np.zeros((self.channel.fluid.n_species, self.n_ele))

        for i in range(len(mole_source)):
            mole_source[i] = current_density * self.active_area_dx \
                * self.n_stoi[i] / (self.n_charge * self.faraday)

        # water cross flow
        mole_source[self.id_h2o] += self.active_area_dx * self.w_cross_flow \
            * self.channel.flow_direction
        mass_source = (mole_source.transpose()
                       / self.channel.fluid.species.mw).transpose()
        return mass_source, mole_source

    def calc_fuel_flow(self, current_density, stoi=None):
        """
        Calculates the reactant molar flow [mol/s]
        """
        if stoi is None:
            stoi = self.stoi
        mol_flow_in = self.target_cd * self.active_area * stoi \
            * abs(self.n_stoi[self.id_fuel]) / (self.n_charge * self.faraday)
        dmol = current_density * self.active_area_dx \
            * self.n_stoi[self.id_fuel] / (self.n_charge * self.faraday)
        # g_func.add_source(self.mol_flow[self.id_fuel], dmol,
        #                   self.flow_direction)
        return mol_flow_in, dmol

    def calc_water_flow(self, current_density, air_flow_in):
        """"
        Calculates the water molar flow [mol/s]
        """
        if not isinstance(self.channel.fluid, fluids.TwoPhaseMixture):
            raise TypeError('Fluid in channel must be of type TwoPhaseMixture')
        id_in = self.channel.id_in
        humidity_in = self.channel.fluid.humidity[id_in]
        sat_p = self.channel.fluid.saturation_pressure[id_in]
        mol_flow_in = air_flow_in * sat_p * humidity_in / \
            (self.channel.p[id_in] - humidity_in * sat_p)
        dmol = np.zeros_like(current_density)
        h2o_prod = self.active_area_dx * self.n_stoi[self.id_h2o] \
            * current_density / (self.n_charge * self.faraday)
        dmol += h2o_prod
        h2o_cross = self.active_area_dx * self.w_cross_flow \
            * self.channel.flow_direction
        dmol += h2o_cross
        return mol_flow_in, dmol

    def update_voltage_loss(self, current_density):
        self.calc_electrode_loss(current_density)

    def calc_activation_loss(self, current_density, reac_conc_ele,
                             reac_conc_in):
        """
        Calculates the activation voltage loss,
        according to (Kulikovsky, 2013).
        """
        try:
            self.act_loss[:] = self.tafel_slope \
                * np.arcsinh((current_density / self.i_sigma) ** 2.
                             / (2. * (reac_conc_ele / reac_conc_in)
                                * (1. - np.exp(-current_density /
                                               (2. * self.i_cd_char)))))
        except FloatingPointError:
            print(current_density)
            raise

    def calc_transport_loss_catalyst_layer(self, current_density, var,
                                           reac_conc_ele):
        """
        Calculates the diffusion voltage loss in the catalyst layer
        according to (Kulikovsky, 2013).
        """
        i_hat = current_density / self.i_cd_char
        short_save = np.sqrt(2. * i_hat)
        beta = short_save / (1. + np.sqrt(1.12 * i_hat) * np.exp(short_save))\
            + np.pi * i_hat / (2. + i_hat)
        self.cl_diff_loss[:] = \
            ((self.prot_con_cl * self.tafel_slope ** 2.)
             / (4. * self.faraday
                * self.diff_coeff_cl * reac_conc_ele)
             * (current_density / self.i_cd_char
                - np.log10(1. + np.square(current_density) /
                           (self.i_cd_char ** 2. * beta ** 2.)))) / var

    def calc_transport_loss_diffusion_layer(self, var):
        """
        Calculates the diffusion voltage loss in the gas diffusion layer
        according to (Kulikovsky, 2013).
        """
        self.gdl_diff_loss[:] = -self.tafel_slope * np.log10(var)
        nan_list = np.isnan(self.gdl_diff_loss)
        if nan_list.any():
            self.gdl_diff_loss[np.argwhere(nan_list)[0, 0]:] = 1.e50

    def calc_electrode_loss(self, current_density):
        """
        Calculates the full voltage losses of the electrode
        """
        reac_conc = self.channel.fluid.concentration[self.id_fuel]
        reac_conc_ele = ip.interpolate_1d(reac_conc)
        if self.channel.flow_direction == 1:
            reac_conc_in = reac_conc[:-1]
        else:
            reac_conc_in = reac_conc[1:]

        i_lim = 4. * self.faraday * reac_conc_in \
            * self.diff_coeff_gdl / self.th_gdl
        var = 1. - current_density / i_lim * reac_conc_in / reac_conc_ele

        self.calc_activation_loss(current_density, reac_conc_ele, reac_conc_in)
        self.calc_transport_loss_catalyst_layer(current_density, var,
                                                reac_conc_ele)
        self.calc_transport_loss_diffusion_layer(var)
        if not self.calc_gdl_diff_loss:
            self.gdl_diff_loss[:] = 0.
        if not self.calc_cl_diff_loss:
            self.cl_diff_loss[:] = 0.
        if not self.calc_act_loss:
            self.act_loss[:] = 0.
        self.v_loss[:] = self.act_loss + self.cl_diff_loss + self.gdl_diff_loss
