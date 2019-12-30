import numpy as np
import data.global_parameters as g_par
import system.fluid as fluids
import system.global_functions as g_func
import system.interpolation as ip
from abc import ABC, abstractmethod
from system.output_object import OutputObject


class Channel(ABC, OutputObject):
    def __new__(cls, channel_dict, fluid=None):
        # Create fluid object
        nx = g_par.dict_case['elements'] + 1
        if fluid is None:
            fluid = \
                fluids.Fluid(nx, channel_dict['fluid_name'],
                             channel_dict.get('fluid_components', None),
                             mole_fractions_init=
                             channel_dict.get('inlet_composition', None),
                             fluid_props=
                             channel_dict.get('liquid_properties', None))

        if type(fluid) is fluids.IncompressibleFluid:
            return super(Channel, cls).__new__(IncompressibleFluidChannel)
        elif type(fluid) is fluids.GasMixture:
            return super(Channel, cls).__new__(GasMixtureChannel)
        elif type(fluid) is fluids.TwoPhaseMixture:
            return super(Channel, cls).__new__(TwoPhaseMixtureChannel)
        else:
            raise NotImplementedError('Only Channel types of '
                                      'IncompressibleFluidChannel, '
                                      'GasMixtureChannel, and '
                                      'TwoPhaseMixtureChannel are '
                                      'implemented')

    def __init__(self, channel_dict, fluid):
        super().__init__()
        self.name = channel_dict['name']
        self.fluid = fluid
        self.length = channel_dict['channel_length']
        # channel length
        self.n_ele = g_par.dict_case['elements']
        self.n_nodes = self.n_ele + 1
        self.x = np.linspace(0.0, self.length, self.n_nodes)
        self.dx = np.diff(self.x)
        # element length
        self.p_out = channel_dict['p_out']
        self.p = np.full(self.n_nodes, self.p_out)
        # outlet and initial pressure
        self.temp_in = channel_dict['temp_in']
        self.temp = np.full_like(self.p, self.temp_in)
        # inlet temperature

        # Geometry
        self.flow_direction = channel_dict['flow_direction']
        if self.flow_direction not in (-1, 1):
            raise ValueError('Member variable flow_direction '
                             'must be either 1 or -1')
        # flow direction
        if self.flow_direction == 1:
            self.id_in = 0
        else:
            self.id_in = -1
        self.width = channel_dict['channel_width']
        # channel width
        self.height = channel_dict['channel_height']
        # channel height
        self.n_bends = channel_dict.get('bend_number', 0)
        # number of channel bends
        self.zeta_bends = channel_dict.get('bend_friction_factor', 0.0)
        self.zeta_other = \
            channel_dict.get('additional_friction_fractor', 0.0)
        # bend friction factor
        self.base_area = self.width * self.length
        # planar area of the channel
        self.base_area_dx = self.width * self.dx
        # planar area of an element of the channel
        self.cross_area = self.width * self.height
        # channel cross-sectional area
        self.d_h = 4. * self.cross_area / (2. * (self.width + self.height))
        # channel hydraulic diameter

        # Flow
        self.velocity = np.zeros(self.n_nodes)
        self.reynolds = np.zeros(self.n_nodes)
        self.mass_flow_total = np.zeros(self.n_nodes)
        self.vol_flow = np.zeros(self.n_nodes)

        # Heat Transfer
        self.k_coeff = np.zeros(self.n_ele)

        self.add_print_data(self.temp, 'Fluid Temperature', 'K')
        self.add_print_data(self.p, 'Fluid Pressure', 'Pa')

    def update(self, *args, **kwargs):
        self.calc_mass_balance(*args, **kwargs)
        self.fluid.update(self.temp, self.p, *args, **kwargs)
        self.calc_flow_velocity()
        self.calc_pressure()
        self.calc_heat_transfer_coeff()

    def calc_flow_velocity(self):
        """
        Calculates the gas phase velocity.
        """
        self.vol_flow[:] = self.mass_flow_total / self.fluid.density
        self.velocity[:] = self.vol_flow / self.cross_area
        self.reynolds[:] = \
            self.velocity * self.d_h * self.fluid.density / self.fluid.viscosity

    @abstractmethod
    def calc_mass_balance(self, *args, **kwargs):
        pass

    def calc_heat_transfer_coeff(self):
        """
        Calculates heat transfer coefficient of channel assuming element-wise
        constant wall temperature
        (Correlations should be reviewed)
        """
        prandtl = self.fluid.density * self.fluid.specific_heat / \
            self.fluid.thermal_conductivity

        d_by_l = self.d_h / self.length
        nu_1 = 3.66
        nu_2 = 1.66 * np.sqrt(self.reynolds * prandtl * d_by_l)
        nu_3 = (2. / (1. + 22. * prandtl)) ** (1. / 6.) \
            * np.sqrt(self.reynolds * prandtl * d_by_l)
        nu_lam = (nu_1 ** 3. + 0.7 ** 3. + (nu_2 - 0.7) ** 3.
                  + nu_3 ** 3.) ** (1. / 3.)
        # laminar nusselt number
        zeta = (1.8 * np.log(self.reynolds) - 1.5) ** -2.
        nu_turb = zeta / 8. * self.reynolds * prandtl \
            / (1. + 12.7 * np.sqrt(zeta / 8.)
                * (prandtl ** 2. / 3.) - 1.) \
            * (1. + d_by_l ** 2. / 3.)
        # if reynolds <= 2300.:
        #     nusselt = nu_lam
        # elif 2300. < reynolds < 1.e4:
        #     gamma = (reynolds - 2300.) / 7700.
        #     nusselt = (1. - (reynolds - 2300.) / 7700.) * nu_lam + (reynolds - 2300.) / 7700. * nu_turb
        # else:
        #     nusselt = nu_turb
        nusselt = \
            np.where(self.reynolds < 2300.0, nu_lam,
                     np.where(self.reynolds < 1e4, nu_lam,
                              (1. - (self.reynolds - 2300.) / 7700.) * nu_lam
                              + (self.reynolds - 2300.) / 7700. * nu_turb))

        ht_coeff = nusselt * self.fluid.thermal_conductivity / self.d_h
        # convection coefficient between the coolant and the channel wall
        conv_area = self.d_h * np.pi * self.length / self.n_ele
        # convection area of the channel wall
        self.k_coeff[:] = ht_coeff * conv_area

    def calc_pressure(self):
        """
        Calculates the static channel pressure
        """
        density_ele = ip.interpolate_1d(self.fluid.density)
        velocity_ele = ip.interpolate_1d(self.velocity)
        reynolds_ele = ip.interpolate_1d(self.reynolds)
        zeta_bends = self.zeta_bends * self.n_bends / self.n_ele
        zeta = zeta_bends + self.zeta_other
        f_ele = g_func.calc_friction_factor(reynolds_ele)
        # friction_factor = 64.0 / reynolds_ele
        dp = g_func.calc_pressure_drop(self.velocity, density_ele, f_ele, zeta,
                                       self.dx, self.d_h)
        # dp = (f_ele * self.dx / self.d_h + zeta_bends) \
        #     * density_ele * 0.5 * velocity_ele ** 2.0
        pressure_direction = -self.flow_direction
        self.p.fill(self.p_out)
        g_func.add_source(self.p, dp, pressure_direction)


class IncompressibleFluidChannel(Channel):
    def __init__(self, channel_dict, fluid):
        super().__init__(channel_dict, fluid)

    def update(self, mass_flow_in, dmass=None):
        super().update(mass_flow_in, dmass)
        # self.calc_mass_balance(mass_flow_in, dmass)
        # self.calc_flow_velocity()

    def calc_mass_balance(self, mass_flow_in, dmass=None):
        self.mass_flow_total[:] = mass_flow_in
        if dmass is not None:
            g_func.add_source(self.mass_flow_total, dmass, self.flow_direction)
        super().calc_mass_balance(self.temp, self.p)


class GasMixtureChannel(Channel):
    def __init__(self, channel_dict, fluid):
        super().__init__(channel_dict, fluid)
        self.inlet_composition = self.mole_flow[0]
        self.mole_flow_total = np.zeros(self.n_nodes)
        arr_shape = (self.fluid.n_species, self.n_nodes)
        self.mole_flow = np.zeros(arr_shape)
        self.mass_flow = np.zeros(arr_shape)

        self.add_print_data(self.mole_flow, 'Mole Flow',
                            'mol/s', self.fluid.species.names)

    def update(self, mol_flow_in, dmol=None):
        super().update(mol_flow_in, dmol)
        # self.calc_mass_balance(mol_flow_in, dmol)
        # self.calc_flow_velocity()

    def calc_mass_balance(self, mol_flow_in, dmol=None):
        """
        Calculate mass balance in 1D channel
        :param mol_flow_in: inlet mol flow
        :param dmol: 2D array (n_species x n_elements) of discretized molar
        source
        :return: None
        """
        self.mole_flow[:] = mol_flow_in
        if dmol is not None:
            if np.shape(dmol) == (self.fluid.n_species, self.n_ele):
                for i in range(self.fluid.n_species):
                    g_func.add_source(self.mole_flow[i], dmol[i],
                                      self.flow_direction)
            else:
                raise ValueError('Shape of dmol does not conform '
                                 'to mole_flow array')
        super().calc_mass_balance(self.temp, self.p, self.mole_flow)
        self.mole_flow_total[:] = np.sum(self.mole_flow, axis=0)
        self.mass_flow[:] = self.mole_flow * self.fluid.species.mw
        self.mass_flow_total[:] = self.mole_flow_total * self.fluid.mw


class TwoPhaseMixtureChannel(GasMixtureChannel):
    def __init__(self, channel_dict, fluid):
        super().__init__(channel_dict, fluid)
        self.mole_flow_gas_total = np.zeros(self.n_nodes)
        self.mass_flow_gas_total = np.zeros(self.n_nodes)
        self.vol_flow_gas = np.zeros(self.n_nodes)
        arr_shape = (self.fluid.n_species, self.n_nodes)
        self.mole_flow_liq = np.zeros(arr_shape)
        self.mole_flow_gas = np.zeros(arr_shape)
        self.mass_flow = np.zeros(arr_shape)
        self.mass_flow_gas = np.zeros(arr_shape)
        self.mass_flow_liq = np.zeros(arr_shape)
        self.cond_rate = np.zeros(self.n_nodes)

        self.add_print_data(self.mole_flow_gas, 'Gas Mole Flow', 'mol/s',
                            self.fluid.species.names)

    def calc_flow_velocity(self):
        """
        Calculates the gas phase velocity.
        """
        self.vol_flow_gas[:] = self.mass_flow_gas_total / self.fluid.density
        self.velocity[:] = self.vol_flow_gas / self.cross_area
        self.reynolds[:] = \
            self.velocity * self.d_h * self.fluid.density / self.fluid.viscosity

    def calc_mass_balance(self, mol_flow_in, dmol=None):
        super().calc_mass_balance(mol_flow_in, dmol)
        self.calc_two_phase_flow()

    def calc_two_phase_flow(self):
        """
        Calculates the condensed phase flow and updates mole and mass fractions
        """
        self.mole_flow_gas[:] = self.mole_flow_total * self.fluid.mole_fraction
        self.mole_flow_liq[:] = self.mole_flow_total - self.mole_flow_gas
        self.mass_flow_gas[:] = self.mass_flow_total * self.fluid.mass_fraction
        self.mass_flow_liq[:] = self.mass_flow_total - self.mass_flow_gas
        self.mole_flow_gas_total[:] = np.sum(self.mole_flow_total, axis=0)
        self.mass_flow_gas_total = np.sum(self.mass_flow_total, axis=0)
        self.calc_cond_rate()

    def calc_cond_rate(self):
        """
        Calculates the molar condensation rate of the phase change species in
        the channel.
        """
        cond_rate_ele = np.ediff1d(self.mole_flow_liq[self.fluid.id_pc])
        self.cond_rate[:] = \
            self.flow_direction * ip.interpolate_1d(cond_rate_ele,
                                                    add_edge_points=True)


