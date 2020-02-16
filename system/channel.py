import numpy as np
import data.global_parameters as g_par
import system.fluid2 as fluids
import system.global_functions as g_func
import system.interpolation as ip
from abc import ABC, abstractmethod
from system.output_object import OutputObject


class Channel(ABC, OutputObject):
    def __new__(cls, channel_dict, fluid):

        if type(fluid) is fluids.IncompressibleFluid \
                or type(fluid) is fluids.ConstantFluid:
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
        self.length = channel_dict['length']
        # channel length
        self.n_nodes = len(self.fluid.density)
        self.n_ele = self.n_nodes - 1

        self.x = np.linspace(0.0, self.length, self.n_nodes)
        self.dx = np.diff(self.x)
        # element length
        self.p_out = channel_dict['p_out']
        self.p = np.zeros(self.n_nodes)
        self.p.fill(self.p_out)
        # outlet and initial pressure
        self.temp_in = channel_dict['temp_in']
        self.temp = np.zeros(self.p.shape)
        self.temp.fill(self.temp_in)
        # inlet temperature

        self.tri_mtx = None
        self.id_in = None
        self.id_out = None
        self.flow_direction = channel_dict['flow_direction']

        self.pressure_recovery = False

        # Geometry
        self.width = channel_dict['channel_width']
        self.cross_shape = channel_dict.get('cross_sectional_shape',
                                            'rectangular')
        if self.cross_shape == 'rectangular':
            self.height = channel_dict['channel_height']
            self.cross_area = self.width * self.height
            self.d_h = 4. * self.cross_area / (2. * (self.width + self.height))
        elif self.cross_shape == 'circular':
            self.d_h = channel_dict['channel_diameter']
            self.cross_area = self.d_h ** 2.0 / 4.0 * np.pi
        else:
            raise NotImplementedError
        self.n_bends = channel_dict.get('bend_number', 0)
        self.zeta_bends = channel_dict.get('bend_friction_factor', 0.0)
        self.zeta_other = \
            channel_dict.get('additional_friction_fractor', 0.0)
        self.friction_factor = np.zeros(self.n_ele)
        self.base_area = self.width * self.length
        self.base_area_dx = self.width * self.dx

        # Flow
        self.velocity = np.zeros(self.n_nodes)
        self.reynolds = np.zeros(self.n_nodes)
        self.mass_flow_total = np.zeros(self.n_nodes)
        self.vol_flow = np.zeros(self.n_nodes)
        self.g_fluid = np.zeros(self.n_nodes)

        # Heat Transfer
        self.k_coeff = np.zeros(self.n_ele)

        self.add_print_data(self.temp, 'Fluid Temperature', 'K')
        self.add_print_data(self.p, 'Fluid Pressure', 'Pa')

    @abstractmethod
    def update(self, *args, **kwargs):
        pass

    def calc_flow_velocity(self):
        """
        Calculates the gas phase velocity.
        """
        self.vol_flow[:] = self.mass_flow_total / self.fluid.density
        self.velocity[:] = np.maximum(self.vol_flow / self.cross_area, 0.0)
        self.reynolds[:] = self.velocity * self.d_h * self.fluid.density \
            / self.fluid.viscosity

    @property
    def flow_direction(self):
        return self._flow_direction

    @flow_direction.setter
    def flow_direction(self, flow_direction):
        if flow_direction not in (-1, 1):
            raise ValueError('Member variable flow_direction '
                             'must be either 1 or -1')
        self._flow_direction = flow_direction
        if self._flow_direction == 1:
            self.id_in = 0
            self.id_out = -1
        else:
            self.id_in = -1
            self.id_out = 0
        ones = np.zeros((self.n_ele, self.n_ele))
        ones.fill(1.0)
        if self._flow_direction == 1:
            self.tri_mtx = np.tril(ones)
        else:
            self.tri_mtx = np.triu(ones)

    @abstractmethod
    def calc_mass_balance(self, *args, **kwargs):
        pass

    def calc_heat_transfer_coeff(self):
        """
        Calculates heat transfer coefficient of channel assuming element-wise
        constant wall temperature
        (Correlations should be reviewed)
        """
        prandtl = self.fluid.viscosity * self.fluid.specific_heat / \
            self.fluid.thermal_conductivity

        d_by_l = self.d_h / self.length
        sqrt_re_pr_dbyl = np.sqrt(self.reynolds * prandtl * d_by_l)
        nu_1 = 3.66
        nu_2 = 1.66 * sqrt_re_pr_dbyl
        nu_3 = (2. / (1. + 22. * prandtl)) ** 0.166667 * sqrt_re_pr_dbyl
        nu_lam = (nu_1 ** 3. + 0.7 ** 3. + (nu_2 - 0.7) ** 3.
                  + nu_3 ** 3.) ** 0.333333

        if np.any(self.reynolds >= 2300.0):
            zeta = \
                (1.8 * np.log(self.reynolds,
                              where=self.reynolds != 0) - 1.5) ** -2.
            nu_turb = zeta / 8. * self.reynolds * prandtl \
                / (1. + 12.7 * np.sqrt(zeta / 8.)
                   * (prandtl ** 0.666667) - 1.) * (1. + d_by_l ** 0.666667)
            nusselt = \
                np.where(self.reynolds < 2300.0, nu_lam,
                         np.where(self.reynolds < 1e4,
                                  (1. - (self.reynolds - 2300.) / 7700.) * nu_lam
                                  + (self.reynolds - 2300.) / 7700. * nu_turb,
                                  nu_turb))
        else:
            nusselt = nu_lam

        ht_coeff = nusselt * self.fluid.thermal_conductivity / self.d_h
        # convection coefficient between the coolant and the channel wall
        conv_area = self.d_h * np.pi * self.length / self.n_ele
        # convection area of the channel wall
        self.k_coeff[:] = ip.interpolate_1d(ht_coeff) * conv_area

    def calc_heat_capacitance(self):
        self.g_fluid[:] = self.mass_flow_total * self.fluid.specific_heat

    def calc_pressure(self):
        """
        Calculates the static channel pressure
        """
        density_ele = ip.interpolate_1d(self.fluid.density)
        # velocity_ele = ip.interpolate_1d(self.velocity)
        reynolds_ele = ip.interpolate_1d(self.reynolds)
        zeta_bends = self.zeta_bends * self.n_bends / self.n_ele
        zeta = zeta_bends + self.zeta_other
        self.friction_factor[:] = g_func.calc_friction_factor(reynolds_ele)
        # friction_factor = 64.0 / reynolds_ele
        dp = g_func.calc_pressure_drop(self.velocity, density_ele,
                                       self.friction_factor, zeta, self.dx,
                                       self.d_h, self.pressure_recovery)
        # dp = (f_ele * self.dx / self.d_h + zeta_bends) \
        #     * density_ele * 0.5 * velocity_ele ** 2.0
        pressure_direction = -self._flow_direction
        self.p.fill(self.p_out)
        g_func.add_source(self.p, dp, pressure_direction)


class IncompressibleFluidChannel(Channel):
    def __init__(self, channel_dict, fluid):
        super().__init__(channel_dict, fluid)
        self.mass_source = np.zeros(self.n_ele)

    def update(self, mass_flow_in=None, mass_source=None):
        self.calc_mass_balance(mass_flow_in, mass_source)
        self.fluid.update(self.temp, self.p)
        self.calc_flow_velocity()
        self.calc_pressure()
        self.calc_heat_transfer_coeff()
        self.calc_heat_capacitance()

    def calc_mass_balance(self, mass_flow_in=None, mass_source=None):
        if mass_flow_in is not None:
            self.mass_flow_total[:] = mass_flow_in
        if mass_source is not None:
            self.mass_source[:] = mass_source
        if mass_source is not None:
            g_func.add_source(self.mass_flow_total, self.mass_source,
                              self.flow_direction, self.tri_mtx)


class GasMixtureChannel(Channel):
    def __init__(self, channel_dict, fluid):
        super().__init__(channel_dict, fluid)
        self.mole_flow_total = np.zeros(self.n_nodes)
        arr_shape = (self.fluid.n_species, self.n_nodes)
        self.mole_flow = np.zeros(arr_shape)
        self.mass_flow = np.zeros(arr_shape)
        self.mass_source = np.zeros((self.fluid.n_species, self.n_ele))
        self.mole_source = np.zeros((self.fluid.n_species, self.n_ele))

        self.add_print_data(self.mole_flow, 'Mole Flow',
                            'mol/s', self.fluid.species.names)

    def update(self, mass_flow_in=None, mass_source=None):
        self.calc_mass_balance(mass_flow_in, mass_source)
        # if mass_flow_in is None:
        #     mass_flow_in = self.mass_flow[:, self.id_in]
        self.fluid.update(self.temp, self.p, self.mole_flow)
        self.calc_flow_velocity()
        self.calc_pressure()
        self.calc_heat_transfer_coeff()
        self.calc_heat_capacitance()
        # self.calc_mass_balance(mol_flow_in, dmol)
        # self.calc_flow_velocity()

    def calc_mass_balance(self, mass_flow_in=None, mass_source=None):
        """
        Calculate mass balance in 1D channel
        :param mass_flow_in: inlet mass flow
        :param mass_source: 2D array (n_species x n_elements) of discretized
                            mass source
        :return: None
        """
        if mass_flow_in is not None:
            mass_flow_in = np.asarray(mass_flow_in)
            if mass_flow_in.shape == self.mass_flow.shape:
                # mass_flow = mass_flow_in
                mass_flow_in = mass_flow_in[:, self.id_in]
                mass_flow = \
                    g_func.fill_transposed(mass_flow_in, self.mass_flow.shape)
            elif mass_flow_in.shape == self.mass_flow_total.shape:
                mass_flow = np.outer(self.fluid.mass_fraction[:, self.id_in],
                                     mass_flow_in)
            elif mass_flow_in.shape == (self.mass_flow.shape[0],):
                mass_flow = \
                    g_func.fill_transposed(mass_flow_in, self.mass_flow.shape)
            elif np.ndim(mass_flow_in) == 0:
                mass_flow_total = np.zeros(self.mass_flow_total.shape)
                mass_flow_total[:] = mass_flow_in
                mass_flow = np.outer(self.fluid.mass_fraction[:, self.id_in],
                                     mass_flow_total)
            else:
                raise ValueError(
                    'provided mass flow cannot be converted to array'
                    ' of shape: ', self.mass_flow.shape)
            self.mass_flow[:] = mass_flow
            # self.mole_flow[:] = \
            #     (mass_flow.transpose() / self.fluid.species.mw).transpose()
        if mass_source is not None:
            if np.shape(mass_source) == (self.fluid.n_species, self.n_ele):
                self.mass_source[:] = mass_source
                self.mole_source[:] = \
                    (mass_source.transpose()
                     / self.fluid.species.mw).transpose()
            else:
                raise ValueError('shape of mass_source does not conform '
                                 'to mole_source array')
        for i in range(self.fluid.n_species):
            g_func.add_source(self.mass_flow[i], self.mass_source[i],
                              self.flow_direction, self.tri_mtx)

        self.mass_flow_total[:] = np.sum(self.mass_flow, axis=0)
        # self.mass_flow[:] = \
        #     (self.mole_flow.transpose() * self.fluid.species.mw).transpose()
        self.mole_flow[:] = \
            (self.mass_flow.transpose() / self.fluid.species.mw).transpose()
        self.mole_flow_total[:] = np.sum(self.mole_flow, axis=0)


class TwoPhaseMixtureChannel(GasMixtureChannel):
    def __init__(self, channel_dict, fluid):
        super().__init__(channel_dict, fluid)
        self.mole_flow_gas_total = np.zeros(self.n_nodes)
        self.mass_flow_gas_total = np.zeros(self.n_nodes)
        self.vol_flow_gas = np.zeros(self.n_nodes)
        arr_shape = (self.fluid.n_species, self.n_nodes)
        self.mole_flow_liq = np.zeros(arr_shape)
        self.mole_flow_gas = np.zeros(arr_shape)
        self.mass_flow_gas = np.zeros(arr_shape)
        self.mass_flow_liq = np.zeros(arr_shape)
        self.cond_rate = np.zeros(self.n_nodes)

        self.add_print_data(self.mole_flow_gas, 'Gas Mole Flow', 'mol/s',
                            self.fluid.species.names)

    def update(self, mass_flow_in=None, mass_source=None):
        self.calc_mass_balance(mass_flow_in, mass_source)
        self.fluid.update(self.temp, self.p, self.mole_flow)
        self.calc_two_phase_flow()
        self.calc_flow_velocity()
        self.calc_pressure()
        self.calc_heat_transfer_coeff()
        self.calc_heat_capacitance()

    def calc_flow_velocity(self):
        """
        Calculates the gas phase velocity.
        """
        self.vol_flow_gas[:] = self.mass_flow_gas_total / self.fluid.density
        self.vol_flow[:] = self.vol_flow_gas
        self.velocity[:] = np.maximum(self.vol_flow_gas / self.cross_area, 0.0)
        self.reynolds[:] = self.velocity * self.d_h * self.fluid.density \
            / self.fluid.viscosity

    def calc_mass_balance(self, mass_flow_in=None, mass_source=None):
        super().calc_mass_balance(mass_flow_in, mass_source)
        # self.calc_two_phase_flow()

    def calc_two_phase_flow(self):
        """
        Calculates the condensed phase flow and updates mole and mass fractions
        """
        mole_flow_liq_total = self.mole_flow_total \
            * self.fluid.liquid_mole_fraction
        mass_flow_liq_total = self.mass_flow_total \
            * self.fluid.liquid_mass_fraction
        self.mole_flow_gas_total[:] = \
            self.mole_flow_total - mole_flow_liq_total
        self.mass_flow_gas_total[:] = \
            self.mass_flow_total - mass_flow_liq_total
        self.mole_flow_gas[:] = \
            self.mole_flow_gas_total * self.fluid.mole_fraction
        self.mass_flow_gas[:] = \
            self.mass_flow_gas_total * self.fluid.mass_fraction
        self.mole_flow_liq[:] = self.mole_flow - self.mole_flow_gas
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

    def calc_heat_capacitance(self):
        self.g_fluid[:] = self.mass_flow_total * self.fluid.specific_heat
