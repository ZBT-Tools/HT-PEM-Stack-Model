import numpy as np
import data.global_parameters as g_par
from abc import ABC, abstractmethod
from system.output_object import OutputObject
import system.species as species


class Fluid(ABC, OutputObject):

    PROPERTY_NAMES = ['Density', 'Specific Heat', 'Viscosity',
                      'Thermal Conductivity']

    def __new__(cls, nx, name, species_dict=None, temp_init=298.15,
                pressure_init=101325.0, **kwargs):
        print("__new__ for Fluid")
        if species_dict is None:
            return super(Fluid, cls).__new__(IncompressibleFluid)
        else:
            species_types = list(species_dict.values())
            species_types_str = ' '.join(species_types)
            if 'gas' in species_types_str and 'liquid' not in species_types_str:
                return super(Fluid, cls).__new__(GasMixture)
            elif 'gas' in species_types_str and 'liquid' in species_types_str:
                return super(Fluid, cls).__new__(TwoPhaseMixture)
            elif 'liquid' in species_types_str and 'gas' \
                    not in species_types_str:
                return super(Fluid, cls).__new__(IncompressibleFluid)
            else:
                raise NotImplementedError('Only fluid types of GasMixture, '
                                          'Liquid, and TwoPhaseMixture are '
                                          'implemented and must be indicated '
                                          'in the species_dict')

    def __init__(self, nx, name, species_dict=None, temp_init=298.15,
                 press_init=101325.0, **kwargs):
        print("__init__ for Fluid")
        super().__init__()
        self.n_species = None
        try:
            if len(temp_init) == nx:
                self.temperature = temp_init
            else:
                raise ValueError('temp_init must be of length nx')
        except TypeError:
            self.temperature = np.full(nx, temp_init)
        try:
            if len(press_init) == nx:
                self.pressure = press_init
            else:
                raise ValueError('press_init must be of length nx')
        except TypeError:
            self.pressure = np.full(nx, press_init)

        self.property = dict()
        for name in self.PROPERTY_NAMES:
            self.property[name] = np.zeros(nx)
        self.density = self.property['Density']
        self.viscosity = self.property['Viscosity']
        self.thermal_conductivity = self.property['Thermal Conductivity']
        self._specific_heat = self.property['Specific Heat']

    @abstractmethod
    def update(self, temperature, pressure, *args, **kwargs):
        self.temperature[:] = temperature
        self.pressure[:] = pressure

    @property
    def specific_heat(self):
        return self._specific_heat


class IncompressibleFluid(Fluid):
    def __init__(self, nx, name, fluid_props, temp_init=298.15,
                 press_init=101325.0, **kwargs):
        print("__init__ for IncompressibleFluid")
        super().__init__(nx, temp_init, press_init, **kwargs)
        self.name = name
        if not isinstance(fluid_props, species.FluidProperties):
            raise TypeError('Argument fluid_props must be of type '
                            'FluidProperties')
        else:
            self.property['Density'][:] = fluid_props.density
            self.property['Specific Heat'][:] = fluid_props.specific_heat
            self.property['Viscosity'][:] = fluid_props.viscosity
            self.property['Thermal Conductivity'][:] = \
                fluid_props.thermal_conductivity

    def update(self, temperature, pressure, *args, **kwargs):
        super().update(temperature, pressure)


class GasMixture(Fluid):
    def __init__(self, nx, name, species_dict, mole_fractions_init,
                 temp_init=298.15, press_init=101325.0, **kwargs):
        print("__init__ for Gas Mixture")
        super().__init__(nx, temp_init, press_init, **kwargs)
        self.name = name
        species_names = list(species_dict.keys())
        self.gas_constant = g_par.constants['R']
        self.species = species.GasProperties(species_names)
        self.species_viscosity = \
            self.species.calc_viscosity(np.full(nx, temp_init))

        self.n_species = len(self.species.names)
        if isinstance(mole_fractions_init, (list, tuple)):
            mole_fractions_init = np.asarray(mole_fractions_init)
        if len(mole_fractions_init) != self.n_species \
                or np.sum(mole_fractions_init) != 1.0:
            raise ValueError('Initial mole fractions must be provided '
                             'for all species and add up to unity')

        array_shape = (nx, self.n_species)
        self._mole_fraction = \
            np.ones(array_shape) * mole_fractions_init
        self.mw = np.sum(self._mole_fraction * self.species.mw, axis=-1)
        mass_fraction_init = \
            mole_fractions_init * self.species.mw / \
            np.sum(mole_fractions_init * self.species.mw)
        self._mass_fraction = \
            np.ones(array_shape) * mass_fraction_init
        self._concentration = np.zeros(array_shape)
        if isinstance(type(self), GasMixture):
            self.calc_properties(temp_init, press_init, method='ideal')
            self._concentration[:] = self.calc_concentration().transpose()

        self.add_print_data(self.mole_fraction, 'Mole Fraction',
                            sub_names=self.species.names)

    @property
    def mole_fraction(self):
        return self._mole_fraction.transpose()

    @property
    def mass_fraction(self):
        return self._mass_fraction.transpose()

    @property
    def concentration(self):
        return self._concentration.transpose()

    def update(self, temperature, pressure, mole_flow=None, method='ideal',
               *args, **kwargs):
        super().update(temperature, pressure)
        if mole_flow is None:
            raise TypeError('Argument mole_flow must be provided')
        elif np.shape(mole_flow)[0] != self.n_species:
            raise ValueError('First dimension of mole_flow must be equal to '
                             'number of species')
        if isinstance(self, GasMixture):
            self.calc_mole_fraction(mole_flow)
            self.calc_molar_mass()
            self._mass_fraction[:] = \
                self.calc_mass_fraction(self._mole_fraction)
            self.calc_properties(self.temperature, self.pressure, method)
            self._concentration[:] = \
                self.calc_concentration().transpose()

    @staticmethod
    def calc_fraction(composition, axis=0):
        """
        Calculates the species mixture fractions based on a multi-dimensional
        array with different species along the provided axis.
        """
        return composition / np.sum(composition, axis)

    def calc_concentration(self):
        """
        Calculates the gas phase molar concentrations.
        """
        total_mol_conc = self.pressure / (self.gas_constant * self.temperature)
        return self._mole_fraction.transpose() * total_mol_conc

    def calc_mole_fraction(self, mole_flow):
        self._mole_fraction[:] = self.calc_fraction(mole_flow).transpose()

    def calc_molar_mass(self, mole_fraction=None):
        if mole_fraction is None:
            self.mw[:] = np.sum(self._mole_fraction * self.species.mw, axis=-1)
        else:
            return np.sum(mole_fraction.transpose() * self.species.mw, axis=-1)

    def calc_mass_fraction(self, mole_fraction, molar_mass=None):
        if molar_mass is None:
            return np.outer(1.0 / self.mw, self.species.mw) * mole_fraction
        else:
            return np.outer(1.0 / molar_mass, self.species.mw) * mole_fraction

    def calc_specific_heat(self, temperature):
        species_cp = self.species.calc_specific_heat(temperature).transpose()
        return np.sum(self._mass_fraction * species_cp, axis=-1)

    def calc_viscosity(self, temperature):
        """
        Calculates the mixture viscosity of a
        gas according to Herning and Zipperer.
        """
        self.species_viscosity[:] = self.species.calc_viscosity(temperature)
        spec_visc = self.species_viscosity.transpose()
        x_sqrt_mw = self._mole_fraction * np.sqrt(self.species.mw)
        return np.sum(spec_visc * x_sqrt_mw, axis=-1)/np.sum(x_sqrt_mw, axis=-1)

    def calc_wilke_coefficients(self):
        """
        Calculates the wilke coefficients for
        each species combination of a gas.
        """
        visc = self.species_viscosity
        mw = self.species.mw
        alpha = []
        for i in range(self.n_species):
            beta = []
            for j in range(self.n_species):
                a = np.power((1. + np.power((visc[i] / visc[j]), 0.5)
                             * np.power(mw[j] / mw[i], 0.25)),  2.0)
                b = np.power(np.sqrt(8.) * (1. + mw[j] / mw[i]), -0.5)
                beta.append(a/b)
            alpha.append(beta)
        return np.asarray(alpha)

    def calc_thermal_conductivity(self, temperature, pressure):
        """
        Calculates the heat conductivity of a gas mixture,
        according to Wilkes equation.
        """
        self._mole_fraction[1:] = np.maximum(1e-16, self._mole_fraction[1:])
        wilke_coeffs = self.calc_wilke_coefficients()
        lambda_species = \
            self.species.calc_thermal_conductivity(temperature, pressure)
        lambda_mix = np.zeros(len(temperature))
        for i in range(self.n_species):
            a = self._mole_fraction[:, i] * lambda_species[i]
            b = np.sum(self._mole_fraction.transpose() * wilke_coeffs[i],
                       axis=0)
            b += 1e-16
            lambda_mix += a / b
        return lambda_mix

    def calc_density(self, temperature, pressure, method="ideal"):
        """
        Calculate gas mixture density
        :param temperature: temperature array
        :param pressure: pressure array
        :param method: string indicating the calculation method
        :return: density of the mixture
        """
        if method == "ideal":
            return pressure * self.mw / (self.gas_constant * temperature)
        else:
            raise NotImplementedError('Method {} to calculate '
                                      'gas mixture density has not '
                                      'been implemented'.format(method))

    def calc_property(self, property_name, temperature, pressure=101325.0,
                      method='ideal'):
        """
        Wrapper function for the native property functions
        :param property_name: name of self.PROPERTY_NAMES to calculate
        :param temperature: 1D temperature array
        :param pressure: 1D pressure array
        :param method: method density calculation (at the moment only ideal)
        :return: the calculated 1D array of the specific property
        """
        if property_name == 'Specific Heat':
            return self.calc_specific_heat(temperature)
        elif property_name == 'Viscosity':
            return self.calc_viscosity(temperature)
        elif property_name == 'Thermal Conductivity':
            return self.calc_thermal_conductivity(temperature, pressure)
        elif property_name == 'Density':
            return self.calc_density(temperature, pressure, method)
        else:
            raise ValueError('property_name '
                             '{} not valid'.format(property_name))

    def calc_properties(self, temperature, pressure=101325.0, method='ideal'):
        """
        Wrapper function to calculate the classes properties
        :param temperature: 1D temperature array
        :param pressure: 1D pressure array
        :param method: method density calculation (at the moment only ideal)
        :return: the calculated 1D array of the specific property
        """
        for prop in self.PROPERTY_NAMES:
            self.property[prop][:] = \
                self.calc_property(prop, temperature, pressure, method)


class TwoPhaseMixture(GasMixture):
    def __init__(self, nx, name, species_dict, mole_fractions_init,
                 liquid_props=None, temp_init=298.15, press_init=101325.0,
                 **kwargs):
        print("__init__ for TwoPhaseMixture")
        if not isinstance(species_dict, dict):
            raise TypeError('Argument species_names must be a dict '
                            'containing all species names and their '
                            'expected aggregate states in terms of "gas", '
                            '"gas-liquid", or "liquid"')
        gas_species_dict = {k: v for k, v in species_dict.items() if 'gas' in v}
        super().__init__(nx, name, gas_species_dict, mole_fractions_init,
                         temp_init, press_init, **kwargs)
        phase_change_species_names = [key for key in species_dict
                                      if 'gas' and 'liquid' in
                                      species_dict[key]]
        if liquid_props is None:
            liquid_props = \
                species.FluidProperties(phase_change_species_names[0])

        ids_pc = [self.species.names.index(name) for name in
                  phase_change_species_names]
        if len(ids_pc) > 1:
            raise NotImplementedError('At the moment only one species '
                                      'undergoing phase change is allowed')
        self.id_pc = ids_pc[0]
        if liquid_props is None:
            liquid = species.FluidProperties(phase_change_species_names[0])
            self.phase_change_species = \
                species.PhaseChangeProperties({liquid.name: liquid})
        elif isinstance(liquid_props, dict):
            self.phase_change_species = species.PhaseChangeProperties(liquid_props)
        else:
            raise TypeError('Data for PhaseChangeSpecies object '
                            'can only be provided as dictionary with species '
                            'name as key and FluidProperties object as value '
                            'for the liquid properties')

        # self._mole_fraction_gas = np.zeros((nx, self.n_species))
        # self._mass_fraction_gas = np.zeros((nx, self.n_species))

        # Total properties (both phases)
        array_shape = (nx, self.n_species)
        self._mole_fraction_total = \
            np.ones(array_shape) * mole_fractions_init
        mass_fraction_init = \
            mole_fractions_init * self.species.mw / \
            np.sum(mole_fractions_init * self.species.mw)
        self._mass_fraction_total = \
            np.ones(array_shape) * mass_fraction_init

        total_specific_heat_name = 'Total Specific Heat'
        self.PROPERTY_NAMES.append(total_specific_heat_name)
        self.property[total_specific_heat_name] = np.zeros(nx)
        self.total_specific_heat = self.property[total_specific_heat_name]

        self.liquid_mass_fraction = np.zeros(nx)
        self.liquid_mole_fraction = np.zeros(nx)
        self.humidity = np.zeros(nx)
        self.saturation_pressure = np.zeros(nx)

        # Print data
        self.add_print_data(self.humidity, 'Humidity')

    # @property
    # def mole_fraction_gas(self):
    #     return self._mole_fraction_gas.transpose()
    #
    # @property
    # def mass_fraction_gas(self):
    #     return self._mass_fraction_gas.transpose()
    @property
    def mole_fraction_total(self):
        return self._mole_fraction_total.transpose()

    @property
    def mass_fraction_total(self):
        return self._mass_fraction_total.transpose()

    @property
    def specific_heat(self):
        return self.total_specific_heat

    @property
    def specific_heat_gas(self):
        return self._specific_heat

    def update(self, temperature, pressure, mole_flow=None, method='ideal',
               *args, **kwargs):
        super().update(temperature, pressure, mole_flow, method)
        if mole_flow is not None:
            if np.sum(mole_flow) > 0.0:
                self._mole_fraction_total[:] = \
                    self.calc_fraction(mole_flow).transpose()
        mw_total = self.calc_molar_mass(self.mole_fraction_total)
        self._mass_fraction_total[:] = \
            self.calc_mass_fraction(self.mole_fraction_total, mw_total)
        self.saturation_pressure[:] = \
            self.phase_change_species.calc_saturation_pressure(temperature)
        gas_conc, total_conc = self.calc_concentration()
        self._concentration[:] = gas_conc.transpose()
        self._mole_fraction[:] = \
            (self.concentration /
             np.sum(self._concentration, axis=-1)).transpose()
        self.calc_molar_mass()
        self._mass_fraction[:] = self.calc_mass_fraction(self._mole_fraction)
        self.liquid_mole_fraction[:] = \
            1.0 - np.sum(gas_conc, axis=0)/np.sum(total_conc, axis=0)
        self.liquid_mass_fraction[:] = \
            1.0 - np.sum(gas_conc * self.mw, axis=0)/np.sum(total_conc *
                                                            mw_total, axis=0)
        self.calc_properties(temperature, pressure, method)
        self.calc_humidity()

    def calc_property(self, property_name, temperature, pressure=101325.0,
                      method='ideal'):
        """
        Calculate a single property listed in PROPERTY_NAMES. Wrapper function
        for the native property functions
        :param property_name: name of self.PROPERTY_NAMES to calculate
        :param temperature: 1D temperature array
        :param pressure: 1D pressure array
        :param method: method density calculation (at the moment only ideal)
        :return: the calculated 1D array of the specific property
        """
        if property_name == 'Total Specific Heat':
            return self.calc_total_specific_heat()
        else:
            return super().calc_property(property_name, temperature, pressure,
                                         method)

    def calc_total_specific_heat(self):
        return self.liquid_mass_fraction \
            * self.phase_change_species.liquid.specific_heat \
            + (1.0 - self.liquid_mass_fraction) * self._specific_heat

    def calc_concentration(self):
        """
        Calculates the gas phase molar concentrations.
        :return: gas phase concentration and total concentration arrays
        (n_species x n_nodes)
        """
        conc = super().calc_concentration()
        all_gas_conc = np.copy(conc)
        total_mol_conc = self.pressure / (self.gas_constant * self.temperature)
        p_sat = self.saturation_pressure
        sat_conc = p_sat / (self.gas_constant * self.temperature)
        dry_mole_fraction = np.copy(self.mole_fraction_total)
        dry_mole_fraction[self.id_pc] = 0.0
        dry_mole_fraction = self.calc_fraction(dry_mole_fraction)
        for i in range(self.n_species):
            if i == self.id_pc:
                conc[self.id_pc] = np.where(conc[self.id_pc] > sat_conc,
                                            sat_conc, conc[self.id_pc])
            else:
                conc[i] = \
                    np.where(conc[self.id_pc] > sat_conc,
                             (total_mol_conc - sat_conc) * dry_mole_fraction[i],
                             conc[i])
        return np.maximum(conc, 1e-6), all_gas_conc

    def calc_humidity(self):
        """
        Calculates the relative humidity of the fluid.
        """
        p_sat = self.saturation_pressure
        self.humidity[:] = \
            self._mole_fraction[:, self.id_pc] * self.pressure / p_sat


# test_species = species.GasSpecies(['O2', 'N2', 'H2'])

# temp = np.array([[300.0, 400.0], [300.0, 400.0]])
# temp = np.array([300.0, 400.0])
# press = np.array([[100000.0, 100000.0], [500000.0, 500000.0]])
# press = 101325.0
# print(species.coeff_dict_dict)
# print(species.coeff_dict_dict2)

# print(species.coeff_dict_arr['Thermal Conductivity'][0][0])
# test = species.calc_thermal_conductivity(temp, press)
# test = species.calc_specific_heat(temp)
# test = species.calc_viscosity(temp)


# print(species.coeff_dict_arr['Thermal Conductivity'][0][0])
# temp = np.linspace(300, 400, 10)
# press = np.linspace(100000, 100000, 10)
# air = Fluid(10, 'air', {'O2': 'gas', 'N2': 'gas'},
#             mole_fractions_init=[0.21, 0.79])
# print(temp)
# print(gas.calc_viscosity(temp))
# print(gas.calc_thermal_conductivity(temp, press))
# print(air.calc_specific_heat(temp))
# print(gas.mw)
# print(gas.calc_density(temp, press))
# liquid_water_props = species.FluidProperties('H2O')
# water = species.PhaseChangeSpecies({'H2O': liquid_water_props})

# print(water.gas_props.calc_specific_heat(temp))
# print(water.calc_saturation_pressure(temp))
# print(water.calc_vaporization_enthalpy(temp))
# print(water.names)
# liquid_water = Fluid(10, 'liquid water', fluid_props=liquid_water_props)
# print(type(liquid_water))
#
# wet_air = Fluid(10, 'wet air', {'O2': 'gas', 'N2': 'gas', 'H2O': 'gas-liquid'},
#                 mole_fractions_init=[0.205, 0.785, 0.01],
#                 liquid_props={'H2O': liquid_water_props})
#
# wet_air.update(temp, press, (1.5, 2.3, 4.0))
#
# print(wet_air.mole_fraction)
