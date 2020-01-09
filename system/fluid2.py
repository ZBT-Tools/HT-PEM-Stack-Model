import numpy as np
import data.global_parameters as g_par
from abc import ABC, abstractmethod
from system.output_object import OutputObject
import system.species as species


class Fluid(ABC, OutputObject):

    PROPERTY_NAMES = ['Density', 'Specific Heat', 'Viscosity',
                      'Thermal Conductivity']

    def __init__(self, nx, name, temperature=298.15, pressure=101325.0,
                 **kwargs):
        print("__init__ for Fluid")
        super().__init__()
        self.name = name
        try:
            if len(temperature) == nx:
                self.temperature = temperature
            else:
                raise ValueError('Argument temperature must be of length nx')
        except TypeError:
            self.temperature = np.full(nx, temperature)
        try:
            if len(pressure) == nx:
                self.pressure = pressure
            else:
                raise ValueError('Argument pressure must be of length nx')
        except TypeError:
            self.pressure = np.full(nx, pressure)

        self.property = dict()
        for name in self.PROPERTY_NAMES:
            self.property[name] = np.zeros(nx)

    @abstractmethod
    def update(self, temperature, pressure, *args, **kwargs):
        self.temperature[:] = temperature
        self.pressure[:] = pressure

    @abstractmethod
    def calc_properties(self, temperature, pressure=101325.0, **kwargs):
        pass

    @staticmethod
    def calc_fraction(composition, axis=0):
        """
        Calculates mixture fractions based on a multi-dimensional
        array with different species along the provided axis.
        """
        return composition / np.sum(composition, axis)

    @property
    def density(self):
        return self.property['Density']

    @property
    def viscosity(self):
        return self.property['Viscosity']

    @property
    def thermal_conductivity(self):
        return self.property['Thermal Conductivity']

    @property
    def specific_heat(self):
        return self.property['Specific Heat']


class IncompressibleFluid(Fluid):
    def __init__(self, nx, name, fluid_props, temperature=298.15,
                 pressure=101325.0, **kwargs):
        print("__init__ for IncompressibleFluid")
        super().__init__(nx, temperature, pressure, **kwargs)
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

    def calc_properties(self, temperature, pressure=101325.0, **kwargs):
        pass


class GasMixture(Fluid):
    def __init__(self, nx, name, species_dict, mole_fractions,
                 temperature=298.15, pressure=101325.0, **kwargs):
        print("__init__ for Gas Mixture")
        super().__init__(nx, temperature, pressure, **kwargs)
        self.name = name
        species_names = list(species_dict.keys())
        self.n_species = len(species_names)
        self.gas_constant = g_par.constants['R']
        self.species = species.GasSpecies(species_names)
        self.species_viscosity = \
            self.species.calc_viscosity(np.full(nx, temperature))

        self.n_species = len(self.species.names)
        if isinstance(mole_fractions, (list, tuple)):
            mole_fractions = np.asarray(mole_fractions)
        if len(mole_fractions) != self.n_species \
                or np.sum(mole_fractions) != 1.0:
            raise ValueError('Initial mole fractions must be provided '
                             'for all species and add up to unity')

        array_shape = (nx, self.n_species)
        self._mole_fraction = \
            np.ones(array_shape) * mole_fractions
        self.mw = np.sum(self._mole_fraction * self.species.mw, axis=-1)
        mass_fraction_init = \
            mole_fractions * self.species.mw / \
            np.sum(mole_fractions * self.species.mw)
        self._mass_fraction = \
            np.ones(array_shape) * mass_fraction_init
        self._concentration = np.zeros(array_shape)
        if isinstance(type(self), GasMixture):
            self.calc_properties(temperature, pressure, method='ideal')
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

    @concentration.setter
    def concentration(self, value):
        self._concentration[:] = value.transpose()

    def update(self, temperature, pressure, composition=None, method='ideal',
               *args, **kwargs):
        super().update(temperature, pressure)
        if composition is None:
            composition = self.mole_fraction
        elif np.shape(composition)[0] != self.n_species:
            raise ValueError('First dimension of composition must be equal to '
                             'number of species')
        self.calc_mole_fraction(composition)
        self.calc_molar_mass()
        self._mass_fraction[:] = \
            self.calc_mass_fraction(self._mole_fraction)
        self.calc_properties(self.temperature, self.pressure, method)
        self._concentration[:] = \
            self.calc_concentration().transpose()

    def calc_concentration(self):
        """
        Calculates the gas phase molar concentrations.
        """
        total_mol_conc = self.pressure / (self.gas_constant * self.temperature)
        return self._mole_fraction.transpose() * total_mol_conc

    def calc_mole_fraction(self, composition):
        self._mole_fraction[:] = self.calc_fraction(composition).transpose()

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


class TwoPhaseMixture(Fluid):
    def __init__(self, nx, name, species_dict, mole_fractions,
                 liquid_props=None, temperature=298.15, pressure=101325.0,
                 **kwargs):
        print("__init__ for TwoPhaseMixture")
        super().__init__(nx, name, temperature, pressure, **kwargs)
        if not isinstance(species_dict, dict):
            raise TypeError('Argument species_dict must be a dict '
                            'containing all species names and their '
                            'expected aggregate states in terms of "gas", '
                            '"gas-liquid", or "liquid"')
        gas_species_dict = {k: v for k, v in species_dict.items() if 'gas' in v}
        phase_change_species_names = \
            [key for key in species_dict
             if 'gas' and 'liquid' in species_dict[key]]

        if liquid_props is None:
            liquid_props = \
                species.FluidProperties(phase_change_species_names[0])
            self.phase_change_species = \
                species.PhaseChangeSpecies({liquid_props.name: liquid_props})
        elif isinstance(liquid_props, dict):
            self.phase_change_species = species.PhaseChangeSpecies(liquid_props)
        else:
            raise TypeError('Data for PhaseChangeSpecies object '
                            'can only be provided as dictionary with species '
                            'name as key and FluidProperties object as value '
                            'for the liquid properties')
        self.liquid = IncompressibleFluid(nx, name, fluid_props=liquid_props,
                                          temperature=self.temperature,
                                          pressure=self.pressure)
        self.gas = GasMixture(nx, name, species_dict=gas_species_dict,
                              mole_fractions=mole_fractions,
                              temperature=self.temperature,
                              pressure=self.pressure)
        ids_pc = [self.species.names.index(name) for name in
                  phase_change_species_names]
        if len(ids_pc) > 1:
            raise NotImplementedError('At the moment only one species '
                                      'undergoing phase change is allowed')
        self.id_pc = ids_pc[0]
        all_ids = np.array(list(range(len(self.species.names))), dtype='int32')
        ids_no_pc = np.delete(all_ids, ids_pc)
        self.ids_no_pc = list(ids_no_pc)

        # Total properties (both phases)
        self.n_species = len(species_dict)
        array_shape = (nx, self.n_species)
        self._mole_fraction = \
            np.ones(array_shape) * mole_fractions
        mass_fraction_init = \
            mole_fractions * self.gas.species.mw / \
            np.sum(mole_fractions * self.gas.species.mw)
        self._mass_fraction = \
            np.ones(array_shape) * mass_fraction_init
        self.mw = np.zeros(nx)

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
    def mole_fraction(self):
        return self._mole_fraction.transpose()

    @property
    def mass_fraction(self):
        return self._mass_fraction.transpose()

    @property
    def species(self):
        return self.gas.species

    def update(self, temperature, pressure, mole_flow=None, method='ideal',
               *args, **kwargs):
        super().update(temperature, pressure)
        if mole_flow is not None:
            if np.sum(mole_flow) > 0.0:
                self._mole_fraction[:] = \
                    self.calc_fraction(mole_flow).transpose()

        self.mw[:] = self.gas.calc_molar_mass(self.mole_fraction)
        self._mass_fraction[:] = \
            self.gas.calc_mass_fraction(self._mole_fraction, self.mw)
        self.saturation_pressure[:] = \
            self.phase_change_species.calc_saturation_pressure(temperature)
        gas_conc = self.calc_concentration()
        self.gas.update(temperature, pressure, gas_conc, method)
        self.liquid.update(temperature, pressure)
        dry_conc = np.copy(gas_conc)
        dry_conc[self.id_pc] = 0.0
        total_conc = dry_conc
        total_conc[self.id_pc] = np.sum(dry_conc, axis=0) \
            * self.mole_fraction[self.id_pc] \
            / (1.0 - self.mole_fraction[self.id_pc])
        self.liquid_mole_fraction[:] = \
            1.0 - np.sum(gas_conc, axis=0)/np.sum(total_conc, axis=0)
        self.liquid_mass_fraction[:] = \
            1.0 - np.sum(gas_conc * self.gas.mw, axis=0)/np.sum(total_conc *
                                                                self.mw, axis=0)

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
        if property_name == 'Specific Heat':
            return self.calc_specific_heat()
        else:
            return self.gas.property[property_name]

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

    def calc_specific_heat(self):
        return self.liquid_mass_fraction \
            * self.liquid.specific_heat \
            + (1.0 - self.liquid_mass_fraction) * self.gas.specific_heat

    def calc_concentration(self):
        """
        Calculates the gas phase molar concentrations.
        :return: gas phase concentration and total concentration arrays
        (n_species x n_nodes)
        """
        r_t = self.gas.gas_constant * self.temperature
        total_gas_conc = self.pressure / r_t
        conc = self.mole_fraction * total_gas_conc
        all_conc = np.copy(conc)
        sat_conc = self.saturation_pressure / r_t
        dry_mole_fraction = np.copy(self.mole_fraction)
        dry_mole_fraction[self.id_pc] = 0.0
        dry_mole_fraction = self.calc_fraction(dry_mole_fraction)
        for i in range(self.n_species):
            if i == self.id_pc:
                conc[self.id_pc] = np.where(conc[self.id_pc] > sat_conc,
                                            sat_conc, conc[self.id_pc])
            else:
                conc[i] = \
                    np.where(conc[self.id_pc] > sat_conc,
                             (total_gas_conc - sat_conc) * dry_mole_fraction[i],
                             conc[i])
        return np.maximum(conc, 0.0)

    # def calc_concentration(self):
    #     """
    #     Calculates the gas phase molar concentrations.
    #     :return: gas phase concentration and total concentration arrays
    #     (n_species x n_nodes)
    #     """
    #     r_t = self.gas.gas_constant * self.temperature
    #     total_gas_conc = self.pressure / r_t
    #     conc = self.mole_fraction * total_gas_conc
    #
    #     all_gas_conc = np.copy(conc)
    #     sat_conc = self.saturation_pressure / r_t
    #     no_pc_total_gas_conc = total_gas_conc - sat_conc
    #     conc[self.id_pc] = \
    #         np.where(conc[self.id_pc] > sat_conc, sat_conc, conc[self.id_pc])
    #     liquid_conc = all_gas_conc - conc
    #     gas_mole_fraction = self.calc_fraction(conc)
    #     conc = gas_mole_fraction * total_gas_conc
    #     print(conc)
    #     print(all_gas_conc)
    #     return np.maximum(conc, 0.0), all_gas_conc + liquid_conc

    def calc_humidity(self):
        """
        Calculates the relative humidity of the fluid.
        """
        p_sat = self.saturation_pressure
        self.humidity[:] = \
            self._mole_fraction[:, self.id_pc] * self.pressure / p_sat


def fluid_factory(nx, name, liquid_props=None, species_dict=None,
                  mole_fractions=None, temperature=298.15,
                  pressure=101325.0, **kwargs):
    if species_dict is None:
        return IncompressibleFluid(nx, name, liquid_props, temperature,
                                   pressure)
    else:
        species_types = list(species_dict.values())
        species_types_str = ' '.join(species_types)
        if 'gas' in species_types_str and 'liquid' not in species_types_str:
            return GasMixture(nx, name, species_dict, mole_fractions,
                              temperature, pressure)
        elif 'gas' in species_types_str and 'liquid' in species_types_str:
            return TwoPhaseMixture(nx, name, species_dict, mole_fractions,
                                   liquid_props, temperature, pressure)
        elif 'liquid' in species_types_str and 'gas' \
                not in species_types_str:
            return IncompressibleFluid(nx, name, liquid_props, temperature,
                                       pressure)
        else:
            raise NotImplementedError('Only fluid types of GasMixture, '
                                      'Liquid, and TwoPhaseMixture are '
                                      'implemented and must be indicated '
                                      'in the species_dict')
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
