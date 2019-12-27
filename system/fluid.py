import numpy as np
from numpy.polynomial.polynomial import polyval
import data.global_parameters as g_par
from abc import ABC, abstractmethod


class FluidProperties(ABC):
    def __new__(cls,  density, viscosity, specific_heat,
                thermal_conductivity, prop_type='constant', **kwargs):
        if prop_type is 'constant':
            return super(FluidProperties, cls).__new__(ConstantProperties)
        else:
            raise NotImplementedError('Provided prop_type is hast not yet '
                                      'been implemented as subclass of '
                                      'FluidProperties')

    def __init__(self, density, viscosity, specific_heat,
                 thermal_conductivity, **kwargs):
        self.density = None
        self.viscosity = None
        self.specific_heat = None
        self.thermal_conductivity = None


class ConstantProperties(FluidProperties):

    def __init__(self, density, viscosity, specific_heat, thermal_conductivity):
        super().__init__(density, viscosity, specific_heat,
                         thermal_conductivity)
        self.density = density
        self.viscosity = viscosity
        self.specific_heat = specific_heat
        self.thermal_conductivity = thermal_conductivity


class GasSpecies:
    MW = \
        {
            'O2': 0.032,
            'N2': 0.028,
            'H2': 0.002,
            'H2O': 0.018
        }
    PROPERTY_NAMES = ['Specific Heat', 'Viscosity', 'Thermal Conductivity']
    COEFFS = \
        {
            'Specific Heat':
            {
                'O2':
                np.array([-8.66817221e-19, 4.73115002e-15, -1.15709215e-11,
                          1.66145697e-08, -1.50422620e-05, 8.23507193e-03,
                          -2.12067742, 1.10964386e+03]),
                'H2':
                np.array([2.75575856e-17, -1.58350769e-13, 3.93319791e-10,
                          -5.48239691e-07, 4.61978745e-04, -2.32058478e-01,
                          6.35361636e+01, 7.25459677e+03]),
                'N2':
                np.array([-5.97839654e-18, 2.58035023e-14, -4.51777701e-11,
                          4.09644416e-08, -2.06285776e-05, 6.06476999e-03,
                          -9.88549011e-01, 1.10768881e+03]),
                'H2O':
                np.array([4.91358141e-18, -1.85174687e-14, 2.53707252e-11,
                          -1.22872163e-08, -4.19918931e-06, 7.45080766e-03,
                          -2.56601734e+00, 2.12709233e+03])
            },
            'Viscosity':
            {
                'O2':
                np.array([1.18758866e-26, -6.48183635e-23, 1.55753837e-19,
                          -2.18596122e-16, 2.02541399e-13, -1.38130567e-10,
                          1.02085148e-07, -1.50063345e-06]),
                'H2':
                np.array([6.62607149e-27, -3.42854972e-23, 7.69171320e-20,
                          -9.86354577e-17, 8.10498717e-14, -4.74743879e-11,
                          3.50618247e-08, 1.14582528e-06]),
                'N2':
                np.array([1.53556927e-26, -8.08312960e-23, 1.85403378e-19,
                          -2.44766231e-16, 2.08881853e-13, -1.27545734e-10,
                          8.53886431e-08, -3.89241700e-07]),
                'H2O':
                np.array([1.45183739e-25, -7.27081451e-22, 1.55686360e-18,
                          -1.85885301e-15, 1.35158418e-12, -6.12002995e-10,
                          2.08858479e-07, -1.90872183e-05])
            },
            'Thermal Conductivity':
            {
                'O2':
                np.array(((-1.58986421e-22, 8.03802084e-19, -1.67882604e-15,
                           1.84325862e-12, -1.11449899e-09, 3.57769046e-07,
                           2.11463976e-05, 8.31514294e-03),
                          (-1.79821722e-22, 9.07145474e-19, -1.89752639e-15,
                           2.10059733e-12, -1.29817069e-09, 4.38946254e-07,
                           -4.88815163e-07, 1.15869249e-02))),
                'H2':
                np.array(((-1.24645824e-21, 6.61764024e-18, -1.52054216e-14,
                           1.99690675e-11, -1.49509678e-08, 5.68819226e-06,
                           -4.82527146e-04, 7.12531055e-02),
                          (-1.27194170e-21, 6.74517329e-18, -1.54782877e-14,
                           2.02945953e-11, -1.51875168e-08, 5.79536974e-06,
                           -5.12222368e-04, 7.61461664e-02))),
                'N2':
                np.array(((-2.70705853e-22, 1.16336874e-18, -1.95546587e-15,
                           1.51768486e-12, -4.03713326e-10, -1.15746366e-07,
                           1.46652557e-04, -3.27592873e-03),
                          (-2.89265541e-22, 1.25583526e-18, -2.15215826e-15,
                           1.75051772e-12, -5.71060853e-10, -4.11730776e-08,
                           1.26579523e-04, -1.99639546e-04))),
                'H2O':
                np.array(((-7.69988150e-22, 3.81045861e-18, -7.88736102e-15,
                           8.77358057e-12, -5.61460795e-09, 2.11880777e-06,
                           -3.22515696e-04, 2.21293426e-02),
                          (-9.94179604e-21, 4.66529326e-17, -9.19773736e-14,
                           9.85634165e-11, -6.18830008e-08, 2.27982875e-05,
                           -4.45126170e-03, 3.69981235e-01)))
            }
        }

    def __init__(self, species_list):
        species_list = list(species_list)
        self.coeff_dict_dict = dict()
        self.coeff_dict_dict2 = dict()
        self.mw = list()
        for prop_name in self.PROPERTY_NAMES:
            self.coeff_dict_dict[prop_name] = dict()
            for name in species_list:
                if name not in self.MW:
                    raise AttributeError('No data available for '
                                         'provided specie: ' + name)
                else:
                    self.coeff_dict_dict[prop_name][name] = \
                        self.COEFFS[prop_name][name]

        for name in species_list:
            if name not in self.MW:
                raise AttributeError('No data available for '
                                     'provided specie: ' + name)
            else:
                self.mw.append(self.MW[name])
                self.coeff_dict_dict2[name] = dict()
                for prop_name in self.PROPERTY_NAMES:
                    self.coeff_dict_dict2[name][prop_name] = \
                        self.COEFFS[prop_name][name]
        self.list = species_list
        self.mw = np.asarray(self.mw)

        self.coeff_dict_arr = dict()
        for prop_name in self.PROPERTY_NAMES:
            self.coeff_dict_arr[prop_name] = \
                np.stack([np.flip(self.coeff_dict_dict[prop_name][item],
                                  axis=-1)
                          for item in self.coeff_dict_dict[prop_name]], axis=-1)

    def calc_specific_heat(self, temperature):
        return polyval(temperature, self.coeff_dict_arr['Specific Heat'])

    def calc_viscosity(self, temperature):
        return polyval(temperature, self.coeff_dict_arr['Viscosity'])

    def calc_thermal_conductivity(self, temperature, pressure):
        lambda_1_bar = \
            polyval(temperature,
                    self.coeff_dict_arr['Thermal Conductivity'][:][0])
        lambda_10_bar = \
            polyval(temperature,
                    self.coeff_dict_arr['Thermal Conductivity'][:][1])
        return lambda_1_bar + \
            (pressure - 1.e5) / 9.e5 * (lambda_10_bar - lambda_1_bar)

    def calc_property(self, property_name, temperature, pressure=101325.0):
        if property_name == 'Specific Heat':
            return self.calc_specific_heat(temperature)
        elif property_name == 'Viscosity':
            return self.calc_viscosity(temperature)
        elif property_name == 'Thermal Conductivity':
            return self.calc_thermal_conductivity(temperature, pressure)
        else:
            raise ValueError('property_name {} not valid'.format(property_name))


class PhaseChangeSpecies:
    PROPERTY_NAMES = ('Saturation Pressure', 'Vaporization Enthalpy')
    COEFFS = \
        {
            'Saturation Pressure':
            {
                'H2O':
                np.asarray([-4.66691122e-18, 2.19146750e-14, -4.56208833e-11,
                            5.54957241e-08, -4.37186346e-05, 2.33207549e-02,
                            -8.53414571e+00, 2.11600925e+03, -3.40245294e+05,
                            3.20415279e+07, -1.34211567e+09]),
            },
            'Vaporization Enthalpy':
            {
                'H2O':
                np.asarray([-2.01107201e-18, 8.92669752e-15, -1.76751771e-11,
                            2.05547260e-08, -1.55445645e-05, 7.98692642e-03,
                            -2.82333561e+00, 6.77951176e+02, -1.05826022e+05,
                            9.69666280e+06, -3.95999179e+08]),
            }
        }

    def __init__(self, liquids_dict):
        # print("Constructor of Two Phase Species")
        if not isinstance(liquids_dict, dict):
            raise TypeError('Input data must be provided as dictionary '
                            'of ConstantProperties with species names as keys')
        self.list = list(liquids_dict.keys())
        for liquid in self.list:
            if liquid not in self.COEFFS[self.PROPERTY_NAMES[0]]:
                raise NotImplementedError('No phase change data available for '
                                          'provided liquid specie: ' + liquid)
        self.gas_props = GasSpecies(self.list)
        if len(liquids_dict) == 1:
            self.liquid_props = next(iter(liquids_dict.values()))
        else:
            density = []
            specific_heat = []
            viscosity = []
            thermal_conductivity = []
            for liquid in liquids_dict:
                density.append(liquid.density)
                viscosity.append(liquid.viscosity)
                specific_heat.append(liquid.specific_heat)
                thermal_conductivity.append(liquid.thermal_conductivity)
            self.liquid_props = \
                FluidProperties(np.asarray(density),
                                np.asarray(viscosity),
                                np.asarray(specific_heat),
                                np.asarray(thermal_conductivity))

        self.coeff_dict_dict = dict()
        for prop_name in self.PROPERTY_NAMES:
            self.coeff_dict_dict[prop_name] = dict()
            for name in self.list:
                    self.coeff_dict_dict[prop_name][name] = \
                        self.COEFFS[prop_name][name]

        self.coeff_dict_arr = dict()
        for prop_name in self.PROPERTY_NAMES:
            self.coeff_dict_arr[prop_name] = \
                np.stack([np.flip(self.coeff_dict_dict[prop_name][item],
                                  axis=-1)
                          for item in self.coeff_dict_dict[prop_name]], axis=-1)

    def calc_saturation_pressure(self, temperature):
        return polyval(temperature,
                       self.coeff_dict_arr['Saturation Pressure'])

    def calc_vaporization_enthalpy(self, temperature):
        return polyval(temperature,
                       self.coeff_dict_arr['Vaporization Enthalpy'])

    def calc_property(self, property_name, temperature, **kwargs):
        if property_name == 'Saturation Pressure':
            return self.calc_saturation_pressure(temperature)
        elif property_name == 'Vaporization Enthalpy':
            return self.calc_vaporization_enthalpy(temperature)
        else:
            raise ValueError('property_name {} not valid'.format(
                             property_name))


class Fluid(ABC):
    PROPERTY_NAMES = ['Density', 'Specific Heat', 'Viscosity',
                      'Thermal Conductivity']

    def __new__(cls, nx, name, species_dict=None, pressure_init=101325.0,
                temp_init=298.15, **kwargs):
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

    def __init__(self, nx, name, species_dict=None, pressure_init=101325.0,
                 temp_init=298.15, **kwargs):
        print("__init__ for Fluid")
        self.temperature = np.full(nx, temp_init)
        self.pressure = np.full(nx, pressure_init)
        self.property = dict()
        self.n_species = None
        for name in self.PROPERTY_NAMES:
            self.property[name] = np.zeros(nx)
        self.density = self.property['Density']
        self.viscosity = self.property['Viscosity']
        self.thermal_conductivity = self.property['Thermal Conductivity']
        self.specific_heat = self.property['Specific Heat']

    def update(self, temperature, pressure, *args, **kwargs):
        self.temperature[:] = temperature
        self.pressure[:] = pressure


class IncompressibleFluid(Fluid):
    def __init__(self, nx, name, fluid_props, pressure_init=101325.0,
                 temp_init=298.15, **kwargs):
        print("__init__ for IncompressibleFluid")
        super().__init__(nx, pressure_init, temp_init, **kwargs)
        self.name = name
        self.properties = fluid_props
        if not isinstance(fluid_props, FluidProperties):
            raise TypeError('Argument fluid_props must be of type '
                            'ConstantProperties')
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
                 pressure_init=101325.0, temp_init=298.15, **kwargs):
        print("__init__ for Gas Mixture")
        super().__init__(nx, pressure_init, temp_init, **kwargs)
        self.name = name
        species_names = list(species_dict.keys())
        self.gas_constant = g_par.constants['R']
        self.species = GasSpecies(species_names)
        self.species_viscosity = \
            self.species.calc_viscosity(np.full(nx, temp_init))

        self.n_species = len(self.species.list)
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
            self.calc_properties(temp_init, pressure_init, method='ideal')
            self._concentration[:] = self.calc_dry_concentration().transpose()

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
                             'n_species')
        if isinstance(self, GasMixture):
            self.calc_mole_fraction(mole_flow)
            self.calc_molar_mass()
            self._mass_fraction[:] = \
                self.calc_mass_fraction(self._mole_fraction)
            self.calc_properties(temperature, pressure, method)
            self._concentration[:] = self.calc_dry_concentration().transpose()

    @staticmethod
    def calc_fraction(species_flow):
        """
        Calculates the species mixture fractions based on a multi-dimensional
        array with different species along the first (0th) axis.
        """
        return species_flow / np.sum(species_flow, axis=0)

    def calc_dry_concentration(self):
        """
        Calculates the gas phase molar concentrations.
        """
        total_mol_conc = self.pressure / (self.gas_constant * self.temperature)
        return self._mole_fraction.transpose() * total_mol_conc

    def calc_mole_fraction(self, mole_flow):
        self._mole_fraction[:] = self.calc_fraction(mole_flow).transpose()

    def calc_molar_mass(self):
        self.mw[:] = np.sum(self._mole_fraction * self.species.mw, axis=-1)

    def calc_mass_fraction(self, mole_fraction):
        return np.outer(1.0 / self.mw, self.species.mw) * mole_fraction

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
            b = np.sum(self._mole_fraction.transpose() * wilke_coeffs[i], axis=0)
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
                 liquid_props, pressure_init=101325.0,
                 temp_init=298.15, **kwargs):
        print("__init__ for TwoPhaseMixture")
        super().__init__(nx, name, species_dict, mole_fractions_init,
                         pressure_init, temp_init, **kwargs)
        if not isinstance(species_dict, dict):
            raise TypeError('Argument species_names must be a dict '
                            'containing all species names and their '
                            'expected aggregate states in terms of "gas", '
                            '"gas-liquid", or "liquid"')
        gas_species_names = [key for key in species_dict
                             if 'gas' in species_dict[key]]
        phase_change_species_names = [key for key in species_dict
                                      if 'gas-liquid' in species_dict[key]]
        ids_pc = [self.species.list.index(name) for name in
                  phase_change_species_names]
        if len(ids_pc) > 1:
            raise NotImplementedError('At the moment only one species '
                                      'undergoing phase change is allowed')
        self.id_pc = ids_pc[0]

        # self._mole_fraction_gas = np.zeros((nx, self.n_species))
        # self._mass_fraction_gas = np.zeros((nx, self.n_species))
        self.liquid_fraction = np.zeros(nx)
        self.humidity = np.zeros(nx)
        self.saturation_pressure = np.zeros(nx)
        self.phase_change_species = PhaseChangeSpecies(liquid_props)

    # @property
    # def mole_fraction_gas(self):
    #     return self._mole_fraction_gas.transpose()
    #
    # @property
    # def mass_fraction_gas(self):
    #     return self._mass_fraction_gas.transpose()

    def update(self, temperature, pressure, mole_flow=None, method='ideal',
               *args, **kwargs):
        super().update(temperature, pressure, mole_flow, method)
        self.saturation_pressure[:] = \
            self.phase_change_species.calc_saturation_pressure(temperature)
        self._concentration[:] = self.calc_concentration(mole_flow).transpose()
        self._mole_fraction[:] = \
            self._concentration / np.sum(self._concentration, axis=-1)
        self.calc_molar_mass()
        self._mass_fraction[:] = self.calc_mass_fraction(self._mole_fraction)
        self.calc_properties(temperature, pressure, method)
        self.calc_humidity()

    def calc_concentration(self, mole_flow):
        """
        Calculates the gas phase molar concentrations.
        :param mole_flow: total molar flow of each species (n_species x
        n_nodes)
        :return: gas phase concentration array (n_species x n_nodes)
        """
        conc = self.calc_dry_concentration()
        total_mol_conc = self.pressure / (self.gas_constant * self.temperature)
        p_sat = self.saturation_pressure
        sat_conc = p_sat / (self.gas_constant * self.temperature)
        dry_mole_flow = np.copy(mole_flow)
        dry_mole_flow[self.id_pc] = 0.0
        dry_mole_fraction = self.calc_fraction(dry_mole_flow)
        for i in range(self.n_species):
            if i == self.id_pc:
                conc[self.id_pc] = np.where(conc[self.id_pc] > sat_conc,
                                            sat_conc, conc[self.id_pc])
            else:
                conc[i] = \
                    np.where(conc[self.id_pc] > sat_conc,
                             (total_mol_conc - sat_conc) * dry_mole_fraction[i],
                             conc[i])
        return np.maximum(conc, 1e-6)

    def calc_humidity(self):
        """
        Calculates the relative humidity of the fluid.
        """
        p_sat = self.saturation_pressure
        self.humidity[:] = \
            self._mole_fraction[:, self.id_pc] * self.pressure / p_sat


species = GasSpecies(['O2', 'N2', 'H2'])

temp = np.array([[300.0, 400.0], [300.0, 400.0]])
#temp = np.array([300.0, 400.0])
#press = np.array([[100000.0, 100000.0], [500000.0, 500000.0]])
press = 101325.0
#print(species.coeff_dict_dict)
#print(species.coeff_dict_dict2)

#print(species.coeff_dict_arr['Thermal Conductivity'][0][0])
#test = species.calc_thermal_conductivity(temp, press)
#test = species.calc_specific_heat(temp)
#test = species.calc_viscosity(temp)


#print(species.coeff_dict_arr['Thermal Conductivity'][0][0])
temp = np.linspace(300, 400, 10)
#press = np.linspace(100000, 100000, 10)
air = Fluid(10, 'air', {'O2': 'gas', 'N2': 'gas'}, [0.21, 0.79])
print(temp)
#print(gas.calc_viscosity(temp))
#print(gas.calc_thermal_conductivity(temp, press))
print(air.calc_specific_heat(temp))
#print(gas.mw)
#print(gas.calc_density(temp, press))
liquid_water_props = FluidProperties(1000.0, 1e-3, 4000.0, 0.2)
water = PhaseChangeSpecies({'H2O': liquid_water_props})

#print(water.gas_props.calc_specific_heat(temp))
print(water.calc_saturation_pressure(temp))
print(water.calc_vaporization_enthalpy(temp))
print(water.list)
liquid_water = Fluid(10, 'liquid water', fluid_props=liquid_water_props)
print(type(liquid_water))

wet_air = Fluid(10, 'wet air', {'O2': 'gas', 'N2': 'gas', 'H2O': 'gas-liquid'},
                mole_fractions_init=[0.205, 0.785, 0.01],
                liquid_props={'H2O': liquid_water_props})

wet_air.update(temp, press, (1.5, 2.3, 4.0))

print(np.sum(wet_air.mole_fraction,axis=0))
