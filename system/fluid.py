import numpy as np
from numpy.polynomial.polynomial import polyval
import data.global_parameters as g_par
from abc import ABC


class LiquidProperties:

    def __init__(self, density, viscosity, specific_heat, thermal_conductivity):
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
                            'of LiquidProperties with species names as keys')
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
            self.liquid_props = LiquidProperties(np.asarray(density),
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


class Fluid(ABC):
    PROPERTY_NAMES = ['Density', 'Specific Heat', 'Viscosity',
                      'Thermal Conductivity']

    def __new__(cls, nx, species_dict, pressure_init=101325.0,
                temp_init=298.15, **kwargs):
        species_types = list(species_dict.values())
        species_types_str = ' '.join(species_types)
        if 'gas' in species_types_str and 'liquid' not in species_types_str:
            return super(Fluid, cls).__new__(GasMixture)
        elif 'gas' in species_types_str and 'liquid' in species_types_str:
            return super(Fluid, cls).__new__(TwoPhaseMixture)
        elif 'liquid' in species_types_str and 'gas' not in species_types_str:
            return super(Fluid, cls).__new__(Liquid)
        else:
            raise NotImplementedError('Only fluid types of GasMixture, '
                                      'Liquid, and TwoPhaseMixture are '
                                      'implemented and must be indicated in '
                                      'the species_dict')

    def __init__(self, nx, species_dict, pressure_init=101325.0,
                 temp_init=298.15, **kwargs):
        self.temperature = np.full(nx, temp_init)
        self.pressure = np.full(nx, pressure_init)
        self.property = dict()
        for name in self.PROPERTY_NAMES:
            self.property[name] = np.zeros(nx)


class Liquid(Fluid):
    def __init__(self, nx, species_dict, liquid_props, pressure_init=101325.0,
                 temp_init=298.15, **kwargs):
        super().__init__(nx, species_dict, pressure_init, temp_init, **kwargs)
        self.properties = liquid_props


class GasMixture(Fluid):
    def __init__(self, nx, species_dict, mole_fractions_init,
                 pressure_init=101325.0, temp_init=298.15, **kwargs):
        super().__init__(nx, species_dict, pressure_init, temp_init, **kwargs)
        print("Constructor for Gas Mixture")
        species_names = list(species_dict.keys())
        self.gas_constant = g_par.constants['R']
        self.species = GasSpecies(species_names)
        self.species_viscosity = \
            self.species.calc_viscosity(np.full(nx, temp_init))

        if len(mole_fractions_init) != len(self.species.list) \
                or np.sum(mole_fractions_init) != 1.0:
            raise ValueError('Initial mole fractions must be provided '
                             'for all species and add up to unity')
        if isinstance(mole_fractions_init, (list, tuple)):
            mole_fractions_init = np.asarray(mole_fractions_init)

        self.mole_fraction = \
            np.ones((nx, len(self.species.list))) * mole_fractions_init
        self.concentration, total_mol_conc = self.calc_dry_concentration()

        mass_fraction_init = \
            mole_fractions_init * self.species.mw / \
            np.sum(mole_fractions_init * self.species.mw)

        self.mass_fraction = \
            np.ones((nx, len(self.species.list))) * mass_fraction_init
        self.mw = np.sum(self.mole_fraction*self.species.mw, axis=-1)

    @staticmethod
    def calc_fraction(species_flow):
        """
        Calculates the species mixture fractions based on a multi-dimensional
        array with different species along the first (0th) axis.
        """
        return species_flow / np.sum(species_flow, axis=-1)

    def calc_dry_concentration(self):
        """
        Calculates the gas phase molar concentrations.
        """
        total_mol_conc = self.pressure / (self.gas_constant * self.temperature)
        return self.mole_fraction.transpose() * total_mol_conc, total_mol_conc

    def calc_mole_fraction(self, mole_flow):
        self.mole_fraction = self.calc_fraction(mole_flow)

    def calc_mass_fraction(self):
        self.mass_fraction = self.mole_fraction * self.species.mw / self.mw

    def calc_molar_mass(self):
        self.mw = np.sum(self.mole_fraction*self.species.mw, axis=-1)

    def calc_specific_heat(self, temperature):
        species_cp = self.species.calc_specific_heat(temperature).transpose()
        return np.sum(self.mass_fraction * species_cp, axis=-1)

    def calc_viscosity(self, temperature):
        """
        Calculates the mixture viscosity of a
        gas according to Herning and Zipperer.
        """
        self.species_viscosity = self.species.calc_viscosity(temperature)
        spec_visc = self.species_viscosity.transpose()
        x_sqrt_mw = self.mole_fraction * np.sqrt(self.species.mw)
        return np.sum(spec_visc * x_sqrt_mw, axis=-1)/np.sum(x_sqrt_mw, axis=-1)

    def calc_wilke_coefficients(self):
        """
        Calculates the wilke coefficients for
        each species combination of a gas.
        """
        visc = self.species_viscosity
        mw = self.species.mw
        n_species = len(self.species.list)
        alpha = []
        for i in range(n_species):
            beta = []
            for j in range(n_species):
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
        self.mole_fraction[1:] = np.maximum(1e-16, self.mole_fraction[1:])
        wilke_coeffs = self.calc_wilke_coefficients()
        lambda_species = \
            self.species.calc_thermal_conductivity(temperature, pressure)
        lambda_mix = np.zeros(len(temperature))
        for i in range(len(self.species.list)):
            a = self.mole_fraction[:, i] * lambda_species[i]
            b = np.sum(self.mole_fraction.transpose() * wilke_coeffs[i], axis=0)
            b += 1e-16
            lambda_mix += a / b
        return lambda_mix

    def calc_density(self, temperature, pressure, method="ideal gas"):
        """
        Calculate gas mixture density
        :param temperature: temperature array
        :param pressure: pressure array
        :param method: string indicating the calculation method
        :return: density of the mixture
        """
        if method == "ideal gas":
            return pressure * self.mw / (self.gas_constant * temperature)
        else:
            raise NotImplementedError('Method "' + method + '" to calculate '
                                      'gas mixture density has not '
                                      'been implemented')


class TwoPhaseMixture(GasMixture):
    def __init__(self, nx, species_dict, mole_fractions_init,
                 liquid_props=None, pressure_init=101325.0,
                 temp_init=298.15, **kwargs):
        super().__init__(nx, species_dict, mole_fractions_init,
                         pressure_init, temp_init, **kwargs)
        if not isinstance(species_dict, dict):
            raise TypeError('Argument "species_names" must be a dict '
                            'containing all species names and their '
                            'expected aggregate states in terms of "gas", '
                            '"gas-liquid", or "liquid"')
        gas_species_names = [key for key in species_dict
                             if 'gas' in species_dict[key]]
        phase_change_species_names = [key for key in species_dict
                                      if 'gas-liquid' in species_dict[key]]
        self.pc_ids = [self.species.list.index(name) for name in
                       phase_change_species_names]
        if len(self.pc_ids) > 1:
            raise NotImplementedError('At the moment only one species '
                                      'undergoing phase change is allowed')
        self.pc_id = self.pc_ids[0]

        self.liquid_fraction = np.full(nx, 0.0)
        self.phase_change_species = PhaseChangeSpecies(liquid_props)

    def calc_concentration(self, mole_flow):
        conc, total_mol_conc = self.calc_dry_concentration()
        p_sat = \
            self.phase_change_species.calc_saturation_pressure(self.temperature)
        sat_conc = p_sat / (self.gas_constant * self.temperature)
        dry_mole_flow = np.copy(mole_flow)
        dry_mole_flow[self.pc_id] = 0.0
        dry_mole_fraction = self.calc_fraction(dry_mole_flow)
        for i in range(len(self.species.list)):
            if i == self.pc_id:
                self.concentration[self.pc_id] = \
                    np.where(self.concentration[self.pc_id] > sat_conc,
                             sat_conc, self.concentration[self.pc_id])
            else:
                self.concentration[i] = \
                    np.where(self.concentration[self.pc_id] > sat_conc,
                             (total_mol_conc - sat_conc) * dry_mole_fraction[i],
                             self.concentration[i])
        self.concentration[:] = np.maximum(self.concentration, 1e-6)

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
air = Fluid(10, {'O2': 'gas', 'N2': 'gas'}, [0.21, 0.79])
print(temp)
#print(gas.calc_viscosity(temp))
#print(gas.calc_thermal_conductivity(temp, press))
print(air.calc_specific_heat(temp))
#print(gas.mw)
#print(gas.calc_density(temp, press))

liquid_water = LiquidProperties(1000.0, 1e-3, 4000.0, 0.2)
water = PhaseChangeSpecies({'H2O': liquid_water})

#print(water.gas_props.calc_specific_heat(temp))
print(water.calc_saturation_pressure(temp))
print(water.calc_vaporization_enthalpy(temp))
print(water.list)

wet_air = Fluid(10, {'O2': 'gas', 'N2': 'gas', 'H2O': 'gas-liquid'},
                mole_fractions_init=[0.205, 0.785, 0.01],
                liquid_props={'H2O': liquid_water})

print(wet_air.mole_fraction)
