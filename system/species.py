import numpy as np
from numpy.polynomial.polynomial import polyval
from abc import ABC


class FluidProperties(ABC):

    PROPERTY_NAMES = ['Density', 'Specific Heat', 'Viscosity',
                      'Thermal Conductivity']
    ATTRIBUTE_NAMES = [name.replace(' ', '_').lower()
                       for name in PROPERTY_NAMES]

    def __new__(cls,  name, density=None, viscosity=None, specific_heat=None,
                thermal_conductivity=None, prop_type='constant', **kwargs):
        if prop_type is 'constant':
            return super(FluidProperties, cls).__new__(ConstantProperties)
        else:
            raise NotImplementedError('Provided prop_type has not yet '
                                      'been implemented')

    def __init__(self, name, **kwargs):
        self.name = name
        for attr_name in self.ATTRIBUTE_NAMES:
            setattr(self, attr_name, None)


class ConstantProperties(FluidProperties):
    COEFFS = \
        {
            'Density':
                {
                    'H2O': 1000.0
                },
            'Specific Heat':
                {
                    'H2O': 4000.0
                },
            'Viscosity':
                {
                    'H2O': 1e-3
                },
            'Thermal Conductivity':
                {
                    'H2O': 0.2
                }
        }

    def __init__(self, name, **kwargs):
        super().__init__(name, **kwargs)
        for i, attr_name in enumerate(self.ATTRIBUTE_NAMES):
            prop = kwargs.get(attr_name, None)
            if name in self.COEFFS[self.PROPERTY_NAMES[i]] and prop is None:
                setattr(self, attr_name,
                        self.COEFFS[self.PROPERTY_NAMES[i]][self.name])
            elif prop is not None:
                setattr(self, attr_name, prop)
            else:
                raise ValueError('Value for {} must be in database or '
                                 'directly provided'.format(attr_name))


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
        self.names = species_list
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
    PROPERTY_NAMES = ['Saturation Pressure', 'Vaporization Enthalpy']
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
                            'with species names as keys and objects'
                            'type FluidProperties as values')
        self.names = list(liquids_dict.keys())
        for name in self.names:
            if name not in self.COEFFS[self.PROPERTY_NAMES[0]]:
                raise NotImplementedError('No phase change data available for '
                                          'provided liquid specie: ' + name)
        self.gas = GasSpecies(self.names)
        if len(liquids_dict) == 1:
            self.liquid = next(iter(liquids_dict.values()))
        else:
            density = [value.density for key, value in liquids_dict]
            specific_heat = [value.specific_heat for key, value in liquids_dict]
            viscosity = [value.viscosity for key, value in liquids_dict]
            thermal_conductivity = \
                [value.thermal_conductivity for key, value in liquids_dict]
            self.liquid = \
                FluidProperties('liquids', density=np.asarray(density),
                                viscosity=np.asarray(viscosity),
                                specific_heat=np.asarray(specific_heat),
                                thermal_conductivity=
                                np.asarray(thermal_conductivity))

        self.coeff_dict_dict = dict()
        for prop_name in self.PROPERTY_NAMES:
            self.coeff_dict_dict[prop_name] = dict()
            for name in self.names:
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
