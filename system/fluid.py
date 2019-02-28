import numpy as np
from numpy.polynomial.polynomial import polyval
import data.global_parameters as g_par


class Species:
    MW = \
        {
            'O2': 32.0,
            'N2': 18.0,
            'H2': 2.0,
            'H2O': 28.0
        }
    PROPERTY_NAMES = {'Specific Heat', 'Viscosity', 'Thermal Conductivity'}
    COEFFS1 = \
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
                        self.COEFFS1[prop_name][name]

        for name in species_list:
            if name not in self.MW:
                raise AttributeError('No data available for '
                                     'provided specie: ' + name)
            else:
                self.mw.append(self.MW[name])
                self.coeff_dict_dict2[name] = dict()
                for prop_name in self.PROPERTY_NAMES:
                    self.coeff_dict_dict2[name][prop_name] = \
                        self.COEFFS1[prop_name][name]
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


class Gas:
    def __init__(self, n_ele, species_names,
                 mole_fraction_init, temp_init=298.15):
        self.species = Species(species_names)
        #self.mole_fraction = np.ones((n_ele, len(self.species.list)))
        self.property = dict()
        self.species_viscosity = \
            self.species.calc_viscosity(np.full(n_ele, temp_init))
        self.temp = np.full(n_ele, temp_init)
        for name in self.species.PROPERTY_NAMES:
            self.property[name] = np.zeros(n_ele)
        if len(mole_fraction_init) != len(self.species.list) \
                or np.sum(mole_fraction_init) != 1.0:
            raise ValueError('Initial mole fractions must be provided '
                             'for all species and add up to unity')
        if isinstance(mole_fraction_init, (list, tuple)):
            mole_fraction_init = np.asarray(mole_fraction_init)

        self.mole_fraction = \
            np.ones((n_ele, len(self.species.list))) * mole_fraction_init
        #self.mole_fraction.transpose()

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
        phi = []
        for i in range(len(self.species.list)):
            for j in range(len(self.species.list)):
                a = (1. + (visc[i] / visc[j]) ** 0.5
                     * (mol_w[j] / mol_w[i]) ** 0.25) ** 2.
                b = np.sqrt(8.) * (1. + mol_w[i] / mol_w[j]) ** 0.5
                phi.append(a / b)
        return np.asarray(phi)

    def calc_lambda_mix(self, lambdax, mol_f, visc, mol_w):
        """
        Calculates the heat conductivity of a gas mixture,
        according to Wilkes equation.
        """
        mol_f[1:] = np.minimum(1.e-20, mol_f[1:])
        psi = self.calc_psi(visc, mol_w)
        counter = 0
        outcome = 0.
        for q in range(len(visc)):
            a = mol_f[q] * lambdax[q]
            b = 1.e-20
            for w in range(len(visc)):
                b = b + mol_f[q] * psi[counter]
                counter = counter + 1
            outcome = outcome + a / b
        return outcome


p_saturation_param = \
    (-4.66691122e-18, 2.19146750e-14, -4.56208833e-11, 5.54957241e-08,
     -4.37186346e-05, 2.33207549e-02, -8.53414571e+00, 2.11600925e+03,
     -3.40245294e+05, 3.20415279e+07, -1.34211567e+09)

vaporization_enthalpy_param =\
    (-2.01107201e-18, 8.92669752e-15, -1.76751771e-11,
     2.05547260e-08, -1.55445645e-05,  7.98692642e-03,
     -2.82333561e+00,  6.77951176e+02, -1.05826022e+05,
     9.69666280e+06, -3.95999179e+08)


class Liquid:

    def __init__(self, p_sat_param, h_vap_param):
        self.p_sat_param = p_sat_param
        self.h_vap_param = h_vap_param

    def calc_p_sat(self, t_in):
        return np.polyval(self.p_sat_param, t_in)

    def calc_h_vap(self, t_in):
        return np.polyval(self.h_vap_param, t_in)


species = Species(['O2', 'N2', 'H2'])

temp = np.array([[300.0, 400.0], [300.0, 400.0]])
#temp = np.array([300.0, 400.0])
press = np.array([[100000.0, 100000.0], [500000.0, 500000.0]])
#press = 100000.0
#print(species.coeff_dict_dict)
#print(species.coeff_dict_dict2)

#print(species.coeff_dict_arr['Thermal Conductivity'][0][0])
#test = species.calc_thermal_conductivity(temp, press)
#test = species.calc_specific_heat(temp)
#test = species.calc_viscosity(temp)


#print(species.coeff_dict_arr['Thermal Conductivity'][0][0])
temp = np.linspace(300, 400, 10)
gas = Gas(10, ['O2', 'N2'], [0.21, 0.79])
print(temp)
print(gas.calc_viscosity(temp))
