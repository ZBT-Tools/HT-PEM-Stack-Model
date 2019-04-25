import numpy as np


p_saturation_param = \
    (-4.66691122e-18, 2.19146750e-14, -4.56208833e-11, 5.54957241e-08,
     -4.37186346e-05, 2.33207549e-02, -8.53414571e+00, 2.11600925e+03,
     -3.40245294e+05, 3.20415279e+07, -1.34211567e+09)

vaporization_enthalpy_param =\
    (-2.01107201e-18, 8.92669752e-15, -1.76751771e-11,
     2.05547260e-08, -1.55445645e-05,  7.98692642e-03,
     -2.82333561e+00,  6.77951176e+02, -1.05826022e+05,
     9.69666280e+06, -3.95999179e+08)


class TwoPhaseSpecies:

    def __init__(self, p_sat_param, h_vap_param):
        self.p_sat_param = p_sat_param
        self.h_vap_param = h_vap_param

    def calc_p_sat(self, t_in):
        return np.polyval(self.p_sat_param, t_in)

    def calc_h_vap(self, t_in):
        return np.polyval(self.h_vap_param, t_in)


water = TwoPhaseSpecies(p_saturation_param, vaporization_enthalpy_param)
