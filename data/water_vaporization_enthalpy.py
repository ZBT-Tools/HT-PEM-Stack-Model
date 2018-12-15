import numpy as np


vaporization_enthalpy_param =\
    (-2.01107201e-18, 8.92669752e-15, -1.76751771e-11,
     2.05547260e-08, -1.55445645e-05,  7.98692642e-03,
     -2.82333561e+00,  6.77951176e+02, -1.05826022e+05,
     9.69666280e+06, -3.95999179e+08)


class Fluid:

    def __init__(self, h_vap_param):
        self.h_vap_param = h_vap_param

    def calc_h_vap(self, t_in):
        return np.polyval(self.h_vap_param, t_in)


water = Fluid(vaporization_enthalpy_param)
