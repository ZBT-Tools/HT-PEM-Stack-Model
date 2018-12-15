import numpy as np

# build a common class structure with the water vap ent and the gas properties
p_saturation_param = \
    (-4.66691122e-18, 2.19146750e-14, -4.56208833e-11, 5.54957241e-08,
     -4.37186346e-05, 2.33207549e-02, -8.53414571e+00, 2.11600925e+03,
     -3.40245294e+05, 3.20415279e+07, -1.34211567e+09)


class Fluid:

    def __init__(self, p_sat_param):
        self.p_sat_param = p_sat_param

    def calc_p_sat(self, t_in):
        return np.polyval(self.p_sat_param, t_in)


water = Fluid(p_saturation_param)
