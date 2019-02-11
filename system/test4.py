import numpy as np
a = [0,1,2,3,4,5]
print(a)
print(a[:-1])

ptotal = 1.e5 + 1250.
T = 298.15
R = 8.31472

def calc_v_flow(visc, length, diameter, density, p, zeta):
    a = 32. * visc * length / ( diameter**2. * zeta * density)
    b = 2. * p / (zeta * density)
    u = -a + np.sqrt(a **2. + b)
    return u, np.pi*.25 * diameter**2. * ptotal / T/ R



print(calc_v_flow(17.e-6, 0.65, 1.e-3, 1., 2500., 0.1))