# Global functions:
import numpy as np


def dw(t):
    return 2.1e-7 * np.exp(-2436. / t)


def d_pos(t):
    return 1.6e-8 * np.exp(-1683. / t)


def toarray(var, m, n):
    return np.reshape(var.flatten(order='C'), (m, n))


def calc_roh(p, r, t):
    return p/(r*t)


def calc_re(roh, v, d, dyn_visc):
    return roh * v *d / dyn_visc


def calc_fanning_friction_factor(re):
    f = np.full(len(re), 0.)
    for q, item in enumerate(re):
        if 0. <= re[q] <= 2100.:
            f[q] = 16./re[q]
        else:
            f[q] = 0.079*re[q]**-0.25
    return f


def calc_header_pressure_drop(roh, v1, v2, f, kf, le, dh):
    a = (v1**2. - v2**2.) / 2.
    b = v2**2. * (2.*f * le / dh + kf/2)
    return roh * (a + b)
