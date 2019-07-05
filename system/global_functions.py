import numpy as np
from matplotlib import pyplot as plt
import os
from numba import jit


def add_source(var, source, direction=1):
    n = len(var) - 1
    if len(source) != n:
        raise ValueError('source variable must be of length (var-1)')
    if direction == 1:
        fwd_mat = np.tril(np.full((n, n), 1.))
        var[1:] += np.matmul(fwd_mat, source)
    elif direction == -1:
        bwd_mat = np.triu(np.full((n, n), 1.))
        var[:-1] += np.matmul(bwd_mat, source)
    return var


def calc_diff(vec):
    """
    Calculates the difference between the i+1 and i position of an 1-d-array.
    """
    return vec[:-1] - vec[1:]


def calc_rho(p, r, t):
    """
    Calculates the density of an ideal gas.
    """
    return p / (r * t)


def calc_reynolds_number(roh, v, d, visc):
    """"
    Calculates the reynolds number of an given fluid.
    """
    return roh * v * d / visc


def calc_fan_fri_fac(reynolds_number):
    """
    Calculates the fanning friction factor between a wall
    and a fluid for the laminar and turbulent case.
    """
    return np.where(reynolds_number < 2100.0, 16. / reynolds_number,
                    0.0791 * reynolds_number ** (-0.25))


def calc_head_p_drop(rho, v1, v2, f, kf, le, dh):
    """
    Calculates the pressure drop at a defined t-junction,
    according to (Koh, 2003).
    """
    a = (np.square(v1) - np.square(v2)) * .5
    b = np.square(v2) * (2. * f * le / dh + kf * .5)
    return rho * (a + b)


def calc_visc_mix(species_viscosity, mol_fraction, mol_mass):
    """
    Calculates the mixture viscosity of a gas acording to Herning ad Zipperer.
    """
    #species_viscosity_t = np.asarray(species_viscosity).transpose()
    # mol_fraction_t = np.array(mol_fraction).transpose()
    #mw_mix = np.sum(mol_fraction * mol_mass, axis=-1)
    # visc_mix = []
    # for q, item in enumerate(visc):
    #     denominator = sum(item * mol_fraction_t[q] * np.sqrt(mw))
    #     divisor = sum(mol_fraction[q] * np.sqrt(mw))
    #     visc_mix.append(denominator / divisor)
    # return visc_mix
    spec_visc = np.asarray(species_viscosity).transpose()
    x_sqrt_mw = (mol_fraction.transpose() * np.sqrt(mol_mass))
    return np.sum(spec_visc * x_sqrt_mw, axis=-1)/np.sum(x_sqrt_mw, axis=-1)


def calc_wilke_coefficients(species_viscosity, mol_mass):
    """
    Calculates the wilke coefficients for each species combination of a gas.
    """
    psi = []
    n = len(mol_mass)
    for i in range(n):
        for j in range(n):
            a = (1. + np.power(species_viscosity[i] / species_viscosity[j],
                               0.5)
                 * np.power(np.power((mol_mass[j] / mol_mass[i]), 0.25), 2.))
            b = np.sqrt(8.) * np.power((1. + mol_mass[i] / mol_mass[j]), 0.5)
            psi.append(a / b)
    return np.asarray(psi)


def calc_lambda_mix(species_lambda, mol_fraction, species_viscosity, mol_mass):
    """
    Calculates the heat conductivity of a gas mixture,
    according to Wilkes equation.
    """
    # mol_f[1:] = np.minimum(1.e-20, mol_f[1:])
    psi = calc_wilke_coefficients(species_viscosity, mol_mass)
    n = len(mol_mass)
    lambda_mix = np.zeros_like(mol_fraction[0])
    for i in range(n):
        a = mol_fraction[i] * species_lambda[i]
        b = 1.e-20
        for j in range(n):
            b += mol_fraction[i] * psi[j + n * i]
        lambda_mix += a / b
    return lambda_mix
