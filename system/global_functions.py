import numpy as np
from matplotlib import pyplot as plt
import os
# from numba import jit


def add_source(var, source, direction=1, tri_mtx=None):
    """
    Add discrete 1d source of length n-1 to var of length n
    :param var: 1d array of quantity variable
    :param source: 1d array of source to add to var
    :param direction: flow_circuit.py direction (1: along array counter, -1: opposite to
    array counter)
    :param tri_mtx: if triangle matrix (2D array, nxn) is not provided,
    it will be created temporarily
    :return:
    """
    n = len(var) - 1
    if len(source) != n:
        raise ValueError('parameter source must be of length (var-1)')
    if direction == 1:
        if tri_mtx is None:
            ones = np.zeros((n, n))
            ones.fill(1.0)
            fwd_mat = np.tril(ones)
        else:
            fwd_mat = tri_mtx
        var[1:] += np.matmul(fwd_mat, source)
    elif direction == -1:
        if tri_mtx is None:
            ones = np.zeros((n, n))
            ones.fill(1.0)
            bwd_mat = np.triu(ones)
        else:
            bwd_mat = tri_mtx
        var[:-1] += np.matmul(bwd_mat, source)
    else:
        raise ValueError('parameter direction must be either 1 or -1')
    return var


def construct_empty_stack_array(cell_array, n_cells):
    """
    Construct zeroed stack array from one- or multidimensional cell variable
    array and number of cells
    :param cell_array: array of variable discretized for a unit cell
    :param n_cells: number of cells in stack
    :return:
    """
    if isinstance(cell_array, np.ndarray):
        cell_array_shape = cell_array.shape
    else:
        cell_array_shape = (len(cell_array))
    stack_array_shape = (n_cells,) + cell_array_shape
    return np.full(stack_array_shape, 0.0)


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
    return np.divide(roh * v * d, visc, where=visc != 0.0)


def calc_friction_factor(reynolds, method='Blasius', type='Darcy'):
    """
    Calculates the fanning friction factor between a wall
    and a fluid for the laminar and turbulent case.
    """
    f = 1.0
    if type == 'Darcy':
        f = 4.0
    elif type == 'Fanning':
        f = 1.0
    else:
        ValueError('friction factor type can only be Darcy or Fanning')
    lam = np.zeros(reynolds.shape)
    turb = np.zeros(reynolds.shape)
    if method == 'Blasius':
        lam = np.divide(f * 16.0, reynolds, out=lam, where=reynolds > 0.0)
        turb = f * 0.0791 * np.power(reynolds, -0.25, out=turb,
                                     where=reynolds > 0.0)
        return np.where(reynolds < 2100.0, lam, turb)
    else:
        raise NotImplementedError


def calc_pressure_drop(velocity, density, f, zeta, length, diameter,
                       pressure_recovery=False):
    """
    Calculates the pressure drop at a defined t-junction,
    according to (Koh, 2003).
    :param density: fluid density array (element-wise)
    :param velocity: velocity array (node-wise)
    :param f: darcy friction factor (element-wise)
    :param zeta: additional loss factors
    :param length: length of element (element-wise)
    :param diameter: hydraulic diameter of pipe
    :param pressure_recovery: for dividing manifolds with momentum effects
    :return: pressure drop (element-wise)
    """
    if np.shape(velocity)[0] != (np.shape(length)[0] + 1):
        raise ValueError('velocity array must be provided as a 1D'
                         'nodal array (n+1), while the other input arrays '
                         'must be element-wise (n)')
    v1 = velocity[:-1]
    v2 = velocity[1:]
    a = v2 ** 2.0 * (f * length / diameter + zeta) * 0.5
    b = (v1 ** 2.0 - v2 ** 2.0) * .5
    return density * (a + b)


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
