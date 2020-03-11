import numpy as np
from scipy import ndimage
import data.global_parameters as g_par
from matplotlib import pyplot as plt
import os
# from numba import jit

SMALL = 1e-12


def ensure_list(variable):
    if isinstance(variable, (list, tuple)):
        return variable
    else:
        return [variable]


def full_like(array):
    """faster than native numpy version"""
    result = np.zeros(array.shape)
    result[:] = array
    return result


def full(shape, value):
    """faster than native numpy version"""
    result = np.zeros(shape)
    result[:] = value
    return result


def zeros_like(array):
    """faster than native numpy version"""
    return np.zeros(array.shape)


def fill_transposed(in_array, shape):
    transposed_array = np.zeros(shape).transpose()
    transposed_array[:] = in_array
    return transposed_array.transpose()


def add_source(var, source, direction=1, tri_mtx=None):
    """
    Add discrete 1d source of length n-1 to var of length n
    :param var: 1d array of quantity variable
    :param source: 1d array of source to add to var
    :param direction: flow direction (1: along array counter, -1: opposite to
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


def fill_last_zeros(array, axis=-1, axis_sum=None):
    if axis == 0:
        array_t = array.transpose()
        return fill_last_zeros(array_t, axis=-1).transpose()
    if axis_sum is None:
        axis_sum = np.abs(np.sum(array, axis=0))
    shape = array.shape
    prev = np.arange(shape[-1])
    prev[axis_sum < SMALL] = 0
    prev = np.maximum.accumulate(prev)
    return array[:, prev]


def fill_first_zeros(array, axis=-1, axis_sum=None):
    array = np.flip(array, axis)
    return np.flip(fill_last_zeros(array, axis), axis, axis_sum)


def fill_zero_sum(array, axis=-1, axis_sum=None):
    if axis == 0:
        array_t = array.transpose()
        return fill_zero_sum(array_t, axis=-1).transpose()
    elif axis in (-1, 1):
        if axis_sum is None:
            axis_sum = np.abs(np.sum(array, axis=0))
        else:
            axis_sum = np.abs(axis_sum)
    else:
        raise ValueError('axis must be 0, 1 or -1. only 2D-arrays allowed.')
    nonzero = np.nonzero(axis_sum)[0]
    if nonzero[-1] != len(array) - 1:
        array = fill_last_zeros(array, axis_sum=axis_sum)
    if nonzero[0] != 0:
        array = fill_first_zeros(array, axis_sum=axis_sum)
    return array


def fill_surrounding_average_1d(array, axis=0):
    footprint = np.zeros((3, 3))
    weights = np.array([1.0, 0.0, 1.0])
    if axis == 0:
        footprint[1, :] = weights
    elif axis in (-1, 1):
        footprint[:, 1] = weights
    else:
        raise ValueError('argument axis can only be 0, 1 or -1')

    mask_array = np.sum(np.abs(array), axis) * np.ones(array.shape)
    averaged = ndimage.generic_filter(array, np.nanmean, footprint=footprint,
                                      mode='constant', cval=np.NaN)
    return np.where(mask_array < SMALL, averaged, array)


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
    return np.zeros(stack_array_shape)


def calc_temp_heat_transfer(wall_temp, fluid_temp, capacity_rate, heat_coeff,
                            flow_direction):
    wall_temp = np.asarray(wall_temp)
    fluid_temp = np.asarray(fluid_temp)
    capacity_rate = np.asarray(capacity_rate)
    heat_coeff = np.asarray(heat_coeff)
    assert capacity_rate.shape == wall_temp.shape
    assert heat_coeff.shape == wall_temp.shape
    fluid_temp_avg = np.asarray(fluid_temp[:-1] + fluid_temp[1:]) * .5
    id_range = range(len(wall_temp))
    if flow_direction == -1:
        id_range = reversed(id_range)
    for i in id_range:
        fluid_avg = fluid_temp_avg[i]
        fluid_out_old = 5e5
        error = 1e3
        iter = 0
        itermax = 10
        while error > 1e-4 and iter <= itermax:
            if flow_direction == -1:
                fluid_in = fluid_temp[i + 1]
            else:
                fluid_in = fluid_temp[i]
            delta_temp = wall_temp[i] - fluid_avg
            q = heat_coeff[i] * delta_temp
            fluid_out = fluid_in + q / capacity_rate[i]
            if fluid_in < wall_temp[i]:
                fluid_out = np.minimum(wall_temp[i] - 1e-3, fluid_out)
            else:
                fluid_out = np.maximum(wall_temp[i] + 1e-3, fluid_out)
            fluid_avg = (fluid_in + fluid_out) * 0.5
            error = np.abs(fluid_out_old - fluid_out) / fluid_out
            fluid_out_old = np.copy(fluid_out)
            iter += 1
        if flow_direction == -1:
            fluid_temp[i] = fluid_out
        else:
            fluid_temp[i + 1] = fluid_out
    fluid_temp_avg = np.asarray(fluid_temp[:-1] + fluid_temp[1:]) * .5
    heat = heat_coeff * (wall_temp - fluid_temp_avg)
    return fluid_temp, heat


def calc_diff(vec):
    """
    Calculates the difference between the i+1 and i position of an 1-d-array.
    """
    return vec[:-1] - vec[1:]


def calc_rho(pressure, gas_constant, temperature):
    """
    Calculates the density of an ideal gas.
    """
    return pressure / (gas_constant * temperature)


def calc_reynolds_number(rho, v, d, visc):
    """"
    Calculates the reynolds number of an given fluid.
    """
    return np.divide(rho * v * d, visc, where=visc != 0.0)


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
        return np.where(reynolds < 2200.0, lam, turb)
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
    a = density * v2 ** 2.0 * (f * length / diameter + zeta) * 0.5
    # b = 0.0
    b = (density * v2 ** 2.0 - density * v1 ** 2.0) * .5
    return a + b


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
