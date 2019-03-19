import numpy as np
from matplotlib import pyplot as plt
import os


def dw(t):
    """
    Calculates the free water content diffusion coefficient in the membrane.
    """

    return 2.1e-7 * np.exp(-2436. / t)


def to_array(var, m, n):
    """
    Changes a sorted 1-d-array to an sorted 2-d-array.
    """

    return np.reshape(var.flatten(order='C'), (m, n))


def calc_dif(vec):
    """
    Calculates the difference between the i+1 and i position of an 1-d-array.
    """

    return vec[:-1] - vec[1:]


def calc_rho(p, r, t):
    """
    Calculates the density of an ideal gas.
    """

    return p / (r * t)


def calc_Re(roh, v, d, visc):
    """"
    Calculates the reynolds number of an given fluid.
    """

    return roh * v * d / visc


def calc_fan_fri_fac(re):
    """
    Calculates the fanning friction factor between a wall
    and a fluid for the laminar and turbulent case.
    """

    f = np.full(len(re), 0.)
    for q, item in enumerate(re):
        if 0. <= re[q] <= 2100.:
            f[q] = 16./re[q]
        else:
            f[q] = 0.079 * re[q]**-0.25
    return f


def calc_head_p_drop(roh, v1, v2, f, kf, le, dh):
    """
    Calculates the pressure drop at an defined t-junction,
    according to (Koh, 2003).
    """

    a = (v1**2. - v2**2.) * .5
    b = v2**2. * (2. * f * le / dh + kf * .5)
    return roh * (a + b)


def calc_visc_mix(visc, mol_f, mol_m):
    """
    Calculates the mixture viscosity of a gas acording to Herning ad Zipperer.
    """

    visc = np.array(visc).transpose()
    mol_f = np.array(mol_f).transpose()
    visc_mix = []
    for q, item in enumerate(visc):
        denominator = sum(item * mol_f[q] * np.sqrt(mol_m))
        divisor = sum(mol_f[q] * np.sqrt(mol_m))
        visc_mix.append(denominator / divisor)
    return visc_mix


def calc_psi(visc, mol_w):
    """
    Calculates the wilke coefficients for each species combination of a gas.
    """
    psi = []
    for q in range(len(visc)):
        for w in range(len(visc)):
            a = (1. + (visc[q] / visc[w])**0.5
                 * (mol_w[w] / mol_w[q])**0.25)**2.
            b = np.sqrt(8.) * (1. + mol_w[q] / mol_w[w])**0.5
            psi.append(a / b)
    return psi


def calc_lambda_mix(lambdax, mol_f, visc, mol_w):
    """
    Calculates the heat conductivity of a gas mixture,
    according to Wilkes equation.
    """

    mol_f[1:] = np.minimum(1.e-20, mol_f[1:])
    psi = calc_psi(visc, mol_w)
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


def calc_elements_1_d(node_vec):
    """
    Calculates an element 1-d-array from a node 1-d-array.
    """

    return np.array((node_vec[:-1] + node_vec[1:])) * .5


def calc_elements_2d(node_mat):
    """
    Calculates an element 2-d-array from a node 2-d-array.
    """

    return np.array((node_mat[:, :-1] + node_mat[:, 1:])) * .5


def calc_nodes_2_d(ele_mat):
    """
    Calculates an node 2-d-array from an element 2-d-array,
    uses the [:, 1], [:, -2] entries of the calculated node 2-d-array
    to fill the first als last row of the node 2-d-array.
    """

    mat = np.array((ele_mat[:, :-1] + ele_mat[:, 1:])) * .5
    return np.hstack([mat[:, [0]], mat, mat[:, [-1]]])


def i_e_polate_nodes_2_d(ele_vec):
    """
    Calculates an node 2-d-array from an element 2-d-array,
    and interpolates the first and last row of the node 2-d-array.
    """

    re_array = np.array((ele_vec[:, :-1] + ele_vec[:, 1:])) * 0.5
    first_node = np.array([2. * re_array[:, 0] - re_array[:, 1]])
    last_node = np.array([2. * re_array[:, -1] - re_array[:, -2]])
    return np.concatenate((first_node.T, re_array, last_node.T), axis=1)


def calc_nodes_1_d(ele_vec):
    """
    Calculates an node 1-d-array from an element 1-d-array,
    uses the [:, 1], [:, -2] entries of the calculated node 1-d-array
    to fill the first als last row of the node 1-d-array.
    """

    vec = np.array((ele_vec[:-1] + ele_vec[1:])) * .5
    return np.hstack([vec[0], vec, vec[-1]])


def calc_fluid_water_enthalpy(t):
    """
    Calculates the enthalpy of fluid water by a given temperature.
    """

    return (t-273.15) * 4182.


def output(y_values, y_label, x_label, y_scale, color,
           title, xlim_low, xlim_up, val_label, path):
    if val_label is not False:
        for l in range(len(y_values)):
            plt.plot(y_values[l], color=color[l],
                     marker='.', label=val_label[l])
    else:
        for l in range(len(y_values)):
            plt.plot(y_values[l], color=color[l], marker='.')

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.yscale(y_scale)
    plt.xlim(xlim_low, xlim_up)
    plt.tight_layout()
    plt.grid()
    if val_label is not False:
        plt.legend()
    plt.savefig(os.path.join(path + title + '.jpg'))
    plt.close()


def output_x(y_values, x_values, y_label, x_label,
             y_scale, title, val_label, lim, path):
    if val_label is not False:
        for l in range(len(y_values)):
            plt.plot(x_values, y_values[l],
                     color=plt.cm.coolwarm(l/len(y_values)),
                     marker='.', label=val_label[l])
    else:
        for l in range(len(y_values)):
            plt.plot(x_values, y_values[l],
                     color=plt.cm.coolwarm(l/len(y_values)), marker='.')

    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    plt.yscale(y_scale)
    plt.tick_params(labelsize=14)
    plt.autoscale(tight=True, axis='both', enable=True)
    plt.xlim(lim[0], lim[1])
    plt.tight_layout()
    plt.grid()
    if val_label is not False:
        plt.legend()
    plt.savefig(os.path.join(path + title + '.jpg'))
    plt.close()


def calc_fluid_temp_out(temp_in, temp_wall, g, k):
    """
    Calculates the linearised fluid outlet temperature
    of an element by the wall temperature and the fluid inlet temperature.
    The function is limited to cases with g > 0.5 * k.
    """

    return (temp_in * (g - .5 * k) + temp_wall * k) / (g + k * .5)
