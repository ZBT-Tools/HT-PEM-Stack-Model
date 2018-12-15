import numpy as np
from matplotlib import pyplot as plt
import os


def dw(t):
    return 2.1e-7 * np.exp(-2436. / t)


def d_pos(t):
    return 1.6e-8 * np.exp(-1683. / t)


def to_array(var, m, n):
    return np.reshape(var.flatten(order='C'), (m, n))


def calc_dif(vec):
    return vec[:-1] - vec[1:]


def calc_rho(p, r, t):
    return p / (r * t)


def calc_Re(roh, v, d, visc):
    return roh * v * d / visc


def calc_fan_fri_fac(re):
    f = np.full(len(re), 0.)
    for q, item in enumerate(re):
        if 0. <= re[q] <= 2100.:
            f[q] = 16./re[q]
        else:
            f[q] = 0.079 * re[q]**-0.25
    return f


def calc_head_p_drop(roh, v1, v2, f, kf, le, dh):
    a = (v1**2. - v2**2.) * .5
    b = v2**2. * (2. * f * le / dh + kf * .5)
    return roh * (a + b)


def calc_visc_mix(visc, mol_f, mol_m):
    visc = np.array(visc).transpose()
    mol_f = np.array(mol_f).transpose()
    visc_mix = []
    for q, item in enumerate(visc):
        denominator = sum(item * mol_f[q] * np.sqrt(mol_m))
        divisor = sum(mol_f[q] * np.sqrt(mol_m))
        visc_mix.append(denominator / divisor)
    return visc_mix


def calc_psi(visc, mol_w):
    psi = []
    for q in range(len(visc)):
        for w in range(len(visc)):
            a = (1. + (visc[q] / visc[w])**0.5
                 * (mol_w[w] / mol_w[q])**0.25)**2.
            b = np.sqrt(8.) * (1. + mol_w[q] / mol_w[w])**0.5
            psi.append(a / b)
    return psi


def calc_lambda_mix(lambdax, mol_f, visc, mol_w):
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

lama = []
mol_f = []
visc = []
mol_w = []


def calc_elements_1_d(node_vec):
    return np.array((node_vec[:-1] + node_vec[1:])) * .5


def calc_elements_2d(node_mat):
    return np.array((node_mat[:, :-1] + node_mat[:, 1:])) * .5


def calc_nodes_2_d(ele_mat):
    mat = np.array((ele_mat[:, :-1] + ele_mat[:, 1:])) * .5
    return np.hstack([mat[:, [0]], mat, mat[:, [-1]]])


def i_e_polate_nodes_2_d(ele_vec):
    re_array = np.array((ele_vec[:, :-1] + ele_vec[:, 1:])) * 0.5
    first_node = np.array([2. * re_array[:, 0] - re_array[:, 1]])
    last_node = np.array([2. * re_array[:, -1] - re_array[:, -2]])
    return np.concatenate((first_node.T, re_array, last_node.T), axis=1)


def calc_nodes_1_d(ele_vec):
    vec = np.array((ele_vec[:-1] + ele_vec[1:])) * .5
    return np.hstack([vec[0], vec, vec[-1]])


def calc_ie_dx(vec1, vec2):
    vec1 = np.array(vec1)
    vec2 = np.array(vec2)
    diff1 = vec1[1:] - vec1[:-1]
    diff2 = vec2[1:] - vec2[:-1]
    res_vec = diff1 / diff2
    res_vec = (res_vec[:-1] + res_vec[1:]) * .5
    return np.concatenate(([2. * res_vec[0] - res_vec[1]],
                           res_vec,
                           [2. * res_vec[-1] - res_vec[-2]]))


def calc_fluid_water_enthalpy(t):
    return (t-273.15) * 4182.


def stack_mat_bc(mat, times):
    array_stack = []
    for q in range(times):
        array_stack.append(mat)
    array_stack[0] = 0.5 * array_stack[0]
    array_stack[-1] = 0.5 * array_stack[-1]
    return np.array(array_stack)


def stack_mat(mat, times):
    array_stack = []
    for q in range(times):
        array_stack.append(mat)
    return np.array(array_stack)


def change_concrete_4d_too_3d(array):
    list = []
    for q in range(len(array)):
        for w, item in enumerate(array[q]):
            list.append(item)
    return np.array(list)


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
    return (temp_in * (g - .5 * k) + temp_wall * k) / (g + k * .5)
