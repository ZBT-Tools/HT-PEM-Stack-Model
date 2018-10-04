from copy import copy
import numpy as np
from matplotlib import pyplot as plt
import os, errno

def dw(t):
    return 2.1e-7 * np.exp(-2436. / t)


def d_pos(t):
    return 1.6e-8 * np.exp(-1683. / t)


def toarray(var, m, n):
    return np.reshape(var.flatten(order='C'), (m, n))


def calc_rho(p, r, t):
    return p/(r*t)


def calc_re(roh, v, d, dyn_visc):
    return roh * v * d / dyn_visc


def calc_fanning_friction_factor(re):
    f = np.full(len(re), 0.)
    for q, item in enumerate(re):
        if 0. <= re[q] <= 2100.:
            f[q] = 16./re[q]
        else:
            f[q] = 0.079*re[q]**-0.25
    return f


def calc_header_pressure_drop(roh, v1, v2, f, kf, le, dh):
    a = (v1**2. - v2**2.) * .5
    b = v2**2. * (2.*f * le / dh + kf *.5)
    return roh * (a + b)


def calc_turbo_nu_numb(re, pr, dh, x):
    if x is 0:
        return 3.66
    f = (1.8 * np.log(re-1.5))**-2.
    a = f/8. * re * pr * (1. + 1./3.*(dh/x)**(2./3.))
    b = 1. + 12.7 * np.sqrt(f/8.) * (pr**(2./3.) -1.)
    return a/b


def calc_visc_mix(visc, mol_f, mol_m):
    visc = np.array(visc).transpose()
    mol_f = np.array(mol_f).transpose()
    visc_mix = []
    for q, item in enumerate(visc):
        denominator = sum(visc[q] * mol_f[q] * np.sqrt(mol_m))
        divisor = sum(mol_f[q] * np.sqrt(mol_m))
        visc_mix.append(denominator / divisor)
    return visc_mix


def calc_psi(visc, mol_w):
    psi = []
    for q, item in enumerate(visc):
        for w, item in enumerate(visc):
            a = (1.+ (visc[q] / visc[w])**0.5 * (mol_w[w] / mol_w[q])**0.25)**2.
            b = np.sqrt(8.) * (1. + mol_w[q] / mol_w[w])**0.5
            psi.append(a / b)
    return psi


def calc_lambda_mix(lambdax, mol_f, visc, mol_w):
    if mol_f[0][0] < 0.999:
        if mol_f[1][0] < 0.01: mol_f[1][0] = 0.01
        psi = calc_psi(visc, mol_w)
        counter = 0
        sum = 0.
        for q, item in enumerate(visc):
            a = mol_f[q] * lambdax[q]
            b = 0.
            for w, item in enumerate(visc):
                b = b + mol_f[q] * psi[counter]
                counter = counter + 1
            sum = sum + a / b
    else: sum = lambdax[0]
    return sum


def calc_medium_temp(t_layer, t_medium, coef, m1, m2):
    return coef/(coef + 1.) * np.matmul(m1, t_layer)\
           + 1./coef * np.matmul(m2,t_medium)


def calc_elements(node_vec):
    node_vec = list(node_vec)
    a = copy(node_vec)
    node_vec.pop(0), a.pop(-1)
    return (np.array(a) + np.array(node_vec)) *.5


def calc_nodes(ele_vec):
    a = copy(ele_vec)
    ele_vec = np.delete(ele_vec, (0), axis =1)
    a = np.delete(a, (-1), axis=1)
    mat = (a + ele_vec) * 0.5
    mat = np.hstack([mat[:, [0]], mat, mat[:, [-1]]])
    return mat


def iepolate_nodes(ele_vec):
    ele_vec = np.array(ele_vec)
    re_array = np.array((ele_vec[:, :-1] + ele_vec[:, 1:]) * .5)
    first_node = np.array([2. * re_array[:, 0] - re_array[:, 1]])
    last_node = np.array([2. * re_array[:, -1] - re_array[:, -2]])
    re_array = np.concatenate((first_node.T, re_array, last_node.T), axis=1)
    return re_array


def calc_nodes_1d (ele_vec):
    ele_vec = list(ele_vec)
    a = copy(ele_vec)
    ele_vec.pop(0), a.pop(-1)
    ele_vec = (np.array(a) + np.array(ele_vec)) *.5
    return np.hstack((ele_vec[0],ele_vec,ele_vec[-1]))


def calc_dif(vec):
    return vec[:-1] - vec[1:]


def calc_mat_vert_dif(mat):
    return -mat[1:, :] + mat[:-1, :]


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

def calc_fw_eh(t):
    return (t-273.15) * 4182.


def output(y_values, y_label, x_label, y_scale, color,
           title, q, xlim_low, xlim_up, val_label):
    try:
        os.makedirs(os.path.join(os.path.dirname(__file__), 'Plots'+q+'/'))
    except OSError as e:
        if e.errno !=errno.EEXIST:
            raise
    if val_label is not False:
        for l, item in enumerate(y_values):
            plt.plot(y_values[l],color=color[l], marker='.', label=val_label[l])
    else:
        for l, item in enumerate(y_values):
            plt.plot(y_values[l],color=color[l], marker='.')

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.yscale(y_scale)
    plt.xlim(xlim_low,xlim_up)
    plt.tight_layout()
    plt.grid()
    if val_label is not False:
        plt.legend()
    plt.savefig(os.path.join(os.path.dirname(__file__),
                             'Plots'+q+'/'+title+'.jpg'))
    plt.close()


def output_x(y_values,x_values, y_label, x_label,
             y_scale, title, q, val_label,lim):
    try:
        os.makedirs(os.path.join(os.path.dirname(__file__), 'Plots'+q+'/'))
    except OSError as e:
        if e.errno !=errno.EEXIST:
            raise
    if val_label is not False:
        for l, item in enumerate(y_values):
            plt.plot(x_values, y_values[l],
                     color=plt.cm.coolwarm(l/len(y_values)),
                     marker='.', label=val_label[l])
    else:
        for l, item in enumerate(y_values):
            plt.plot(x_values, y_values[l],
                     color=plt.cm.coolwarm(l/len(y_values)), marker='.')

    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    plt.yscale(y_scale)
    plt.tick_params(labelsize=14)
    plt.autoscale(tight=True, axis='both', enable=True)
    plt.xlim(lim[0],lim[1])
    plt.tight_layout()
    plt.grid()
    if val_label is not False:
        plt.legend()
    plt.savefig(os.path.join(os.path.dirname(__file__),
                             'Plots'+q+'/'+title+'.jpg'))
    plt.close()
