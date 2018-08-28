# Global functions:
from copy import copy
from numpy import reshape, full, sqrt, log, array, exp, matmul, array, hstack, delete, vstack, zeros
from matplotlib import pyplot as plt
import os, errno

def dw(t):
    return 2.1e-7 * exp(-2436. / t)


def d_pos(t):
    return 1.6e-8 * exp(-1683. / t)


def toarray(var, m, n):
    return reshape(var.flatten(order='C'), (m, n))


def calc_rho(p, r, t):
    return p/(r*t)


def calc_re(roh, v, d, dyn_visc):
    return roh * v * d / dyn_visc


def calc_fanning_friction_factor(re):
    f = full(len(re), 0.)
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


def calc_turbo_nu_numb(re, pr, dh, x):
    if x is 0:
        return 3.66
    f = (1.8 * log(re-1.5))**-2.
    a = f/8. * re * pr * (1. + 1./3.*(dh/x)**(2./3.))
    b = 1. + 12.7 * sqrt(f/8.) * (pr**(2./3.) -1.)
    print(a/b)
    return a/b


def calc_visc_mix(visc, mol_f, mol_m):## Herning Zipper
    visc = array(visc).transpose()
    mol_f = array(mol_f).transpose()
    visc_mix = []
    for q, item in enumerate(visc):
        nenner = sum(visc[q] * mol_f[q] * sqrt(mol_m))
        teiler = sum(mol_f[q] * sqrt(mol_m))
        visc_mix.append(nenner/teiler)
    return visc_mix


def calc_psi(visc, mol_w):
    psi = []
    for q, item in enumerate(visc):
        for w, item in enumerate(visc):
            a = (1.+ (visc[q]/visc[w])**0.5 * (mol_w[w]/mol_w[q])**0.25)**2.
            b = sqrt(8.) * (1 + mol_w[q]/mol_w[w])**0.5
            psi.append(a/b)
    return psi


def calc_lambda_mix(lambdax, mol_f, visc, mol_w):
    if mol_f[0][0] < 0.999:
        if mol_f[1][0] < 0.01: mol_f[1][0] = 0.01
        psi = calc_psi(visc, mol_w)
        counter = 0
        summe = 0.
        for q, item in enumerate(visc):
            a = mol_f[q] * lambdax[q]
            b = 0.
            for w, item in enumerate(visc):
                b = b + mol_f[q] * psi[counter]
                counter = counter + 1
            summe = summe + a / b
    else: summe = lambdax[0]
    return summe


def calc_medium_temp(t_layer, t_medium, coef, m1, m2):
    return coef/(coef + 1.) * matmul(m1,t_layer) + 1./coef * matmul(m2,t_medium)

def calc_elements(node_vec):
    node_vec = list(node_vec)
    a = copy(node_vec)
    node_vec.pop(0), a.pop(-1)
    return (array(a) + array(node_vec)) *.5

def calc_nodes(ele_vec):
    a = copy(ele_vec)
    ele_vec = delete(ele_vec, (0), axis =1)
    a = delete(a, (-1), axis=1)
    mat = (a + ele_vec) * 0.5
    mat = hstack([mat[:, [0]], mat, mat[:, [-1]]])
    return mat

def calc_nodes_1d (ele_vec):
    ele_vec = list(ele_vec)
    a = copy(ele_vec)
    ele_vec.pop(0), a.pop(-1)
    ele_vec = (array(a) + array(ele_vec)) *.5
    return hstack((ele_vec[0],ele_vec,ele_vec[-1]))

def output(y_values, y_label, x_label, y_scale, color, title, q, xlim_low, xlim_up, val_label):
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
    #plt.autoscale(tight=True, axis='both', enable=True)
    plt.xlim(xlim_low,xlim_up)
    plt.tight_layout()
    plt.grid()
    if val_label is not False:
        plt.legend()
    plt.savefig(os.path.join(os.path.dirname(__file__), 'Plots'+q+'/'+title+'.jpg'))
    plt.close()


def output_x(y_values,x_values, y_label, x_label, y_scale, color, title, q, val_label,lim):
    try:
        os.makedirs(os.path.join(os.path.dirname(__file__), 'Plots'+q+'/'))
    except OSError as e:
        if e.errno !=errno.EEXIST:
            raise
    if val_label is not False:
        for l, item in enumerate(y_values):
            plt.plot(x_values, y_values[l],color=color[l], marker='.', label=val_label[l])
    else:
        for l, item in enumerate(y_values):
            plt.plot(x_values, y_values[l], color=color[l], marker='.')

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.yscale(y_scale)
    plt.autoscale(tight=True, axis='both', enable=True)
    plt.xlim(lim[0],lim[1])
    plt.tight_layout()
    plt.grid()
    if val_label is not False:
        plt.legend()
    plt.savefig(os.path.join(os.path.dirname(__file__), 'Plots'+q+'/'+title+'.jpg'))
    plt.close()

