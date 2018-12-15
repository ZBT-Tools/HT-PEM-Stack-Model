import numpy as np
import scipy.optimize as sp_op
import understoi_dict as u_dict
import data.global_parameter as g_par



const_param = None
dyn_param = None
self.i_theta = np.sqrt(2. * self.vol_ex_cd * self.prot_con
                       * self.tafel_slope)
i_star = prot_con * tafel_slope / th_cat


def calc_support_param(i_ca, gas_con_ele, gas_con, diff_coef_gdl, th_gdl):
    i_lim = 4. * g_par.dict_uni['F'] * gas_con * diff_coef_gdl / th_gdl
    i_hat = i_ca / i_star
    short_save = np.sqrt(2. * i_hat)
    beta = short_save / (1. + np.sqrt(1.12 * i_hat) * np.exp(short_save))\
        + np.pi * i_hat / (2. + i_hat)
    var = 1. - i_ca / (i_lim * gas_con_ele / gas_con)
    return i_lim, i_hat, beta, var


def calc_activation_losses(i_ca, tafel_slope, i_theta,
                           gas_con_ele, gas_con, i_star):
    act_ov = tafel_slope\
             * np.arcsinh((i_ca / i_theta) ** 2.
                          / (2. * (gas_con_ele / gas_con[0, :-1])
                          * (1. - np.exp(i_ca / (2. * i_star)))))
    return act_ov


def calc_transport_losses_catalyst_layer(i_ca, tafel_slope, prot_con,
                                         diff_coef_cat, gas_con_ele,
                                         i_star, beta, var):
    cat_diff_los = ((prot_con * tafel_slope ** 2.)
                    / (4. * g_par.dict_uni['F'] * diff_coef_cat
                    * gas_con_ele) * (i_ca / i_star
                                      - np.log10(1. + i_ca**2.
                                                 / (i_star ** 2.
                                                    * beta ** 2.)))) / var
    return cat_diff_los


def calc_transport_losses_diffusion_layer(tafel_slope, var):
    v_gdl_diff = -tafel_slope * np.log10(var)
    if np.isnan(v_gdl_diff) is True:
        return 1.e20
    return -tafel_slope * np.log10(var)


def calc_membrane_losses(i_ca, omega_ca):
    return omega_ca * i_ca


def calc_sum_voltage_losses(v_act_cat, v_act_ano,
                            v_diff_el_cat, v_diff_el_ano,
                            v_diff_gdl_cat, v_diff_dgl_ano, v_mem):
    v_losses = v_act_cat + v_act_ano + v_diff_el_cat + v_diff_el_ano\
               + v_diff_gdl_cat + v_diff_dgl_ano + v_mem
    return v_losses


def calc_over_voltages(i_ca):
    sup_cat = calc_support_param(i_ca, dyn_param[0, 0], dyn_param[0, 1],
                                 const_param[0, 3], const_param[0, 4])
    v_act_cat =\
        calc_activation_losses(i_ca, const_param[0, 0], const_param[0, 1],
                               dyn_param[0, 0], const_param[0, 2])
    v_diff_el_cat =\
        calc_transport_losses_catalyst_layer(i_ca, const_param[0, 0],
                                             const_param[0, 5],
                                             const_param[0, 6],
                                             dyn_param[0, 0],
                                             const_param[0, 2],
                                             sup_cat[2], sup_cat[3])
    v_diff_gdl_cat =\
        calc_transport_losses_diffusion_layer(i_ca, dyn_param[0, 2])
    sup_ano = calc_support_param(i_ca, dyn_param[1, 0], dyn_param[1, 1],
                                 const_param[1, 3], const_param[1, 4])
    v_act_ano =\
        calc_activation_losses(i_ca, const_param[1, 0], const_param[1, 1],
                               dyn_param[1, 0], const_param[1, 2])
    v_diff_el_ano = \
        calc_transport_losses_catalyst_layer(i_ca, const_param[1, 0],
                                             const_param[1, 5],
                                             const_param[1, 6],
                                             dyn_param[1, 0],
                                             const_param[1, 2],
                                             sup_ano[2], sup_ano[3])

    v_diff_gdl_ano = calc_transport_losses_diffusion_layer()
    v_mem = calc_membrane_losses()
    return calc_sum_voltage_losses(v_act_cat, v_act_ano,
                                   v_diff_el_cat, v_diff_el_ano,
                                   v_diff_gdl_cat, v_diff_gdl_ano, v_mem)


def solve_current_density(i_start):
    i_ca = sp_op.fsolve(const_param[2]['e_0'] - calc_over_voltages, i_start)

