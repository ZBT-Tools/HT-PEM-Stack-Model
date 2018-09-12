##global parameter
## universal parameter
import input_parameter as i_p
dict_uni = {'r': 8.3144598, 'f': 96485.33289, 'h_vap': i_p.h_vap,'cp_liq': 4182.}
##
dict_case = {'tu': i_p.t_u, 'tar_cd': i_p.tar_cd[0], 'nodes': i_p.nodes,
             'elements': i_p.nodes-1, 'vap_m_t_coef': i_p.gamma,
             'mol_con_m': i_p.alpha, 'e_0': i_p.e_o, 't_u':i_p.t_u,
             'vtn': i_p.vtn, 'pem_type': i_p.pem_type,
             'conv_coef': i_p.alpha_conv, 'header_p_in_cat': i_p.p_cat_in,
             'header_p_in_ano':i_p.p_ano_in}