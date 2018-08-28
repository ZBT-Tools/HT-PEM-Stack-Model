from scipy.linalg import block_diag
import numpy as np
from scipy import sparse


def t_mat_no_bc_col(nodes, cells, r_g, r_m, r_p, r_gp, r_gm, r_pp,
                    r_col, r_cat, r_ano, r_conv_gm, r_conv_gp, r_conv_pp,
                    r_conv_egm, r_conv_e_gp, r_conv_e_pp, bc_ch):
    zero = np.array([0., 0., 0., 0., 0.])
    m_mat_stack = []
    t5_vec_stack = []
    t3_vec_stack = []
    node_r_stack_l_side = []
    node_r_stack_r_side = []
    for q in range(cells):
        m = np.full((5, 5), 0.)
        m[0, 0] = 0.
        m[0, 1] = -1. / r_g[q]
        m[1, 0] = -1. / r_g[q]
        m[1, 1] = 0.
        m[1, 2] = -1. / r_m[q]
        m[2, 1] = -1. / r_m[q]
        m[2, 2] = 0.
        m[2, 3] = -1. / r_g[q]
        m[3, 2] = -1. / r_g[q]
        m[3, 3] = 0.
        m[3, 4] = -1. / r_p[q]
        m[4, 3] = 1. / r_p[q]
        m[4, 4] = 0.
        m_cell = block_diag(m * 0.5, np.kron(np.eye(nodes - 2), m), m * 0.5)
        node_r = -1. / np.array([r_gp[q], r_gm[q], r_gm[q], r_gp[q], -r_pp[q]])
        node_r_cell_l_side = np.hstack((np.tile(node_r, nodes - 1), zero))
        node_r_cell_mid = -block_diag(*np.hstack((node_r, 2. * np.tile(node_r, nodes - 2), node_r)))
        node_r_cell_r_side = np.hstack((np.tile(node_r, nodes - 1), zero))
        if q is 0:
            diag_vec_mid = np.array([1. / r_p[q] + 1. / r_g[q] + 1./r_ano[q] + 1./r_conv_gp[q],
                                     1. / r_g[q] + 1. / r_m[q] + 1./r_conv_gm[q],
                                     1. / r_g[q] + 1. / r_m[q] + 1./r_conv_gm[q],
                                     1. / r_g[q] + 1. / r_p[q] + 1./r_cat[q] + 1./r_conv_gp[q],
                                     -1. / r_p[q] - 1./r_conv_pp[q]])
            diag_vec_bc = np.array([1. / r_p[q] + 1. / r_g[q] + 1./r_ano[q] + 1./r_conv_gp[q] + 2./r_conv_e_gp[q],
                                    1. / r_g[q] + 1. / r_m[q] + 1./r_conv_gm[q] + 2./r_conv_egm[q],
                                    1. / r_g[q] + 1. / r_m[q] + 1./r_conv_gm[q] + 2./r_conv_egm[q],
                                   1. / r_g[q] + 1. / r_p[q] + 1./r_cat[q] + 1./r_conv_gp[q] + 2./r_conv_e_gp[q],
                                    -1. / r_p[q] - 1./r_conv_pp[q] -2./r_conv_e_pp[q]])
            if bc_ch is True:
                diag_vec_mid[-1] = diag_vec_mid[-1] - 1./r_col[q]
                diag_vec_bc[-1] = diag_vec_bc[-1] -1./r_col[q]
            t3_vec = np.array([-1. / r_p[q], 0., 0., 0., 0.])
            t5_vec = np.array([0., 0., 0., 0., 0.])
        elif q is cells-1:
            diag_vec_mid = np.array([1. / r_g[q] + 1./r_ano[q] + 1./r_conv_gp[q],
                                 1. / r_g[q] + 1. / r_m[q] + 1./r_conv_gm[q],
                                 1. / r_g[q] + 1. / r_m[q] + 1./r_conv_gm[q],
                                 1. / r_g[q] + 1. / r_p[q] + 1./r_cat[q] + 1./r_conv_gp[q],
                                -2. / r_p[q] -1. / r_col[q] - 1./r_conv_pp[q]])
            diag_vec_bc = np.array([1. / r_g[q] + 1./r_ano[q] + 1./r_conv_gp[q] + 2. / r_conv_e_gp[q],
                                 1. / r_g[q] + 1. / r_m[q] + 1./r_conv_gm[q] + 2. / r_conv_egm[q],
                                 1. / r_g[q] + 1. / r_m[q] + 1./r_conv_gm[q] + 2. / r_conv_egm[q],
                                 1. / r_g[q] + 1. / r_p[q] + 1./r_cat[q] + 1./r_conv_gp[q] + 2./r_conv_e_gp[q],
                                -2. / r_p[q] -1. / r_col[q] - 1./r_conv_pp[q] - 2./r_conv_e_pp[q]])
            if bc_ch is True:
                diag_vec_mid[0] = diag_vec_mid[0] + 1./r_col[q]
                diag_vec_bc[0] = diag_vec_bc[0] + 1./r_col[q]
            t_3_vec_cell = np.full(nodes*5, 0.)
            t5_vec = np.array([1. / r_p[q], 0., 0., 0., 0.])
        else:
            diag_vec_mid = np.array([1. / r_p[q] + 1. / r_g[q] +1./r_ano[q] + 1./r_conv_gp[q],
                                 1. / r_g[q] + 1. / r_m[q] + 1./r_conv_gm[q],
                                 1. / r_g[q] + 1. / r_m[q] + 1./r_conv_gm[q],
                                 1. / r_g[q] + 1. / r_p[q] + 1./r_cat[q] + 1./r_conv_gp[q],
                                 -2. / r_p[q]  -1. / r_col[q] -1./r_conv_pp[q]])
            diag_vec_bc = np.array([1. / r_p[q] + 1. / r_g[q] +1./r_ano[q] + 1./r_conv_gp[q] + 2./r_conv_e_gp[q],
                                 1. / r_g[q] + 1. / r_m[q] + 1./r_conv_gm[q] + 2./r_conv_egm[q],
                                 1. / r_g[q] + 1. / r_m[q] + 1./r_conv_gm[q] +2./r_conv_egm[q],
                                 1. / r_g[q] + 1. / r_p[q] + 1./r_cat[q] + 1./r_conv_gp[q] + 2./r_conv_e_gp[q],
                                 -2. / r_p[q] -1. / r_col[q] -1./r_conv_pp[q] - 2./r_conv_e_pp[q]])
            t3_vec = np.array([-1. / r_p[q], 0., 0., 0., 0.])
            t5_vec = np.array([1. / r_p[q], 0., 0., 0., 0.])
        t_5_vec_cell = np.hstack((t5_vec / 2., np.tile(t5_vec, nodes - 2), t5_vec / 2.))
        t_3_vec_cell = np.hstack((t3_vec / 2., np.tile(t3_vec, nodes - 2), t3_vec / 2.))
        diag_vec_cell = block_diag(*np.hstack((diag_vec_bc * 0.5, np.tile(diag_vec_mid, nodes - 2), diag_vec_bc * .5)))
        m_mat_stack.append(m_cell + diag_vec_cell + node_r_cell_mid)
        node_r_stack_l_side.append(node_r_cell_l_side)
        node_r_stack_r_side.append(node_r_cell_r_side)
        t3_vec_stack.append(t_3_vec_cell)
        t5_vec_stack.append(t_5_vec_cell)
    t5_vec_stack.pop(0)
    t3_vec_stack = np.hstack((np.full(nodes*5+4, 0.),np.block(t3_vec_stack)))
    t5_vec_stack = np.block(t5_vec_stack)
    node_r_stack_r_side = np.hstack((np.full(5,0.),np.block(node_r_stack_r_side)))
    node_r_stack_l_side = np.block(node_r_stack_l_side)
    m_final = block_diag(*m_mat_stack)\
              + sparse.spdiags(t3_vec_stack, [5 * nodes + 4], cells * nodes * 5, cells * nodes * 5)\
              + sparse.spdiags(t5_vec_stack, [ - 5 * nodes-4], cells * nodes * 5, cells * nodes * 5)\
              + sparse.spdiags(node_r_stack_r_side, [5], cells*nodes*5, cells*nodes*5)\
              + sparse.spdiags(node_r_stack_l_side, [-5], cells*nodes*5, cells*nodes*5)
    #print(m_final)
    return m_final


#f = t_mat_no_bc_col(3,3,1./500,1./5000,1./2000,1./50.,1./100.,1./200.,1./4000.)
#print(f)