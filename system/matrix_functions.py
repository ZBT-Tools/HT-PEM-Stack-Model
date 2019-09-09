import numpy as np
from scipy import linalg as sp_la


def tile_add_overlap(array, n, m=1):
    if n <= 1:
        return array
    else:
        add_on = np.copy(array[:-m])
        add_on[:m] += array[-m:]
        last = np.copy(array)
        last[:m] += array[-m:]
        return np.hstack((array[:-m], np.tile(add_on, n-2), last))


def overlapping_vector(vector, reps, overlap_size):
    n_vector = len(vector)
    n_result = reps*n_vector - (reps-1)*overlap_size
    # n_result = reps * (n_vector - overlap_size) + overlap_size
    result = np.zeros(n_result)
    non_overlap_size = int(n_vector-overlap_size)
    for i in range(reps):
        start_id = i * non_overlap_size
        end_id = start_id + n_vector
        result[start_id:end_id] += vector
    return result


def build_1d_conductance_matrix(cond_vector, offset=1):
    n_layer = len(cond_vector) + 1
    center_diag = overlapping_vector(cond_vector, 2, n_layer-2)
    center_diag *= -1.0
    off_diag = np.copy(cond_vector)
    return np.diag(center_diag, k=0) \
        + np.diag(off_diag, k=-offset) \
        + np.diag(off_diag, k=offset)


def build_stack_1d_conductance_matrix(cond_array, n_ele, n_cells):
    n_layer = len(cond_array) + 1
    mat_base = build_1d_conductance_matrix(cond_array[:-1])
    # heat conductance matrix in z-direction for the cells 0-(n-1)
    mat_n = np.full((n_layer, n_layer), 0.)
    mat_n[:-1, :-1] = mat_base.copy()
    mat_n[-2, -2] -= cond_array[-1]
    mat_n[-2, -1] += cond_array[-1]
    mat_n[-1, -2] += cond_array[-1]
    mat_n[-1, -1] = -cond_array[-1]
    # heat conductance matrix in z-direction for the cell n

    list_mat = []
    for i in range(n_cells - 1):
        for j in range(n_ele):
            list_mat.append(mat_base)
    for i in range(n_ele):
        list_mat.append(mat_n)
    # list of all the heat conductance matrix in z-direction
    # for all cells and all elements
    # for i in range(len(list_mat)):
    #     print(list_mat[i])
    return sp_la.block_diag(*list_mat)


def build_z_cell_conductance_matrix(cond_vector, n_ele):
    list_mat = []
    for i in range(n_ele):
        list_mat.append(build_1d_conductance_matrix(cond_vector[:, i]))
    return sp_la.block_diag(*list_mat)


def build_x_cell_conductance_matrix(cond_vector, n_ele):
    n_layer = len(cond_vector)
    center_diag = np.concatenate([cond_vector[:, i] for i in range(n_ele)])
    center_diag[n_layer:-n_layer] *= 2.0
    center_diag *= -1.0
    off_diag = np.concatenate([cond_vector[:, i] for i in range(n_ele-1)])
    return np.diag(center_diag, k=0) \
        + np.diag(off_diag, k=-n_layer) \
        + np.diag(off_diag, k=n_layer)


def build_cell_conductance_matrix(x_cond_vector, z_cond_vector, n_ele):
    return build_x_cell_conductance_matrix(x_cond_vector, n_ele) \
        + build_z_cell_conductance_matrix(z_cond_vector, n_ele)


def build_heat_conductance_matrix(k_layer, k_cool, k_alpha_env,
                                  n_layer, n_ele, n_cells,
                                  cool_ch_bc, cells):
    mat_base = np.full((n_layer, n_layer), 0.)
    mat_base[0, 0] = - k_layer[0, 2, 0]
    mat_base[0, 1] = k_layer[0, 2, 0]
    mat_base[1, 0] = k_layer[0, 2, 0]
    mat_base[1, 1] = - k_layer[0, 2, 0] - k_layer[0, 1, 0]
    mat_base[1, 2] = + k_layer[0, 1, 0]
    mat_base[2, 1] = + k_layer[0, 1, 0]
    mat_base[2, 2] = - k_layer[0, 1, 0] - k_layer[0, 0, 0]
    mat_base[2, 3] = + k_layer[0, 0, 0]
    mat_base[3, 2] = + k_layer[0, 0, 0]
    mat_base[3, 3] = - k_layer[0, 0, 0] - k_layer[0, 1, 0]
    mat_base[3, 4] = + k_layer[0, 1, 0]
    mat_base[4, 3] = + k_layer[0, 1, 0]
    mat_base[4, 4] = - k_layer[0, 1, 0]
    # heat conductance matrix in z-direction for the cells 0-(n-1)
    mat_n = np.full((n_layer + 1, n_layer + 1), 0.)
    mat_n[0:5, 0:5] = mat_base
    mat_n[4, 4] -= k_layer[0, 2, 0]
    mat_n[4, 5] = k_layer[0, 2, 0]
    mat_n[5, 4] = k_layer[0, 2, 0]
    mat_n[5, 5] = -k_layer[0, 2, 0]
    # heat conductance matrix in z-direction for the cell n

    list_mat = []
    for i in range(n_cells - 1):
        for j in range(n_ele):
            list_mat.append(mat_base)
    for i in range(n_ele):
        list_mat.append(mat_n)
    # list of all the heat conductance matrix in z-direction
    # for all cells and all elements
    # for i in range(len(list_mat)):
    #     print(list_mat[i])
    mat_const = sp_la.block_diag(*list_mat)
    # uncoupled heat conductance matrix in z-direction
    #mat_const_1 = mat_const.copy()
    print(mat_const)

    cond_array = [k_layer[0, 2, 0],
                  k_layer[0, 1, 0],
                  k_layer[0, 0, 0],
                  k_layer[0, 1, 0],
                  k_layer[0, 2, 0]]

    mat_const_1 = \
        build_stack_1d_conductance_matrix(cond_array, n_ele, n_cells)
    print(np.sum(np.abs(mat_const-mat_const_1)))
    """Setting the x-axis heat conductance"""
    x_con_base = np.array([k_layer[1, 2, 0],
                           k_layer[1, 1, 0],

                           k_layer[1, 0, 0],
                           k_layer[1, 0, 0],
                           k_layer[1, 1, 0]])
    # heat conductance vec for one element of the 1-(n-1) cell
    x_con_n = np.hstack((x_con_base, 0.5 * k_layer[1, 2, 0]))
    # heat conductance vec for one element of the n cell
    x_con_zero = np.copy(x_con_base)
    x_con_zero[0] *= 0.5

    x_con_base_new = np.array([k_layer[1, 2, 0],
                               k_layer[1, 1, 0],
                               k_layer[1, 0, 0],
                               k_layer[1, 0, 0],
                               k_layer[1, 1, 0],
                               k_layer[1, 2, 0]])
    x_con_base_new[0] *= 0.5
    x_con_base_new[-1] *= 0.5

    # heat conductance vec for one element of the 0th cell
    # sub_array = np.hstack((np.tile(x_con_base, n_ele - 1),
    #                        np.zeros(n_layer)))
    sub_array_base = np.tile(x_con_base, n_ele)
    sub_array_base[-n_layer:] *= 0.0
    sub_array_zero = np.tile(x_con_zero, n_ele)
    sub_array_zero[-n_layer:] *= 0.0
    sub_array_base_new = np.tile(x_con_base_new, n_ele)
    sub_array_base_new[-(n_layer + 1):] *= 0.0

    # x_con_side_base = np.hstack((
    #     np.hstack((np.tile(x_con_zero, n_ele - 1),
    #                np.zeros(n_layer))),
    #     np.tile(np.hstack((np.tile(x_con_base, n_ele - 1),
    #                        np.zeros(n_layer))), n_cells - 2),
    #     np.zeros((n_layer + 1) * (n_ele - 1) + 1)))
    x_con_side_base = np.hstack((
        sub_array_zero, np.tile(sub_array_base, n_cells - 2),
        np.zeros((n_layer + 1) * (n_ele - 1) + 1)))
    # heat conductance vec for the right and left side
    # of the matrix main diagonal for the cells 0-(n-1)
    x_con_side_n = \
        np.hstack((np.zeros((n_cells - 1) * n_ele
                            * n_layer),
                   np.tile(x_con_n, n_ele - 1)))
    # heat conductance vec for the right and left side
    # of the matrix main diagonal for the cell n
    x_con_mid = \
        np.hstack((np.hstack((x_con_zero,
                              *np.tile(2. * x_con_zero, n_ele - 2),
                              x_con_zero)),
                   np.tile(np.hstack((x_con_base,
                                      *np.tile(2. * x_con_base,
                                               n_ele - 2),
                                      x_con_base)),
                           n_cells - 2),
                   np.hstack((x_con_n, *np.tile(2. * x_con_n,
                                                n_ele - 2),
                              x_con_n))))
    # heat conductance vec for
    # the diagonal of the main heat conductance matrix for the cells 0-n
    mat_const = mat_const \
        - np.diag(x_con_mid) \
        + np.diag(x_con_side_base, n_layer) \
        + np.diag(x_con_side_base, -n_layer) \
        + np.diag(x_con_side_n, n_layer + 1) \
        + np.diag(x_con_side_n, -(n_layer + 1))
    mat_const_2 = mat_const.copy()
    print(mat_const_2)
    print(mat_const_2-mat_const_1)

    list_mat = [cell.heat_cond_mtx for cell in cells]
    mat_const_2_2 = sp_la.block_diag(*list_mat)
    print(mat_const_2_2 - mat_const_1)
    print(mat_const_2_2 - mat_const_2)
    print('mat_const_2_2')
    print(np.sum(np.abs(mat_const_2-mat_const_2_2)))

    """Setting the coolant channel heat conductance"""
    cool_pos_n_up = \
        np.arange(n_ele * (n_cells - 1) * n_layer,
                  n_ele * (n_layer * n_cells + 1),
                  n_layer + 1)
    # upper cool ch pos for the n cell

    if cool_ch_bc:
        cool_pos_base = \
            np.arange(0, n_ele * (n_cells - 1) * n_layer,
                      n_layer)
        # cool ch pos for the 0-(n-1) cell
        cool_pos_n_down = \
            np.arange(n_ele * (n_cells - 1) * n_layer +
                      n_layer,
                      n_ele * (n_layer * n_cells + 1),
                      n_layer + 1)
        # lower cool ch pos for the n cell
        for pos in cool_pos_n_up:
            mat_const[pos, pos] -= k_cool

    else:
        cool_pos_base = \
            np.arange(n_ele,
                      n_ele * (n_cells - 1) * n_layer,
                      n_layer)
        # cool ch pos for the 1-(n-1) cell
    for pos in cool_pos_base:
        mat_const[pos, pos] -= k_cool
    if cool_ch_bc:
        for pos in cool_pos_n_down:
            mat_const[pos, pos] -= k_cool

    mat_const_3 = mat_const.copy()
    print(mat_const_3)
    print(mat_const_3-mat_const_2)

    """Setting the cell connecting heat conductance, z-direction"""
    pos_r, pos_c = [], []
    for i in range(n_ele * (n_cells - 1)):
        pos_r.append(4 + n_layer * i)
        if i <= n_ele * (n_cells - 2):
            pos_c.append(n_layer * n_ele + n_layer * i)
        else:
            pos_c.append(pos_c[-1] + n_layer + 1)
    for i in range(len(pos_c)):
        mat_const[pos_c[i], pos_r[i]] += k_layer[0, 2, 0]
        mat_const[pos_r[i], pos_c[i]] += k_layer[0, 2, 0]
    # heat conductance outside the main diagonal

    pos_base_out = \
        np.arange(4, n_ele * (n_cells - 1) * n_layer,
                  n_layer)
    # coordinates of the heat conductance
    # of the last layer of the elements for the cells 0-(n-1)
    pos_base_in = \
        np.arange(n_ele * n_layer,
                  n_ele * (n_cells - 1) * n_layer,
                  n_layer)
    # coordinates of the heat conductance
    # of the first layer of the cells 1-(n-1)
    pos_n_in = np.arange((n_cells - 1) * n_ele * 5,
                         n_ele * (n_layer * n_cells + 1),
                         n_layer + 1)
    # coordinates of the heat conductance
    # of first layer of the elements for the last cell
    pos_in = np.hstack((pos_base_in, pos_n_in))
    # coordinates inside the main diagonal
    for i in range(len(pos_in)):
        mat_const[pos_in[i], pos_in[i]] -= k_layer[0, 2, 0]
        mat_const[pos_base_out[i], pos_base_out[i]] -= k_layer[0, 2, 0]

    mat_const_4 = mat_const.copy()
    print(mat_const_4)
    print(mat_const_4-mat_const_3)
    """Adding the environment heat conductance"""
    env_con_base = np.array([-k_alpha_env[0, 2, 0],
                             -k_alpha_env[0, 1, 0],
                             -k_alpha_env[0, 0, 0],
                             -k_alpha_env[0, 0, 0],
                             -k_alpha_env[0, 1, 0]])
    # environment heat conductance for the cells 1-(n-1)
    env_con_n = np.hstack((env_con_base,
                           - .5 * k_alpha_env[0, 2, 0]))
    # environment heat conductance for the cell n
    env_con_zero = np.hstack((-.5 * k_alpha_env[0, 2, 0],
                              env_con_base[1:]))
    # environment heat conductance for the zeroth cell
    env_con_vec = \
        np.hstack((np.tile(env_con_zero, n_ele),
                   np.tile(np.tile(env_con_base, n_cells - 2),
                           n_ele),
                   np.tile(env_con_n, n_ele)))
    # vector of the main diagonal of the heat conductance matrix
    mat_const = mat_const + np.diag(env_con_vec)
    mat_const_5 = mat_const.copy()
    print(mat_const_5)
    print(mat_const_5-mat_const_4)
    return mat_const
    # mat_dyn = np.copy(mat_const)