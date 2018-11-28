import numpy as np
import scipy.linalg as sp_l
import data.global_parameter as g_par
import system.global_functions as g_func
np.set_printoptions(linewidth=10000, threshold=None, precision=2)


class TemperatureSystem:

    def __init__(self, t_sys_const_dict):
        # Handover
        self.cell_num = t_sys_const_dict['cell_num']
        self.nodes = t_sys_const_dict['nodes']
        self.elements = self.nodes - 1
        self.cool_ch_bc = t_sys_const_dict['cool_ch_bc']
        self.t_gas_in = t_sys_const_dict['t_gas_in']
        self.t_cool_in = t_sys_const_dict['t_cool_in']
        self.t_layer_init = t_sys_const_dict['t_layer_init']
        self.g_cool = t_sys_const_dict['g_cool']
        self.heat_pow = t_sys_const_dict['heat_pow']
        self.k_cool_ch = t_sys_const_dict['k_cool_ch']
        self.k_layer = t_sys_const_dict['k_layer']
        self.k_alpha_env = t_sys_const_dict['k_alpha_env']
        # Variable
        self.t_env = g_par.dict_case['t_u']
        self.h_vap = g_par.dict_uni['h_vap']
        self.v_tn = g_par.dict_case['vtn']
        self.const_var = - self.k_cool_ch / self.g_cool
        self.layer = 5
        # Array
        self.mat = None
        self.zero = np.full(5, .0)
        self.k_gas_ch = np.full((2, self.cell_num, self.elements), 0.)
        self.t_layer_1d =\
            np.full(self.elements * (5 * (self.cell_num - 1) + 6), 0.)
        self.r_side = np.full(self.elements * (5 * (self.cell_num - 1) + 6), 0.)
        t_layer_else = np.full((5, self.elements), self.t_layer_init)
        t_layer_n = np.full((6, self.elements), self.t_layer_init)
        self.t_layer = []
        for q in range(self.cell_num):
            if q is not self.cell_num-1:
                self.t_layer.append(t_layer_else)
            else:
                self.t_layer.append(t_layer_n)
        self.t_layer = np.asarray(self.t_layer)
        self.t_gas = np.full((2, self.cell_num, self.nodes), self.t_gas_in[0])
        self.q_cat_mem = np.full((self.cell_num, self.elements), 0.)
        self.q_ano_mem = np.full((self.cell_num, self.elements), 0.)
        self.q_cool = np.full(self.cell_num, 0.)
        self.q_sum = np.full((self.cell_num, self.elements), 0.)
        if self.cool_ch_bc is True:
            self.t_cool = np.full((self.cell_num + 1, self.nodes),
                                  self.t_cool_in)
            self.t_cool_ele = np.full((self.cell_num + 1, self.elements), 0.)
        else:
            self.t_cool = np.full((self.cell_num, self.nodes),
                                  self.t_cool_in)
            self.t_cool_ele = np.full((self.cell_num, self.elements), 0.)
        self.g_gas = np.full((2, self.cell_num, self.elements), 0.)
        self.k_gas_ch = np.full((2, self.cell_num, self.elements), 0.)
        self.m_reac_flow_delta = np.full((self.cell_num, self.nodes), 0.)
        self.cp_h2 = np.full((self.cell_num, self.nodes), 0.)
        self.gamma = np.full((2, self.cell_num, self.nodes), 0.)
        self.i = np.full((self.cell_num, self.elements), 0.)
        self.v_los = np.full((2, self.cell_num, self.elements), 0.)
        self.omega = np.full((self.cell_num, self.elements), 0.)
        self.k_alpha_env = self.k_alpha_env * 0.
        # Initialize temperature matrix arrays
        stack_mat_cell = []
        stack_k_node_side = []
        stack_k_up_cell = []
        stack_k_down_cell_n = []
        stack_k_down_cell_else = []
        stack_k_mid_cell = []
        for q in range(self.cell_num):
            k_node = np.array([-self.k_layer[1, 2, q],
                               self.k_layer[1, 1, q],
                               self.k_layer[1, 0, q],
                               self.k_layer[1, 0, q],
                               + self.k_layer[1, 1, q]])
            if q is self.cell_num - 1:
                layer = 6
                mat = np.full((layer, layer), 0.)
                k_node = np.hstack((*k_node, [-self.k_layer[1, 2, q]]))
            else:
                layer = 5
                mat = np.full((layer, layer), 0.)
            mat[0, 0] = -self.k_layer[0, 2, q]\
                - self.k_alpha_env[0, 2, q]
            mat[0, 1] = self.k_layer[0, 2, q]
            mat[1, 0] = -self.k_layer[0, 2, q]
            mat[1, 1] = self.k_layer[0, 2, q] + self.k_layer[0, 1, q]\
                + self.k_alpha_env[0, 1, q]
            mat[1, 2] = -self.k_layer[0, 1, q]
            mat[2, 1] = -self.k_layer[0, 1, q]
            mat[2, 2] = self.k_layer[0, 1, q] + self.k_layer[0, 0, q]\
                + self.k_alpha_env[0, 0, q]
            mat[2, 3] = -self.k_layer[0, 0, q]
            mat[3, 2] = -self.k_layer[0, 0, q]
            mat[3, 3] = \
                self.k_layer[0, 0, q] + self.k_layer[0, 1, q]\
                + self.k_alpha_env[0, 0, q]
            mat[3, 4] = -self.k_layer[0, 1, q]
            mat[4, 3] = -self.k_layer[0, 1, q]
            mat[4, 4] = self.k_layer[0, 1, q] + self.k_alpha_env[0, 1, q]
            if q is 0:
                mat[0, 0] += .5 * self.k_alpha_env[0, 2, 0]
                mat[4, 4] += self.k_layer[0, 2, q]
                if self.cool_ch_bc is True:
                    mat[0, 0] -= self.k_cool_ch
            elif 0 < q < self.cell_num-1:
                mat[0, 0] -= self.k_layer[0, 2, q] + self.k_cool_ch
                mat[4, 4] += self.k_layer[0, 2, q]
            else:
                mat[0, 0] -= self.k_layer[0, 2, q] + self.k_cool_ch
                mat[4, 4] += self.k_layer[0, 2, q]
                mat[4, 5] -= self.k_layer[0, 2, q]
                mat[5, 5] -=\
                    self.k_layer[0, 2, q] + .5 * self.k_alpha_env[0, 2, q]
                mat[5, 4] += self.k_layer[0, 2, q]
                if self.cool_ch_bc is True:
                    mat[5, 5] -= self.k_cool_ch
            k_node_mid = np.hstack((k_node,
                                    *np.tile(2. * k_node, self.elements-2),
                                    k_node))
            k_node_side = np.hstack((*np.tile(-k_node, self.elements-1),
                                     np.full(layer, 0.)))



            stack_mat_cell.append(g_func.stack_mat(mat, self.elements))
            stack_k_mid_cell.append(k_node_mid)
            stack_k_node_side.append(k_node_side)
        k_node_else = np.array([-self.k_layer[1, 2, q],
                                self.k_layer[1, 1, q],
                                self.k_layer[1, 0, q],
                                self.k_layer[1, 0, q],
                                + self.k_layer[1, 1, q]])
        k_node_n = np.array([-self.k_layer[1, 2, q],
                             self.k_layer[1, 1, q],
                             self.k_layer[1, 0, q],
                             self.k_layer[1, 0, q],
                             + self.k_layer[1, 1, q],
                             -self.k_layer[1, 2, q]])
        k_node_side_else = np.hstack((np.tile(
            np.hstack((np.tile(k_node_else, self.elements-1), np.zeros(5))),
            self.cell_num-1), np.zeros(6*(self.elements-1)+1)))
        k_node_side_n = np.hstack((np.zeros((self.cell_num-1)
                                            * self.elements * 5),
                                   np.tile(k_node_n, self.elements - 1)))
        stack_mat_cell = g_func.change_concrete_4d_too_3d(stack_mat_cell)
        stack_k_node_side = g_func.change_concrete_4d_too_3d(stack_k_node_side)
        stack_k_node_mid = g_func.change_concrete_4d_too_3d(stack_k_mid_cell)
        self.stack_mat = sp_l.block_diag(*stack_mat_cell)\
            + np.diag(stack_k_node_mid)\
            - np.diag(k_node_side_else, 5)\
            - np.diag(k_node_side_else, -5)\
            - np.diag(k_node_side_n, 6)\
            - np.diag(k_node_side_n, -6)

        down_cell_r, down_cell_c = [], []
        for ct in range(self.elements * (self.cell_num - 1)):
            down_cell_r.append(4 + 5 * ct)
            if ct <= self.elements * (self.cell_num - 2):
                down_cell_c.append(5 * self.elements + 5 * ct)
            else:
                down_cell_c.append(down_cell_c[-1] + 6)
        for ct, item in enumerate(down_cell_c):
            self.stack_mat[down_cell_c[ct],
                           down_cell_r[ct]] = self.k_layer[0, 2, -1]
            self.stack_mat[down_cell_r[ct],
                           down_cell_c[ct]] = -self.k_layer[0, 2, -1]
        """print(self.stack_mat)
        for q in range(len(self.stack_mat)):
             print('Zeile:', q, np.sum(self.stack_mat[q]))"""

    def update_values(self, t_sys_dyn_dict):
        self.g_gas =\
            np.array([g_func.calc_elements_2d(t_sys_dyn_dict['g_gas'][0]),
                      g_func.calc_elements_2d(t_sys_dyn_dict['g_gas'][1])])
        self.k_gas_ch = t_sys_dyn_dict['k_gas_ch']
        self.m_reac_flow_delta = t_sys_dyn_dict['m_reac_flow_delta']
        self.cp_h2 = t_sys_dyn_dict['cp_h2']
        self.gamma = t_sys_dyn_dict['gamma']
        self.i = t_sys_dyn_dict['i']
        self.v_los = t_sys_dyn_dict['v_los']
        self.omega = t_sys_dyn_dict['omega']

    def change_value_shape(self):
        self.v_los = np.array(self.v_los)
        self.k_gas_ch = np.array([g_func.calc_elements_2d(self.k_gas_ch[0]),
                                 g_func.calc_elements_2d(self.k_gas_ch[1])])
        self.cp_h2 = g_func.calc_elements_2d(self.cp_h2)

    def update(self):
        self.change_value_shape()
        self.update_gas_channel()
        self.update_coolant_channel()
        self.update_t_layer()

    def update_t_layer(self):
        self.update_matrix()
        self.update_right_side()
        self.solve_system()
        self.sort_results()

    def update_gas_channel(self):
        for q in range(self.cell_num):
            for w in range(1, self.nodes):
                self.t_gas[0, q, w] = \
                    (self.t_gas[0, q, w - 1] - self.t_layer[q][1, w-1])\
                    * np.exp(-self.k_gas_ch[0, q, w-1]
                             / self.g_gas[0, q, w - 1])\
                    + self.t_layer[q][1, w-1]
            for w in range(self.elements - 1, -1, -1):
                self.t_gas[1, q, w] = \
                    (self.t_gas[1, q, w + 1] - self.t_layer[q][ -1, w]) \
                    * np.exp(-self.k_gas_ch[1, q, w] / (self.g_gas[1, q, w])) \
                    + self.t_layer[q][-1, w]

    def update_coolant_channel(self):
        for q in range(self.cell_num):
            for w in range(1, self.nodes):
                self.t_cool[q, w] =\
                    (self.t_cool[q, w - 1] - self.t_layer[q][0, w-1])\
                    * np.exp(self.const_var) + self.t_layer[q][0, w-1]
        if self.cool_ch_bc is True:
            for w in range(1, self.nodes):
                self.t_cool[-1, w] = \
                    (self.t_cool[-1, w - 1] - self.t_layer[-1][-1, w - 1])\
                    * np.exp(self.const_var) + self.t_layer[-1][-1, w-1]
        self.t_cool_ele = g_func.calc_elements_2d(self.t_cool)

    def update_matrix(self):
        dyn_vec = np.full(self.elements*(5 * (self.cell_num-1) + 6), 0.)
        ct = 0
        for q in range(self.cell_num):
            if q is not self.cell_num-1:
                cr = 5
            else:
                cr = 6
            for w in range(self.elements):
                dyn_vec[ct + 1] = self.k_gas_ch[0, q, w]
                dyn_vec[ct + 4] = self.k_gas_ch[1, q, w]
                ct = ct + cr
        self.mat = self.stack_mat + np.diag(dyn_vec)

    def update_right_side(self):
        self.r_side = np.full(self.elements * (5 * (self.cell_num - 1) + 6), 0.)
        s = self
        s.r_s = self.r_side
        ct = 0
        for q in range(self.cell_num):
            for w in range(self.elements):
                s.r_s[ct] -= s.t_env * s.k_alpha_env[0, 2, q]
                s.r_s[ct + 1] += \
                    s.t_env * s.k_alpha_env[0, 1, q]\
                    + s.t_gas[0, q, w] * s.k_gas_ch[0, q, w]\
                    + s.h_vap * s.gamma[0, q, w]
                s.r_s[ct + 2] += s.t_env * s.k_alpha_env[0, 0, q]\
                    + (s.v_tn - g_par.dict_case['e_0'] + self.v_los[0, q, w]
                       + .5 * s.omega[q, w] * s.i[q, w]) * s.i[q, w]
                s.r_s[ct + 3] += s.t_env * s.k_alpha_env[0, 0, q]\
                    + (self.v_los[1, q, w] + s.omega[q, w] * s.i[q, w] * .5)\
                    * self.i[q, w]
                s.r_s[ct + 4] += s.t_env * s.k_alpha_env[0, 1, q]\
                    + s.t_gas[1, q, w] * s.k_gas_ch[1, q, w]\
                    + s.h_vap * s.gamma[1, q, w]
                if q is 0:
                    s.r_s[ct] -=\
                        s.heat_pow + 0.5 * s.t_env * s.k_alpha_env[0, 2, q]
                    if s.cool_ch_bc is True:
                        s.r_s[ct] -= s.k_cool_ch * s.t_cool_ele[0, w]
                    cr = 5
                elif 0 < q < self.cell_num-1:
                    s.r_s[ct] -= s.k_cool_ch * s.t_cool_ele[q, w]
                    cr = 5
                else:
                    s.r_s[ct] -= s.k_cool_ch * s.t_cool_ele[q, w]
                    s.r_s[ct + 5] -=\
                        s.heat_pow + .5 * s.t_env * s.k_alpha_env[0, 2, q]
                    if s.cool_ch_bc is True:
                        s.r_s[ct + 5] -= s.k_cool_ch * s.t_cool_ele[-1, w]
                    cr = 6
                ct += cr

    def solve_system(self):
        self.t_layer_1d = np.linalg.tensorsolve(self.mat, self.r_side)

    def sort_results(self):
        ct = 0
        for q in range(self.cell_num):
            if q is not self.cell_num - 1:
                cr = 5
            else:
                cr = 6
            for w in range(self.elements):
                self.t_layer[q][:, w] = self.t_layer_1d[ct: ct + cr]
                ct = ct + cr
