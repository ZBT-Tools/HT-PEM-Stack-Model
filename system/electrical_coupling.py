import numpy as np
import data.global_parameter as g_par
import system.global_functions as g_func


class ElectricalCoupling:

    def __init__(self, dict_electrical_coupling_const):
        # Handover
        self.cell_num = dict_electrical_coupling_const['cell_num']
        self.d_x = dict_electrical_coupling_const['d_x']
        self.th_plate = dict_electrical_coupling_const['th_plate']
        self.w_ch = dict_electrical_coupling_const['w_ch']
        # Variables
        self.nodes = g_par.dict_case['nodes']
        self.elements = self.nodes - 1
        self.v_end_plate = 0.
        c_x = self.w_ch * self.th_plate\
            / (self.d_x * g_par.dict_case['plate_resistivity'])
        # Arrays
        self.v_los = []
        self.cell_c = []
        self.cell_c_mid = []
        self.cell_r = np.full((self.cell_num, self.elements), 0.)
        self.mat = np.full((self.elements, self.cell_num+1), 0.)
        self.r_side = np.full((self.cell_num + 1) * self.elements, 0.)
        self.i_ca = np.full((self.cell_num, self.elements), 0.)
        c_x_cell = np.hstack(([c_x],
                              np.full(self.elements - 2, 2. * c_x), [c_x]))
        c_x_stack = np.tile(c_x_cell * 2, self.cell_num-1)
        c_x_cell_sr = np.hstack((np.full(self.elements - 1, c_x), 0.))
        c_x_stack_sr = np.tile(2. * c_x_cell_sr, self.cell_num-1)
        self.mat_const = - np.diag(c_x_stack)\
            + np.diag(c_x_stack_sr[:-1], 1)\
            + np.diag(c_x_stack_sr[:-1], -1)

    def update_values(self, dict_electrical_coupling_dyn):
        self.cell_r = dict_electrical_coupling_dyn['r_cell']
        self.v_los = dict_electrical_coupling_dyn['v_los']
        self.cell_c = self.w_ch * self.th_plate / self.cell_r
        self.cell_c_mid = np.hstack((self.cell_c[:-self.elements]
                                     + self.cell_c[self.elements:]))

    def update(self):
        self.update_mat()
        self.update_right_side()
        self.calc_i()

    def update_mat(self):
        self.mat = self.mat_const\
            - np.diag(self.cell_c_mid, 0)\
            + np.diag(self.cell_c[:-self.elements][self.elements:],
                      self.elements)\
            + np.diag(self.cell_c[:-self.elements][self.elements:],
                      -self.elements)

    def update_right_side(self):
        self.v_end_plate = np.sum(self.v_los) / self.elements
        i_end = self.v_end_plate * self.cell_c[:self.elements]
        self.r_side = np.hstack((-i_end,
                                 np.full((self.cell_num-2)
                                         * self.elements, 0.)))

    def calc_i(self):
        v_new = np.linalg.tensorsolve(self.mat, self.r_side)
        v_new = np.hstack((np.full(self.elements, self.v_end_plate),
                           v_new, np.full(self.elements, 0.)))
        v_dif = v_new[:-self.elements] - v_new[self.elements:]
        i_ca_vec = v_dif / self.cell_r
        self.i_ca = g_func.to_array(i_ca_vec, self.cell_num, self.nodes - 1)
        self.i_ca = self.i_ca / np.average(self.i_ca)\
            * g_par.dict_case['tar_cd']
