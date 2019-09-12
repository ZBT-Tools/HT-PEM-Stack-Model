import numpy as np
import data.global_parameters as g_par
import system.global_functions as g_func
from scipy import sparse
from scipy.sparse.linalg import spsolve


class ElectricalCoupling:

    def __init__(self, electrical_dict):
        # Handover
        self.n_cells = electrical_dict['cell_number']
        # number of the stack cells
        self.dx = electrical_dict['dx']
        # length of an element
        self.th_plate = electrical_dict['th_bpp']
        # thickness of the bipolar plate
        self.c_width = electrical_dict['conducting_width']
        # width of the channel
        # Variables
        self.nodes = g_par.dict_case['nodes']
        # number of the nodes along the channel
        self.n_ele = self.nodes - 1
        # number of the elements along the channel
        self.v_end_plate = 0.
        # accumulated voltage loss over the stack at the lower end plate
        c_x = self.c_width * self.th_plate \
              / (self.dx * g_par.dict_case['bpp_resistivity'])
        # electrical conductance of a single bipolar plate in x-direction
        self.v_loss = []
        # 2-d-array of of the voltage loss over the stack in z-direction
        self.cell_c = []
        # 1-d-array of the  combined cell & bipolar plate
        # conductance in z-direction
        self.cell_c_mid = []
        # cell_c array for the main diagonal of the matrix self.mat_const
        self.cell_r = np.full((self.n_cells, self.n_ele), 0.)
        # 2-d-array of the combined cell & bipolar plate
        # resistance in z-direction
        self.mat = np.full((self.n_ele, self.n_cells + 1), 0.)
        # electrical conductance matrix
        self.rhs = np.full((self.n_cells + 1) * self.n_ele, 0.)
        # right hand side terms, here the current
        self.i_cd = np.full((self.n_cells, self.n_ele), 0.)
        # current density of the elements in z-direction
        c_x_cell = np.full(self.n_ele, c_x)
        c_x_cell[1:-1] *= 2.0
        # bipolar conductance 1-d-array in x-direction over one cell
        c_x_stack = np.tile(c_x_cell * 2, self.n_cells - 1)
        # bipolar conductance 1-d-array  in x-direction over the stack
        c_x_cell_sr = np.hstack((np.full(self.n_ele - 1, c_x), 0.))
        # bipolar conductance side 1-d-array of the over one cell
        c_x_stack_sr = np.tile(2. * c_x_cell_sr, self.n_cells - 1)
        # bipolar conductance side 1-d-array of the over the stack
        self.mat_const = \
            - np.diag(c_x_stack) \
            + np.diag(c_x_stack_sr[:-1], 1) \
            + np.diag(c_x_stack_sr[:-1], -1)
        self.mat_const = sparse.csr_matrix(self.mat_const)

    def update_values(self, cell_r, v_loss):
        """
        Updates the dynamic parameters

            Access to:
            -dict_electrical_coupling_dyn

            Manipulate:
            -self.cell_r
            -self.v_loss
            -self.cell_c
            -self.cell_c_mid
        """
        self.cell_r = cell_r
        self.v_loss = v_loss
        self.cell_c = self.c_width * self.th_plate / self.cell_r
        self.cell_c_mid = np.hstack((self.cell_c[:-self.n_ele]
                                     + self.cell_c[self.n_ele:]))
        #print(self.cell_r)
        #print(self.v_loss)

    def update(self):
        """
        This function coordinates the program sequence
        """
        self.update_mat()
        self.update_right_side()
        self.calc_i()

    def update_mat(self):
        """
        This function updates the conductance matrix

            Access to:
            -self.mat_const
            -self.cell_c_mid
            -self.cell_c
            -self.elements

            Manipulate:
            -self.mat
        """
        mat_dyn = \
            - np.diag(self.cell_c_mid, 0) \
            + np.diag(self.cell_c[:-self.n_ele][self.n_ele:], self.n_ele) \
            + np.diag(self.cell_c[:-self.n_ele][self.n_ele:], -self.n_ele)
        mat_dyn = sparse.csr_matrix(mat_dyn)
        self.mat = self.mat_const + mat_dyn

    def update_right_side(self):
        """
        This function sets the right hand side.

            Access to:
            -self.v_loss
            -self.elements
            -self.cell_c
            -self.cell_numb

            Manipulate:
            -self.v_end_plate
            -self.rhs
        """

        self.v_end_plate = np.sum(self.v_loss) / self.n_ele
        i_end = self.v_end_plate * self.cell_c[:self.n_ele]
        self.rhs = \
            np.hstack((-i_end, np.full((self.n_cells - 2) * self.n_ele, 0.)))

    def calc_i(self):
        """
        This function calculates the current density
        of the elements in z-direction.

            Access to:
            -self.mat
            -self.rhs
            -self.elements
            -self.v_end_plate
            -self.cell_r
            -self.cell_numb
            -self.nodes
            -g_par.dict_case['tar_cd']

            Manipulate:
            -self.i_ca
        """

        # v_new = np.linalg.tensorsolve(self.mat, self.rhs)
        v_new = spsolve(self.mat, self.rhs)
        v_new = np.hstack((np.full(self.n_ele, self.v_end_plate),
                           v_new, np.full(self.n_ele, 0.)))
        v_diff = v_new[:-self.n_ele] - v_new[self.n_ele:]
        try:
            i_ca_vec = v_diff / self.cell_r
        except ValueError:
            print('test')
        i_cd = np.reshape(i_ca_vec.flatten(order='C'),
                          (self.n_cells, self.n_ele))
        self.i_cd = i_cd / np.average(i_cd) * g_par.dict_case['tar_cd']
