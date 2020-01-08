import numpy as np
import data.global_parameters as g_par
import system.global_functions as g_func
from scipy import sparse
from scipy.sparse.linalg import spsolve
import system.matrix_functions as mtx


class ElectricalCoupling:

    def __init__(self, electrical_dict, stack, cells):
        self.stack = stack
        self.cells = cells
        # Handover
        self.n_cells = self.stack.n_cells
        # number of the stack cells
        self.dx = self.cells[0].dx
        # length of an element
        self.th_plate = self.cells[0].cathode.bpp.thickness
        # thickness of the bipolar plate

        # Variables
        self.i_cd_tar = g_par.dict_case['target_current_density']
        self.nodes = g_par.dict_case['nodes']
        # number of the nodes along the channel
        self.n_ele = self.nodes - 1
        # number of the elements along the channel
        self.v_end_plate = 0.
        # accumulated voltage loss over the stack at the lower end plate
        self.v_loss = []
        # 2-d-array of of the voltage loss over the stack in z-direction
        self.cell_c = []
        # 1-d-array of the  combined cell & bipolar plate
        # conductance in z-direction
        self.cell_c_mid = []
        # cell_c array for the main diagonal of the matrix self.mat_const
        self.resistance = np.zeros((self.n_cells, self.n_ele))
        # 2-d-array of the combined cell & bipolar plate
        # resistance in z-direction
        self.mat = np.zeros((self.n_ele, self.n_cells + 1))
        # electrical conductance matrix
        self.rhs = np.zeros((self.n_cells + 1) * self.n_ele)
        # right hand side terms, here the current
        self.i_cd = np.zeros((self.n_cells, self.n_ele))
        # current density of the elements in z-direction
        self.c_width = \
            self.cells[0].cathode.rib_width \
            * (self.cells[0].cathode.n_channel + 1)
        # width of the channel
        # c_x = self.c_width * self.th_plate / self.dx \
        #     * self.cells[0].cathode.bpp.electrical_conductivity
        # c_x = self.cells[0].cathode.bpp.electrical_conductance[1]
        # print('c_x', c_x)
        # # electrical conductance of a single bipolar plate in x-direction
        # c_x_cell = np.full(self.n_ele, c_x)
        # print('c_x_cell', c_x_cell)
        # c_x_cell[1:-1] *= 2.0
        # print('c_x_cell', c_x_cell)
        # # bipolar conductance 1-d-array in x-direction over one cell
        # c_x_stack = np.tile(c_x_cell * 2, self.n_cells - 1)
        # # bipolar conductance 1-d-array  in x-direction over the stack
        # c_x_cell_sr = np.hstack((np.full(self.n_ele - 1, c_x[0]), 0.))
        # print('c_x_cell_sr', c_x_cell_sr)
        # # bipolar conductance side 1-d-array of the over one cell
        # c_x_stack_sr = np.tile(2. * c_x_cell_sr, self.n_cells - 1)
        # print('c_x_stack_sr', c_x_stack_sr)
        # # bipolar conductance side 1-d-array of the over the stack
        # self.mat_const_2 = \
        #     - np.diag(c_x_stack) \
        #     + np.diag(c_x_stack_sr[:-1], 1) \
        #     + np.diag(c_x_stack_sr[:-1], -1)

        cell_mat_x_list = [cell.elec_x_mat_const for cell in self.cells]
        self.mat_const = mtx.block_diag_overlap(cell_mat_x_list,
                                                (self.n_ele, self.n_ele))
        self.mat_const = \
            self.mat_const[self.n_ele:-self.n_ele, self.n_ele:-self.n_ele]
        self.mat_const = sparse.csr_matrix(self.mat_const)

    def update(self):
        """
        Coordinates the program sequence
        """
        resistance = \
            np.asarray([cell.resistance for cell in self.cells])
        self.resistance = resistance.flatten()
        conductance = (self.c_width * self.dx / resistance).flatten()
        self.update_mat(conductance)
        self.update_right_side(conductance)
        self.calc_i()

    def update_mat(self, conductance):
        """
        Updates the conductance matrix
        """
        cell_c_mid = np.hstack((conductance[:-self.n_ele]
                                + conductance[self.n_ele:]))
        mat_dyn = \
            - np.diag(cell_c_mid, 0) \
            + np.diag(conductance[:-self.n_ele][self.n_ele:], self.n_ele) \
            + np.diag(conductance[:-self.n_ele][self.n_ele:], -self.n_ele)
        mat_dyn = sparse.csr_matrix(mat_dyn)
        self.mat = self.mat_const + mat_dyn

    def update_right_side(self, conductance):
        """
        Updates the right hand side of the linear system.
        """
        v_loss = \
            np.asarray([cell.v_loss for cell in self.cells]).flatten()
        self.v_end_plate = np.sum(v_loss) / self.n_ele
        i_end = self.v_end_plate * conductance[:self.n_ele]
        self.rhs = \
            np.hstack((-i_end, np.full((self.n_cells - 2) * self.n_ele, 0.)))

    def calc_i(self):
        """
        Calculates the current density
        of the elements in z-direction.
        """

        # v_new = np.linalg.tensorsolve(self.mat, self.rhs)
        v_new = spsolve(self.mat, self.rhs)
        v_new = np.hstack((np.full(self.n_ele, self.v_end_plate),
                           v_new, np.full(self.n_ele, 0.)))
        v_diff = v_new[:-self.n_ele] - v_new[self.n_ele:]
        try:
            i_ca_vec = v_diff / self.resistance
        except ValueError:
            print('test')
        i_cd = np.reshape(i_ca_vec.flatten(order='C'),
                          (self.n_cells, self.n_ele))
        self.i_cd = i_cd / np.average(i_cd) * self.i_cd_tar
