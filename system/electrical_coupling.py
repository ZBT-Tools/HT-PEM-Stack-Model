import numpy as np
import data.global_parameters as g_par
import system.global_functions as g_func


class ElectricalCoupling:

    def __init__(self, dict_electrical_coupling_const):
        # Handover
        self.cell_numb = dict_electrical_coupling_const['cell_numb']
        # number of the stack cells
        self.dx = dict_electrical_coupling_const['dx']
        # length of an element
        self.th_plate = dict_electrical_coupling_const['th_bpp']
        # thickness of the bipolar plate
        self.width_channels = dict_electrical_coupling_const['width_channels']
        # width of the channel
        # Variables
        self.nodes = g_par.dict_case['nodes']
        # number of the nodes along the channel
        self.elements = self.nodes - 1
        # number of the elements along the channel
        self.v_end_plate = 0.
        # accumulated voltage loss over the stack at the lower end plate
        c_x = self.width_channels * self.th_plate \
            / (self.dx * g_par.dict_case['bpp_resistivity'])
        # electrical conductance of the bipolar plate in x-direction
        self.v_loss = []
        # 2-d-array of of the voltage loss over the stack in z-direction
        self.cell_c = []
        # 1-d-array of the  combined cell & bipolar plate
        # conductance in z-direction
        self.cell_c_mid = []
        # cell_c array for the main diagonal of the matrix self.mat_const
        self.cell_r = np.full((self.cell_numb, self.elements), 0.)
        # 2-d-array of the combined cell & bipolar plate
        # resistance in z-direction
        self.mat = np.full((self.elements, self.cell_numb + 1), 0.)
        # electrical conductance matrix
        self.rhs = np.full((self.cell_numb + 1) * self.elements, 0.)
        # right hand side terms, here the current
        self.i_cd = np.full((self.cell_numb, self.elements), 0.)
        # current density of the elements in z-direction
        c_x_cell = np.hstack(([c_x],
                              np.full(self.elements - 2, 2. * c_x), [c_x]))
        # bipolar conductance 1-d-array in x-direction over one cell
        c_x_stack = np.tile(c_x_cell * 2, self.cell_numb - 1)
        # bipolar conductance 1-d-array  in x-direction over the stack
        c_x_cell_sr = np.hstack((np.full(self.elements - 1, c_x), 0.))
        # bipolar conductance side 1-d-array of the over one cell
        c_x_stack_sr = np.tile(2. * c_x_cell_sr, self.cell_numb - 1)
        # bipolar conductance side 1-d-array of the over the stack
        self.mat_const = - np.diag(c_x_stack)\
            + np.diag(c_x_stack_sr[:-1], 1)\
            + np.diag(c_x_stack_sr[:-1], -1)

    def update_values(self, dict_electrical_coupling_dyn):
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
        self.cell_r = dict_electrical_coupling_dyn['r_cell']
        self.v_loss = dict_electrical_coupling_dyn['v_loss']
        self.cell_c = self.width_channels * self.th_plate / self.cell_r
        self.cell_c_mid = np.hstack((self.cell_c[:-self.elements]
                                     + self.cell_c[self.elements:]))

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
        self.mat = self.mat_const\
            - np.diag(self.cell_c_mid, 0)\
            + np.diag(self.cell_c[:-self.elements][self.elements:],
                      self.elements)\
            + np.diag(self.cell_c[:-self.elements][self.elements:],
                      -self.elements)

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

        self.v_end_plate = np.sum(self.v_loss) / self.elements
        i_end = self.v_end_plate * self.cell_c[:self.elements]
        self.rhs = np.hstack((-i_end,
                              np.full((self.cell_numb - 2) * self.elements, 0.)))

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

        v_new = np.linalg.tensorsolve(self.mat, self.rhs)
        v_new = np.hstack((np.full(self.elements, self.v_end_plate),
                           v_new, np.full(self.elements, 0.)))
        v_diff = v_new[:-self.elements] - v_new[self.elements:]
        i_ca_vec = v_diff / self.cell_r
        i_cd = g_func.to_array(i_ca_vec, self.cell_numb, self.elements)
        self.i_cd = i_cd / np.average(i_cd) * g_par.dict_case['tar_cd']
        print(self.i_cd)
