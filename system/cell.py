import numpy as np
import data.global_parameters as g_par
import system.global_functions as g_func
import system.half_cell as h_c
import system.matrix_functions as mtx
import system.membrane as membrane
import system.interpolation as ip
from system.output_object import OutputObject
import copy


class Cell(OutputObject):

    def __init__(self, number, cell_dict, membrane_dict, half_cell_dicts,
                 channels):
        name = 'Cell ' + str(number)
        super().__init__(name)
        self.cell_dict = cell_dict
        print('Initializing: ', self.name)
        self.n_layer = 5
        self.n_electrodes = 2

        self.first_cell = cell_dict['first_cell']
        self.last_cell = cell_dict['last_cell']

        if self.last_cell:
            self.n_layer += 1
        n_nodes = g_par.dict_case['nodes']
        # number of nodes along the channel
        self.n_ele = n_nodes - 1
        n_ele = self.n_ele

        self.e_0 = g_par.dict_case['e_0']

        self.width = self.cell_dict['width']
        self.length = self.cell_dict['length']
        self.active_area = self.width * self.length

        half_cell_dicts = copy.deepcopy(half_cell_dicts)
        # Create half cell objects
        for i in range(len(half_cell_dicts)):
            name = self.name + ': ' + half_cell_dicts[i]['name']
            half_cell_dicts[i]['name'] = name

        self.half_cells = [h_c.HalfCell(half_cell_dicts[i], cell_dict,
                                        channels[i])
                           for i in range(2)]
        self.cathode = self.half_cells[0]
        self.anode = self.half_cells[1]

        self.dx = self.cathode.channel.dx
        # cathode - object of the class HalfCell
        # self.th_mem = cell_dict['th_mem']
        # thickness membrane
        # self.lambda_bpp = [cell_dict['thermal conductivity bpp'][0],
        #                    cell_dict['thermal conductivity bpp'][1]]
        # # heat conductivity of the bipolar plate
        # self.lambda_gde = [cell_dict['thermal conductivity gde'][0],
        #                    cell_dict['thermal conductivity gde'][1]]
        # # heat conductivity of the gas diffusion layer
        # self.lambda_mem = [membrane_dict['thermal conductivity'][0],
        #                    membrane_dict['thermal conductivity'][1]]
        # heat conductivity of the membrane

        # Setup membrane
        membrane_dict['width'] = self.cathode.width_straight_channels
        membrane_dict['length'] = self.cathode.length_straight_channels

        self.membrane = membrane.Membrane(membrane_dict, self.dx)

        self.thickness = self.cathode.thickness + self.membrane.thickness \
            + self.anode.thickness

        # Cell coordinates in z-direction (stacking/current direction)
        # will be initialized correctly through stack class
        self.coords = [0.0, 0.0]

        self.mem_base_r = cell_dict['mem_base_r']
        # basic electrical resistance of the membrane
        self.mem_acl_r = cell_dict['mem_acl_r']
        # thermal related electrical resistance gain of the membrane

        self.is_ht_pem = self.cell_dict['is_ht_pem']

        """heat conductivity along and through the cell layers"""
        self.width_straight_channels = self.cathode.width_straight_channels
        self.active_area_dx = self.width_straight_channels * self.dx

        # heat conductivity along the gas diffusion electrode and membrane
        self.th_layer = \
            np.asarray([self.cathode.th_bpp,
                        self.cathode.th_gde,
                        self.membrane.thickness,
                        self.anode.th_gde,
                        self.anode.th_bpp])

        self.thermal_conductance_z = \
            np.asarray([self.cathode.bpp.thermal_conductance[0],
                        self.cathode.gde.thermal_conductance[0],
                        self.membrane.thermal_conductance[0],
                        self.anode.gde.thermal_conductance[0],
                        self.anode.bpp.thermal_conductance[0]])

        self.thermal_conductance_x = \
            np.asarray([self.cathode.bpp.thermal_conductance[1],
                        self.cathode.gde.thermal_conductance[1],
                        self.membrane.thermal_conductance[1],
                        self.anode.gde.thermal_conductance[1],
                        self.anode.bpp.thermal_conductance[1]])

        self.thermal_conductance_x = \
            (self.thermal_conductance_x
             + np.roll(self.thermal_conductance_x, 1, axis=0)) * 0.5
        self.thermal_conductance_x = \
            np.vstack((self.thermal_conductance_x,
                       self.thermal_conductance_x[0]))
        if self.first_cell:
            self.thermal_conductance_x[0] *= 0.5
        if self.last_cell:
            self.thermal_conductance_x[-1] *= 0.5

        if self.last_cell:
            self.heat_cond_mtx = \
                mtx.build_cell_conductance_matrix(self.thermal_conductance_x,
                                                  self.thermal_conductance_z,
                                                  n_ele)
        else:
            self.heat_cond_mtx = \
                mtx.build_cell_conductance_matrix(self.thermal_conductance_x[:-1],
                                                  self.thermal_conductance_z[:-1],
                                                  n_ele)
        self.heat_mtx_const = self.heat_cond_mtx.copy()
        self.heat_mtx_dyn = np.zeros(self.heat_mtx_const.shape)
        self.heat_mtx = self.heat_mtx_dyn.copy()

        self.heat_rhs_const = np.full(self.n_layer * n_ele, 0.0)
        self.heat_rhs_dyn = np.zeros_like(self.heat_rhs_const)
        self.heat_rhs = self.heat_rhs_dyn.copy()

        # Create array for each thermal layer with indices according to
        # corresponding position in center diagonal of conductance matrix and
        # right hand side vector
        index_list = []
        for i in range(self.n_layer):
            index_list.append([(j * self.n_layer) + i for j in range(n_ele)])
        self.index_array = np.asarray(index_list)

        # Set constant thermal boundary conditions
        if self.first_cell:
            heat_bc = cell_dict['heat_pow']
            self.add_explicit_layer_source(self.heat_rhs_const, heat_bc, 0)
        if self.last_cell:
            heat_bc = cell_dict['heat_pow']
            self.add_explicit_layer_source(self.heat_rhs_const, heat_bc, -1)

        # Create electric conductance matrix
        self.elec_cond = \
            np.asarray([self.cathode.bpp.electrical_conductance[1],
                        self.anode.bpp.electrical_conductance[1]])
        self.elec_cond = \
            (self.elec_cond + np.roll(self.elec_cond, 1, axis=1)) * 0.5
        self.elec_cond = self.elec_cond[:, :-1]
        self.elec_x_mat_const = \
            mtx.build_z_cell_conductance_matrix(self.elec_cond.transpose(),
                                                len(self.elec_cond))
        # print(self.elec_x_mat_const)

        """boolean alarms"""
        self.v_alarm = False
        # True if :voltage loss > cell voltage
        self.break_program = False
        # True if the program aborts because of some critical impact

        """general parameter"""
        self.height = self.membrane.thickness \
            + self.cathode.th_bpp \
            + self.cathode.th_gde \
            + self.anode.th_bpp \
            + self.anode.th_gde

        self.v_loss = np.zeros(n_ele)
        # voltage loss
        self.temp_layer = \
            g_func.full((self.n_layer, n_ele), cell_dict['temp_init'])
        # layer temperature
        # coolant inlet temperature
        self.temp_names = ['Cathode BPP-BPP',
                           'Cathode BPP-GDE',
                           'Cathode GDE-MEM',
                           'Anode MEM-GDE',
                           'Anode GDE-BPP',
                           'Anode BPP-BPP']
        # interface names according to temperature array
        self.temp_mem = np.zeros(n_ele)
        # membrane temperature
        self.i_cd = np.zeros(n_ele)
        # current density
        self.v = np.zeros(n_ele)
        # cell voltage
        # self.resistance_z = np.zeros(n_ele)
        self.conductance_z = np.zeros(n_ele)
        # cell resistance

        self.add_print_data(self.i_cd, 'Current Density', 'A/m²')
        self.add_print_data(self.temp_layer, 'Temperature', 'K',
                            self.temp_names[:self.n_layer-1])

    def calc_ambient_conductance(self, alpha_amb):
        """
        :param alpha_amb: heat transfer coefficient for free or forced
        convection to the ambient
        :return: discretized conductance for each layer of cell based on
        external surface
        """
        # geometry model
        ext_surface_factor = (self.width + self.length) \
            / (self.cathode.channel.length * self.width_straight_channels)
        # k_amb = np.full((self.n_layer, self.n_ele), 0.)
        # convection conductance to the environment
        th_layer_amb = (self.th_layer + np.roll(self.th_layer, 1)) * 0.5
        # if self.last_cell:
        th_layer_amb = np.hstack((th_layer_amb, th_layer_amb[0]))
        k_amb = np.outer(th_layer_amb, self.dx) * alpha_amb / ext_surface_factor
        if self.first_cell:
            k_amb[0] *= 0.5
        if self.last_cell:
            k_amb[-1] *= 0.5
        return k_amb

    def add_explicit_layer_source(self, rhs_vector, source_term,
                                  layer_id=None):
        if layer_id is None:
            if np.isscalar(source_term):
                source_vector = np.full_like(rhs_vector, -source_term)
            else:
                source_vector = np.asarray(-source_term)
        else:
            source_vector = np.zeros_like(rhs_vector)
            np.put(source_vector, self.index_array[layer_id], -source_term)
        rhs_vector += source_vector
        return rhs_vector, source_vector

    def add_implicit_layer_source(self, matrix, coefficient, layer_id=None):
        matrix_size = matrix.shape[0]
        if layer_id is None:
            if np.isscalar(coefficient):
                source_vector = np.full(matrix_size, coefficient)
            else:
                source_vector = np.asarray(coefficient)
        else:
            source_vector = np.zeros(matrix_size)
            np.put(source_vector, self.index_array[layer_id], coefficient)
        matrix += np.diag(source_vector)
        return matrix, source_vector

    def update(self, current_density, channel_update=False):
        """
        This function coordinates the program sequence
        """
        self.i_cd[:] = current_density
        # self.temp_mem[:] = .5 * (self.temp_layer[2] + self.temp_layer[3])
        self.membrane.temp = .5 * (self.temp_layer[2] + self.temp_layer[3])
        if isinstance(self.membrane, membrane.WaterTransportMembrane):
            self.cathode.w_cross_flow[:] = self.membrane.w_cross_flow
            self.anode.w_cross_flow[:] = self.membrane.w_cross_flow
        # self.cathode.set_layer_temperature([self.temp[2], self.temp[3],
        #                                     self.temp[4]])
        # self.anode.set_layer_temperature([self.temp[0], self.temp[1]])
        self.cathode.update(self.i_cd, channel_update)
        self.anode.update(self.i_cd, channel_update)
        if self.anode.break_program or self.cathode.break_program:
            self.break_program = True
        else:
            # if not self.is_ht_pem:
            #     self.membrane.update()
            #     self.calc_mem_resistivity_springer()
            #     # self.calc_mem_resistivity_gossling()
            # else:
            #     self.calc_mem_resistivity_kvesic()
            # self.membrane.update()
            # self.calc_membrane_loss()
            humidity = np.asarray([self.cathode.channel.fluid.humidity,
                                   self.anode.channel.fluid.humidity])
            humidity_ele = \
                np.array([ip.interpolate_1d(humidity[0]),
                          ip.interpolate_1d(humidity[1])])
            self.membrane.update(self.i_cd, humidity_ele)
            self.calc_voltage_loss()
            self.calc_conductance()

    def calc_voltage_loss(self):
        """
        Calculates the cell voltage loss. If the cell voltage loss is larger
        than the open circuit cell voltage, the cell voltage is set to zero.
        """
        self.v_loss[:] = \
            self.membrane.v_loss + self.cathode.v_loss + self.anode.v_loss
        if np.any(self.v_loss) >= self.e_0:
            self.v_alarm = True
        self.v_loss[:] = np.minimum(self.v_loss, self.e_0)
        self.v[:] = self.e_0 - self.v_loss

    def calc_conductance(self):
        """
        Calculates the area-specific electrical resistance of the element in
        z-direction
        """
        current = self.i_cd * self.active_area_dx
        resistance_z = self.v_loss / current
        self.conductance_z[:] = 1 / resistance_z
