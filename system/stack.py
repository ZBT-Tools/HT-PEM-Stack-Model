import numpy as np
import copy as copy
import data.global_parameters as g_par
import system.cell as cl
import system.manifold as mfd
import system.electrical_coupling as el_cpl
import data.electrical_coupling_dict as el_cpl_dict
import system.temperature_system as therm_cpl
import system.flow_circuit2 as flow_circuit
import system.channel as chl
import system.fluid2 as fluid
import data.temperature_system_dict as therm_dict
import data.input_dicts as in_dicts


class Stack:

    def __init__(self):

        # Read input dictionaries
        stack_dict = in_dicts.dict_stack

        # Handover
        self.n_cells = stack_dict['cell_number']
        # number of cells of the stack

        # global stoichiometry for cathode and anode
        # self.stoi = [stack_dict['stoi_cat'], stack_dict['stoi_ano']]
        # # inlet stoichiometry of the cathode header
        # # inlet stoichiometry of the anode header
        n_nodes = g_par.dict_case['nodes']
        # node points along the x-axis
        self.calc_temp = stack_dict['calc_temperature']
        # switch to calculate the temperature distribution
        self.calc_cd = stack_dict['calc_current_density']
        # switch to calculate the current density distribution
        self.calc_flow_dis = stack_dict['calc_flow_distribution']
        # switch to calculate the flow distribution

        cell_dict = in_dicts.dict_cell
        membrane_dict = in_dicts.dict_membrane
        anode_dict = in_dicts.dict_anode
        cathode_dict = in_dicts.dict_cathode
        ano_channel_dict = in_dicts.dict_anode_channel
        cat_channel_dict = in_dicts.dict_cathode_channel
        ano_fluid_dict = in_dicts.dict_anode_fluid
        cat_fluid_dict = in_dicts.dict_cathode_fluid
        ano_in_manifold_dict = in_dicts.dict_anode_in_manifold
        cat_in_manifold_dict = in_dicts.dict_cathode_in_manifold
        ano_out_manifold_dict = in_dicts.dict_anode_out_manifold
        cat_out_manifold_dict = in_dicts.dict_cathode_out_manifold
        electrical_dict = in_dicts.dict_electrical_coupling
        temperature_dict = in_dicts.dict_temp_sys

        half_cell_dicts = [cathode_dict, anode_dict]
        channel_dicts = [cat_channel_dict, ano_channel_dict]
        fluid_dicts = [cat_fluid_dict, ano_fluid_dict]
        manifold_in_dicts = [cat_in_manifold_dict, ano_in_manifold_dict]
        manifold_out_dicts = [cat_out_manifold_dict, ano_out_manifold_dict]
        flow_circuit_dicts = [in_dicts.dict_cathode_flow_circuit,
                              in_dicts.dict_anode_flow_circuit]

        # Initialize fluid channels
        fluids, channels = [], []
        for i in range(len(half_cell_dicts)):
            fluid_dicts[i]['temp_in'] = manifold_in_dicts[i]['temp_in']
            fluid_dicts[i]['p_out'] = manifold_out_dicts[i]['p_out']
            fluids.append(fluid.factory2(fluid_dicts[i]))
            channels.append([chl.Channel(channel_dicts[i],
                                         copy.deepcopy(fluids[i]))
                             for j in range(self.n_cells)])

        # Initialize fuel cells
        self.cells = []
        for i in range(self.n_cells):
            if self.n_cells == 1:
                cell_dict['first_cell'] = True
                cell_dict['last_cell'] = True
                cell_dict['heat_pow'] = temperature_dict['heat_pow']
            elif i == 0:
                cell_dict['first_cell'] = True
                cell_dict['last_cell'] = False
                cell_dict['heat_pow'] = temperature_dict['heat_pow']
            elif i == self.n_cells-1:
                cell_dict['first_cell'] = False
                cell_dict['last_cell'] = True
            else:
                cell_dict['first_cell'] = False
                cell_dict['last_cell'] = False
                cell_dict['heat_pow'] = temperature_dict['heat_pow']

            cell_channels = [channels[0][i], channels[1][i]]
            # Cell constructor
            cell = cl.Cell(i, cell_dict, membrane_dict, half_cell_dicts,
                           cell_channels)
            if i == 0:
                cell.coords[0] = 0.0
                cell.coords[1] = cell.thickness
            else:
                cell.coords[0] = self.cells[i - 1].coords[1]
                cell.coords[1] = cell.coords[0] + cell.thickness
            self.cells.append(cell)

        # self.set_stoichiometry(np.full(self.n_cells, self.stoi_cat),
        #                        np.full(self.n_cells, self.stoi_ano))

        # Initialize flow circuits
        manifold_length = \
            self.cells[-1].coords[-1] - self.cells[0].coords[0]
        self.fluid_circuits = []
        for i in range(len(half_cell_dicts)):
            manifold_in_dicts[i]['length'] = manifold_length
            manifold_out_dicts[i]['length'] = manifold_length
            sub_channel_number = self.cells[0].half_cells[i].n_channel
            self.fluid_circuits.append(
                flow_circuit.factory2(flow_circuit_dicts[i],
                                      manifold_in_dicts[i],
                                      manifold_out_dicts[i],
                                      channels[i], sub_channel_number))

        coolant_dict = in_dicts.dict_coolant_fluid
        in_dicts.dict_coolant_in_manifold['length'] = manifold_length
        in_dicts.dict_coolant_out_manifold['length'] = manifold_length
        coolant_dict['temp_in'] = in_dicts.dict_coolant_in_manifold['temp_in']
        coolant_dict['p_out'] = in_dicts.dict_coolant_out_manifold['p_out']
        coolant = fluid.factory2(coolant_dict)

        if temperature_dict['cool_ch_bc']:
            n_cool = self.n_cells + 1
        else:
            n_cool = self.n_cells - 1

        n_cool_cell = temperature_dict['cool_ch_numb']

        cool_channels = [chl.Channel(in_dicts.dict_coolant_channel,
                                     copy.deepcopy(coolant))
                         for i in range(n_cool)]
        self.coolant_circuit = \
            flow_circuit.factory2(in_dicts.dict_coolant_flow_circuit,
                                  in_dicts.dict_coolant_in_manifold,
                                  in_dicts.dict_coolant_out_manifold,
                                  cool_channels, n_cool_cell)

        # Initialize the electrical coupling
        # if self.n_cells > 1:
        self.elec_sys = \
            el_cpl.ElectricalCoupling(electrical_dict, self, self.cells)

        # Initialize temperature system
        self.temp_sys = therm_cpl.TemperatureSystem(temperature_dict,
                                                    self.cells, cool_channels)

        """boolean alarms"""
        self.v_alarm = False
        # True if :voltage loss > cell voltage
        self.break_program = False
        # True if the program aborts because of some critical impact

        """General data"""
        n_ele = n_nodes - 1

        # target current density
        self.i_target = g_par.dict_case['target_current_density']

        # current density array
        self.i_cd = np.zeros((self.n_cells, n_ele))
        self.i_cd[:] = self.i_target

        # current density array of previous iteration step
        self.i_cd_old = np.copy(self.i_cd)

        self.temp_old = np.zeros(self.temp_sys.temp_layer_vec.shape)
        self.temp_old[:] = self.temp_sys.temp_layer_vec

    def update(self, current_density=None):
        """
        This function coordinates the program sequence
        """
        update_inflows = False
        if current_density is not None:
            self.i_cd[:] = current_density
            update_inflows = True
        self.update_flows(update_inflows, self.calc_flow_dis)
        for i, cell in enumerate(self.cells):
            cell.update(self.i_cd[i, :])
            if cell.break_program:
                self.break_program = True
                break
        self.i_cd_old[:] = self.i_cd
        self.temp_old[:] = self.temp_sys.temp_layer_vec
        if not self.break_program:
            if self.calc_temp:
                self.update_temperature_coupling()
            if self.calc_cd:
                self.update_electrical_coupling()

    def update_flows(self, update_inflows=False, calc_flow_dist=False):
        """
        This function updates the flow distribution of gas over the stack cells
        """
        mass_flows_in = [None, None]
        if update_inflows:
            mass_flows_in[:] = self.calc_mass_flows()
        for i in range(len(self.fluid_circuits)):
            self.fluid_circuits[i].update(mass_flows_in[i], calc_flow_dist)
        dtemp_target = 10.0
        if update_inflows:
            n_cool_cell = self.coolant_circuit.n_subchannels
            over_potential = 0.5
            heat = self.i_target * over_potential  \
                * self.cells[0].active_area * self.n_cells * n_cool_cell
            cp_cool = np.average([np.average(channel.fluid.specific_heat)
                                  for channel in self.coolant_circuit.channels])
            cool_mass_flow = heat / (cp_cool * dtemp_target)
        else:
            # id_in = self.coolant_circuit.manifolds[0].id_in
            # temp_in = self.coolant_circuit.manifolds[0].temp[id_in]
            # temp_out = np.sum([channel.temp[channel.id_out]
            #                    * channel.g_fluid[channel.id_out]
            #                    for channel in self.coolant_circuit.channels])
            # temp_out /= np.sum([channel.g_fluid[channel.id_out]
            #                     for channel in self.coolant_circuit.channels])
            # temp_ratio = abs(temp_out - temp_in) / dtemp_target
            # cool_mass_flow = self.coolant_circuit.mass_flow_in * temp_ratio
            cool_mass_flow = None
        self.coolant_circuit.update(cool_mass_flow, calc_flow_dist)

    def update_electrical_coupling(self):
        """
        This function updates current distribution over the stack cells
        """
        # self.el_cpl_stack.update_values(self.stack_cell_r, self.v_loss)
        # if self.n_cells > 1:
        self.elec_sys.update()
        self.i_cd[:] = self.elec_sys.i_cd[:]
        # else:
        #     self.i_cd[0] = self.cells[0].i_cd

    def update_temperature_coupling(self):
        """
        This function updates the layer and fluid temperatures of the stack
        """
        self.temp_sys.update()

    def calc_mass_flows(self):
        mass_flows_in = []
        for i in range(len(self.cells[0].half_cells)):
            cell_mass_flow, cell_mole_flow = \
                self.cells[0].half_cells[i].calc_inlet_flow()
            cell_mass_flow = np.sum(cell_mass_flow, axis=0)
            mass_flow = cell_mass_flow \
                * self.cells[0].half_cells[i].n_channel * self.n_cells

            mass_flows_in.append(mass_flow)

        return mass_flows_in

