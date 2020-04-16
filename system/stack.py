import numpy as np
import data.global_parameters as g_par
import system.cell as cl
import system.electrical_coupling as el_cpl
import system.temperature_system as therm_cpl
import system.flow_circuit as flow_circuit
import system.channel as chl
import system.fluid as fluid
import data.input_dicts as in_dicts
import system.global_functions as g_func


class Stack:

    def __init__(self, current_control=False):

        # Read settings dictionaries
        stack_dict = in_dicts.dict_stack

        self.n_cells = stack_dict['cell_number']
        # number of cells of the stack
        n_nodes = g_par.dict_case['nodes']
        n_ele = n_nodes - 1
        # node points/elements along the x-axis
        self.calc_temp = stack_dict['calc_temperature']
        # switch to calculate the temperature distribution
        self.calc_electric = stack_dict['calc_current_density']
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
        # electrical_dict = in_dicts.dict_electrical_coupling
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
            # fluids.append(fluid.factory2(fluid_dicts[i]))
            channels.append([chl.Channel(channel_dicts[i],
                                         fluid.dict_factory(fluid_dicts[i]))
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
            cell = cl.Cell(cell_dict, membrane_dict, half_cell_dicts,
                           cell_channels, number=i)
            if i == 0:
                cell.coords[0] = 0.0
                cell.coords[1] = cell.thickness
            else:
                cell.coords[0] = self.cells[i - 1].coords[1]
                cell.coords[1] = cell.coords[0] + cell.thickness
            self.cells.append(cell)

        # Initialize flow circuits
        manifold_length = \
            self.cells[-1].coords[-1] - self.cells[0].coords[0]
        self.fuel_circuits = []
        for i in range(len(half_cell_dicts)):
            manifold_in_dicts[i]['length'] = manifold_length
            manifold_out_dicts[i]['length'] = manifold_length
            sub_channel_number = self.cells[0].half_cells[i].n_channel
            self.fuel_circuits.append(
                flow_circuit.factory2(flow_circuit_dicts[i],
                                      manifold_in_dicts[i],
                                      manifold_out_dicts[i],
                                      channels[i], sub_channel_number))

        cool_flow = stack_dict['cool_flow']
        if cool_flow:
            coolant_dict = in_dicts.dict_coolant_fluid
            in_dicts.dict_coolant_in_manifold['length'] = manifold_length
            in_dicts.dict_coolant_out_manifold['length'] = manifold_length
            coolant_dict['temp_in'] = \
                in_dicts.dict_coolant_in_manifold['temp_in']
            coolant_dict['p_out'] = in_dicts.dict_coolant_out_manifold['p_out']

            if temperature_dict['cool_ch_bc']:
                n_cool = self.n_cells + 1
            else:
                n_cool = self.n_cells - 1

            n_cool_cell = temperature_dict['cool_ch_numb']
            cool_channels = []
            for i in range(n_cool):
                cool_channels.append(
                    chl.Channel(in_dicts.dict_coolant_channel,
                                fluid.dict_factory(coolant_dict)))
                cool_channels[i].name += ' ' + str(i)
                cool_channels[i].fluid.name = \
                    cool_channels[i].name + ': ' \
                    + cool_channels[i].fluid.TYPE_NAME

            if n_cool > 0:
                self.coolant_circuit = \
                    flow_circuit.factory2(in_dicts.dict_coolant_flow_circuit,
                                          in_dicts.dict_coolant_in_manifold,
                                          in_dicts.dict_coolant_out_manifold,
                                          cool_channels, n_cool_cell)
            else:
                self.coolant_circuit = None
        else:
            self.coolant_circuit = None

        self.flow_circuits = \
            [self.fuel_circuits[0], self.fuel_circuits[1], self.coolant_circuit]

        self.current_control = current_control
        self.i_target = g_par.dict_case['target_current_density']
        self.target_cell_voltage = g_par.dict_case['average_cell_voltage']
        self.v_target = self.n_cells * self.target_cell_voltage

        # Initialize the electrical coupling
        self.elec_sys = el_cpl.ElectricalCoupling(self)
        # Initialize temperature system
        self.temp_sys = therm_cpl.TemperatureSystem(self, temperature_dict)

        """boolean alarms"""
        self.v_alarm = False
        # True if :voltage loss > cell voltage
        self.break_program = False
        # True if the program aborts because of some critical impact

        # target current density

        # current density array
        self.i_cd = np.zeros((self.n_cells, n_ele))
        # self.i_cd[:] = self.i_target
        for i in range(self.n_cells):
            self.i_cd[i, :] = \
                g_func.exponential_distribution(self.i_target, n_ele,
                                                a=2.0, b=0.0)
        # current density array of previous iteration step
        self.i_cd_old = np.copy(self.i_cd)
        # voltage array
        self.v = np.zeros(self.n_cells)
        self.v_stack = 0.0

        # old temperature for convergence calculation
        self.temp_old = np.zeros(self.temp_sys.temp_layer_vec.shape)
        self.temp_old[:] = self.temp_sys.temp_layer_vec

    def update(self, current_density=None, voltage=None):
        """
        This function coordinates the program sequence
        """
        update_inflows = False
        if current_density is not None:
            self.i_cd[:] = current_density
            update_inflows = True
        elif voltage is not None:
            self.v_stack = voltage
            update_inflows = True
        self.update_flows(update_inflows, self.calc_flow_dis)
        for i, cell in enumerate(self.cells):
            cell.update(self.i_cd[i, :], check_stoi=self.current_control)
            if cell.break_program:
                self.break_program = True
                break
        self.i_cd_old[:] = self.elec_sys.i_cd
        self.temp_old[:] = self.temp_sys.temp_layer_vec
        if not self.break_program:
            if self.calc_temp:
                self.temp_sys.update()
            if self.calc_electric:
                self.elec_sys.update(current_density=current_density,
                                     voltage=voltage)
                self.i_cd[:] = self.elec_sys.i_cd
            self.v[:] = \
                np.asarray([np.average(cell.v, weights=cell.active_area_dx)
                            for cell in self.cells])
            if self.current_control:
                self.v_stack = np.sum(self.v)

    def update_flows(self, update_inflows=False, calc_flow_dist=False):
        """
        This function updates the flow distribution of gas over the stack cells
        """
        mass_flows_in = [None, None]
        if update_inflows:
            mass_flows_in[:] = self.calc_mass_flows()
        for i in range(len(self.fuel_circuits)):
            self.fuel_circuits[i].update(mass_flows_in[i], calc_flow_dist)
        dtemp_target = 10.0
        if self.coolant_circuit is not None:
            if update_inflows:
                n_cool_cell = self.coolant_circuit.n_subchannels
                over_potential = 0.5
                heat = self.i_target * over_potential  \
                    * self.cells[0].active_area * self.n_cells * n_cool_cell
                cp_cool = \
                    np.average([np.average(channel.fluid.specific_heat)
                                for channel in self.coolant_circuit.channels])
                cool_mass_flow = heat / (cp_cool * dtemp_target)
            else:
                # id_in = self.coolant_circuit.manifolds[0].id_in
                # temp_in = self.coolant_circuit.manifolds[0].temp[id_in]
                # temp_out = \
                #   np.sum([channel.temp[channel.id_out]
                #           * channel.g_fluid[channel.id_out]
                #           for channel in self.coolant_circuit.channels])
                # temp_out /= \
                #   np.sum([channel.g_fluid[channel.id_out]
                #           for channel in self.coolant_circuit.channels])
                # temp_ratio = abs(temp_out - temp_in) / dtemp_target
                # cool_mass_flow = \
                #   self.coolant_circuit.mass_flow_in * temp_ratio
                cool_mass_flow = None
            self.coolant_circuit.update(cool_mass_flow, calc_flow_dist)

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

