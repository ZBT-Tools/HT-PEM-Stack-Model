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
        # switch to calculate the flow_circuit.py distribution

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
        # Initialize the electrical coupling
        if self.n_cells > 1:
            self.elec_sys = \
                el_cpl.ElectricalCoupling(electrical_dict, self, self.cells)

        # Initialize temperature system
        self.temp_sys = therm_cpl.TemperatureSystem(temperature_dict,
                                                    self.cells)

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

        cool_channels = [chl.Channel(in_dicts.dict_coolant_channel,
                                     copy.deepcopy(coolant))
                         for i in range(self.temp_sys.n_cool)]
        self.coolant_circuit = \
            flow_circuit.factory2(in_dicts.dict_coolant_flow_circuit,
                                  in_dicts.dict_coolant_in_manifold,
                                  in_dicts.dict_coolant_out_manifold,
                                  cool_channels,
                                  self.temp_sys.n_cool_cell)

        # # Initialize the manifolds
        # cathode_channels = \
        #     [cell.half_cells[0].channel for cell in self.cells]
        # anode_channels = \
        #     [cell.half_cells[1].channel for cell in self.cells]
        # cathode_channel_multiplier = self.cells[0].half_cells[0].n_channel
        # anode_channel_multiplier = self.cells[0].half_cells[1].n_channel

        # self.manifolds = \
        #     [mfd.Manifold(mfd_dicts[0],
        #                   self.cells, cathode_channels,
        #                   cathode_channel_multiplier),
        #      mfd.Manifold(mfd_dicts[1],
        #                   self.cells, anode_channels,
        #                   anode_channel_multiplier)]
        # self.manifolds[0].head_stoi = self.cells[0].cathode.stoi
        # self.manifolds[1].head_stoi = self.cells[0].anode.stoi

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
        self.i_cd = np.full((self.n_cells, n_ele), self.i_target)

        # current density array of previous iteration step
        self.i_cd_old = np.copy(self.i_cd)

        # self.v_cell = []
        # # cell voltage
        # self.v_loss = []
        # # cell voltage loss
        # self.v_loss_cat = []
        # # cathode voltage loss
        # self.v_loss_ano = []
        # # anode voltage loss
        # self.stack_cell_r = []
        # # cell resistance in z-direction
        # self.vol_flow_gas_cat = []
        # # molar flow at the inlet and outlet of the cathode channels
        # self.vol_flow_gas_ano = []
        # # molar flow  at the inlet and outlet of the anode channels
        # self.m_sum_cat = []
        # # mass flow at the inlet and outlet of the cathode channels
        # self.m_sum_ano = []
        # # mass flow at the inlet and outlet of the anode channels
        # self.cp_cat = []
        # # heat capacity in the cathode channels
        # self.cp_ano = []
        # # heat capacity in the anode channels
        # self.visc_cat = []
        # # viscosity at the inlet and outlet of the cathode channels
        # self.visc_ano = []
        # # viscosity at the inlet and outlet of the anode channels
        # self.p_cat = []
        # # pressure at the inlet and outlet of the cathode channels
        # self.p_ano = []
        # # pressure at the inlet and outlet of the anode channels
        # self.r_cat = []
        # # gas constant of air fluid at the inlet
        # # and outlet of the cathode channels
        # self.r_ano = []
        # # gas constant of the hydrogen fluid at the
        # # inlet and outlet of the anode channels
        # self.temp_fluid_cat = []
        # # inlet and outlet temperature of the cathode channel fluid
        # self.temp_fluid_ano = np.zeros(self.n_cells)
        # # inlet and outlet temperature of the anode channel fluid

    def update(self):
        """
        This function coordinates the program sequence
        """
        for i, cell in enumerate(self.cells):
            # self.cells[j].set_current_density(self.i_cd[j, :])
            cell.i_cd[:] = self.i_cd[i, :]
            cell.update()
            if cell.break_program:
                self.break_program = True
                break
        if not self.break_program:
            # if self.n_cells > 1:
            #     if self.calc_flow_dis:
            self.update_flows()
            # self.stack_dynamic_properties()
            if self.calc_temp:
                self.update_temperature_coupling()

            self.i_cd_old[:] = self.i_cd
            if self.calc_cd:
                self.update_electrical_coupling()
        # print('Current density', self.i_cd)

    def update_flows(self):
        """
        This function updates the flow distribution of gas over the stack cells
        """
        mass_flows_in = self.calc_mass_flows()
        for i in range(len(self.fluid_circuits)):
            self.fluid_circuits[i].update(mass_flows_in[i])

        self.coolant_circuit.update()

        # n_ch = self.cells[0].cathode.n_channel
        # self.manifolds[0].update_values(
        #     self.vol_flow_gas_cat * n_ch, self.temp_fluid_cat,
        #     self.cp_cat, self.visc_cat, self.p_cat, self.r_cat,
        #     self.m_sum_f_cat * n_ch,
        #     self.m_sum_g_cat * n_ch)
        # n_ch = self.cells[0].anode.n_channel
        # self.manifolds[1].update_values(
        #     self.vol_flow_gas_ano[::-1] * n_ch, self.temp_fluid_ano[::-1],
        #     self.cp_ano[::-1], self.visc_ano[::-1], self.p_ano[::-1],
        #     self.r_ano[::-1], self.m_sum_f_ano[::-1] * n_ch,
        #     self.m_sum_g_ano[::-1] * n_ch)
        # self.manifolds[0].update()
        # self.manifolds[1].update()
        # self.set_stoichiometry(self.manifolds[0].cell_stoi,
        #                        self.manifolds[1].cell_stoi)
        # self.set_channel_outlet_pressure(self.manifolds[0].head_p[-1],
        #                                  self.manifolds[1].head_p[-1])

    def update_electrical_coupling(self):
        """
        This function updates current distribution over the stack cells
        """
        # self.el_cpl_stack.update_values(self.stack_cell_r, self.v_loss)
        if self.n_cells > 1:
            self.elec_sys.update()
            self.i_cd[:] = self.elec_sys.i_cd[:]
        else:
            self.i_cd[0] = self.cells[0].i_cd

    def update_temperature_coupling(self):
        """
        This function updates the layer and fluid temperatures of the stack
        """
        # current = self.i_cd * self.cells[0].active_area_dx
        # n_ch = self.cells[0].cathode.n_channel
        # self.temp_sys.update_values(self.k_alpha_ch,
        #                             self.cond_rate,
        #                             self.omega,
        #                             np.array([self.v_loss_cat,
        #                                       self.v_loss_ano]),
        #                             self.g_fluid, current)
        self.temp_sys.update()
        # self.set_temperature()

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

    # def stack_dynamic_properties(self):
    #     """
    #     This function sums up the dynamic values inside the cells
    #     necessary to calculate the flow_circuit.py distribution,
    #     the electrical coupling or the temperature coupling
    #     """
    #     v_alarm = []
    #     k_alpha_cat, k_alpha_ano = [], []
    #     g_fluid_cat, g_fluid_ano = [], []
    #     v_cell, v_loss, resistance = [], [], []
    #     v_loss_cat, v_loss_ano = [], []
    #     vol_flow_gas_cat_in, vol_flow_gas_cat_out = [], []
    #     vol_flow_gas_ano_in, vol_flow_gas_ano_out = [], []
    #     cp_cat_in, cp_cat_out = [], []
    #     cp_ano_in, cp_ano_out = [], []
    #     visc_cat_in, visc_cat_out = [], []
    #     visc_ano_in, visc_ano_out = [], []
    #     p_cat_in, p_cat_out = [], []
    #     p_ano_in, p_ano_out = [], []
    #     r_cat_in, r_cat_out = [], []
    #     r_ano_in, r_ano_out = [], []
    #     temp_fluid_cat_in, temp_fluid_cat_out = [], []
    #     temp_fluid_ano_in, temp_fluid_ano_out = [], []
    #     cond_rate_cat, cond_rate_ano = [], []
    #     m_sum_f_cat_in, m_sum_f_cat_out = [], []
    #     m_sum_f_ano_in, m_sum_f_ano_out = [], []
    #     m_sum_g_cat_in, m_sum_g_cat_out = [], []
    #     m_sum_g_ano_in, m_sum_g_ano_out = [], []
    #     omega = []
    #     for cell in self.cells:
    #         v_alarm.append(cell.v_alarm)
    #         cond_rate_cat.append(cell.cathode.channel.cond_rate)
    #         cond_rate_ano.append(cell.anode.channel.cond_rate)
    #         k_alpha_cat.append(cell.cathode.channel.k_coeff)
    #         k_alpha_ano.append(cell.anode.channel.k_coeff)
    #         g_fluid_cat.append(cell.cathode.channel.g_fluid)
    #         g_fluid_ano.append(cell.anode.channel.g_fluid)
    #         v_cell.append(cell.v)
    #         v_loss = np.hstack((v_loss, cell.v_loss))
    #         v_loss_cat.append(cell.cathode.v_loss)
    #         v_loss_ano.append(cell.anode.v_loss)
    #         # omega.append(cell.omega)
    #         resistance = np.hstack((resistance, cell.resistance))
    #         vol_flow_gas_cat_in = \
    #             np.hstack((vol_flow_gas_cat_in,
    #                        cell.cathode.channel.vol_flow_gas[0]))
    #         vol_flow_gas_cat_out = \
    #             np.hstack((vol_flow_gas_cat_out,
    #                        cell.cathode.channel.vol_flow_gas[-1]))
    #         vol_flow_gas_ano_in = \
    #             np.hstack((vol_flow_gas_ano_in,
    #                        cell.anode.channel.vol_flow_gas[0]))
    #         vol_flow_gas_ano_out = \
    #             np.hstack((vol_flow_gas_ano_out,
    #                        cell.anode.channel.vol_flow_gas[-1]))
    #         m_sum_f_cat_in = np.hstack((m_sum_f_cat_in,
    #                                     cell.cathode.channel.mass_flow_total[
    #                                         0]))
    #         m_sum_f_cat_out = np.hstack((m_sum_f_cat_out,
    #                                      cell.cathode.channel.mass_flow_total[-1]))
    #         m_sum_f_ano_in = np.hstack((m_sum_f_ano_in,
    #                                     cell.anode.channel.mass_flow_total[0]))
    #         m_sum_f_ano_out = np.hstack((m_sum_f_ano_out,
    #                                      cell.anode.channel.mass_flow_total[
    #                                          -1]))
    #         m_sum_g_cat_in = np.hstack((m_sum_g_cat_in,
    #                                     cell.cathode.channel.mass_flow_gas_total[0]))
    #         m_sum_g_cat_out = np.hstack((m_sum_g_cat_out,
    #                                      cell.cathode.channel.mass_flow_gas_total[-1]))
    #         m_sum_g_ano_in = np.hstack((m_sum_g_ano_in,
    #                                     cell.anode.channel.mass_flow_gas_total[0]))
    #         m_sum_g_ano_out = np.hstack((m_sum_g_ano_out,
    #                                      cell.anode.channel.mass_flow_gas_total[-1]))
    #         cp_cat_in = np.hstack((cp_cat_in,
    #                                cell.cathode.channel.fluid.cp_fluid[0]))
    #         cp_cat_out = np.hstack((cp_cat_out,
    #                                 cell.cathode.channel.fluid.cp_fluid[-1]))
    #         cp_ano_in = np.hstack((cp_ano_in,
    #                                cell.anode.channel.fluid.cp_fluid[0]))
    #         cp_ano_out = np.hstack((cp_ano_out,
    #                                 cell.anode.channel.fluid.cp_fluid[-1]))
    #         p_cat_in = np.hstack((p_cat_in, cell.cathode.channel.p[0]))
    #         p_cat_out = np.hstack((p_cat_out, cell.cathode.channel.p[-1]))
    #         p_ano_in = np.hstack((p_ano_in, cell.anode.channel.p[0]))
    #         p_ano_out = np.hstack((p_ano_out, cell.anode.channel.p[-1]))
    #         r_cat_in = np.hstack((r_cat_in, cell.cathode.r_gas[0]))
    #         r_cat_out = np.hstack((r_cat_out, cell.cathode.r_gas[-1]))
    #         r_ano_in = np.hstack((r_ano_in, cell.anode.r_gas[0]))
    #         r_ano_out = np.hstack((r_ano_out, cell.anode.r_gas[-1]))
    #         visc_cat_in = np.hstack((visc_cat_in, cell.cathode.visc_gas[0]))
    #         visc_cat_out =
    #             np.hstack((visc_cat_out, cell.cathode.visc_gas[-1]))
    #         visc_ano_in = np.hstack((visc_ano_in, cell.anode.visc_gas[0]))
    #         visc_ano_out = np.hstack((visc_ano_out, cell.anode.visc_gas[-1]))
    #         temp_fluid_cat_in = np.hstack((temp_fluid_cat_in,
    #                                        cell.cathode.temp_fluid[0]))
    #         temp_fluid_cat_out = np.hstack((temp_fluid_cat_out,
    #                                         cell.cathode.temp_fluid[-1]))
    #         temp_fluid_ano_in = np.hstack((temp_fluid_ano_in,
    #                                        cell.anode.temp_fluid[0]))
    #         temp_fluid_ano_out = np.hstack((temp_fluid_ano_out,
    #                                         cell.anode.temp_fluid[-1]))
    #     self.v_cell, self.v_loss, self.stack_cell_r =
    #         v_cell, v_loss, resistance
    #     self.v_loss_cat, self.v_loss_ano = v_loss_cat, v_loss_ano
    #     self.vol_flow_gas_cat = \
    #         np.array([vol_flow_gas_cat_in, vol_flow_gas_cat_out])
    #     self.vol_flow_gas_ano = \
    #         np.array([vol_flow_gas_ano_in, vol_flow_gas_ano_out])
    #     self.m_sum_f_cat = np.array([m_sum_f_cat_in, m_sum_f_cat_out])
    #     self.m_sum_f_ano = np.array([m_sum_f_ano_in, m_sum_f_ano_out])
    #     self.m_sum_g_cat = np.array([m_sum_g_cat_in, m_sum_g_cat_out])
    #     self.m_sum_g_ano = np.array([m_sum_g_ano_in, m_sum_g_ano_out])
    #     self.cp_cat = np.array([cp_cat_in, cp_cat_out])
    #     self.cp_ano = np.array([cp_ano_in, cp_ano_out])
    #     self.visc_cat = np.array([visc_cat_in, visc_cat_out])
    #     self.visc_ano = np.array([visc_ano_in, visc_ano_out])
    #     self.p_cat = np.array([p_cat_in, p_cat_out])
    #     self.p_ano = np.array([p_ano_in, p_ano_out])
    #     self.r_cat = np.array([r_cat_in, r_cat_out])
    #     self.r_ano = np.array([r_ano_in, r_ano_out])
    #     self.temp_fluid_cat =
    #         np.array([temp_fluid_cat_in, temp_fluid_cat_out])
    #     self.temp_fluid_ano =
    #         np.array([temp_fluid_ano_in, temp_fluid_ano_out])
    #     self.v_alarm = np.array(v_alarm)
    #
    # def set_stoichiometry(self, stoi_cat, stoi_ano):
    #     """
    #     This function sets up the inlet stoichiometry
    #     of the cathode and anode channels.
    #     """
    #     for i, cell in enumerate(self.cells):
    #         cell.cathode.stoi = stoi_cat[i]
    #         cell.anode.stoi = stoi_ano[i]
    #
    # def set_channel_outlet_pressure(self, p_cat, p_ano):
    #     """
    #     This function sets up the inlet pressure
    #     of the cathode and the anode channels.
    #     """
    #     for i, cell in enumerate(self.cells):
    #         cell.cathode.channel.p_out = p_cat[i]
    #         cell.anode.channel.p_out = p_ano[i]

    # def set_temperature(self):
    #     """
    #     This function sets up the layer and fluid temperatures in the cells.
    #     """
    #     for i, cell in enumerate(self.cells):
    #         cell.temp_layer[:] = self.temp_sys.temp_layer[i][:, :]
    #         cell.cathode.temp_fluid[:] = self.temp_sys.temp_fluid[0, i]
    #         cell.anode.temp_fluid[:] = self.temp_sys.temp_fluid[1, i]
    #         cell.temp_cool[:] = self.temp_sys.temp_cool_ele[i]
