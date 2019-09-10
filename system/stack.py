import numpy as np
import copy as copy
import data.global_parameters as g_par
import system.cell as cl
import system.manifold as m_fold
import system.electrical_coupling as el_cpl
import data.electrical_coupling_dict as el_cpl_dict
import system.temperature_system as therm_cpl
import data.temperature_system_dict as therm_dict


class Stack:

    def __init__(self, stack_dict, cell_dict, anode_dict, cathode_dict,
                 ano_channel_dict, cat_channel_dict, ano_manifold_dict,
                 cat_manifold_dict, electrical_dict, temperature_dict):
        # Handover
        self.n_cells = stack_dict['cell_number']
        # number of cells of the stack
        # self.stoi_cat = stack_dict['stoi_cat']
        # # inlet stoichiometry of the cathode header
        # self.stoi_ano = stack_dict['stoi_ano']
        # # inlet stoichiometry of the anode header
        n_nodes = g_par.dict_case['nodes']
        # node points along the x-axis
        self.calc_temp = stack_dict['calc_temperature']
        # switch to calculate the temperature distribution
        self.calc_cd = stack_dict['calc_current_density']
        # switch to calculate the current density distribution
        self.calc_flow_dis = stack_dict['calc_flow_distribution']
        # switch to calculate the flow distribution
        self.cells = []
        # Initialize individual cells
        for i in range(self.n_cells):
            if i == 0:
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
            self.cells.append(cl.Cell(i, cell_dict, anode_dict, cathode_dict,
                                      ano_channel_dict, cat_channel_dict))

        # self.set_stoichiometry(np.full(self.n_cells, self.stoi_cat),
        #                        np.full(self.n_cells, self.stoi_ano))

        # Initialize the manifolds
        self.manifold = [m_fold.Manifold(cat_manifold_dict),
                         m_fold.Manifold(ano_manifold_dict)]
        self.manifold[0].head_stoi = self.cells[0].cathode.stoi
        self.manifold[1].head_stoi = self.cells[0].anode.stoi

        # Initialize the electrical coupling
        self.el_cpl_stack = el_cpl.ElectricalCoupling(electrical_dict)

        """boolean alarms"""
        self.v_alarm = False
        # True if :voltage loss > cell voltage
        self.break_program = False
        # True if the program aborts because of some critical impact

        """General data"""
        self.i_cd = np.full((self.n_cells, n_nodes - 1),
                            g_par.dict_case['tar_cd'])
        # current density
        self.i_cd_old = np.copy(self.i_cd)
        # current density of the last iteration
        #self.i = np.full((self.n_cells, n_nodes - 1), 20.)
        # current
        self.v_cell = []
        # cell voltage
        self.v_loss = []
        # cell voltage loss
        self.v_loss_cat = []
        # cathode voltage loss
        self.v_loss_ano = []
        # anode voltage loss
        self.stack_cell_r = []
        # cell resistance in z-direction
        self.q_sum_cat = []
        # molar flow at the inlet and outlet of the cathode channels
        self.q_sum_ano = []
        # molar flow  at the inlet and outlet of the anode channels
        self.m_sum_cat = []
        # mass flow at the inlet and outlet of the cathode channels
        self.m_sum_ano = []
        # mass flow at the inlet and outlet of the anode channels
        self.cp_cat = []
        # heat capacity in the cathode channels
        self.cp_ano = []
        # heat capacity in the anode channels
        self.visc_cat = []
        # viscosity at the inlet and outlet of the cathode channels
        self.visc_ano = []
        # viscosity at the inlet and outlet of the anode channels
        self.p_cat = []
        # pressure at the inlet and outlet of the cathode channels
        self.p_ano = []
        # pressure at the inlet and outlet of the anode channels
        self.r_cat = []
        # gas constant of air fluid at the inlet
        # and outlet of the cathode channels
        self.r_ano = []
        # gas constant of the hydrogen fluid at the
        # inlet and outlet of the anode channels
        self.temp_fluid_cat = []
        # inlet and outlet temperature of the cathode channel fluid
        self.temp_fluid_ano = np.zeros(self.n_cells)
        # inlet and outlet temperature of the anode channel fluid

        self.k_alpha_ch = None
        # convection conductance between the channel and the fluid
        self.cond_rate = np.full((2, self.n_cells, n_nodes), 0.)
        # molar condensation rate
        self.omega = np.full((self.n_cells, n_nodes), 0.)
        # electrical resistance of the membrane
        self.m_reac_flow_delta = np.full((self.n_cells, n_nodes), 0.)
        # mass flow of the consumed oxygen in the cathode channels
        self.g_fluid = []
        # heat capacity flow of the channel fluids
        self.cp_h2 = np.full((self.n_cells, n_nodes), 0.)
        k_p, k_g, k_m = [], [], []
        k_pp, k_gp, k_gm = [], [], []
        for cell in self.cells:
            k_m = np.hstack((k_m, cell.k_mem_z))
            k_g = np.hstack((k_g, cell.k_gde_z))
            k_p = np.hstack((k_p, cell.k_bpp_z))
            k_gm = np.hstack((k_gm, cell.k_gm))
            k_gp = np.hstack((k_gp, cell.k_gp))
            k_pp = np.hstack((k_pp, cell.k_bpp_x))
        k_layer = np.array([[k_m, k_g, k_p], [k_gm, k_gp, k_pp]])
        # heat conductivity of the cell layer

        """"Calculation of the environment heat conductivity"""
        # free convection geometry model
        cell_width = self.cells[0].width
        cell_length = self.cells[0].length
        fac = (cell_width + cell_length)\
            / (self.cells[0].cathode.channel.length
               * self.cells[0].width_channels)
        k_alpha_amb = np.full((2, 3, self.n_cells), 0.)
        # convection conductance to the environment
        alpha_amb = temperature_dict['alpha_amb']
        for i, cell in enumerate(self.cells):
            avg_dx = np.average(cell.cathode.channel.dx)
            k_alpha_amb[0, 1, i] =\
                .5 * alpha_amb * avg_dx\
                * (cell.cathode.th_bpp + cell.cathode.th_gde) / fac
            k_alpha_amb[0, 0, i] =\
                .5 * (alpha_amb * avg_dx
                      * (cell.cathode.th_gde + cell.th_mem)) / fac
            k_alpha_amb[0, 2, i] = \
                alpha_amb * avg_dx * cell.cathode.th_bpp / fac
        # Initialize the thermal coupling
        temperature_dict['k_layer'] = k_layer
        temperature_dict['k_alpha_amb'] = k_alpha_amb
        self.temp_sys = therm_cpl.TemperatureSystem(temperature_dict,
                                                    self.cells)

    def update(self):
        """
        This function coordinates the program sequence
        """
        for i, cell in enumerate(self.cells):
            #self.cells[j].set_current_density(self.i_cd[j, :])
            cell.i_cd[:] = self.i_cd[i, :]
            cell.update()
            if cell.break_program:
                self.break_program = True
                break
        if not self.break_program:
            self.stack_dynamic_properties()
            if self.calc_temp:
                self.update_temperature_coupling()
            if self.n_cells > 1:
                if self.calc_flow_dis:
                    self.update_flows()
            self.i_cd_old[:] = self.i_cd
            if self.calc_cd:
                self.update_electrical_coupling()
        #print('Current density', self.i_cd)

    def update_flows(self):
        """
        This function updates the flow distribution of gas over the stack cells
        """
        n_ch = self.cells[0].cathode.n_chl
        self.manifold[0].update_values(
            self.q_sum_cat * n_ch, self.temp_fluid_cat,
            self.cp_cat, self.visc_cat, self.p_cat, self.r_cat,
            self.m_sum_f_cat * n_ch,
            self.m_sum_g_cat * n_ch)
        self.manifold[1].update_values(
            self.q_sum_ano[::-1] * n_ch, self.temp_fluid_ano[::-1],
            self.cp_ano[::-1], self.visc_ano[::-1], self.p_ano[::-1],
            self.r_ano[::-1], self.m_sum_f_ano[::-1] * n_ch,
            self.m_sum_g_ano[::-1] * n_ch)
        self.manifold[0].update()
        self.manifold[1].update()
        self.set_stoichiometry(self.manifold[0].cell_stoi,
                               self.manifold[1].cell_stoi)
        self.set_channel_outlet_pressure(self.manifold[0].head_p[-1],
                                         self.manifold[1].head_p[-1])

    def update_electrical_coupling(self):
        """
        This function updates current distribution over the stack cells
        """
        self.el_cpl_stack.update_values(self.v_loss, self.stack_cell_r)
        self.el_cpl_stack.update()
        self.i_cd = self.el_cpl_stack.i_cd

    def update_temperature_coupling(self):
        """
        This function updates the layer and fluid temperatures of the stack
        """
        current = self.i_cd * self.cells[0].active_area_dx
        n_ch = self.cells[0].cathode.n_chl
        self.temp_sys.update_values(self.k_alpha_ch,
                                    self.cond_rate,
                                    self.omega,
                                    np.array([self.v_loss_cat,
                                              self.v_loss_ano]),
                                    self.g_fluid, current)
        self.temp_sys.update()
        self.set_temperature()

    def stack_dynamic_properties(self):
        """
        This function sums up the dynamic values inside the cells
        necessary to calculate the flow distribution,
        the electrical coupling or the temperature coupling
        """
        v_alarm = []
        k_alpha_cat, k_alpha_ano = [], []
        g_fluid_cat, g_fluid_ano = [], []
        v_cell, v_loss, resistance = [], [], []
        v_loss_cat, v_loss_ano = [], []
        q_sum_cat_in, q_sum_cat_out = [], []
        q_sum_ano_in, q_sum_ano_out = [], []
        cp_cat_in, cp_cat_out = [], []
        cp_ano_in, cp_ano_out = [], []
        visc_cat_in, visc_cat_out = [], []
        visc_ano_in, visc_ano_out = [], []
        p_cat_in, p_cat_out = [], []
        p_ano_in, p_ano_out = [], []
        r_cat_in, r_cat_out = [], []
        r_ano_in, r_ano_out = [], []
        temp_fluid_cat_in, temp_fluid_cat_out = [], []
        temp_fluid_ano_in, temp_fluid_ano_out = [], []
        cond_rate_cat, cond_rate_ano = [], []
        m_sum_f_cat_in, m_sum_f_cat_out = [], []
        m_sum_f_ano_in, m_sum_f_ano_out = [], []
        m_sum_g_cat_in, m_sum_g_cat_out = [], []
        m_sum_g_ano_in, m_sum_g_ano_out = [], []
        omega = []
        for cell in self.cells:
            v_alarm.append(cell.v_alarm)
            cond_rate_cat.append(cell.cathode.cond_rate)
            cond_rate_ano.append(cell.anode.cond_rate)
            k_alpha_cat.append(cell.cathode.k_ht_coeff_ca)
            k_alpha_ano.append(cell.anode.k_ht_coeff_ca)
            g_fluid_cat.append(cell.cathode.g_fluid)
            g_fluid_ano.append(cell.anode.g_fluid)
            v_cell.append(cell.v)
            v_loss = np.hstack((v_loss, cell.v_loss))
            v_loss_cat.append(cell.cathode.v_loss)
            v_loss_ano.append(cell.anode.v_loss)
            omega.append(cell.omega)
            resistance = np.hstack((resistance, cell.resistance))
            q_sum_cat_in = \
                np.hstack((q_sum_cat_in, cell.cathode.vol_flow_gas[0]))
            q_sum_cat_out = \
                np.hstack((q_sum_cat_out, cell.cathode.vol_flow_gas[-1]))
            q_sum_ano_in = \
                np.hstack((q_sum_ano_in, cell.anode.vol_flow_gas[0]))
            q_sum_ano_out = \
                np.hstack((q_sum_ano_out, cell.anode.vol_flow_gas[-1]))
            m_sum_f_cat_in = np.hstack((m_sum_f_cat_in,
                                        cell.cathode.mass_flow_total[0]))
            m_sum_f_cat_out = np.hstack((m_sum_f_cat_out,
                                         cell.cathode.mass_flow_total[-1]))
            m_sum_f_ano_in = np.hstack((m_sum_f_ano_in,
                                        cell.anode.mass_flow_total[0]))
            m_sum_f_ano_out = np.hstack((m_sum_f_ano_out,
                                         cell.anode.mass_flow_total[-1]))
            m_sum_g_cat_in = np.hstack((m_sum_g_cat_in,
                                        cell.cathode.mass_flow_gas_total[0]))
            m_sum_g_cat_out = np.hstack((m_sum_g_cat_out,
                                         cell.cathode.mass_flow_gas_total[-1]))
            m_sum_g_ano_in = np.hstack((m_sum_g_ano_in,
                                        cell.anode.mass_flow_gas_total[0]))
            m_sum_g_ano_out = np.hstack((m_sum_g_ano_out,
                                         cell.anode.mass_flow_gas_total[-1]))
            cp_cat_in = np.hstack((cp_cat_in, cell.cathode.cp_fluid[0]))
            cp_cat_out = np.hstack((cp_cat_out, cell.cathode.cp_fluid[-1]))
            cp_ano_in = np.hstack((cp_ano_in, cell.anode.cp_fluid[0]))
            cp_ano_out = np.hstack((cp_ano_out, cell.anode.cp_fluid[-1]))
            p_cat_in = np.hstack((p_cat_in, cell.cathode.p[0]))
            p_cat_out = np.hstack((p_cat_out, cell.cathode.p[-1]))
            p_ano_in = np.hstack((p_ano_in, cell.anode.p[0]))
            p_ano_out = np.hstack((p_ano_out, cell.anode.p[-1]))
            r_cat_in = np.hstack((r_cat_in, cell.cathode.r_gas[0]))
            r_cat_out = np.hstack((r_cat_out, cell.cathode.r_gas[-1]))
            r_ano_in = np.hstack((r_ano_in, cell.anode.r_gas[0]))
            r_ano_out = np.hstack((r_ano_out, cell.anode.r_gas[-1]))
            visc_cat_in = np.hstack((visc_cat_in, cell.cathode.visc_gas[0]))
            visc_cat_out = np.hstack((visc_cat_out, cell.cathode.visc_gas[-1]))
            visc_ano_in = np.hstack((visc_ano_in, cell.anode.visc_gas[0]))
            visc_ano_out = np.hstack((visc_ano_out, cell.anode.visc_gas[-1]))
            temp_fluid_cat_in = np.hstack((temp_fluid_cat_in,
                                           cell.cathode.temp_fluid[0]))
            temp_fluid_cat_out = np.hstack((temp_fluid_cat_out,
                                            cell.cathode.temp_fluid[-1]))
            temp_fluid_ano_in = np.hstack((temp_fluid_ano_in,
                                           cell.anode.temp_fluid[0]))
            temp_fluid_ano_out = np.hstack((temp_fluid_ano_out,
                                            cell.anode.temp_fluid[-1]))
        self.k_alpha_ch = np.array([k_alpha_cat, k_alpha_ano])
        self.omega = np.array(omega)
        self.cond_rate = np.array([cond_rate_cat, cond_rate_ano])
        self.g_fluid = np.array([g_fluid_cat, g_fluid_ano])
        self.v_cell, self.v_loss, self.stack_cell_r = v_cell, v_loss, resistance
        self.v_loss_cat, self.v_loss_ano = v_loss_cat, v_loss_ano
        self.q_sum_cat = np.array([q_sum_cat_in, q_sum_cat_out])
        self.q_sum_ano = np.array([q_sum_ano_in, q_sum_ano_out])
        self.m_sum_f_cat = np.array([m_sum_f_cat_in, m_sum_f_cat_out])
        self.m_sum_f_ano = np.array([m_sum_f_ano_in, m_sum_f_ano_out])
        self.m_sum_g_cat = np.array([m_sum_g_cat_in, m_sum_g_cat_out])
        self.m_sum_g_ano = np.array([m_sum_g_ano_in, m_sum_g_ano_out])
        self.cp_cat = np.array([cp_cat_in, cp_cat_out])
        self.cp_ano = np.array([cp_ano_in, cp_ano_out])
        self.visc_cat = np.array([visc_cat_in, visc_cat_out])
        self.visc_ano = np.array([visc_ano_in, visc_ano_out])
        self.p_cat = np.array([p_cat_in, p_cat_out])
        self.p_ano = np.array([p_ano_in, p_ano_out])
        self.r_cat = np.array([r_cat_in, r_cat_out])
        self.r_ano = np.array([r_ano_in, r_ano_out])
        self.temp_fluid_cat = np.array([temp_fluid_cat_in, temp_fluid_cat_out])
        self.temp_fluid_ano = np.array([temp_fluid_ano_in, temp_fluid_ano_out])
        self.v_alarm = np.array(v_alarm)

    def set_stoichiometry(self, stoi_cat, stoi_ano):
        """
        This function sets up the inlet stoichiometry
        of the cathode and anode channels.
        """
        for i, cell in enumerate(self.cells):
            cell.cathode.stoi = stoi_cat[i]
            cell.anode.stoi = stoi_ano[i]

    def set_channel_outlet_pressure(self, p_cat, p_ano):
        """
        This function sets up the inlet pressure
        of the cathode and the anode channels.
        """
        for i, cell in enumerate(self.cells):
            cell.cathode.channel.p_out = p_cat[i]
            cell.anode.channel.p_out = p_ano[i]

    def set_temperature(self):
        """
        This function sets up the layer and fluid temperatures in the cells.
        """
        for i, cell in enumerate(self.cells):
            cell.temp_layer[:] = self.temp_sys.temp_layer[i][:, :]
            cell.cathode.temp_fluid[:] = self.temp_sys.temp_fluid[0, i]
            cell.anode.temp_fluid[:] = self.temp_sys.temp_fluid[1, i]
            cell.temp_cool[:] = self.temp_sys.temp_cool_ele[i]
