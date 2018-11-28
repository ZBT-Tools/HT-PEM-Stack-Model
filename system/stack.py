import numpy as np
import copy as copy
import data.global_parameter as g_par
import system.cell as cl
import data.cell_dict as c_dict
import system.manifold as m_fold
import data.manifold_dict as m_fold_dict
import system.electrical_coupling as el_cpl
import data.electrical_coupling_dict as el_cpl_dict
import system.temperature_system as therm_cpl
import data.temperature_system_dict as therm_dict


class Stack:

    def __init__(self, dict_stack):
        # Handover
        self.cell_num = dict_stack['cell_numb']
        self.heat_pow = dict_stack['heat_power']
        self.stoi_cat = dict_stack['stoi_cat']
        self.stoi_ano = dict_stack['stoi_ano']
        self.alpha_env = g_par.dict_case['conv_coef']
        self.ch_bc_cool = dict_stack['cool_ch_bc']
        self.height_cool = dict_stack['h_cool']
        self.m_cool = dict_stack['m_flow_cool']
        self.cp_cool = dict_stack['cp_cool']
        self.alpha_cool = dict_stack['alpha_cool']
        nodes = g_par.dict_case['nodes']
        # Initialize the stack cells
        self.cells = []
        for w in range(self.cell_num):
            x = cl.Cell(c_dict.cell)
            self.cells.append(x)
        self.set_stoi(np.full(self.cell_num, self.stoi_cat),
                      np.full(self.cell_num, self.stoi_ano))
        # Initialize the manifolds
        self.manifold = [m_fold.Manifold(m_fold_dict.m_fold_const_cathode),
                         m_fold.Manifold(m_fold_dict.m_fold_const_anode)]
        self.manifold[0].head_stoi = self.stoi_cat
        self.manifold[1].head_stoi = self.stoi_ano
        # Initialize the electrical coupling
        self.el_cpl_stack = el_cpl\
            .ElectricalCoupling(el_cpl_dict.electrical_coupling_dict_const)
        # Prove the coolant channel settings
        if self.ch_bc_cool is False:
            self.t_cool = np.full((self.cell_num, nodes),
                                  c_dict.cell['t_cool_in'])
            self.a_cool = np.full(self.cell_num, self.alpha_cool)
            self.a_cool[0] = 1.e-20
        else:
            self.a_cool = np.full(self.cell_num + 1, self.alpha_cool)
            self.t_cool = np.full((self.cell_num + 1, nodes),
                                  c_dict.cell['t_cool_in'])
        # Scalar variables
        self.cathode_mfd_criteria = 0.
        self.anode_mfd_criteria = 0.
        self.g_cool = self.m_cool * self.cp_cool
        self.d_cool = 4. * (self.cells[0].cathode.channel.width
                            * self.height_cool) / (2. * (self.height_cool
                                                         + self.cells[0].
                                                         cathode.channel.width))
        self.break_program = False
        # Arrays
        self.v_alarm = np.full(self.cell_num, False)
        self.i_ca = np.full((self.cell_num, nodes - 1),
                            g_par.dict_case['tar_cd'])
        self.i_old = copy.deepcopy(self.i_ca)
        self.i = np.full((self.cell_num, nodes-1), 20.)
        self.zero = np.full(self.cell_num, 0.)
        self.inf = np.full(self.cell_num, 1.e50)
        self.k_alpha_ch = np.full((2, self.cell_num), 0.)
        self.v = copy.deepcopy(self.zero)
        self.v_los = copy.deepcopy(self.zero)
        self.v_los_cat = copy.deepcopy(self.zero)
        self.v_los_ano = copy.deepcopy(self.zero)
        self.stack_cell_r = copy.deepcopy(self.inf)
        self.q_sum_cat = copy.deepcopy(self.zero)
        self.q_sum_ano = copy.deepcopy(self.zero)
        self.m_sum_cat = copy.deepcopy(self.zero)
        self.m_sum_ano = copy.deepcopy(self.zero)
        self.cp_cat = copy.deepcopy(self.zero)
        self.cp_ano = copy.deepcopy(self.zero)
        self.visc_cat = copy.deepcopy(self.zero)
        self.visc_ano = copy.deepcopy(self.zero)
        self.p_cat = copy.deepcopy(self.zero)
        self.p_ano = copy.deepcopy(self.zero)
        self.r_cat = copy.deepcopy(self.zero)
        self.r_ano = copy.deepcopy(self.zero)
        self.t_gas_cat = copy.deepcopy(self.zero)
        self.t_gas_ano = copy.deepcopy(self.zero)
        self.k_alpha_env = np.full((2, 3, self.cell_num), -1.e50)
        self.gamma = np.full((2, self.cell_num, nodes), 0.)
        self.omega = np.full((self.cell_num, nodes), 0.)
        self.m_reac_flow_delta = np.full((self.cell_num, nodes), 0.)
        self.g_gas = np.full((self.cell_num, nodes), 0.)
        self.cp_h2 = np.full((self.cell_num, nodes), 0.)
        # Initialize the constant properties
        self.stack_constant_properties()
        fac = 1.
        # free convection geometry model
        n_ch_x = int(self.cells[0].cathode.channel.bend_num * .5 + 1)
        n_ch_y = n_ch_x - 1
        l_ch_y = 2. * self.cells[0].cathode.channel.width\
            + self.cells[0].cathode.channel.rack_width
        l_ch_x = (self.cells[0].cathode.channel.length
                  - n_ch_y * l_ch_y
                  - 2. * self.cells[0].cathode.channel.rack_width)\
            / n_ch_x
        h_flow_y = l_ch_y * n_ch_y\
            + 2. * self.cells[0].cathode.channel.rack_width
        h_flow_x = l_ch_x + 2. * self.cells[0].cathode.channel.rack_width
        extent_flow = 2. * (h_flow_x + h_flow_y)
        fac = extent_flow / self.cells[0].cathode.channel.extent
        for q, item in enumerate(self.cells):
            self.k_alpha_env[0, 1, q] =\
                .5 * self.alpha_env * item.cathode.channel.d_x\
                * (item.cathode.th_plate + item.cathode.th_gde) / fac
            self.k_alpha_env[0, 0, q] =\
                .5 * (self.alpha_env * item.cathode.channel.d_x
                      * (item.cathode.th_plate + item.th_mem)) / fac
            self.k_alpha_env[0, 2, q] = \
                self.alpha_env * item.cathode.channel.d_x\
                * item.cathode.th_plate / fac
        # Initialize the thermal coupling
        therm_dict.t_sys_const_dict['k_layer'] = self.k_layer
        therm_dict.t_sys_const_dict['k_alpha_env'] = self.k_alpha_env
        self.temp_cpl_stack = therm_cpl.\
            TemperatureSystem(therm_dict.t_sys_const_dict)

    def update(self):
        for j in range(self.cell_num):
            self.cells[j].set_i(self.i_ca[j, :])
            self.cells[j].update()
            if self.cells[j].break_program is True:
                self.break_program = True
                break
        if self.break_program is False:
            self.stack_dynamic_properties()
            self.update_temperature_coupling()
            if self.cell_num > 1:
                self.update_flows()
            self.i_old = copy.deepcopy(self.i_ca)
            self.update_electrical_coupling()

    def update_flows(self):
        self.manifold[0].update_values(
            m_fold_dict.dict_manifold_dyn(self.q_sum_cat, self.t_gas_cat,
                                          self.cp_cat, self.visc_cat,
                                          self.p_cat, self.r_cat,
                                          self.m_sum_cat))
        self.manifold[1].update_values(
            m_fold_dict.dict_manifold_dyn(self.q_sum_ano[::-1],
                                          self.t_gas_ano[::-1],
                                          self.cp_ano[::-1],
                                          self.visc_ano[::-1],
                                          self.p_ano[::-1],
                                          self.r_ano[::-1],
                                          self.m_sum_ano[::-1]))
        self.manifold[0].update()
        self.manifold[1].update()
        self.set_stoi(self.manifold[0].cell_stoi, self.manifold[1].cell_stoi)
        self.set_p(self.manifold[0].head_p[0], self.manifold[1].head_p[0])
        self.cathode_mfd_criteria = self.manifold[0].criteria
        self.anode_mfd_criteria = self.manifold[1].criteria

    def update_electrical_coupling(self):

        self.el_cpl_stack.update_values(
            el_cpl_dict.electrical_coupling_dict_dyn(self.v_los,
                                                     self.stack_cell_r))
        self.el_cpl_stack.update()
        self.i_ca = self.el_cpl_stack.i_ca

    def update_temperature_coupling(self):
        self.i = self.i_ca * self.cells[0].cathode.channel.plane_dx
        self.temp_cpl_stack.update_values(
            therm_dict.t_sys_dyn_dict(self.k_alpha_ch, self.gamma, self.omega,
                                      [self.v_los_cat, self.v_los_ano],
                                      self.m_reac_flow_delta, self.g_gas,
                                      self.cp_h2, self.i))
        self.temp_cpl_stack.update()
        self.set_t()

    def stack_constant_properties(self):
        k_p, k_g, k_m = [], [], []
        k_pp, k_gp, k_gm = [], [], []
        for i, item in enumerate(self.cells):
            k_p = np.hstack((k_p, self.cells[i].k_p))
            k_g = np.hstack((k_g, self.cells[i].k_g))
            k_m = np.hstack((k_m, self.cells[i].k_m))
            k_pp = np.hstack((k_pp, self.cells[i].k_pp))
            k_gp = np.hstack((k_gp, self.cells[i].k_gp))
            k_gm = np.hstack((k_gm, self.cells[i].k_gm))
        self.k_layer = np.array([[k_m, k_g, k_p], [k_gm, k_gp, k_pp]])

    def stack_dynamic_properties(self):
        v_alarm = []
        k_alpha_cat, k_alpha_ano = [], []
        g_gas_cat, g_gas_ano = [], []
        v, v_los, resistance = [], [], []
        v_los_cat, v_los_ano = [],[]
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
        t_gas_cat_in, t_gas_cat_out = [], []
        t_gas_ano_in, t_gas_ano_out = [], []
        gamma_cat, gamma_ano = [], []
        m_sum_cat_in, m_sum_cat_out = [], []
        m_sum_ano_in, m_sum_ano_out = [], []
        omega = []
        m_reac_flow = []
        cp_h2 = []
        for w, item in enumerate(self.cells):
            v_alarm.append(item.v_alarm)
            cp_h2.append(item.anode.cp[0])
            gamma_cat.append(item.cathode.gamma)
            gamma_ano.append(item.anode.gamma)
            k_alpha_cat.append(item.cathode.k_ht_coef_ca)
            k_alpha_ano.append(item.anode.k_ht_coef_ca)
            g_gas_cat.append(item.cathode.g_full)
            g_gas_ano.append(item.anode.g_full)
            m_reac_flow.append(item.anode.m_reac_flow_delta)
            v.append(item.v)
            v_los = np.hstack((v_los, item.v_los))
            v_los_cat.append(item.cathode.v_los)
            v_los_ano.append(item.anode.v_los)
            omega.append(item.omega)
            resistance = np.hstack((resistance, item.resistance))
            q_sum_cat_in = np.hstack((q_sum_cat_in, item.cathode.q_gas[0]))
            q_sum_cat_out = np.hstack((q_sum_cat_out, item.cathode.q_gas[-1]))
            q_sum_ano_in = np.hstack((q_sum_ano_in, item.anode.q_gas[0]))
            q_sum_ano_out = np.hstack((q_sum_ano_out, item.anode.q_gas[-1]))
            m_sum_cat_in = np.hstack((m_sum_cat_in,
                                      item.cathode.m_full_flow[0]))
            m_sum_cat_out = np.hstack((m_sum_cat_out,
                                       item.cathode.m_full_flow[-1]))
            m_sum_ano_in = np.hstack((m_sum_ano_in,
                                      item.anode.m_full_flow[0]))
            m_sum_ano_out = np.hstack((m_sum_ano_out,
                                       item.anode.m_full_flow[-1]))
            cp_cat_in = np.hstack((cp_cat_in, item.cathode.cp_mix[0]))
            cp_cat_out = np.hstack((cp_cat_out, item.cathode.cp_mix[-1]))
            cp_ano_in = np.hstack((cp_ano_in, item.anode.cp_mix[0]))
            cp_ano_out = np.hstack((cp_ano_out, item.anode.cp_mix[-1]))
            p_cat_in = np.hstack((p_cat_in, item.cathode.p[0]))
            p_cat_out = np.hstack((p_cat_out, item.cathode.p[-1]))
            p_ano_in = np.hstack((p_ano_in, item.anode.p[0]))
            p_ano_out = np.hstack((p_ano_out, item.anode.p[-1]))
            r_cat_in = np.hstack((r_cat_in, item.cathode.r_mix[0]))
            r_cat_out = np.hstack((r_cat_out, item.cathode.r_mix[-1]))
            r_ano_in = np.hstack((r_ano_in, item.anode.r_mix[0]))
            r_ano_out = np.hstack((r_ano_out, item.anode.r_mix[-1]))
            visc_cat_in = np.hstack((visc_cat_in, item.cathode.visc_mix[0]))
            visc_cat_out = np.hstack((visc_cat_out, item.cathode.visc_mix[-1]))
            visc_ano_in = np.hstack((visc_ano_in, item.anode.visc_mix[0]))
            visc_ano_out = np.hstack((visc_ano_out, item.anode.visc_mix[-1]))
            t_gas_cat_in = np.hstack((t_gas_cat_in, item.cathode.t_gas[0]))
            t_gas_cat_out = np.hstack((t_gas_cat_out, item.cathode.t_gas[-1]))
            t_gas_ano_in = np.hstack((t_gas_ano_in, item.anode.t_gas[0]))
            t_gas_ano_out = np.hstack((t_gas_ano_out, item.anode.t_gas[-1]))
        self.k_alpha_ch = np.array([k_alpha_cat, k_alpha_ano])
        self.omega = np.array(omega)
        self.gamma = np.array([gamma_cat, gamma_ano])
        self.g_gas = np.array([g_gas_cat, g_gas_ano])
        self.m_reac_flow_delta = np.array(m_reac_flow)
        self.v, self.v_los, self.stack_cell_r = v, v_los, resistance
        self.v_los_cat, self.v_los_ano = v_los_cat, v_los_ano
        self.q_sum_cat = np.array([q_sum_cat_in, q_sum_cat_out])
        self.q_sum_ano = np.array([q_sum_ano_in, q_sum_ano_out])
        self.m_sum_cat = np.array([m_sum_cat_in, m_sum_cat_out])
        self.m_sum_ano = np.array([m_sum_ano_in, m_sum_ano_out])
        self.cp_cat = np.array([cp_cat_in, cp_cat_out])
        self.cp_h2 = np.array(cp_h2)
        self.cp_ano = np.array([cp_ano_in, cp_ano_out])
        self.visc_cat = np.array([visc_cat_in, visc_cat_out])
        self.visc_ano = np.array([visc_ano_in, visc_ano_out])
        self.p_cat = np.array([p_cat_in, p_cat_out])
        self.p_ano = np.array([p_ano_in, p_ano_out])
        self.r_cat = np.array([r_cat_in, r_cat_out])
        self.r_ano = np.array([r_ano_in, r_ano_out])
        self.t_gas_cat = np.array([t_gas_cat_in, t_gas_cat_out])
        self.t_gas_ano = np.array([t_gas_ano_in, t_gas_ano_out])
        self.v_alarm = np.array(v_alarm)

    def set_stoi(self, stoi_cat, stoi_ano):
        for w, item in enumerate(self.cells):
            item.cathode.stoi = stoi_cat[w]
            item.anode.stoi = stoi_ano[w]

    def set_p(self, p_cat, p_ano):
        for w, item in enumerate(self.cells):
            item.cathode.channel.p_in = p_cat[w]
            item.anode.channel.p_in = p_ano[w]

    def set_t(self):
        for w, item in enumerate(self.cells):
            item.t = self.temp_cpl_stack.t_layer[w][0:5, :]
            item.cathode.t_gas = self.temp_cpl_stack.t_gas[0, w]
            item.anode.t_gas = self.temp_cpl_stack.t_gas[1, w]
