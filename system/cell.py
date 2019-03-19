import numpy as np
import system.global_functions as g_func
import data.global_parameters as g_par
import system.half_cell as h_c
import data.half_cell_dict as hc_dict


class Cell:

    def __init__(self, dict_cell):
        # Handover
        # anode - object of the class HalfCell
        self.anode = h_c.HalfCell(hc_dict.dict_anode)
        # cathode - object of the class HalfCell
        self.cathode = h_c.HalfCell(hc_dict.dict_cathode)
        # thickness membrane [m]
        self.th_mem = dict_cell['th_mem']
        # heat conductivity of the bipolar plate [W/m/K]
        self.lambda_bpp = [dict_cell['lambda_z_bpp'], dict_cell['lambda_x_bpp']]
        # heat conductivity of the gas diffusion layer [W/m/K]
        self.lambda_gde = [dict_cell['lambda_z_gde'], dict_cell['lambda_x_gde']]
        # heat conductivity of the membrane [W/m/K]
        self.lambda_mem = [dict_cell['lambda_z_mem'], dict_cell['lambda_x_mem']]
        # inlet temperature of the coolant [K]
        self.temp_cool_in = dict_cell['temp_cool_in']
        # basic electrical resistance of the membrane(Kvesic, 2013) [Ohm/m²]
        self.mem_base_r = dict_cell['mem_base_r']
        # thermal related electrical resistance
        # of the membrane (Kvesic, 2013) [Ohm/m²/K]
        self.mem_acl_r = dict_cell['mem_acl_r']

        """heat conductivity along and through the cell layers"""
        # width of the the meander including n+1 ribs and n channels
        self.width_channels = self.cathode.channel.width\
            * self.cathode.channel_numb\
            + self.cathode.channel.rib_width\
            * (self.cathode.channel_numb + 1)
        # active area of an element
        self.active_area_dx = self.width_channels * self.cathode.channel.dx
        # heat conductivity through the bipolar plate [W/m/K]
        self.k_bpp_z = self.lambda_bpp[0] * self.active_area_dx \
            / self.cathode.th_bpp
        # heat conductivity through the gas diffusion electrode [W/m/K]
        self.k_gde_z = self.lambda_gde[0]\
            * self.active_area_dx \
            / self.cathode.th_gde
        # heat conductivity through the membrane [W/m/K]
        self.k_mem_z = self.lambda_mem[0]\
            * self.active_area_dx \
            / self.th_mem
        # heat conductivity along the bipolar plate [W/m/K]
        self.k_bpp_x = self.width_channels * self.lambda_bpp[1] \
            * self.cathode.th_bpp / self.cathode.channel.dx
        # heat conductivity along the gas diffusion
        # electrode and the bipolar plate [W/m/K]
        self.k_gp = (self.width_channels
                     * (self.lambda_bpp[1] * self.cathode.th_bpp
                        + self.lambda_gde[1] * self.cathode.th_gde))\
            / (2. * self.cathode.channel.dx)
        # heat conductivity along the bipolar plate and
        # the gas diffusion electrode [W/m/K]
        self.k_gm = (self.width_channels
                     * (self.lambda_mem[1] * self.th_mem
                        + self.lambda_gde[1] * self.cathode.th_gde))\
            / (2. * self.cathode.channel.dx)

        """boolean alarms"""
        # True if :voltage loss > cell voltage
        self.v_alarm = False
        # True if the program aborts because of some critical impact
        self.break_program = False

        """general parameter"""
        # height of the cell
        self.height = self.th_mem\
            + 2. * self.cathode.th_bpp\
            + 2. * self.cathode.th_gde
        # number of nodes along the channel
        nodes = g_par.dict_case['nodes']
        # water cross flux through the membrane [mol / s]
        self.w_cross_flow = np.zeros(nodes-1)
        # area specific membrane resistance [ohm m²]
        self.omega_ca = np.zeros(nodes)
        # cell voltage loss [V]
        self.v_loss = np.full(nodes - 1, 0.)
        # layer temperature [K]
        self.temp = np.full((5, nodes - 1), dict_cell['temp_init'])
        # membrane temperature [K]
        self.temp_mem = np.zeros(nodes)
        # current density [A/m²]
        self.i_ca = np.full((nodes - 1), 0.)
        # resistance of the membrane [Ohm]
        self.omega = np.full(nodes-1, 0.)
        # voltage loss at the membrane [V]
        self.mem_loss = np.full((nodes - 1), 0.)
        # cell voltage [V]
        self.v = np.full((nodes - 1), 0.)
        # cell resistance [Ohm]
        self.resistance = np.full((nodes - 1), 0.)

    def update(self):
        """
        This function coordinates the program sequence
        """

        self.temp_mem = .5 * (self.temp[2] + self.temp[3])
        if g_par.dict_case['pem_type'] is False:
            self.cathode.set_pem_type(False)
            self.anode.set_pem_type(False)
            self.cathode.set_water_cross_flux(self.w_cross_flow)
            self.anode.set_water_cross_flux(self.w_cross_flow)
        self.cathode.set_current_density(self.i_ca)
        self.anode.set_current_density(self.i_ca)
        self.cathode.set_layer_temperature([self.temp[2],
                                            self.temp[3],
                                            self.temp[4]])
        self.anode.set_layer_temperature([self.temp[0], self.temp[1]])
        self.cathode.update()
        self.anode.update()
        if self.anode.break_program is True\
                or self.cathode.break_program is True:
            self.break_program = True
        else:
            if g_par.dict_case['pem_type'] is False:
                self.calc_cross_water_flux()
                self.calc_mem_resistivity_springer()
            else:
                self.calc_mem_resistivity_kvesic()
            self.calc_membrane_loss()
            self.calc_voltage()
            self.calc_resistance()

    def set_current_density(self, i_ca):
        """
        This function sets the current density.
        The current density can be obtained
        from the electrical coupling.

            Manipulate:
            - self.i_ca, scalar
        """

        self.i_ca = i_ca

    def calc_cross_water_flux(self):
        """
        Calculates the water cross flux through the membrane
        according to (Springer, 1991).

            Access to:
            -g_par.dict_case['vap_m_t_coef']
            -self.cathode.free_w
            -self.anode.free_w
            -g_par.dict_case['mol_con_m']
            -g_par.dict_uni['F']
            -self.i_ca
            -self.temp_mem
            -self.th_mem

            Manipulate:
            -self.w_cross_flow
        """

        vap_coeff = g_par.dict_case['vap_m_temp_coeff']
        humidity = [self.cathode.humidity, self.anode.humidity]
        humidity_ele = np.array([g_func.calc_elements_1_d(humidity[0]),
                                 g_func.calc_elements_1_d(humidity[1])])
        a = 0.043 + 17.81 * humidity_ele
        b = -39.85 * humidity_ele ** 2. + 36. * humidity_ele ** 3.
        free_w_content = a + b
        zeta_plus = free_w_content[0] + free_w_content[1] \
            + self.i_ca / (2. * vap_coeff
                           * g_par.dict_case['mol_con_m']
                           * g_par.dict_uni['F'])
        zeta_negative =\
            (free_w_content[0]
             - free_w_content[1]
             + 5. * self.i_ca / (2. * vap_coeff * g_par.dict_case['mol_con_m']
                                 * g_par.dict_uni['F'])) \
            / (1. + g_func.dw(self.temp_mem) * zeta_plus
               / (self.th_mem * vap_coeff))
        m_c = 0.5 * (zeta_plus + zeta_negative)
        m_a = 0.5 * (zeta_plus - zeta_negative)
        self.w_cross_flow =\
            self.i_ca / g_par.dict_uni['F'] + g_par.dict_case['mol_con_m']\
            * g_func.dw(self.temp_mem) * (m_a ** 2. - m_c ** 2.)\
            / (2. * self.th_mem)

    def calc_mem_resistivity_kvesic(self):
        """
        Calculates the membrane resistivity and resistance
        according to (Kvesic, 2013).

            Access to:
            -self.mem_base_r
            -self.mem_acl_r
            -self.temp_mem
            -self.omega_ca
            -self.cathode.cathode.active_area_dx_ch

            Manipulate:
            -self.omega_ca
            -self.omega
        """

        self.omega_ca = (self.mem_base_r
                         - self.mem_acl_r * self.temp_mem) * 1.e-4
        self.omega = self.omega_ca / self.active_area_dx

    def calc_mem_resistivity_springer(self):
        """
        Calculates the membrane resistivity
        for NT-PEMFC according to (Springer, 1991).

            Access to:
            -self.cathode.humidity
            -self.anode.humidity
            -self.temp_mem
            -self.th_mem
            -self.cathode.cathode.active_area_dx_ch

            Manipulate:
            -self.omega_ca
            -self.omega
        """
        humidity = (self.cathode.humidity + self.anode.humidity) * 0.5
        humidity_ele = g_func.calc_elements_1_d(humidity)
        a = 0.043 + 17.81 * humidity_ele
        b = -39.85 * humidity_ele ** 2. + 36. * humidity_ele ** 3.
        free_water_content = a + b
        mem_el_con = 0.005139 * free_water_content - 0.00326
        mem_el_con_temp =\
            np.exp(1268 * (0.0033 - 1. / self.temp_mem)) * mem_el_con
        self.omega_ca = self.th_mem / mem_el_con_temp * 1.e-4

    def calc_membrane_loss(self):
        """
        Calculates the voltage loss at the membrane.

            Access to:
            -self.omega_ca
            -self.i_ca

            Manipulate:
            -self.mem_los
        """

        self.mem_loss = self.omega_ca * self.i_ca
        if g_par.dict_case['calc_mem_loss'] is False:
            self.mem_loss = 0.

    def calc_voltage(self):
        """
            Calculates the cell voltage loss.
            Calculates the cell voltage
            If the cell voltage loss is larger than the
            open circuit cell voltage, the cell voltage is set to zero.


            Access to:
            -self.mem_loss
            -self.cathode.v_loss
            -self.anode.v_loss
            -g_par.dict_case['e_0']

            Manipulate:
            -self.v_loss
            -self.v
            -self.v_alarm
        """

        self.v_loss = self.mem_loss + self.cathode.v_loss + self.anode.v_loss
        if any(self.v_loss) >= g_par.dict_case['e_0']:
            self.v_alarm = True
        self.v_loss = np.minimum(self.v_loss, g_par.dict_case['e_0'])
        self.v = g_par.dict_case['e_0'] - self.v_loss

    def calc_resistance(self):
        """
            Calculates the electrical resistance of the element in z-direction

            Access to:
            -self.v_loss
            -self.i_ca
            -self.cathode.th_bpp
            -g_par.dict_case['plate_resistivity']

            Manipulate:
            -self.resistance
        """

        self.resistance = self.v_loss / self.i_ca + 2. \
            * g_par.dict_case['bpp_resistivity'] * self.cathode.th_bpp
