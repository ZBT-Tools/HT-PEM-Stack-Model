from global_parameter import dict_case


class Channel:

    def __init__(self, dict):
        self.length = dict['length']
        self.d_x = self.length / dict_case['elements']
        self.p_in = dict['p_in']
        self.t_in = dict['t_in']
        self.phi = dict['hum_in']
        self.flow_dir = dict['flow_dir']
        self.width = dict['width']
        self.heigth = dict['heigth']
        self.bend_numb = dict['numb_bends']
        self.bend_fri_fac = dict['bend_fri_fac']
        self.cross_area = self.width * self.heigth
        self.extent = 2 * (self.width + self.heigth)
        self.plane = self.width * self.length
        self.plane_dx = self.width * self.d_x

    def set_p_in(self, p_in):
        self.p_in = p_in