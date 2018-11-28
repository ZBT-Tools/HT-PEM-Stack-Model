from data.global_parameter import dict_case


class Channel:

    def __init__(self, dict_ch):
        self.length = dict_ch['length']
        self.d_x = self.length / float(dict_case['elements'])
        self.p_in = dict_ch['p_in']
        self.t_in = dict_ch['t_in']
        self.humidity_in = dict_ch['hum_in']
        self.flow_dir = dict_ch['flow_dir']
        self.width = dict_ch['width']
        self.height = dict_ch['height']
        self.bend_num = dict_ch['num_bends']
        self.bend_fri_fac = dict_ch['bend_fri_fac']
        self.rack_width = dict_ch['rack_width']
        self.cross_area = self.width * self.height
        self.extent = 2 * (self.width + self.height)
        self.plane = self.width * self.length
        self.plane_dx = self.width * self.d_x

    def set_p_in(self, p_in):
        self.p_in = p_in
