from data.global_parameters import dict_case


class Channel:

    def __init__(self, dict_ch):
        """geometry"""
        # channel height [m]
        height = dict_ch['channel_height']
        # channel length [m]
        self.length = dict_ch['channel_length']
        # self.width [m]
        self.width = dict_ch['channel_width']
        # element length [m]
        self.dx = self.length / float(dict_case['elements'])
        # rib self.width [m]
        self.rib_width = dict_ch['rib_width']
        # planar active area of the channel [m²]
        self.active_area = self.width * self.length
        # channel cross area [m²]
        self.cross_area = self.width * height
        # hydraulic diameter of the channel [m]
        self.d_h = 2. * self.cross_area / (self.width + height)
        # number of channel bends
        self.bend_numb = dict_ch['bend_numb']
        """physical properties"""
        # channel outlet pressure [Pa]
        self.p_out = dict_ch['p_out']
        # channel inlet temperature [K]
        self.temp_in = dict_ch['temp_in']
        # channel inlet humidity
        self.humidity_in = dict_ch['hum_in']
        # channel bend factor
        self.bend_fri_fac = dict_ch['bend_fri_fac']
        # flow direction
        self.flow_dir = dict_ch['flow_dir']

    def set_inlet_pressure(self, p_in):
        """
        This function sets the inlet pressure for each channel.
        The inlet pressure can be obtained
        from the pressure distribution of the manifold.

            Manipulate:
            - self.p_in, scalar
        """

        self.p_out = p_in
