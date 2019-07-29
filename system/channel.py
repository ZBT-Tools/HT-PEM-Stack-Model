from data.global_parameters import dict_case


class Channel:

    def __init__(self, dict_ch):
        self.length = dict_ch['channel_length']
        # channel length
        self.dx = self.length / float(dict_case['elements'])
        # element length
        self.p_out = dict_ch['p_in']
        # inlet pressure
        self.temp_in = dict_ch['temp_in']
        # inlet temperature
        self.humidity_in = dict_ch['hum_in']
        # inlet humidity
        self.flow_dir = dict_ch['flow_dir']
        # flow direction
        self.width = dict_ch['channel_width']
        # channel width
        self.height = dict_ch['channel_height']
        # channel height
        self.n_bends = dict_ch['bend_numb']
        # number of channel bends
        self.bend_fri_fac = dict_ch['bend_fri_fac']
        # bend friction factor
        self.rib_width = dict_ch['rib_width']
        # rack width
        self.active_area = self.width * self.length
        # planar area of the channel
        self.active_area_dx = self.width * self.dx
        # planar area of an element of the channel
        self.cross_area = self.width * self.height
        # channel cross area
        self.circum = 2. * (self.width + self.height)
        # channel circumference
        self.d_h = 4. * self.cross_area / self.circum
        # channel hydraulic diameter

    def set_inlet_pressure(self, p_in):
        """
        This function sets the inlet pressure for each channel.
        The inlet pressure can be obtained
        from the pressure distribution of the manifold.

            Manipulate:
            - self.p_in, scalar
        """

        self.p_out = p_in
