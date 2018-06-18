class Channel:

    def __init__(self, length, nodes, k, p_in, t_in, phi, flow_dir, width, heigth):
        self.length = length
        self.nodes = nodes
        self.elements = nodes - 1
        self.d_x = length / self.elements
        self.k = k
        self.p_in = p_in
        self.t_in = t_in
        self.phi = phi
        self.flow_dir = flow_dir
        self.width = width
        self.heigth = heigth
        self.cross_area = self.width * self.heigth
        self.extent = 2 * (self.width + self.heigth)

    def set_p_in(self, p_in):
        self.p_in = p_in
