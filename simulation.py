import input_parameter as i_p
import channel
import halfcell
import cell
import stack
import numpy as np
from scipy.optimize import fsolve, root
import matrix_database as m_d
import global_functions as g_func
import time


class Simulation:

    def __init__(self, stack, k_it, max_it):
        self.stack = stack
        self.k_it = k_it
        self.max_it = max_it
        self.b = m_d.b(i_p.nodes, i_p.cell_numb, self.stack.cell.cathode.channel.d_x)
        self.c = m_d.c(i_p.nodes, i_p.cell_numb, self.stack.cell.cathode.channel.d_x)

    def update(self):
        start_time = time.time()
        for w in range(i_p.n_loop):
            self.stack.update()
            self.calc_initial_current_density()
            self.stack.i = self.stack.i * i_p.tar_cd * i_p.cell_numb\
                           * i_p.nodes /np.sum(self.stack.i)
        #self.calc_i()

        print("--- %s seconds ---" % (time.time() - start_time))

    def calc_initial_current_density_fsolve(self, x, i, w, v):
        a = self.stack.cell.e0 - x * self.stack.cell_list[i].omega[w]
        b = self.stack.cell.cathode.r * self.stack.cell_list[i].t1[w] \
            / self.stack.cell.cathode.f_const
        c = np.log(x * self.stack.cell.c_ref /
                   (self.stack.cell.i_ref
                    * (self.stack.cell_list[i].cathode.c1[w] - self.stack.cell.delta * x)))
        return a - b * c - v

    def calc_initial_current_density(self):
        for i in range(self.stack.cell_numb):
            avv = sum(self.stack.cell_list[i].v) \
                  / self.stack.cell.cathode.channel.nodes
            for q in range(self.stack.cell.cathode.channel.nodes):
                self.stack.i[i, q] = fsolve(self.calc_initial_current_density_fsolve,
                                               self.stack.cell.cathode.tar_cd, args=(i, q, avv))

    #def calc_i(self):




channel_anode = channel.Channel(i_p.channel_length, i_p.nodes, i_p.k_ano,
                                i_p.p_ano_in, i_p.t_ano_in, i_p.phi_ano,
                                i_p.flow_dir_ano)
channel_cathode = channel.Channel(i_p.channel_length, i_p.nodes, i_p.k_cat,
                                  i_p.p_cat_in, i_p.t_cat_in, i_p.phi_cat,
                                  i_p.flow_dir_cat)
anode = halfcell.Halfcell(channel_anode, i_p.stoi_ano, i_p.spec_numb_ano,
                          i_p.val_ano, i_p.tar_cd)
cathode = halfcell.Halfcell(channel_cathode, i_p.stoi_cat, i_p.spec_numb_cat,
                            i_p.val_cat, i_p.tar_cd)
cell1 = cell.Cell(anode, cathode, i_p.gamma, i_p.alpha, i_p.mem_thick, i_p.mu_p,
                  i_p.mu_g, i_p.mu_m, i_p.e_o, i_p.c_ref, i_p.i_ref, i_p.delta,
                  i_p.pem_type, i_p.coolant_channel, i_p.t_cool_in,
                  i_p.channel_width, i_p.gde_thick, i_p.plate_thick, i_p.t_u)
stack1 = stack.Stack(cell1, i_p.cell_numb, i_p.end_plate_t, i_p.alpha_conv, i_p.resistance, i_p.endplate_adiabat)
