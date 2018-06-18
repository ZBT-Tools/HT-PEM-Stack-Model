import input_parameter as i_p
import channel
import halfcell
import cell
import stack
import numpy as np
import copy
from scipy.optimize import fsolve, least_squares, root
import matrix_database as m_d
import global_functions as g_func
import time


class Simulation:

    def __init__(self, stack, k_it, max_it):
        self.stack = stack
        self.k_it = k_it
        self.max_it = max_it
        # self.b = m_d.b(i_p.nodes, i_p.cell_numb, self.stack.cell.cathode.channel.d_x)
        # self.c = m_d.c(i_p.nodes, i_p.cell_numb, self.stack.cell.cathode.channel.d_x)
        self.mdf_criteria_cat_process = []
        self.mdf_criteria_ano_process = []
        self.i_criteria_process = []
        self.b = m_d.b(self.stack.cell.cathode.channel.nodes,
                       self.stack.cell_numb,
                       self.stack.cell.cathode.channel.length
                       / self.stack.cell.cathode.channel.nodes)
        self.c = m_d.c(self.stack.cell.cathode.channel.nodes,
                       self.stack.cell_numb,
                       self.stack.cell.cathode.channel.length
                       / self.stack.cell.cathode.channel.nodes)

    def update(self):
        start_time = time.time()
        self.stack.update()
        self.calc_initial_current_density()
        statement = True
        counter = 0
        while statement is True:
            self.save_old_value()
            self.calc_i()
            self.correct_i()
            self.stack.update()
            self.calc_convergence_criteria()
            counter = counter + 1
            if ((self.i_criteria < self.k_it and counter > 5)
                and (self.mdf_criteria_ano < self.k_it
                     and self.mfd_criteria_cat < self.k_it))\
                    or counter > self.max_it:
                statement = False
        print("--- %s seconds ---" % (time.time() - start_time))

    def calc_initial_current_density_fsolve(self, x, i, w, v):
        # print(x, v)
        a = self.stack.cell.e0 - x * self.stack.cell_list[i].omega[w]
        # print(self.stack.cell_list[i].omega[w])
        b = self.stack.cell.cathode.r * self.stack.cell_list[i].t1[w] \
            / self.stack.cell.cathode.f_const
        # print(self.stack.cell_list[i].t1[w])
        c = np.log(x * self.stack.cell.c_ref /
                   (self.stack.cell.i_ref
                    * (self.stack.cell_list[i].cathode.gas_con[0][w] - self.stack.cell.delta * x)))
        # print(self.stack.cell_list[i].cathode.c1[w])
        return a - b * c - v

    def calc_initial_current_density(self):
        for i in range(self.stack.cell_numb):
            avv = sum(self.stack.cell_list[i].v) \
                  / self.stack.cell.cathode.channel.nodes
            for q in range(self.stack.cell.cathode.channel.nodes):
                self.stack.i[i, q] = fsolve(self.calc_initial_current_density_fsolve,
                                            self.stack.cell.cathode.tar_cd, args=(i, q, avv))

    def calc_initial_current_density_bd(self):
        for i in range(self.stack.cell_numb):
            avv = sum(self.stack.cell_list[i].v) \
                  / self.stack.cell.cathode.channel.nodes
            for q in range(self.stack.cell.cathode.channel.nodes):
                self.stack.i[i, q] = (least_squares(self.calc_initial_current_density_fsolve,
                                                    i_p.tar_cd,
                                                    bounds=(i_p.tar_cd / 3., i_p.tar_cd * 1.5),
                                                    args=(i, q, avv))).x

    def calc_convergence_criteria(self):
        self.mfd_criteria_cat = np.abs(sum(((self.stack.q_x_cat - self.q_x_cat_old)
                                            / self.stack.q_x_cat) ** 2))
        self.mdf_criteria_ano = np.abs(sum(((self.stack.q_x_ano - self.q_x_ano_old)
                                            / self.stack.q_x_ano ** 2)))
        self.i_criteria = np.abs(sum(((self.stack.i.flatten()
                                       - self.i_old.flatten())
                                      / self.stack.i.flatten() ** 2)))
        self.mdf_criteria_cat_process.append(self.mfd_criteria_cat)
        self.mdf_criteria_ano_process.append(self.mdf_criteria_ano)
        self.i_criteria_process.append(self.i_criteria)

    def save_old_value(self):
        self.q_x_cat_old = copy.deepcopy(self.stack.q_x_cat)
        self.q_x_ano_old = copy.deepcopy(self.stack.q_x_ano)
        self.i_old = copy.deepcopy(self.stack.i)

    def calc_n(self):  # electrical coupling on
        self.n = np.matmul(self.b, self.stack.v) - self.stack.resistance\
                 * np.matmul(self.c, self.stack.i.flatten(order='C'))

    def calc_g(self):
        self.s = np.diag(self.stack.dv)
        self.g = self.b*self.s - self.stack.resistance * self.c
        self.g_inv = np.linalg.inv(self.g)

    def correct_i(self):
        self.stack.i = self.stack.i / (np.median(self.stack.i) / i_p.tar_cd)

    def calc_i(self):
        self.calc_n()
        self.calc_g()
        i_pre_cor = self.i_old.flatten() - np.matmul(self.g_inv, self.n)
        self.stack.i = g_func.toarray(i_pre_cor,
                                      self.stack.cell_numb,
                                      self.stack.cell.cathode.channel.nodes)


channel_anode = channel.Channel(i_p.channel_length, i_p.nodes, i_p.k_ano,
                                i_p.p_ano_in, i_p.t_ano_in, i_p.phi_ano,
                                i_p.flow_dir_ano, i_p.channel_width,
                                i_p.channel_heigth)
channel_cathode = channel.Channel(i_p.channel_length, i_p.nodes, i_p.k_cat,
                                  i_p.p_cat_in, i_p.t_cat_in, i_p.phi_cat,
                                  i_p.flow_dir_cat, i_p.channel_width,
                                  i_p.channel_heigth)
anode = halfcell.Halfcell(channel_anode, i_p.spec_numb_ano,
                          i_p.val_ano, i_p.tar_cd)
cathode = halfcell.Halfcell(channel_cathode, i_p.spec_numb_cat,
                            i_p.val_cat, i_p.tar_cd)
cell1 = cell.Cell(anode, cathode, i_p.gamma, i_p.alpha, i_p.mem_thick, i_p.mu_p,
                  i_p.mu_g, i_p.mu_m, i_p.e_o, i_p.c_ref, i_p.i_ref, i_p.delta,
                  i_p.pem_type, i_p.coolant_channel, i_p.t_cool_in,
                  i_p.channel_width, i_p.gde_thick, i_p.plate_thick, i_p.t_u)
stack1 = stack.Stack(cell1, i_p.cell_numb, i_p.end_plate_t, i_p.alpha_conv,
                     i_p.resistance, i_p.endplate_adiabat, i_p.manifold_heigth,
                     i_p.manifold_width, i_p.manifold_kf, i_p.stoi_cat, i_p.stoi_ano)
