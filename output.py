import matplotlib.pyplot as plt
import numpy as np
import simulation as sim
import input_parameter as i_p


class Output:

    def __init__(self, simulation):
        self.simulation = simulation
        self.x = np.linspace(0., self.simulation.stack.cell.cathode.channel.length,
                             self.simulation.stack.cell.cathode.channel.nodes)

    def display_files(self):
        plt.plot(self.simulation.mdf_criteria_ano_process, label='Anode')
        plt.plot(self.simulation.mdf_criteria_cat_process, label='Cathode')
        plt.plot(self.simulation.i_criteria_process, label='Current density')
        plt.title('Convergence')
        plt.xlabel('Iteration')
        plt.yscale('log')
        plt.ylabel('ERR')
        plt.legend()
        plt.show()

        plt.plot(self.simulation.stack.q_x_cat
                 / (self.simulation.stack.q_h_in_cat[-1]
                    / self.simulation.stack.cell_numb), label='Cathode')
        plt.plot(self.simulation.stack.q_x_ano
                 / (self.simulation.stack.q_h_in_ano[-1]
                    / self.simulation.stack.cell_numb), label='Anode')
        plt.xlim(self.simulation.stack.cell_numb, 0)
        plt.legend()
        plt.title('Flow distribution')
        plt.xlabel('Cell numb')
        plt.ylabel(r'$\frac{q(i)}{q(avg)}$', fontsize=14)
        plt.show()

        # for q in self.simulation.stack.cell_list:
        #     plt.plot(self.x, q.cathode.u)
        # plt.title('Cathode channel velocity')
        # plt.show()

        # for q in self.simulation.stack.cell_list:
        #     plt.plot(self.x, q.anode.u)
        # plt.title('Anode channel velocity')
        # plt.show()

        # for q in self.simulation.stack.cell_list:
        #     plt.plot(self.x, q.cathode.re)
        # plt.title('Cathode channel reynolds number')
        # plt.show()

        # for q in self.simulation.stack.cell_list:
        #     plt.plot(self.x, q.anode.re)
        # plt.title('Anode channel reynolds number')
        # plt.show()

        # for q in self.simulation.stack.cell_list:
        #    plt.plot(self.x, q.cathode.roh)
        # plt.title('Cathode channel density')
        # plt.show()

        # for q in self.simulation.stack.cell_list:
        #     plt.plot(self.x, q.anode.roh)
        # plt.title('Anode channel density')
        # plt.show()

        a = 0
        for q in self.simulation.stack.cell_list:
            a = sum(q.i)+a
            plt.plot(self.x, q.i)
        plt.ylim(i_p.tar_cd*0.25, i_p.tar_cd*1.5)
        plt.xlim(0, i_p.channel_length)
        plt.title('Current density')
        # print('Average current density:', a/(i_p.nodes*i_p.cell_numb))
        plt.show()

        for q in self.simulation.stack.cell_list:
            plt.plot(self.x, q.v)
        plt.ylim(0., 1.28)
        plt.xlim(0, i_p.channel_length)
        plt.title('Voltage')
        plt.show()

        for i in self.simulation.stack.cell_list:
            plt.plot(self.x, i.dv)
        plt.title('Dv/di')
        plt.show()

        for i in self.simulation.stack.cell_list:
            plt.plot(self.x, i.psi)
        plt.title('Psi')
        plt.show()

        for q in self.simulation.stack.cell_list:
            plt.plot(self.x, q.omega)
            #plt.plot(self.x, q.omega_goessling)
        #plt.ylim(1.e-4, 1.e-5)
        plt.xlim(0, i_p.channel_length)
        plt.title('Omega')
        plt.show()

        for q in self.simulation.stack.cell_list:
            plt.plot(self.x, q.t)
        plt.xlim(0, i_p.channel_length)
        plt.title('Coolant temperature')
        plt.show()

        for q in self.simulation.stack.cell_list:
            plt.plot(self.x, q.t0)
        plt.xlim(0, i_p.channel_length)
        plt.title('Membran temperature')
        plt.show()

        for q in self.simulation.stack.cell_list:
            plt.plot(self.x, q.t1)
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode temperature')
        plt.show()

        for q in self.simulation.stack.cell_list:
            plt.plot(self.x, q.t2)
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode gas channel temperature')
        plt.show()

        for q in self.simulation.stack.cell_list:
            plt.plot(self.x, q.t3)
        plt.xlim(0, i_p.channel_length)
        plt.title('Coolant plate temperature')
        plt.show()

        for q in self.simulation.stack.cell_list:
            plt.plot(self.x, q.t4)
        plt.xlim(0, i_p.channel_length)
        plt.title('Anode temperature')
        plt.show()

        for q in self.simulation.stack.cell_list:
            plt.plot(self.x, q.t5)
        plt.xlim(0, i_p.channel_length)
        plt.title('Anode gas channel temperature')
        plt.show()

        for q in self.simulation.stack.cell_list:
            plt.plot(self.x, q.cathode.p)
        #plt.ylim(i_p.p_cat_in-10, i_p.p_cat_in)
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode channel pressure')
        plt.show()

        for q in self.simulation.stack.cell_list:
            plt.plot(self.x, q.anode.p)
        #plt.ylim(i_p.p_ano_in-10, i_p.p_ano_in)
        plt.xlim(0, i_p.channel_length)
        plt.title('Anode channel pressure')
        plt.show()

        for q in self.simulation.stack.cell_list:
            plt.plot(self.x, q.cathode.gas_flow[0])
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode oxygen molar flow')
        plt.show()

        for q in self.simulation.stack.cell_list:
            plt.plot(self.x, q.cathode.gas_flow[1])
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode water molar flow')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.cathode.gas_flow[2)
            plt.plot(self.x, q.cathode.gas_flow[2])
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode nitrogen molar flow')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.anode.gas_flow[0)
            plt.plot(self.x, q.anode.gas_flow[0])
        plt.xlim(0, i_p.channel_length)
        plt.title('Anode hydrogen molar flow')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.anode.gas_flow[2)
            plt.plot(self.x, q.anode.gas_flow[1])
        plt.xlim(0, i_p.channel_length)
        plt.title('Anode water molar flow')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.cathode.gas_con[0)
            plt.plot(self.x, q.cathode.gas_con[0])
        plt.xlim(0, i_p.channel_length)
        plt.ylim(0, 100)
        plt.title('Cathode oxygen concentration')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.cathode.gas_con[1)
            plt.plot(self.x, q.cathode.gas_con[1])
        plt.xlim(0, i_p.channel_length)
        plt.ylim(0, 100)
        plt.title('Cathode water vapour concentration')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.anode.gas_con[0)
            plt.plot(self.x, q.anode.gas_con[0])
        plt.xlim(0, i_p.channel_length)
        #plt.ylim(0, 100)
        plt.title('Anode hydrogen concentration')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.anode.gas_con[1)
            plt.plot(self.x, q.anode.gas_con[1])
        plt.xlim(0, i_p.channel_length)
        plt.ylim(0, 100)
        plt.title('Anode water vapour concentration')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.cathode.w)
            plt.plot(self.x, q.cathode.w)
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode liquid water flow')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.anode.w)
            plt.plot(self.x, q.anode.w)
        plt.xlim(0, i_p.channel_length)
        plt.title('Anode liquid water flow')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.cathode.gamma)
            plt.plot(self.x, q.cathode.gamma)
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode condensation rate')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.anode.gamma)
            plt.plot(self.x, q.anode.gamma)
        plt.xlim(0, i_p.channel_length)
        plt.title('Anode condensation rate')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.cathode.humidity)
            plt.plot(self.x, q.cathode.humidity)
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode humidity')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.anode.humidity)
            plt.plot(self.x, q.anode.humidity)
        plt.xlim(0, i_p.channel_length)
        plt.title('Anode humidity')
        plt.show()

        if i_p.pem_type is False:

            for q in self.simulation.stack.cell_list:
                # print(i.cathode.free_water)
                plt.plot(self.x, q.cathode.free_water)
            plt.xlim(0, i_p.channel_length)
            plt.title('Cathode free water')
            plt.show()

            for q in self.simulation.stack.cell_list:
                # print(i.anode.free_water)
                plt.plot(self.x, q.anode.free_water)
            plt.xlim(0, i_p.channel_length)
            plt.title('Cathode free water')
            plt.show()

            for q in self.simulation.stack.cell_list:
                plt.plot(self.x, q.j)
            plt.xlim(0, i_p.channel_length)
            plt.title('Water crossover flux')
            plt.show()


simulation1 = sim.Simulation(sim.stack1, i_p.k_tol, i_p.max_it)
simulation1.update()
output = Output(simulation1)
output.display_files()
