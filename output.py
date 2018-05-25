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

        a = 0
        for q in self.simulation.stack.cell_list:
            a = sum(q.i)+a
            plt.plot(self.x, q.i)
        plt.ylim(i_p.tar_cd*0.25, i_p.tar_cd*1.5)
        plt.xlim(0, i_p.channel_length)
        plt.title('Current density')
        print('Average current density:', a/(i_p.nodes*i_p.cell_numb))
        plt.show()

        for q in self.simulation.stack.cell_list:
            plt.plot(self.x, q.v)
        plt.ylim(0., 1.28)
        plt.xlim(0, i_p.channel_length)
        plt.title('Voltage')
        plt.show()

        # print('dvdi')
        # for i in self.simulation.stack.cell_list:
        #   #print(i.dv)
        #    plt.plot(self.x,i.dv)
        # plt.show()

        # print('psi')
        # for i in self.simulation.stack.cell_list:
        #    #print(i.psi)
        #    plt.plot(self.x,i.psi)
        # plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.omega)
            plt.plot(self.x, q.omega)
        plt.ylim(1.e-4, 1.e-5)
        plt.xlim(0, i_p.channel_length)
        plt.title('Omega')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.t)
            plt.plot(self.x, q.t)
        plt.ylim(i_p.t_cool_in, 450.)
        plt.xlim(0, i_p.channel_length)
        plt.title('Coolant temperature')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.T0)
            plt.plot(self.x, q.t0)
        plt.ylim(280., 450.)
        plt.xlim(0, i_p.channel_length)
        plt.title('Membran temperature')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.T1)
            plt.plot(self.x, q.t1)
        plt.ylim(280., 450.)
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode temperature')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.T2)
            plt.plot(self.x, q.t2)
        plt.ylim(280., 450.)
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode gas channel temperature')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.T3)
            plt.plot(self.x, q.t3)
        plt.ylim(310., 320.)
        plt.xlim(0, i_p.channel_length)
        plt.title('Coolant plate temperature')
        plt.show()

        for q in self.simulation.stack.cell_list:
            #print(q.t4)
            plt.plot(self.x, q.t4)
        plt.ylim(280., 450.)
        plt.xlim(0, i_p.channel_length)
        plt.title('Anode temperature')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.T5)
            plt.plot(self.x, q.t5)
        plt.ylim(280., 450.)
        plt.xlim(0, i_p.channel_length)
        plt.title('Anode gas channel temperature')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.cathode.p)
            plt.plot(self.x, q.cathode.p)
        plt.ylim(i_p.p_cat_in-10, i_p.p_cat_in)
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode channel pressure')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.anode.p)
            plt.plot(self.x, q.anode.p)
        plt.ylim(i_p.p_ano_in-10, i_p.p_ano_in)
        plt.xlim(0, i_p.channel_length)
        plt.title('Anode channel pressure')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.cathode.q1)
            plt.plot(self.x, q.cathode.q1)
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode oxygen molar flux')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.cathode.q2)
            plt.plot(self.x, q.cathode.q2)
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode water molar flux')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.cathode.q3)
            plt.plot(self.x, q.cathode.q3)
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode nitrogen molar flux')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.anode.q1)
            plt.plot(self.x, q.anode.q1)
        plt.xlim(0, i_p.channel_length)
        plt.title('Anode hydrogen molar flux')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.anode.q2)
            plt.plot(self.x, q.anode.q2)
        plt.xlim(0, i_p.channel_length)
        plt.title('Anode water molar flux')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.cathode.c1)
            plt.plot(self.x, q.cathode.c1)
        plt.xlim(0, i_p.channel_length)
        plt.ylim(0, 100)
        plt.title('Cathode oxygen concentration')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.cathode.c2)
            plt.plot(self.x, q.cathode.c2)
        plt.xlim(0, i_p.channel_length)
        plt.ylim(0, 100)
        plt.title('Cathode water vapour concentration')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.anode.c1)
            plt.plot(self.x, q.anode.c1)
        plt.xlim(0, i_p.channel_length)
        #plt.ylim(0, 100)
        plt.title('Anode hydrogen concentration')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.anode.c2)
            plt.plot(self.x, q.anode.c2)
        plt.xlim(0, i_p.channel_length)
        plt.ylim(0, 100)
        plt.title('Anode water vapour concentration')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.cathode.w)
            plt.plot(self.x, q.cathode.w)
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode liquid water flux')
        plt.show()

        for q in self.simulation.stack.cell_list:
            # print(i.cathode.gamma)
            plt.plot(self.x, q.cathode.gamma)
        plt.xlim(0, i_p.channel_length)
        plt.title('Cathode condensation rate')
        plt.show()

        if i_p.pem_type is False:
            print('wc')
            for q in self.simulation.stack.cell_list:
                # print(i.cathode.w)
                plt.plot(self.x, q.cathode.w)
            plt.xlim(0, i_p.channel_length)
            plt.show()

            print('wa')
            for q in self.simulation.stack.cell_list:
                # print(i.anode.w)
                plt.plot(self.x, q.anode.w)
            plt.xlim(0, i_p.channel_length)
            plt.show()

            print('gamma_c')
            for q in self.simulation.stack.cell_list:
                # print(i.cathode.gamma)
                plt.plot(self.x, q.cathode.gamma)
            plt.xlim(0, i_p.channel_length)
            plt.show()

            print('gamma_a')
            for q in self.simulation.stack.cell_list:
                # print(i.anode.gamma)
                plt.plot(self.x, q.anode.gamma)
            plt.xlim(0, i_p.channel_length)
            plt.show()

            print('humidity_c')
            for q in self.simulation.stack.cell_list:
                # print(i.cathode.humidity)
                plt.plot(self.x, q.cathode.humidity)
            plt.xlim(0, i_p.channel_length)
            plt.show()

            print('humidity_a')
            for q in self.simulation.stack.cell_list:
                # print(i.anode.humidity)
                plt.plot(self.x, q.anode.humidity)
            plt.xlim(0, i_p.channel_length)
            plt.show()

            print('free_water_c')
            for q in self.simulation.stack.cell_list:
                # print(i.cathode.free_water)
                plt.plot(self.x, q.cathode.free_water)
            plt.xlim(0, i_p.channel_length)
            plt.show()

            print('free_water_a')
            for q in self.simulation.stack.cell_list:
                # print(i.anode.free_water)
                plt.plot(self.x, q.anode.free_water)
            plt.xlim(0, i_p.channel_length)
            plt.show()


simulation1 = sim.Simulation(sim.stack1, i_p.k_tol, i_p.max_it)
simulation1.update()
output = Output(simulation1)
output.display_files()
