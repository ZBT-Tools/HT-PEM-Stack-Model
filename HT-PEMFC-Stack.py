from pylab import*
from scipy.optimize import fsolve
import copy
import saturation_pressure_vapour
import matrix_database


# # Globabal Variables and Matrix
pem_type = True      #True=HT False=NT
I_t = 6000.          #Target currentdensity                      A/m²
Fday = 96485.        #Faraday's constant                         C/mol
R = 8.3143           #ideal gas constant                        J/(kmol)
N = 50                #knots, elements = N-1
M = 4              #Number of cell
Tu = 298.15 
node_backward = matrix_database.backward_matrix(N)
element_backward = matrix_database.backward_matrix(N-1)
node_forward = matrix_database.forward_matrix(N)
element_forward = matrix_database.forward_matrix(N-1)
node_to_element = matrix_database.node_to_element(N)
element_to_node = matrix_database.element_to_node(N-1)


# # Classes and Functions

# Global functions:

def dw (T):
    return 2.1*10.**-7.*exp(-2436./T)
    
def d_pos(T):
    return 1.6*10.**-8.*exp(-1683./T)

def toarray(var):
    return reshape(var.flatten(order='C'),(M,N))


# 0.Layer

class Channel:
    
    def __init__(self,length,division,k,p_in,T_in,flow_dir):
        self.length = length
        self.division = division
        self.d_x = length/division
        self.k = k
        self.p_in = p_in
        self.T_in = T_in
        self.flow_dir = flow_dir


# 1.Layer

class Halfcell:
    def __init__(self,channel,stoi,spec_numb,val):
        #object channel, stoichiometry
        #number of species and layer, number of valence electrones
        self.channel = channel
        self.stoi = stoi
        self.spec_numb = spec_numb
        self.val = val
        self.pem_type = True #True for HT-PEM  False for NT-PEM
        self.we = full((self.channel.division),0.)
        self.w = full((self.channel.division+1),0.)
        self.pe = full((self.channel.division),self.channel.p_in)
        self.p = full((self.channel.division+1),self.channel.p_in)
        self.gamma = full((self.channel.division+1),0.)
        self.humidity = full((self.channel.division+1),0.)
        self.free_water = full((self.channel.division+1),0.)
        self.i = full((self.channel.division+1),I_t)
        self.ie = full((self.channel.division),I_t)
        for l in range(self.spec_numb):#Init flow and con arrays (1-n)
            exec("self.q%d = full((self.channel.division+1),0.)" % (l + 1))
            #1 reactant(hydrogen/oxygen), 2 water, 3 (if oxygen) nitrogen
            exec("self.c%d = full((self.channel.division+1),0.)" % (l + 1))
            #1 reactant(hydrogen/oxygen), 2 water, 3 (if oxygen) nitrogen
            exec('self.T%d = full((self.channel.division+1),self.channel.T_in)' %(l + 1))
            #1 catalyst layer, 2 channel layer, 3 coolant plate layer
        
    def update(self):
        if self.pem_type is True:#HT
            self.node_to_element()
            self.calc_reac_flow()
            self.calc_water_flow()#stabel with j=0?
            self.calc_pressure()
            self.calc_con_ic()
            self.calc_fluid_water()
            self.calc_cond_rates()
        else:#NT
            self.node_to_element()
            self.calc_reac_flow()
            self.calc_water_flow()
            self.calc_pressure()
            self.calc_con()
            self.calc_fluid_water()
            self.calc_cond_rates()
            self.calc_rel_humidity()
            self.calc_free_water()
            
    
    def set_i(self,i):
        self.i = i
    
    def set_j(self,j):
        self.j = j
    
    def set_pem_type(self,pem_type):
        self.pem_type = pem_type
    
    def set_T(self,var):#übergabe der Layertemperaturen als Array(T1-T2/3)
        for l, item in enumerate(var):
            exec('self.T%d = var[l]'%(l+1))
            
    def node_to_element(self):
        self.ie = matmul(node_to_element,self.i)
        if self.pem_type is False:self.je = matmul(node_to_element,self.j)
        self.T1e = matmul(node_to_element,self.T1)
        self.T2e = matmul(node_to_element,self.T2)
        if self.spec_numb is 3: self.T3e = matmul(node_to_element,self.T3)
        
    def calc_reac_flow(self):#elementwise
        if self.channel.flow_dir ==True:
            matrix = element_forward
            var1 = self.stoi*I_t*self.channel.length/(self.val*Fday)
            self.q1e = var1-matmul(matrix,self.ie)*self.channel.d_x/(self.val*Fday)
            var2 = self.q1e[-1]-0.25*(self.i[-1]+self.ie[-1])*self.channel.d_x/(self.val*Fday)
            self.q1 = matrix_database.element_to_node_func(var1,var2,self.q1e)
        
        else:
            matrix = element_backward
            var1 = self.stoi*I_t*self.channel.length/(self.val*Fday)
            self.q1e = var1-matmul(matrix,self.ie)*self.channel.d_x/(self.val*Fday)
            var2 = self.q1e[0]-0.25*(self.i[0]+self.ie[0])*self.channel.d_x/(self.val*Fday)
            self.q1 = matrix_database.element_to_node_func(var2,var1,self.q1e)
            ###working do not change!
                    
    def calc_water_flow(self):#water molar flux in the channel#elementwise
        if self.pem_type is True:
            if self.channel.flow_dir == True:
                start_cell = 0
                matrix = element_forward
                a = reshape(self.channel.d_x/(self.val*Fday*0.5)*matmul(matrix,self.ie),self.channel.division)
                #production
                b = saturation_pressure_vapour.water.calc_psat(self.channel.T_in)
                c = (0.79/0.21+1)*self.q1[start_cell]*b/(self.channel.p_in-b)#Value at q2[0]
                self.q2e = a+c
                var1 = 0.5*(self.ie[-1]+self.i[-1])/(self.val*0.5*Fday)
                q2last_node = self.q2e[-1] + self.channel.d_x/2.*(var1) 
                self.q2 = matrix_database.element_to_node_func(c,q2last_node,self.q2e)
                self.q3 = full((self.channel.division+1),self.q1[0]*0.79/0.21)
                self.q3e = full((self.channel.division),self.q1[0]*0.79/0.21)

            elif self.channel.flow_dir == False:
                matrix = element_backward
                b = saturation_pressure_vapour.water.calc_psat(self.channel.T_in)
                self.q2e = full(self.channel.division,self.q1[-1]*b/(self.channel.p_in-b))
                self.q2 = full(self.channel.division+1,self.q1[-1]*b/(self.channel.p_in-b))
        else:
            if self.channel.flow_dir == True:
                start_cell = 0
                matrix = element_forward
                a = reshape(self.channel.d_x/(self.val*Fday*0.5)*matmul(matrix,self.ie),self.channel.division)
                #production
                b = saturation_pressure_vapour.water.calc_psat(self.channel.T_in)
                c = (0.79/0.21+1)*self.q1[start_cell]*b/(self.channel.p_in-b)#Value at q2[0]
                d = reshape(self.channel.d_x * matmul(matrix,self.je),self.channel.division)
                #crossover
                e = a + c + d
                self.q2e = e
                var1 = 0.5*(self.ie[-1]+self.i[-1])/(self.val*0.5*Fday)
                var2 = 0.5*(self.je[-1]+self.j[-1])
                q2last_node = self.q2e[-1] + self.channel.d_x/2.*(var1+var2) 
                self.q2 = matrix_database.element_to_node_func(c,q2last_node,self.q2e)
                self.q3 = full((self.channel.division+1),self.q1[0]*0.79/0.21)
                self.q3e = full((self.channel.division),self.q1[0]*0.79/0.21)

            elif self.channel.flow_dir == False:
                last_cell = -1
                matrix = element_backward
                b = saturation_pressure_vapour.water.calc_psat(self.channel.T_in)
                c = self.q1[last_cell]*b/(self.channel.p_in-b)
                e = self.channel.d_x*matmul(-matrix,self.je)
                self.q2e = e + c
                var1 = 0.5*(self.je[0]+self.j[0])
                q2first_node = self.q2e[0] -self.channel.d_x/2.*var1
                self.q2 = matrix_database.element_to_node_func(q2first_node,c,self.q2e)    
    
    def calc_pressure(self):#channel pressure #qo qc qn / qh qa#elementwise  
        a = self.channel.p_in
        b = self.channel.d_x*self.channel.k
        if self.spec_numb == 3:
            c = self.q1e+self.q2e+self.q3e-self.we
            d = R*self.T1e/self.pe 
            self.pe = a-b*matmul(element_forward,c*d)
            var1 = -self.channel.k*self.channel.d_x*0.25
            var2 = self.q1e[-1] + self.q2e[-1] + self.q3e[-1] - self.we[-1]
            var3 = self.q1[-1] + self.q2[-1] + self.q3[-1] - self.w[-1]
            var4 = R*(self.T2e[-1]+self.T2[-1])/(self.pe[-1]+self.p[-1])
            plast_node = self.pe[-1]+var1*(var2+var3)*var4
            self.p = matrix_database.element_to_node_func(self.channel.p_in,plast_node,self.pe)
            
        else:
            c = self.q1e+self.q2e-self.we
            d = R*self.T1e/self.pe 
            self.pe = a-b*matmul(element_backward,c*d) 
            var1 = -self.channel.k*self.channel.d_x*0.25
            var2 = self.q1e[0] + self.q2e[0]  - self.we[0]
            var3 = self.q1[0] + self.q2[0] +  - self.w[0]
            var4 = R*(self.T2e[0]+self.T2[0])/(self.pe[0]+self.p[0])
            pfirst_node = self.pe[0]+var1*(var2+var3)*var4
            ###working do just change for improvement
                
    def calc_con(self):# Gas cocentrations in the channel [moles/m³]#nodewise
        if self.spec_numb==3:###cathode
            for w in range (self.channel.division+1):
                var1 = self.p[w]/(R*self.T2[w])
                var2 = self.q2[w]/(self.q1[w] + self.q2[w] + self.q3[w])
                self.c2[w] = var1*var2
                if saturation_pressure_vapour.water.calc_psat(self.T2[w])/(R*self.T2[w])<=self.c2[w]:
                    #print(self.T2[w])
                    #print(saturation_pressure_vapour.water.calc_psat(self.T2[w]))
                    a = (self.p[w]-saturation_pressure_vapour.water.calc_psat(self.T2[w]))/(R*self.T2[w])
                    b = self.q1[w]/(self.q1[w] + self.q3[w])
                    self.c1[w] = a*b
                    self.c2[w] = saturation_pressure_vapour.water.calc_psat(self.T2[w])/(R*self.T2[w])
                else:
                    var3 = self.p[w]/(R*self.T2[w])
                    var4 = self.q1[w]/(self.q1[w] + self.q2[w] + self.q3[w])
                    self.c1[w] = var3*var4
        else:
            for w in range (self.channel.division+1):
                var1 = self.p[w]/(R*self.T2[w])
                var2 = self.q2[w]/(self.q1[w] + self.q2[w])
                self.c2[w] = var1*var2
                if saturation_pressure_vapour.water.calc_psat(self.T2[w])/(R*self.T2[w])<=self.c2[w]:
                    self.c1[w] = (self.p[w]-saturation_pressure_vapour.water.calc_psat(self.T2[w]))/(R*self.T2[w])
                    self.c2[w] = saturation_pressure_vapour.water.calc_psat(self.T2[w])/(R*self.T2[w])
                else:
                    var3 = self.p[w]/(R*self.T2[w])
                    var4 = self.q1[w]/(self.q1[w] + self.q2[w])
                    self.c1[w] = var3*var4
        self.c1e = matmul(node_to_element,self.c1)
        self.c2e = matmul(node_to_element,self.c2)
        
    def calc_con_ic(self):# Gas cocentrations in the channel [moles/m³]#nodewise
        #no liquid water
        if self.spec_numb==3:###cathode
            for w in range (self.channel.division+1):
                    a = (self.p[w]-saturation_pressure_vapour.water.calc_psat(self.T2[w]))/(R*self.T2[w])
                    b = self.q1[w]/(self.q1[w] + self.q3[w])
                    self.c1[w] = a*b
                    self.c2[w] = saturation_pressure_vapour.water.calc_psat(self.T2[w])/(R*self.T2[w])
        else:
            for w in range (self.channel.division+1):
                    self.c1[w] = (self.p[w]-saturation_pressure_vapour.water.calc_psat(self.T2[w]))/(R*self.T2[w])
                    self.c2[w] = saturation_pressure_vapour.water.calc_psat(self.T2[w])/(R*self.T2[w])
        self.c1e = matmul(node_to_element,self.c1)
        self.c2e = matmul(node_to_element,self.c2)
        
    
    def calc_fluid_water(self):#fluid water in the channel#nodewise
        self.w = self.q2 - self.c2/self.c1*self.q1
        self.we = matmul(node_to_element,self.w)
        #ok
    
    def calc_cond_rates(self):#condensation rates of the vapour in the channel#nodewise
        if self.spec_numb is 3:self.gamma = gradient(self.w,self.channel.d_x)
        else:self.gamma = -gradient(self.w,self.channel.d_x) 
        self.gammae = matmul(node_to_element,self.gamma)
        #element first -> element to node ?
    
    def calc_rel_humidity(self):#relative humidity in the channel#nodewise
        self.humidity = self.c2*R*self.T2/saturation_pressure_vapour.water.calc_psat(self.T2)
        self.humiditye = matmul(node_to_element,self.humidity)
        #ok
        
    def calc_free_water(self):#membrane free water content#nodewise
        a = 0.043 + 17.81*self.humidity
        b = -39.85*self.humidity**2 + 36.*self.humidity**3
        self.free_water = a + b
        self.free_watere = matmul(node_to_element,self.free_water)
        #ok


# 2.Layer

class Cell:
    
    def __init__(self,anode,cathode,Gamma,alpha,mem_thic,mu_p,mu_g,mu_m,e0,c_ref,i_ref,delta,T_cool_in,pem_type):
        self.anode = anode
        self.cathode = cathode
        self.Gamma = Gamma
        self.alpha = alpha
        self.mem_thic = mem_thic
        self.mu_p = mu_p
        self.mu_g = mu_g
        self.mu_m = mu_m
        self.e0 = e0
        self.c_ref = c_ref
        self.i_ref = i_ref
        self.delta = delta
        self.T_cool_in = T_cool_in
        self.pem_type = pem_type
        self.j = full(self.cathode.channel.division+1,0.)
        self.Ke = 0.
        self.Ho = -52300.
        self.Ko = 6.2
        self.m_a = full((self.cathode.channel.division+1),0.)
        self.m_c = full((self.cathode.channel.division+1),0.)
        self.m0 = full((self.cathode.channel.division+1),0.)
        self.m1 = full((self.cathode.channel.division+1),0.)
        self.m_pos0 = full((self.cathode.channel.division+1),0.)
        self.m_pos1 = full((self.cathode.channel.division+1),0.)
        self.omega = full((self.cathode.channel.division+1),0.)
        self.psi = full((self.cathode.channel.division+1),0.)
        self.T = full(self.cathode.channel.division+1,self.T_cool_in)
        self.Te = full(self.cathode.channel.division,self.T_cool_in)
        self.T0 = (self.anode.T1 + self.cathode.T1)/2. 
        self.T1 = self.cathode.T1
        self.T2 = self.cathode.T2
        self.T3 = self.cathode.T3
        self.T4 = self.anode.T1
        self.T5 = self.anode.T2
        self.vtn = 1.28
                
    def update(self):
        if self.pem_type is False:
            self.cathode.set_pem_type(False)
            self.cathode.set_i(self.i)
            self.cathode.set_j(self.j)
            self.cathode.set_T([self.T1,self.T2,self.T3])
            self.cathode.update()
            self.anode.set_pem_type(False)
            self.anode.set_i(self.i)
            self.anode.set_j(self.j)
            self.anode.set_T([self.T4,self.T5])
            self.anode.update()
            self.calc_mem_block_1()
            self.calc_cross_over_water_flux()
            self.calc_mem_block_2()
            self.calc_con_overpotential()
            self.calc_mem_resistivity()
            self.calc_voltage()
            self.calc_dvdi()
            self.T0 =  0.5*self.anode.T1 + 0.5*self.cathode.T1
            self.node_to_element()
        else:
            self.cathode.set_pem_type(True)
            self.cathode.set_i(self.i)
            self.cathode.set_T([self.T1,self.T2,self.T3])
            self.cathode.update()
            self.cathode.set_pem_type(True)
            self.anode.set_i(self.i)
            self.anode.set_T([self.T4,self.T5])
            self.anode.update()
            #self.calc_con_overpotential() need a different function here
            self.psi = full((self.cathode.channel.division+1),0.)
            self.calc_mem_resistivity()
            self.calc_voltage()
            self.calc_dvdi()
            self.T0 =  0.5*self.anode.T1 + 0.5*self.cathode.T1
            self.node_to_element()


    def node_to_element(self):
        self.T0e = matmul(node_to_element,self.T0)
        self.T1e = matmul(node_to_element,self.T1)
        self.T2e = matmul(node_to_element,self.T2)
        self.T3e = matmul(node_to_element,self.T3)
        self.T4e = matmul(node_to_element,self.T4)
        self.T5e = matmul(node_to_element,self.T5)
        
    def set_i(self,i):
        self.i = i
                    
    def calc_mem_block_1(self):#nodewise
        a = self.cathode.free_water + self.anode.free_water
        b = reshape(self.i/(2.*self.Gamma*self.alpha*Fday),self.cathode.channel.division+1)
        self.zeta_plus = a + b
        c = self.cathode.free_water - self.anode.free_water
        d = reshape(5.*self.i/(2.*self.Gamma*self.alpha*Fday),self.cathode.channel.division+1)
        e = 1. + dw(self.T0) * self.zeta_plus / (self.mem_thic*self.Gamma)
        self.zeta_negative = (c+d)/e
        self.m_c = 0.5*(self.zeta_plus + self.zeta_negative)
        self.m_a = 0.5*(self.zeta_plus - self.zeta_negative)
        #ok
        
        
    def calc_cross_over_water_flux(self):# nodewise
        a = reshape(self.i/Fday,self.cathode.channel.division+1)
        aa = self.alpha*dw(self.T0)
        b = (self.m_a**2.-self.m_c**2.)/(2.*self.mem_thic)
        self.j = a+aa*b
        #ok
        
    def calc_mem_block_2(self):#nodewise
        self.Ke = self.Ko*exp(-self.Ho/R*(1./self.T0-1./Tu))
        a = self.m_c**2.
        b = self.m_a**2.-self.m_c**2.
        self.m0 = sqrt(a-b)#z = 0
        self.m1 = sqrt(a)#z = 1
        self.m_pos0 = -self.Ke*self.m0*0.5 + sqrt((self.Ke*self.m0*0.5)**2 + self.Ke*self.m0)
        self.m_pos1 = -self.Ke*self.m1*0.5 + sqrt((self.Ke*self.m1*0.5)**2 + self.Ke*self.m1)
        #ok
        
    def calc_con_overpotential(self):#nodewise
        self.psi = R*self.T0/Fday * log(self.m_pos1/self.m_pos0)
        #ok
    
    def calc_mem_resistivity(self):#nodewise
        self.omega = (0.4025-0.0007*self.T0)/10000.
        #ok
        
    def calc_voltage(self):#nodewise
        a = self.e0 - self.i*self.omega
        #print(a)
        b = R*self.T1/(Fday)
        #print(b)
        c = log(self.i*self.c_ref/(self.i_ref*(self.cathode.c1-self.delta*self.i)))
        #print(c)
        self.v = a+self.psi-b*c
        #ok
        
    def calc_dvdi(self):#nodewise
        a = -self.omega
        b = R*self.T1*self.cathode.c1
        c = Fday*self.i*(-self.cathode.c1-self.delta*self.i)
        self.dv = a+b/c
        #ok


# 3.Layer

class Stack:
    
    def __init__(self,cell,cell_numb):
        self.cell_numb = cell_numb
        self.cell = cell    
        self.cell_list = []
        for i in range(self.cell_numb):self.cell_list.append(copy.deepcopy(cell))
        self.resistance = 4.*10**(-1)
        self.i = full((cell_numb,cell.cathode.channel.division+1),I_t)
        self.T = full((cell_numb,cell.cathode.channel.division+1),self.cell.T_cool_in)
        self.Te = full((cell_numb,cell.cathode.channel.division),self.cell.T_cool_in)
        self.g = 529.# coolant flo
        self.h_vap = 45400.# vaporization heat
        self.a = 4000. # scaled heat transfer factor to the coolant
        #print(len(self.mat))
        
    def update(self):
        for j in range(self.cell_numb):
            self.cell_list[j].set_i(self.i[j,:])
            self.cell_list[j].update()

        #self.calc_coolant_T()
        self.calc_coolant_T_fit()
        self.calc_layer_T()
        self.stack_v()
        self.stack_dv()       
    
    def set_i(self,i):
        self.i = i
    
    def stack_v(self):
        var = []
        for i, item in enumerate(self.cell_list):var = hstack((var,self.cell_list[i].v))
        self.v = var
        #running
    
    def stack_dv(self):
        var = []
        for i, item in enumerate(self.cell_list):var = hstack((var,self.cell_list[i].dv))
        self.dv = var
        #running

    def calc_coolant_T(self):
        for q, item in enumerate(self.cell_list):
            if q==0:
                for w in range(self.cell.cathode.channel.division+1):
                    if w ==0:
                        self.cell_list[q].T[w] = self.cell_list[q].T_cool_in
                    else:
                        var1 = self.cell.mu_p/self.g*self.cell.cathode.channel.d_x
                        var2 =  self.cell_list[q].T2[w] -1 *self.cell_list[q].T3[w]
                        self.cell_list[q].T[w] = var1*var2+self.cell_list[q].T[w-1]
                        self.cell_list[q].T[w] = self.cell.cathode.channel.T_in + w / (
                                    self.cell.cathode.channel.division + 1)
            else:
                for w in range(self.cell.cathode.channel.division+1):
                    if w ==0:
                        self.cell_list[q].T[w] = self.cell_list[q].T_cool_in
                    else:
                        var1 = self.cell.mu_p / self.g * self.cell.cathode.channel.d_x
                        var2 = self.cell_list[q].T2[w]+ self.cell_list[q-1].T5[w] - 2 * self.cell_list[q].T3[w]
                        self.cell_list[q].T[w] = var1 * var2+self.cell_list[q].T[w-1]


    def calc_coolant_T_fit(self):
        for q, item in enumerate(self.cell_list):
            self.cell_list[q].T = linspace(self.cell.T_cool_in, self.cell.cathode.channel.T_in,
                                       self.cell.cathode.channel.division + 1)

        
    def calc_T1(self):
        for q, item in enumerate(self.cell_list):
            f2 = (self.cell.vtn - self.cell_list[q].v - self.cell_list[q].omega / 2. * self.cell_list[q].i) * \
                 self.cell_list[q].i
            f3 = 0.5 * self.cell_list[q].omega * self.cell_list[q].i ** 2.
            f4 = self.h_vap*self.cell_list[q].anode.gamma
            if q >= 1:
                f5 = -self.cell.mu_p*self.cell_list[q-1].T5-self.a*self.cell_list[q].T
            else:
                f5 = -self.cell.mu_p*self.cell_list[q].T5-self.a*self.cell_list[q].T
            var1 = self.cell.mu_p*(self.cell_list[q].T2*self.cell.mu_p-f5)/(2.*self.cell.mu_p+self.a)
            var2 = self.cell.mu_g*(f4+var1)/(self.cell.mu_p+self.cell.mu_g)
            var3 = self.cell.mu_g+self.cell.mu_m-self.cell.mu_g**2./(self.cell.mu_p+self.cell.mu_g)
            var4 = f2+self.cell.mu_g*self.cell_list[q].T2+self.cell.mu_m*(f3+var2)/var3
            var5 = self.cell.mu_m + self.cell.mu_g - self.cell.mu_m**2./var3
            self.cell_list[q].T1 = var4/var5
            #print(self.cell_list[q].T1,'T1')

    def calc_T2(self):
        for q, item in enumerate(self.cell_list):
            f1 = self.h_vap * self.cell_list[q].cathode.gamma
            f2 = (self.cell.vtn - self.cell_list[q].v - self.cell_list[q].omega / 2. * self.cell_list[q].i) * \
                 self.cell_list[q].i
            f3 = 0.5 * self.cell_list[q].omega * self.cell_list[q].i ** 2.
            f4 = self.h_vap * self.cell_list[q].anode.gamma
            if q >= 1:
                f5 = -self.cell.mu_p * self.cell_list[q - 1].T5 - self.a * self.cell_list[q].T
            else:
                f5 = -self.cell.mu_p * self.cell_list[q].T5 - self.a * self.cell_list[q].T
            var1 = (-f5*self.cell.mu_p)/(2.*self.cell.mu_p + self.a)#k
            var2 = self.cell.mu_g*(f4+var1)/(self.cell.mu_p+self.cell.mu_g)#k
            var3 = self.cell.mu_g + self.cell.mu_m-self.cell.mu_g**2./(self.cell.mu_p+self.cell.mu_g)#k
            var4 = f2+self.cell.mu_m*(f3+var2)/var3#k
            var5 = self.cell.mu_g+self.cell.mu_m-self.cell.mu_m**2./var3#k
            var6 = f1 + self.cell.mu_p*self.cell_list[q].T3 + self.cell.mu_g*var4/var5#k
            var7 = self.cell.mu_p**2./(2.*self.cell.mu_p + self.a)#k
            var8 = self.cell.mu_g**2*var7/(self.cell.mu_g+self.cell.mu_p)#k
            var9 = self.cell.mu_m*var8/var3#k
            var10 = self.cell.mu_p+self.cell.mu_g-(self.cell.mu_g**2+var9)/var5
            self.cell_list[q].T2 = var6/var10
            #print(self.cell_list[q].T2,'T2')

    def calc_T3(self):
        for q, item in enumerate(self.cell_list):
            #print(self.cell_list[q].cathode.gamma,'gammac')
            f1 = self.h_vap * self.cell_list[q].cathode.gamma  # f1
            f2 = (self.cell.vtn - self.cell_list[q].v - self.cell_list[q].omega/2.*self.cell_list[q].i)*self.cell_list[q].i
            f3 = 0.5 * self.cell_list[q].omega * self.cell_list[q].i ** 2.
            f4 = self.h_vap * self.cell_list[q].anode.gamma
            if q >= 1:
                f5 = -self.cell.mu_p*self.cell_list[q-1].T5-self.a*self.cell_list[q].T
            else:
                f5 = -self.cell.mu_p*self.cell_list[q].T5-self.a*self.cell_list[q].T
            var1 = self.cell.mu_p*f5/(2.*self.cell.mu_p + self.a)
            var2 = self.cell.mu_g*(f4-var1)/(self.cell.mu_p+self.cell.mu_g)
            var3 = self.cell.mu_g + self.cell.mu_m - self.cell.mu_g**2/(self.cell.mu_p+self.cell.mu_g)
            var4 = self.cell.mu_m*(f3+var2)/var3
            var5 = self.cell.mu_g + self.cell.mu_m - self.cell.mu_m**2./(self.cell.mu_g+self.cell.mu_m-self.cell.mu_g**2./(self.cell.mu_p+self.cell.mu_g))
            var6 = self.cell.mu_g*(f2+var4)/var5
            var7 = self.cell.mu_p**2./(2.*self.cell.mu_p+self.a)
            var8 = self.cell.mu_g*var7/(self.cell.mu_p+self.cell.mu_g)
            var9 = self.cell.mu_m*var8/var3
            var10 = (self.cell.mu_g**2. + var9*self.cell.mu_g)/(var5)
            var11 = f5-self.cell.mu_p*(f1+var6)/(self.cell.mu_p + self.cell.mu_g -var10)
            var12 = -2.*self.cell.mu_p - self.a + self.cell.mu_p**2/(-var10+self.cell.mu_p+self.cell.mu_g)
            self.cell_list[q].T3 = var11/var12
            #print(self.cell_list[q].T3,'T3')
    def calc_T4(self):
        for q, item in enumerate(self.cell_list):
            f3 = 0.5 * self.cell_list[q].omega * self.cell_list[q].i ** 2.
            f4 = self.h_vap * self.cell_list[q].anode.gamma
            if q >= 1:
                f5 = -self.cell.mu_p * self.cell_list[q - 1].T5 - self.a * self.cell_list[q].T
            else:
                f5 = -self.cell.mu_p * self.cell_list[q].T5 - self.a * self.cell_list[q].T
            var1 = (self.cell.mu_p**2.*self.cell_list[q].T2-self.cell.mu_p*f5)/(2.*self.cell.mu_p+self.a)
            var2 = self.cell.mu_g*(f4 + var1)/(self.cell.mu_p + self.cell.mu_g)
            var3 = self.cell.mu_g + self.cell.mu_m - self.cell.mu_g**2./(self.cell.mu_p + self.cell.mu_g)
            self.cell_list[q].T4 = (f3 + self.cell.mu_m*self.cell_list[q].T1 + var2)/var3
            #print(self.cell_list[q].T4,'T4')
    def calc_T5(self):
        for q, item in enumerate(self.cell_list):
            f4 = self.h_vap * self.cell_list[q].anode.gamma
            if q >= 1:
                f5 = -self.cell.mu_p * self.cell_list[q - 1].T5 - self.a * self.cell_list[q].T
            else:
                f5 = -self.cell.mu_p * self.cell_list[q].T5 - self.a * self.cell_list[q].T
            var1 = self.cell.mu_p*(self.cell.mu_p*self.cell_list[q].T2-f5) / (2. * self.cell.mu_p + self.a)
            self.cell_list[q].T5 =  (f4 + var1 + self.cell.mu_g*self.cell_list[q].T4) / (self.cell.mu_p + self.cell.mu_g)
            #print(self.cell_list[q].T5,'T5')

    def calc_layer_T(self):
        #for q, item in enumerate(self.cell_list):
            #self.cell_list[q].T1 = full((self.cell.cathode.channel.division+1), 300.)
            #self.cell_list[q].T2 = full((self.cell.cathode.channel.division+1), 300.)
            #self.cell_list[q].T3 = full((self.cell.cathode.channel.division+1), 300.)
            #self.cell_list[q].T4 = full((self.cell.cathode.channel.division+1), 300.)
            #self.cell_list[q].T5 = full((self.cell.cathode.channel.division+1), 300.)
        self.calc_T3()
        self.calc_T2()
        self.calc_T1()
        self.calc_T4()
        self.calc_T5()

# 4.Layer

class Simulation():
    
    def __init__(self,stack,k_it,max_it):
        self.stack = stack
        self.k_it = k_it
        self.max_it = max_it
        self.b = matrix_database.b(self.stack.cell.cathode.channel.division+1,self.stack.cell_numb,self.stack.cell.cathode.channel.length/(self.stack.cell.cathode.channel.division+1.))
        self.c = matrix_database.c(self.stack.cell.cathode.channel.division+1,self.stack.cell_numb,self.stack.cell.cathode.channel.length/(self.stack.cell.cathode.channel.division+1.))
    
    def update(self):
        for i in range(10):
            self.stack.update()
            self.calc_initial_current_density()
        #for i in range(2):
         #   self.stack.update()
          #  self.calc_sensitivity()
           # self.calc_g()
           # self.calc_n()
            #self.calc_delta()
            #self.calc_i()

    def calc_initial_current_density_fsolve(self,x,i,w,v):
            a = self.stack.cell.e0 - x*self.stack.cell_list[i].omega[w]
            b = R*self.stack.cell_list[i].T1[w]/(Fday)
            c = log(x*self.stack.cell.c_ref/(self.stack.cell.i_ref*(self.stack.cell_list[i].cathode.c1[w]-self.stack.cell.delta*x)))
            return a-b*c-v

    def calc_initial_current_density(self):
        a =[]
        for i in range (self.stack.cell_numb):
            avv  = sum(self.stack.cell_list[i].v)/(self.stack.cell.cathode.channel.division+1)
            b = []
            for q in range (self.stack.cell.cathode.channel.division+1):
                #avv = self.stack.cell_list[i].v[q]
                self.stack.i[i,q] = fsolve(self.calc_initial_current_density_fsolve,I_t,args=(i,q,avv))



    def calc_sensitivity(self):
        self.s = diag(self.stack.dv)
        
    def calc_g(self):
        self.g = matmul(self.b,self.s) - self.stack.resistance*self.c
    
    def calc_n(self):
        var0 = I_t*(self.stack.cell.cathode.channel.division+1)
        self.stack.i[0,0] = var0 - (sum(self.stack.i[0])-self.stack.i[0,0])
        var1 = matmul(self.b,self.stack.v)
        var2 = self.stack.resistance*matmul(self.c,self.stack.i.flatten(order='C'))
        self.n = var1-var2
    
    def calc_delta(self):
        var1 = linalg.pinv(self.g)
        self.delta = matmul(var1,self.n)
        
    def calc_i(self):
        #print(len(self.stack.i.flatten(order='F')),len(self.delta))
        var1 = self.stack.i.flatten(order='C')-self.delta
        #print(var1)
        self.stack.i = toarray(var1)
        #print(self.stack.i)


channel_anode = Channel(0.67,N-1,20.*10.**3.,3.*10.**5.,350.,False)
anode = Halfcell(channel_anode,1.8,2,2.)
channel_cathode = Channel(0.67,N-1,20.*10.**3.,3.2*10.**5.,350.,True)
cathode = Halfcell(channel_cathode,1.8,3,4.)
cell = Cell(anode,cathode,0.62*10.**-5.,1200.,50.*10**-6,2300.,5000.,11200.,0.944,40.9,64.,0.8*10.**-3.,340.,pem_type)
stack = Stack(cell,M)
simulation = Simulation(stack,1.*10.**-3.,10.)
simulation.update()

x = linspace(0.,0.67,N)
xe = linspace(0.,0.67,N-1)

print('current density')
for i in simulation.stack.cell_list:
    #print(i.i)
    plot(x,i.i)
show()
print('voltage')

for i in simulation.stack.cell_list:
    #print(i.v)
    plot(x,i.v)
ylim(0.,1.28)
show()
 
    
print('dvdi')
for i in simulation.stack.cell_list:
    #print(i.dv)
    plot(x,i.dv)
show()

print('psi')
for i in simulation.stack.cell_list:
    #print(i.psi)
    plot(x,i.psi)
show()

print('omega')
for i in simulation.stack.cell_list:
    #print(i.omega)
    plot(x,i.omega)
show()

print('Coolant Temperature')
for i in simulation.stack.cell_list:
    #print(i.T)
    plot(x,i.T)
show()

print('Membran Temperature')
for i in simulation.stack.cell_list:
    #print(i.T0)
    plot(x,i.T0)
show()

print('Cathode Temperature')
for i in simulation.stack.cell_list:
    #print(i.T1)
    plot(x,i.T1)
show()

print('Cathode Gas Channel Temperature')
for i in simulation.stack.cell_list:
    #print(i.T2)
    plot(x,i.T2)
show()

print('Coolant Plate')
for i in simulation.stack.cell_list:
    #print(i.T3)
    plot(x,i.T3)
show()

print('Anode Temperature')
for i in simulation.stack.cell_list:
    #print(i.T4)
    plot(x,i.T4)

show()

print('Anode Gas Channel Temperature')
for i in simulation.stack.cell_list:
    #print(i.T5)
    plot(x,i.T5)
show()

print('pc')
for i in simulation.stack.cell_list:
    #print(i.cathode.p)
    plot(x,i.cathode.p)
#ylim(3*10**5,3.2*10**5)
show()

print('pa')
for i in simulation.stack.cell_list:
    #print(i.anode.p)
    plot(x,i.anode.p)
show()

print('q0')
for i in simulation.stack.cell_list:
    #print(i.cathode.q1)
    plot(x,i.cathode.q1)
show()

print('qc')
for i in simulation.stack.cell_list:
    #print(i.cathode.q2)
    plot(x,i.cathode.q2)
show()

print('qn')
for i in simulation.stack.cell_list:
    #print(i.cathode.q3)
    plot(x,i.cathode.q3)
show()

print('qh')
for i in simulation.stack.cell_list:
    #print(i.anode.q1)
    plot(x,i.anode.q1)
show()

print('qa')
for i in simulation.stack.cell_list:
    #print(i.anode.q2)
    plot(x,i.anode.q2)
show()

print('co')
for i in simulation.stack.cell_list:
    #print(i.cathode.c1)
    plot(x,i.cathode.c1)
show()

print('cc')
for i in simulation.stack.cell_list:
    #print(i.cathode.c2)
    plot(x,i.cathode.c2)
show()

print('ch')
for i in simulation.stack.cell_list:
    #print(i.anode.c1)
    plot(x,i.anode.c1)
show()

print('ca')
for i in simulation.stack.cell_list:
    #print(i.anode.c2)
    plot(x,i.anode.c2)
show()

print('wc')
for i in simulation.stack.cell_list:
    # print(i.cathode.w)
    plot(x, i.cathode.w)
show()

print('gamma_c')
for i in simulation.stack.cell_list:
    #print(i.cathode.gamma)
    plot(x,i.cathode.gamma)
show()

if pem_type is False:
    print('wc')
    for i in simulation.stack.cell_list:
        #print(i.cathode.w)
        plot(x,i.cathode.w)
    show()

    print('wa')
    for i in simulation.stack.cell_list:
        #print(i.anode.w)
        plot(x,i.anode.w)
    show()

    print('gamma_c')
    for i in simulation.stack.cell_list:
        #print(i.cathode.gamma)
        plot(x,i.cathode.gamma)
    show()

    print('gamma_a')
    for i in simulation.stack.cell_list:
        #print(i.anode.gamma)
        plot(x,i.anode.gamma)
    show()

    print('humidity_c')
    for i in simulation.stack.cell_list:
        #print(i.cathode.humidity)
        plot(x,i.cathode.humidity)
    show()

    print('humidity_a')
    for i in simulation.stack.cell_list:
        #print(i.anode.humidity)
        plot(x,i.anode.humidity)
    show()

    print('free_water_c')
    for i in simulation.stack.cell_list:
        #print(i.cathode.free_water)
        plot(x,i.cathode.free_water)
    show()

    print('free_water_a')
    for i in simulation.stack.cell_list:
        #print(i.anode.free_water)
        plot(x,anode.free_water)
    show()

