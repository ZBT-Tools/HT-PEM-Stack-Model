################################################################################
############################Channel#############################################
# channellength m
channel_length = 6.7
# number of Nodes
nodes = 50
# cathode gas channel scaled pressure drop coefficient Pa s /m³
k_cat = 2.e4
# cathode gas channel inlet pressure Pa
p_cat_in = 3.2e5
# cathode gas channel inlet temperature K
t_cat_in = 2.98e2
# cathode phi
phi_cat = 0.5
# cathode flow direction
flow_dir_cat = True
# anode gas channel scaled pressure drop coefficient Pa s /m³
k_ano = 2.e4
# anode gas channel inlet pressure Pa
p_ano_in = 3.e5
# anode gas channel inlet temperature K
t_ano_in = 2.98e2
# ano phi
phi_ano = 0.5
# anode flow direction
flow_dir_ano = False
################################################################################
############################Halfcell############################################
# cathode stoichiometry
stoi_cat = 1.5
# cathode number of species (3 (N2,O2,H2O))
spec_numb_cat = 3
# cathode side reaction, number of valence electrones
val_cat = 4
# cathode stoichiometry
stoi_ano = 1.5
# anode number of species (2 (H2,H2O))
spec_numb_ano = 2
#  side reaction, number of valence electrones
val_ano = 2
# target current density A/m²
tar_cd = 9000.
################################################################################
############################Cell################################################
# fited vapour mass transport coefficient m/s
gamma = 0.62e-5
# molar concentration of membrane acid groups moles/m³
alpha = 1.2e3
# thermal conductivity of the bipolar plate (z adjusted) W/(m²K)
mu_p = 2.3e3
# thermal conductivity of the gde (z adjusted) W/(m²K)
mu_g = 5.e3
# thermal conductivity of the membran (z adjusted) W/(m²K)
mu_m = 1.1e4
# fitted open circuit voltage V
e_o = 0.944
# reference oxygen concentration moles/m³
c_ref = 40.9
# fitted oxygen reduction reaction exchange current density A/m²
i_ref = 64.
# fitted oxygen mass transport coefficient mole/(m*C)
delta = 0.8993e-3
# pemtype (True = HTPEM, False = NTPEM)
pem_type = True
# coolant channel (True = coolant channel, False =  no coolant channel)
coolant_channel = False
# inlet temperature of the coolant K
t_cool_in = 0
# thickness of the membran m
mem_thick = 5.e-3
# thickness of the gde m
gde_thick = 1.e-3
# thickness of the plate m
plate_thick = 0.5e-3
# widht of the channel m
channel_width = 1.e-3
# environment temperature K
t_u = 298.15
################################################################################
############################Stack###############################################
# number of stack cells (min =!3)
cell_numb = 10
# heated endplates ?( True = y or False = no)
endplate_adiabat = False
# Endplate_temperature K
end_plate_t = 310.
# convection coefficient of the stack walls W/(m²K)
alpha_conv = 0.
#resistance of the bipolarplates
resistance = 4.e-3
################################################################################
############################Simulation##########################################
#number of initial parameter calculation loops
n_loop = 20
# tolerance
k_tol = 1.e-3
# max number of iterations
max_it = 1
