from numpy import array, linspace
################################################################################
############################Channel#############################################

# channellength m
channel_length = 0.67
# number of Nodes
nodes = 10
# cathode gas channel scaled pressure drop coefficient Pa s /mÂ³
k_cat = 2.e4
# cathode gas channel inlet pressure Pa
p_cat_in = 3.2e5
# cathode gas channel inlet temperature K
t_cat_in = 350.15
# cathode phi
phi_cat = 0.
# cathode flow direction
flow_dir_cat = True
# anode gas channel scaled pressure drop coefficient Pa s /mÂ³
k_ano = 2.e4
# anode gas channel inlet pressure Pa
p_ano_in = 3.e5
# anode gas channel inlet temperature K
t_ano_in = t_cat_in
# ano phi
phi_ano = 0.
# anode flow direction
flow_dir_ano = False
# width of the channel m
channel_width = 1.e-2
# heigth of the channel m
channel_heigth = 1.e-3
################################################################################
############################Halfcell############################################
# oxygen concentration at the inlet
o2_con_in = 0.21
# cathode stoichiometry
stoi_cat = 1.2
# minimum stoichiometry (do not set it smaller 1.1)
stoi_min = 1.1
# cathode number of species (3 (N2,O2,H2O))
spec_numb_cat = 3
# cathode side reaction, number of valence electrones
val_cat = 4.
# anode stoichiometry
stoi_ano = 5.
# anode number of species (2 (H2,H2O))
spec_numb_ano = 2
#  side reaction, number of valence electrones
val_ano = 2.
# target current density A/mÂ²
#tar_cd = linspace(1000.,16000., 10)
tar_cd = [6000.]
# gas reference concentration cathode [mol/mÂ³]
gas_con_ref_cat = 7.36e0
# gas reference concentration anode [mol/mÂ³}
gas_con_ref_ano = 3.68e1
# activation energy cathode
act_energy_cat = 100.e3
# activation energy anode
act_energy_ano = 10.e3
# reference current density cathode
i_ref_cat = 0.817e3
# reference current density anode
i_ref_ano = 5.e5
# symmetry factor cathode
sym_cat = 1.
# symmetry factor anode
sym_ano = 1.
# molar mass oxygen
m_o2 = 32.
# molar mass hydrogen
m_h2 = 2.
# molar mass water
m_h2o = 18.
#molar mass nitrogen
m_n2 = 28.
#hcell init temperature
t_hcell_init = 353.15
#cathode layer proton conductivity Ohm^-1m^-1
cat_prot_con = 3.e0
#anode layer proton conductivity Ohm^-1m^-1
ano_prot_con = 3.e0
#cathode volumetric exchange current density [A/mÂ³]
cat_vol_ex_cd = 0.817e3
#anode volumetric exchange current density [A/mÂ²]
ano_vol_ex_cd = 0.817e9
#cathode layer diffusion coefficient [mÂ²/s]
cat_layer_dif_coef = 1.36e-8
#anode layer diffuision coefficient [mÂ²/s]
ano_layer_dif_coef = 1.36e-8
#cathode gdl diffusion coefficient [mÂ²/s]
cat_gdl_dif_coef = 2.59e-6
#anode gdl diffusion coefficient [mÂ²/s]
ano_gdl_dif_coef = 2.59e-6
#cathode layer thickness [m]
cat_layer_thick = 1.e-5
#anode layer thickness [m]
ano_layer_thick = 1.e-5
#cathode tafel slope [V]
cat_tafel_slope = 0.03
#anode tafel slope [V]
ano_tafel_slope = 0.03
################################################################################
############################Cell################################################
#cell width
cell_width = 0.01
#cell length
cell_length = 1.e-1
# fited vapour mass transport coefficient m/s
gamma = 0.62e-5
# molar concentration of membrane acid groups moles/mÂ³
alpha = 1.2e3
# thermal conductivity of the bipolar plate through plane  W/(mK)
k_p = 1.e2
ki_p = 1.e2
# thermal conductivity of the bipolar plate in plane  W/(mK)
# thermal conductivity of the gde through plane  W/(mK)
k_g = 1.e0
ki_g = 1.e0
# thermal conductivity of the membran through plane  W/(mK)
k_m = .26e0
ki_m = .26e0
# fitted open circuit voltage V
e_o = 1.18
# thermoneutral voltage for ORR V
vtn = 1.28
# pemtype (True = HTPEM, False = NTPEM)
pem_type = True
# inlet temperature of the coolant K
t_cool_in = t_cat_in
# thickness of the membran m
mem_thick = 50.e-6
# thickness of the gde m
gde_thick = 250.e-6
# thickness of the plate m
plate_thick = 5.e-3
# environment temperature K
t_u = 298.15
# diameter coolant channel m
d_col = plate_thick/2.
# cp value coolant J/(kgK)
cp_col = 4000.
# mass flow coolant kg/s
m_col = 1.e-3
# convection coef coolant W/(mÂ²K)
a_col = 4000.
################################################################################
############################Stack###############################################
# number of stack cells (min =!2)
cell_numb = 10
#boundry conditoon cooling channel (no cooling channel at the endplates= False)
cooling_bc = True
# Endplate heating power W
heat_power = 0.#5.e0
# convection coefficient of the stack walls W/(mÂ²K)
alpha_conv = 0.5e-10
# resitivity of the bipolarplates Ohm/m
resistivity = 8.e-4
# heigth of the manifold channel
manifold_heigth = 10.e-3
# width of the manifold channel
manifold_width = 10.e-3
# geometrical pressure loss coefficient
manifold_kf = 0.27
# latent heat of vaporization
h_vap = 45400.
################################################################################
############################Simulation##########################################
# tolerance
k_tol = 1.e-10
# max number of iterations
max_it = 4000
#############################Dicts##############################################
################################################################################
channel_cat = {'length': channel_length, 'p_in': p_cat_in, 't_in': t_cat_in,
               'hum_in': phi_cat, 'flow_dir': flow_dir_cat, 'width': channel_width,
               'heigth': channel_heigth}
channel_ano = {'length': channel_length, 'p_in': p_ano_in, 't_in': t_ano_in,
               'hum_in': phi_ano, 'flow_dir': flow_dir_ano, 'width': channel_width,
               'heigth': channel_heigth}
cathode = {'spec_numb': spec_numb_cat, 'val': val_cat, 'type': True,
           'gas_con_ref': gas_con_ref_cat, 'act_energy': act_energy_cat,
           'sym_fac': sym_cat, 'thick_gde': gde_thick,
           'plate_thick': plate_thick, 'm_mass': array([m_o2, m_h2o, m_n2]),
           't_init': t_hcell_init, 'tafel_slope': cat_tafel_slope,
           'prot_con':cat_prot_con, 'vol_ex_cd': cat_vol_ex_cd,
           'dif_coef_cat':cat_layer_dif_coef, 'dif_coef_gdl': cat_gdl_dif_coef,
           'cat_thick':cat_layer_thick}
anode = {'spec_numb': spec_numb_ano, 'val': val_ano, 'type': False,
         'gas_con_ref': gas_con_ref_ano, 'act_energy': act_energy_ano,
         'sym_fac': sym_ano, 'thick_gde': gde_thick,
         'plate_thick': plate_thick, 'm_mass': array([m_h2, m_h2o]),
         't_init': t_hcell_init, 'tafel_slope': ano_tafel_slope,
         'prot_con':ano_prot_con, 'vol_ex_cd': ano_vol_ex_cd,
         'dif_coef_cat':ano_layer_dif_coef, 'dif_coef_gdl': ano_gdl_dif_coef,
         'cat_thick':ano_layer_thick}
cell = {'mem_thick': mem_thick, 'plate_h_con': k_p, 'gde_h_con': k_g,
        'mem_h_con': k_m,'t_coolant_in': t_cool_in
        ,'width': cell_width,'length': cell_length, 'plate_hi_con': ki_p,
        'gde_hi_con': ki_g, 'mem_hi_con': ki_m}
stack = {'cell_numb': cell_numb, 'heat_power': heat_power,
         'plate_res': resistivity, 'heigth': manifold_heigth,
         'width': manifold_width, 'dis_dis_fac': manifold_kf,
         'stoi_cat': stoi_cat, 'stoi_ano': stoi_ano,'cool_ch_bc':cooling_bc,
         'd_col': d_col, 'm_flow_col': m_col, 'cp_col': cp_col,'alpha_cool':a_col}
simulation = {'max_it': max_it, 'k_it': k_tol}