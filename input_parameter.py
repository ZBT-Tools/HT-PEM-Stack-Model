from numpy import array, linspace
################################################################################
############################Channel#############################################
# channel length [m]
channel_length = 0.2322
# number of Nodes
nodes = 5
# cathode gas channel inlet pressure [Pa]
p_cat_in = 1.e5
# cathode gas channel inlet temperature [K]
t_cat_in = 433.15
# cathode phi
phi_cat = 0.
# cathode flow direction
flow_dir_cat = True
# anode gas channel inlet pressure [Pa]
p_ano_in = 1.e5
# anode gas channel inlet temperature [K]
t_ano_in = t_cat_in
# ano phi
phi_ano = 0.
# anode flow direction
flow_dir_ano = False
# width of the channel [m]
channel_width = 1.e-3
# height of the channel [m]
channel_height = 1.e-3
# number of channel bends
channel_bends = 48
#bend friction factor [0-1]
channel_fri_fac = 0.1
################################################################################
############################Halfcell############################################
# oxygen concentration at the inlet [mol/m³]
o2_con_in = 0.21
# cathode stoichiometry
stoi_cat = 2.5
# minimum stoichiometry (do not set it smaller 1.1)
stoi_min = 1.1
# cathode number of species (3 [N2,O2,H2O])
spec_numb_cat = 3
# cathode side reaction, number of valence electrons
val_cat = 4.
# anode stoichiometry
stoi_ano = 1.34
# anode number of species (2 (H2,H2O))
spec_numb_ano = 2
# side reaction, number of valence electrons
val_ano = 2.
# target current density [A/m^2]
#tar_cd = linspace(1.e-3, 10., 10)
#tar_cd = array([1111.11,2222.22, 3333.33, 4444.44, 5555.55, 6666.66,8888.88,9999.99,11111.111,12222.22])
#tar_cd = array([1111.11,2222.22, 3333.33, 4444.44, 5555.55, 6666.66,7000.,8000.,8500.,9000.])
tar_cd = array([6000.])
# gas reference concentration cathode [mol/m^3]
gas_con_ref_cat = 4.e0
# gas reference concentration anode [mol/m^3}
gas_con_ref_ano = 3.68e1
# symmetry factor cathode
sym_cat = 1.
# symmetry factor anode
sym_ano = 1.
# molar mass oxygen [g/mol]
m_o2 = 32.
# molar mass hydrogen [g/mol]
m_h2 = 2.
# molar mass water [g/mol]
m_h2o = 18.
# molar mass nitrogen [g/mol]
m_n2 = 28.
# cell init temperature [K]
t_hcell_init = t_cat_in
# cathode layer proton conductivity [Ohm^-1m^-1]
cat_prot_con = 3.e0
# anode layer proton conductivity [Ohm^-1m^-1]
ano_prot_con = 3.e0
# cathode volumetric exchange current density [A/m^3]
cat_vol_ex_cd = 2.3e3
# anode volumetric exchange current density [A/m^3]
ano_vol_ex_cd = 0.817e9
# cathode layer diffusion coefficient [m^2/s]
cat_layer_dif_coef = 1.36e-8
# anode layer diffusion coefficient [m^2/s]
ano_layer_dif_coef = 1.36e-8
# cathode gdl diffusion coefficient [m^2/s]
cat_gdl_dif_coef = 2.59e-6
# anode gdl diffusion coefficient [m^2/s]
ano_gdl_dif_coef = 2.59e-6
# cathode layer thickness [m]
cat_layer_thick = 1.e-5
# anode layer thickness [m]
ano_layer_thick = 1.e-5
# cathode tafel slope [V]
cat_tafel_slope = 0.03
# anode tafel slope [V]
ano_tafel_slope = 0.03
################################################################################
############################Cell################################################
# membrane basic res
r_mem0 = 0.33
#membran temp res acl
r_memm = 0.0007
# cell width [m]
cell_width = 0.01
# cell length [m]
cell_length = 1.e-1
# fitted vapour mass transport coefficient [m/s]
gamma = 0.62e-5
# molar concentration of membrane acid groups [mol/m³]
alpha = 1.2e3
# thermal conductivity of the bipolar plate through plane  [W/(mK)]
k_p = 1.e2
# thermal conductivity of the bipolar plate in plane  [W/(mK)]
ki_p = 1.e2
# thermal conductivity of the gde through plane  [W/(mK)]
k_g = 1.e0
# thermal conductivity of the gde in plane  [W/(mK)]
ki_g = 1.e0
# thermal conductivity of the membrane through plane  [W/(mK)]
k_m = .26e0
# thermal conductivity of the membrane in plane  [W/(mK)]
ki_m = .26e0
# fitted open circuit voltage [V]
e_o = 1.
# thermoneutral voltage for ORR [V]
vtn = 1.28
# pem-type (True = HT-PEM, False = NT-PEM)
pem_type = True
# inlet temperature of the coolant [K]
t_cool_in = t_cat_in
# thickness of the membrane [m]
mem_thick = 50.e-6
# thickness of the gde [m]
gde_thick = 250.e-6
# thickness of the plate [m]
plate_thick = 5.e-3
# environment temperature [K]
t_u = 298.15
# height coolant channel [m]
h_col_ch = plate_thick/ 2.
# cp value coolant [J/(kgK)]
cp_col = 4000.
# mass flow coolant [kg/s]
m_col = 1.e-6
# convection coefficient coolant [W/(m^2K)]
a_col = 4000.
################################################################################
############################Stack###############################################
# number of stack cells (min =!2)
cell_numb = 3
# boundary condition cooling channel (no cooling channel at the endplates = False)
cooling_bc = True
# endplate heating power [W]
heat_power = 0.#5.e0
# convection coefficient of the stack walls [W/(m^2K)]
alpha_conv = 0.5e-1
# resistivity of the bipolarplates [Ohm/m]
resistivity = 2.e-6
# height of the manifold channel [m]
manifold_height = 2.01e-2
# width of the manifold channel [m]
manifold_width = 2.01e-2
# geometrical pressure loss coefficient
manifold_kf = 1.e6
# latent heat of vaporization [J/mol]
h_vap = 45400.
################################################################################
############################Simulation##########################################
# tolerance
k_tol = 1.e-12
# max number of iterations
max_it = 100
# gas channel relaxations factor
channel_fac = 0.5
#############################Dicts##############################################
################################################################################
channel_cat = {'length': channel_length, 'p_in': p_cat_in, 't_in': t_cat_in,
               'hum_in': phi_cat, 'flow_dir': flow_dir_cat, 'width': channel_width,
               'heigth': channel_height, 'numb_bends': channel_bends,
               'bend_fri_fac': channel_fri_fac}
channel_ano = {'length': channel_length, 'p_in': p_ano_in, 't_in': t_ano_in,
               'hum_in': phi_ano, 'flow_dir': flow_dir_ano, 'width': channel_width,
               'heigth': channel_height, 'numb_bends': channel_bends,
               'bend_fri_fac': channel_fri_fac}
cathode = {'spec_numb': spec_numb_cat, 'val': val_cat, 'type': True,
           'gas_con_ref': gas_con_ref_cat,
           'sym_fac': sym_cat, 'thick_gde': gde_thick,
           'plate_thick': plate_thick, 'm_mass': array([m_o2, m_h2o, m_n2]),
           't_init': t_hcell_init, 'tafel_slope': cat_tafel_slope,
           'prot_con': cat_prot_con, 'vol_ex_cd': cat_vol_ex_cd,
           'dif_coef_cat':cat_layer_dif_coef, 'dif_coef_gdl': cat_gdl_dif_coef,
           'cat_thick':cat_layer_thick}
anode = {'spec_numb': spec_numb_ano, 'val': val_ano, 'type': False,
         'gas_con_ref': gas_con_ref_ano,
         'sym_fac': sym_ano, 'thick_gde': gde_thick,
         'plate_thick': plate_thick, 'm_mass': array([m_h2, m_h2o]),
         't_init': t_hcell_init, 'tafel_slope': ano_tafel_slope,
         'prot_con':ano_prot_con, 'vol_ex_cd': ano_vol_ex_cd,
         'dif_coef_cat':ano_layer_dif_coef, 'dif_coef_gdl': ano_gdl_dif_coef,
         'cat_thick':ano_layer_thick}
cell = {'mem_thick': mem_thick, 'plate_h_con': k_p, 'gde_h_con': k_g,
        'mem_h_con': k_m,'t_coolant_in': t_cool_in
        ,'width': cell_width,'length': cell_length, 'plate_hi_con': ki_p,
        'gde_hi_con': ki_g, 'mem_hi_con': ki_m, 'mem_bas_r': r_mem0,
        'mem_acl_r': r_memm}
stack = {'cell_numb': cell_numb, 'heat_power': heat_power,
         'heigth': manifold_height, 'width': manifold_width,
         'dis_dis_fac': manifold_kf,
         'stoi_cat': stoi_cat, 'stoi_ano': stoi_ano,'cool_ch_bc': cooling_bc,
         'h_col': h_col_ch, 'm_flow_col': m_col, 'cp_col': cp_col, 'alpha_cool': a_col}
simulation = {'max_it': max_it, 'k_it': k_tol}