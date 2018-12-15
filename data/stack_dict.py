import input.operating_conditions as oper_con
import input.geometry as geom
import input.physical_property as phy_prop

dict_stack = {
    'cell_numb': oper_con.cell_number,
    'heat_power': oper_con.endplates_heat_power,
    'header_height': geom.manifold_height,
    'header_width': geom.manifold_width,
    'dis_dis_fac': phy_prop.manifold_pressure_loss_coefficient,
    'stoi_cat': oper_con.stoichiometry_cathode,
    'stoi_ano': oper_con.stoichiometry_anode,
    'cool_ch_bc': oper_con.cooling_bc,
    'alpha_env': phy_prop.convection_coefficient_stack_environment
             }