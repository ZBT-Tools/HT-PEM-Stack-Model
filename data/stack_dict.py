import input.operating_conditions as oper_con
import input.geometry as geo
import input.physical_property as phy_prop

stack = {'cell_numb': oper_con.cell_number,
         'heat_power': oper_con.endplates_heat_power,
         'height': geo.manifold_height, 'width': geo.manifold_width,
         'dis_dis_fac': phy_prop.manifold_pressure_loss_coefficient,
         'stoi_cat': oper_con.stoichiometry_cathode,
         'stoi_ano': oper_con.stoichiometry_anode,
         'cool_ch_bc': oper_con.cooling_bc,
         'h_cool': geo.bipolar_plate_thickness * .5,
         'm_flow_cool': oper_con.mass_flow_coolant,
         'cp_cool': phy_prop.heat_capacity_coolant,
         'alpha_cool': phy_prop.convection_coefficient_coolant_channel}