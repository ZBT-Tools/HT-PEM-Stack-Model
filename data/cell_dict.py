import input.geometry as geo
import input.physical_property as phy_prop
import input.operating_conditions as oper_con
cell = {'th_mem': geo.membrane_thickness, 'plate_h_con': phy_prop.k_p,
        'gde_h_con': phy_prop.k_g, 'mem_h_con': phy_prop.k_m,
        't_cool_in': oper_con.t_coolant_in,'plate_hi_con': phy_prop.ki_p,
        'gde_hi_con': phy_prop.ki_g, 'mem_hi_con': phy_prop.ki_m,
        'mem_bas_r': phy_prop.membrane_basic_resistance,
        'mem_acl_r': phy_prop.membrane_temperature_resistance,
        't_init': oper_con.t_initial}
