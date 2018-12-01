import input.geometry as geo
import input.physical_property as phy_prop
import input.operating_conditions as oper_con
cell = {'th_mem': geo.membrane_thickness,
        'plate_h_con': phy_prop.thermal_conductivity_bipolar_plate_z,
        'gde_h_con': phy_prop.thermal_conductivity_gas_diffusion_electrode_z, 'mem_h_con': phy_prop.thermal_conductivity_membrane_z,
        't_cool_in': oper_con.temp_coolant_in, 'plate_hi_con': phy_prop.thermal_conductivity_bipolar_plate_x,
        'gde_hi_con': phy_prop.thermal_conductivity_gas_diffusion_electrode_x, 'mem_hi_con': phy_prop.thermal_conductivity_membrane_x,
        'mem_bas_r': phy_prop.membrane_basic_resistance,
        'mem_acl_r': phy_prop.membrane_temperature_resistance,
        't_init': oper_con.temp_initial}
