import input.geometry as geo
import input.physical_property as phy_prop
import input.operating_conditions as oper_con
import input.simulation as sim
# still in search for shortings to replace mem_bas_r and mem_acl_r
dict_cell = {
    'th_mem': geo.membrane_thickness,
    'lambda_z_bpp': phy_prop.thermal_conductivity_bipolar_plate_z,
    'lambda_z_gde': phy_prop.thermal_conductivity_gas_diffusion_electrode_z,
    'lambda_z_mem': phy_prop.thermal_conductivity_membrane_z,
    'lambda_x_bpp': phy_prop.thermal_conductivity_bipolar_plate_x,
    'lambda_x_gde': phy_prop.thermal_conductivity_gas_diffusion_electrode_x,
    'lambda_x_mem': phy_prop.thermal_conductivity_membrane_x,
    'temp_cool_in': oper_con.temp_coolant_in,
    'mem_base_r': phy_prop.membrane_basic_resistance,
    'mem_acl_r': phy_prop.membrane_temperature_resistance,
    'temp_init': oper_con.temp_initial,
    'calc_mem_loss': sim.calc_membrane_loss
            }
