import input.physical_property as phy_prop
import input.geometry as geom
import input.simulation as sim


dict_cathode = {
 'cl_type': True,
 'th_gdl': geom.gas_diffusion_layer_thickness,
 'th_bpp': geom.bipolar_plate_thickness,
 'tafel_slope': phy_prop.tafel_slope_cathode,
 'prot_con_cl': phy_prop.catalyst_layer_proton_conductivity_cathode,
 'vol_ex_cd': phy_prop.exchange_current_density_cathode,
 'diff_coeff_cl': phy_prop.oxygen_catalyst_layer_diffusion_coefficient,
 'diff_coeff_gdl': phy_prop.oygen_gas_diffusion_layer_diffusion_coefficient,
 'th_cl': geom.catalyst_layer_thickness,
 'calc_act_loss': sim.calc_activation_loss,
 'calc_cl_diff_loss': sim.calc_cl_loss,
 'calc_gdl_diff_loss': sim.calc_gdl_loss
                }

dict_anode = {
    'cl_type': False,
    'th_gdl': geom.gas_diffusion_layer_thickness,
    'th_bpp': geom.bipolar_plate_thickness,
    'tafel_slope': phy_prop.tafel_slope_anode,
    'prot_con_cl': phy_prop.catalyst_layer_proton_conductivity_anode,
    'vol_ex_cd': phy_prop.exchange_current_density_anode,
    'diff_coeff_cl': phy_prop.hydrogen_catalyst_layer_diffusion_coefficient,
    'diff_coeff_gdl': phy_prop.hydrogen_diffusion_layer_diffusion_coefficient,
    'th_cl': geom.catalyst_layer_thickness,
    'calc_act_loss': sim.calc_activation_loss,
    'calc_cl_diff_loss': sim.calc_cl_loss,
    'calc_gdl_diff_loss': sim.calc_gdl_loss
            }
