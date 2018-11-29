import input.physical_property as phy_prop
import input.geometry as geo
import numpy as np

cathode = {'type': True, 'th_gdl': geo.gas_diffusion_layer_thickness,
           'th_plate': geo.bipolar_plate_thickness,
           'tafel_slope': phy_prop.tafel_slope_cathode,
           'prot_con': phy_prop.catalyst_layer_proton_conductivity_cathode,
           'vol_ex_cd': phy_prop.exchange_current_density_cathode,
           'diff_coef_cat':phy_prop.oxygen_catalyst_layer_diffusion_coefficient,
           'diff_coef_gdl': phy_prop.oygen_gas_diffusion_layer_diffusion_coefficient,
           'th_electrode': geo.catalyst_layer_thickness}
anode = {'type': False,
         'th_gdl': geo.gas_diffusion_layer_thickness, 'th_plate': geo.bipolar_plate_thickness,
         'tafel_slope': phy_prop.tafel_slope_anode,
         'prot_con': phy_prop.catalyst_layer_proton_conductivity_anode,
         'vol_ex_cd': phy_prop.exchange_current_density_anode,
         'diff_coef_cat':phy_prop.hydrogen_catalyst_layer_diffusion_coefficient,
         'diff_coef_gdl': phy_prop.hydrogen_diffusion_layer_diffusion_coefficient,
         'th_electrode': geo.catalyst_layer_thickness}
