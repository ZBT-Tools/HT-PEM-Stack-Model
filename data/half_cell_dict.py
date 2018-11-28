import input.physical_property as phy_prop
import input.geometry as geo
import numpy as np

cathode = {'type': True, 'th_gdl': geo.gdl_thickness,
           'th_plate': geo.plate_thickness,
           'tafel_slope': phy_prop.cathode_tafel_slope,
           'prot_con': phy_prop.cathode_proton_conductivity,
           'vol_ex_cd': phy_prop.cathode_exchange_current_density,
           'diff_coef_cat':phy_prop.cathode_layer_diffusion_coefficient,
           'diff_coef_gdl': phy_prop.cathode_gdl_diffusion_coefficient,
           'th_electrode': geo.electrode_thickness}
anode = {'type': False,
         'th_gdl': geo.gdl_thickness, 'th_plate': geo.plate_thickness,
         'tafel_slope': phy_prop.anode_tafel_slope,
         'prot_con': phy_prop.anode_proton_conductivity,
         'vol_ex_cd': phy_prop.anode_exchange_current_density,
         'diff_coef_cat':phy_prop.anode_layer_diffusion_coefficient,
         'diff_coef_gdl': phy_prop.anode_gdl_diffusion_coefficient,
         'th_electrode': geo.electrode_thickness}
