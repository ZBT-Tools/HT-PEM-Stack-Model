import input.physical_properties as phy_prop
import input.geometry as geom
import input.simulation as sim
import input.operating_conditions as op_con


dict_cathode = {
 'species': ['O2', 'H2O', 'N2'],
 'cell_width': geom.cell_width,
 'cell_length': geom.cell_length,
 'channel_numb': op_con.gas_channel_number,
 'cl_type': True,
 'th_gdl': geom.gas_diffusion_layer_thickness,
 'th_bpp': geom.bipolar_plate_thickness,
 'tafel_slope': phy_prop.tafel_slope_cathode,
 'prot_con_cl': phy_prop.catalyst_layer_proton_conductivity_cathode,
 'vol_ex_cd': phy_prop.exchange_current_density_cathode,
 'diff_coeff_cl': phy_prop.oxygen_catalyst_layer_diffusion_coefficient,
 'diff_coeff_gdl': phy_prop.oygen_gas_diffusion_layer_diffusion_coefficient,
 'th_cl': geom.catalyst_layer_thickness
                }

dict_anode = {
    'species': ['H2', 'H2O', 'N2'],
    'cell_width': geom.cell_width,
    'cell_length': geom.cell_length,
    'channel_numb': op_con.gas_channel_number,
    'cl_type': False,
    'th_gdl': geom.gas_diffusion_layer_thickness,
    'th_bpp': geom.bipolar_plate_thickness,
    'tafel_slope': phy_prop.tafel_slope_anode,
    'prot_con_cl': phy_prop.catalyst_layer_proton_conductivity_anode,
    'vol_ex_cd': phy_prop.exchange_current_density_anode,
    'diff_coeff_cl': phy_prop.hydrogen_catalyst_layer_diffusion_coefficient,
    'diff_coeff_gdl': phy_prop.hydrogen_diffusion_layer_diffusion_coefficient,
    'th_cl': geom.catalyst_layer_thickness
            }
