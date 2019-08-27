import input.physical_properties as phy_prop
import input.geometry as geom
import input.simulation as sim
import input.operating_conditions as op_con


dict_cathode = {
    'name': 'Cathode',
    'flow_direction': geom.cathode_flow_direction,
    'cell_width': geom.cell_width,
    'cell_length': geom.cell_length,
    'channel_numb': geom.gas_channel_number,
    'is_cathode': True,
    'species_names': op_con.cathode_species,
    'molar_mass': op_con.cathode_molar_mass,
    'composition': op_con.cathode_inlet_composition,
    'charge_number': op_con.cathode_charge_number,
    'reaction_stoichiometry': op_con.cathode_reaction_stoich,
    'th_gdl': geom.gas_diffusion_layer_thickness,
    'th_bpp': geom.bipolar_plate_thickness,
    'tafel_slope': phy_prop.tafel_slope_cathode,
    'prot_con_cl': phy_prop.catalyst_layer_proton_conductivity_cathode,
    'vol_ex_cd': phy_prop.exchange_current_density_cathode,
    'diff_coeff_cl': phy_prop.oxygen_catalyst_layer_diffusion_coefficient,
    'diff_coeff_gdl': phy_prop.oxygen_gas_diffusion_layer_diffusion_coefficient,
    'th_cl': geom.catalyst_layer_thickness,
    'calc_act_loss': sim.calc_activation_loss,
    'calc_cl_diff_loss': sim.calc_cl_loss,
    'calc_gdl_diff_loss': sim.calc_gdl_loss
    }

dict_anode = {
    'name': 'Anode',
    'flow_direction': geom.anode_flow_direction,
    'cell_width': geom.cell_width,
    'cell_length': geom.cell_length,
    'channel_numb': geom.gas_channel_number,
    'is_cathode': False,
    'species_names': op_con.anode_species,
    'molar_mass': op_con.anode_molar_mass,
    'composition': op_con.anode_inlet_composition,
    'charge_number': op_con.anode_charge_number,
    'reaction_stoichiometry': op_con.cathode_reaction_stoich,
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
