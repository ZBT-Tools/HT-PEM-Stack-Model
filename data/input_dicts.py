import input.physical_properties as phy_prop
import input.geometry as geom
import input.simulation as sim
import input.operating_conditions as op_con
import numpy as np
import copy


dict_stack = {
    'cell_number': geom.cell_number,
    'heat_power': op_con.endplates_heat_power,
    'header_height': geom.manifold_height,
    'header_width': geom.manifold_width,
    'dis_dis_fac': geom.manifold_pressure_loss_coefficient,
    'cool_ch_bc': geom.cooling_bc,
    'calc_temperature': sim.calc_temperature,
    'calc_current_density': sim.calc_current_density,
    'calc_flow_distribution': sim.calc_flow_distribution
    }

simulation_dict = {
    'maximum_iteration': sim.maximum_iteration_number,
    'minimum_iteration': sim.minimum_iteration_number,
    'iteration_criteria': sim.convergence_criteria
    }

dict_cell = {
    'is_ht_pem': op_con.is_ht_pem,
    'th_mem': geom.membrane_thickness,
    'width': geom.cell_width,
    'length': geom.cell_length,
    'thermal conductivity bpp':
        (phy_prop.thermal_conductivity_bipolar_plate_z,
         phy_prop.thermal_conductivity_bipolar_plate_x),
    'thermal conductivity gde':
        (phy_prop.thermal_conductivity_gas_diffusion_electrode_z,
         phy_prop.thermal_conductivity_gas_diffusion_electrode_x),
    'electrical conductivity bpp':
        phy_prop.electrical_conductivity_bipolar_plate,
    'electrical conductivity gde':
        phy_prop.electrical_conductivity_gde,
    'temp_cool_in': op_con.temp_coolant_in,
    'mem_base_r': phy_prop.membrane_basic_resistance,
    'mem_acl_r': phy_prop.membrane_temperature_coefficient,
    'temp_init': op_con.temp_initial
    }

dict_membrane = {
    'type': phy_prop.membrane_type,
    'thickness': geom.membrane_thickness,
    'acid group concentration':
        phy_prop.molar_membrane_acid_group_concentration,
    'vapour transport coefficient': phy_prop.vapour_mass_transport_coefficient,
    'ionic conductivity': phy_prop.membrane_basic_conductivity,
    'basic resistance': phy_prop.membrane_basic_resistance,
    'temperature coefficient':
        phy_prop.membrane_temperature_coefficient,
    'thermal conductivity':
        (phy_prop.thermal_conductivity_membrane_z,
         phy_prop.thermal_conductivity_membrane_x),
    'calc_loss': sim.calc_membrane_loss
}

dict_cathode = {
    'name': 'Cathode',
    'flow_direction': geom.cathode_flow_direction,
    'channel_number': geom.gas_channel_number,
    'stoichiometry': op_con.stoichiometry_cathode,
    'is_cathode': True,
    'species_names': op_con.cathode_species,
    #'molar_mass': op_con.cathode_molar_mass,
    'inlet_composition': op_con.cathode_inlet_composition,
    'charge_number': op_con.cathode_charge_number,
    'reaction_stoichiometry': op_con.cathode_reaction_stoich,
    'th_cl': geom.catalyst_layer_thickness,
    'th_gdl': geom.gas_diffusion_layer_thickness,
    'th_bpp': geom.bipolar_plate_thickness,
    'porosity cl': geom.catalyst_layer_porosity,
    'porosity gdl': geom.gas_diffusion_layer_porosity,
    'tafel_slope': phy_prop.tafel_slope_cathode,
    'prot_con_cl': phy_prop.catalyst_layer_proton_conductivity_cathode,
    'vol_ex_cd': phy_prop.exchange_current_density_cathode,
    'diff_coeff_cl': phy_prop.oxygen_catalyst_layer_diffusion_coefficient,
    'diff_coeff_gdl': phy_prop.oxygen_gas_diffusion_layer_diffusion_coefficient,
    'calc_act_loss': sim.calc_activation_loss,
    'calc_cl_diff_loss': sim.calc_cl_loss,
    'calc_gdl_diff_loss': sim.calc_gdl_loss
    }

dict_anode = {
    'name': 'Anode',
    'flow_direction': geom.anode_flow_direction,
    'channel_number': geom.gas_channel_number,
    'stoichiometry': op_con.stoichiometry_anode,
    'is_cathode': False,
    'species_names': op_con.anode_species,
    #'molar_mass': op_con.anode_molar_mass,
    'inlet_composition': op_con.anode_inlet_composition,
    'charge_number': op_con.anode_charge_number,
    'reaction_stoichiometry': op_con.anode_reaction_stoich,
    'th_cl': geom.catalyst_layer_thickness,
    'th_gdl': geom.gas_diffusion_layer_thickness,
    'th_bpp': geom.bipolar_plate_thickness,
    'porosity cl': geom.catalyst_layer_porosity,
    'porosity gdl': geom.gas_diffusion_layer_porosity,
    'tafel_slope': phy_prop.tafel_slope_anode,
    'prot_con_cl': phy_prop.catalyst_layer_proton_conductivity_anode,
    'vol_ex_cd': phy_prop.exchange_current_density_anode,
    'diff_coeff_cl': phy_prop.hydrogen_catalyst_layer_diffusion_coefficient,
    'diff_coeff_gdl': phy_prop.hydrogen_diffusion_layer_diffusion_coefficient,
    'calc_act_loss': sim.calc_activation_loss,
    'calc_cl_diff_loss': sim.calc_cl_loss,
    'calc_gdl_diff_loss': sim.calc_gdl_loss
    }

dict_cathode_channel = {
    'name': 'Cathode Channel',
    'channel_length': geom.channel_length,
    'p_out': op_con.p_manifold_cathode_out,
    'temp_in': op_con.temp_air_in,
    #'hum_in': op_con.inlet_humidity_cathode,
    'flow_direction': geom.cathode_flow_direction,
    'channel_width': geom.channel_width,
    'channel_height': geom.channel_height,
    'bend_number': geom.channel_bends,
    'bend_friction_factor': geom.bend_pressure_loss_coefficient
    }

dict_cathode_fluid = {
    'fluid_name': 'Cathode Gas',
    'fluid_components': op_con.cathode_species,
    'inlet_composition': op_con.cathode_inlet_composition,
    'temp_init': op_con.temp_air_in,
    'press_init': op_con.p_manifold_cathode_out
}

dict_anode_channel = {
    'name': 'Anode Channel',
    'channel_length': geom.channel_length,
    'p_out': op_con.p_manifold_anode_out,
    'temp_in': op_con.temp_anode_gas_in,
    #'hum_in': op_con.inlet_humidity_anode,
    'flow_direction': geom.anode_flow_direction,
    'channel_width': geom.channel_width,
    'channel_height': geom.channel_height,
    'bend_number': geom.channel_bends,
    'bend_friction_factor': geom.bend_pressure_loss_coefficient
    }

dict_anode_fluid = {
    'fluid_name': 'Anode Gas',
    'fluid_components': op_con.anode_species,
    'inlet_composition': op_con.anode_inlet_composition,
    'temp_init': op_con.temp_anode_gas_in,
    'press_init': op_con.p_manifold_anode_out
}

dict_cathode_manifold = {
    'name': 'Cathode Manifold',
    'configuration': geom.cathode_manifold_configuration,
    'p_out': op_con.p_manifold_cathode_out,
    'temp_in': op_con.temp_air_in,
    # 'hum_in': op_con.inlet_humidity_anode,
    'channel_width': geom.manifold_width,
    'channel_height': geom.manifold_height,
    'additional_friction_factor': geom.manifold_pressure_loss_coefficient
    }

dict_anode_manifold = {
    'name': 'Anode Manifold',
    'configuration': geom.anode_manifold_configuration,
    'p_out': op_con.p_manifold_anode_out,
    'temp_in': op_con.temp_anode_gas_in,
    'channel_width': geom.manifold_width,
    'channel_height': geom.manifold_height,
    'additional_friction_factor': geom.manifold_pressure_loss_coefficient
    }

# dict_anode_manifold = copy.deepcopy(dict_cathode_manifold)

dict_mfold_cat = {
    'name': 'Cathode Manifold',
    'cell_number': geom.cell_number,
    'header_width': geom.manifold_width,
    'header_height': geom.manifold_height,
    'kf': geom.manifold_pressure_loss_coefficient,
    'cell_height': np.full(geom.cell_number,
                           2. * (geom.bipolar_plate_thickness
                                 + geom.gas_diffusion_layer_thickness
                                 + geom.catalyst_layer_thickness)
                           + geom.membrane_thickness),
    # 'cell_channel_length': np.full(geom.cell_number, geom.channel_length),
    # 'cell_channel_cross_area': np.full(geom.cell_number,
    #                                    geom.channel_width
    #                                    * geom.channel_height),
    'p_out': op_con.p_manifold_cathode_out
    }

dict_mfold_ano = copy.deepcopy(dict_mfold_cat)
dict_mfold_ano['name'] = 'Anode Manifold'
dict_mfold_ano['p_out'] = op_con.p_manifold_anode_out

dict_electrical_coupling =\
    {
        'cell_number': geom.cell_number,
        'dx': geom.channel_length / float(sim.elements),
        'th_bpp': geom.bipolar_plate_thickness
        #'conducting_width': geom.rib_width * (geom.gas_channel_number + 1)
    }

dict_temp_sys = {
    'cell_number': geom.cell_number,
    'nodes': sim.elements + 1,
    'channel_length': geom.channel_length,
    'channel_width': geom.coolant_channel_width,
    'channel_height': geom.coolant_channel_height,
    'gas_ch_numb': geom.gas_channel_number,
    'cool_ch_numb': geom.coolant_channel_number,
    'cool_ch_bc': geom.cooling_bc,
    'temp_gas_in': [op_con.temp_air_in, op_con.temp_anode_gas_in],
    'cool_cp': phy_prop.heat_capacity_coolant,
    'cool_m_flow': op_con.mass_flow_coolant,
    'cool_density': phy_prop.density_coolant,
    'cool_visc': phy_prop.dynamic_viscosity_coolant,
    'cool_lambda': phy_prop.thermal_conductivity_coolant,
    'cool_th': geom.bipolar_plate_thickness * .5,
    'heat_pow': op_con.endplates_heat_power / float(sim.elements),
    'temp_layer_init': op_con.temp_initial,
    'temp_amb': op_con.temp_environment,
    'alpha_amb': op_con.convection_coefficient_stack_environment,
    'cool_temp_in': op_con.temp_coolant_in
    }

dict_output = {
    'save_csv': sim.save_csv_data,
    'save_plot': sim.save_plot_data,
    'show_loss': sim.show_voltage_loss
    }
