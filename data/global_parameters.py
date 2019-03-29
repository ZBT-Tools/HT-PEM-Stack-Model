import input.physical_properties as phy_prop
import input.operating_conditions as op_cond
import input.simulation as sim
import numpy as np


dict_uni = {
    'R': 8.314459848, 'F': 96485.3328959,
    'cp_liq': phy_prop.heat_capacity_coolant
            }

dict_case = {
    'tar_cd': np.array(op_cond.target_current_density),
    'nodes': sim.elements + 1, 'elements': sim.elements,
    'mol_con_m': phy_prop.molar_membrane_acid_group_concentration,
    'e_0': op_cond.open_circuit_voltage,
    'temp_env': op_cond.temp_environment, 'v_tn': phy_prop.v_thermo_neutral,
    'pem_type': op_cond.pem_type,
    'conv_coeff': phy_prop.convection_coefficient_stack_environment,
    'header_p_in_cat': op_cond.p_manifold_cathode_out,
    'header_p_in_ano': op_cond.p_manifold_anode_out,
    'bpp_resistivity': phy_prop.bipolar_plate_resistivity,
    'vap_m_temp_coeff': phy_prop.fitted_vapour_vapour_mass_transport_coefficient,
    'calc_act_loss': sim.calc_activation_loss,
    'calc_gdl_diff_loss': sim.calc_gdl_loss,
    'calc_cl_diff_loss': sim.calc_cl_loss,
    'calc_mem_loss': sim.calc_membrane_loss,
    'calc_temperature': sim.calc_temperature,
    'calc_current_density': sim.calc_current_density,
    'calc_flow_distribution': sim.calc_flow_distribution,
    'save_csv': sim.save_csv_data,
    'save_plt': sim.save_plot_data,
    'show_pol_curve': sim.show_voltage_loss
     }