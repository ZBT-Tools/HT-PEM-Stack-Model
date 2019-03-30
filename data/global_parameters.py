import input.physical_properties as phy_prop
import input.operating_conditions as op_con
import input.simulation as sim
import numpy as np


dict_uni = {
    'R': 8.314459848,
    'F': 96485.3328959
    }

dict_case = {
    'tar_cd': np.array(op_con.target_current_density),
    'nodes': sim.elements + 1,
    'elements': sim.elements,
    'mol_con_m': phy_prop.molar_membrane_acid_group_concentration,
    'e_0': op_con.open_circuit_voltage,
    'temp_env': op_con.temp_environment,
    'v_tn': phy_prop.v_thermo_neutral,
    'conv_coeff': phy_prop.convection_coefficient_stack_environment,
    'header_p_in_cat': op_con.p_manifold_cathode_out,
    'header_p_in_ano': op_con.p_manifold_anode_out,
    'bpp_resistivity': phy_prop.bipolar_plate_resistivity,
    'cp_liq': phy_prop.heat_capacity_coolant,
    'vap_m_temp_coeff': phy_prop.fitted_vapour_vapour_mass_transport_coefficient
    }
