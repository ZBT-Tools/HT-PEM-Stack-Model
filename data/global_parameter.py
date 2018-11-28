import input.physical_property as phy_prop
import input.operating_conditions as op_cond
import input.simulation as sim
import numpy as np

dict_uni = {'R': 8.314459848, 'F': 96485.3328959,
            'h_vap': phy_prop.enthalpy_vaporization,
            'cp_liq': phy_prop.cp_coolant}
dict_case = {'tar_cd': np.array(op_cond.target_current_density),
             'nodes': sim.nodes, 'elements': sim.nodes-1,
             'vap_m_t_coef': phy_prop.vapour_mass_transport_coefficient,
             'mol_con_m': phy_prop.molar_membrane_acid_group_concentration,
             'e_0': op_cond.open_circuit_voltage,
             't_u': op_cond.t_environment, 'vtn': phy_prop.v_thermo_neutral,
             'pem_type': op_cond.pem_type,
             'conv_coef': phy_prop.convection_coefficient_stack,
             'header_p_in_cat': op_cond.p_cathode_basic,
             'header_p_in_ano': op_cond.p_anode_basic,
             'plate_resistivity': phy_prop.bipolar_plate_resistivity}
