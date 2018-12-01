import input.physical_property as phy_prop
import input.operating_conditions as op_cond
import input.simulation as sim
import numpy as np

dict_uni = {'R': 8.314459848, 'F': 96485.3328959,
            'cp_liq': phy_prop.heat_capacity_coolant}
dict_case = {'tar_cd': np.array(op_cond.target_current_density),
             'nodes': sim.elements + 1, 'elements': sim.elements,
             'mol_con_m': phy_prop.molar_membrane_acid_group_concentration,
             'e_0': op_cond.open_circuit_voltage,
             't_u': op_cond.temp_environment, 'vtn': phy_prop.v_thermo_neutral,
             'pem_type': op_cond.pem_type,
             'conv_coef': phy_prop.convection_coefficient_stack_environment,
             'header_p_in_cat': op_cond.p_manifold_cathode_out,
             'header_p_in_ano': op_cond.p_manifold_anode_out,
             'plate_resistivity': phy_prop.bipolar_plate_resistivity}
