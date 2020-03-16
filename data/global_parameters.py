import settings.physical_properties as phy_prop
import settings.operating_conditions as op_con
import settings.simulation as sim
import numpy as np

SMALL = 1e-10

GAS_CONSTANT = 8.314459848
FARADAY = 96485.3328959

constants = {
    'R': 8.314459848,
    'F': 96485.3328959
    }

dict_case = {
    'target_current_density': op_con.target_current_density,
    'nodes': sim.elements + 1,
    'elements': sim.elements,
    'mol_con_m': phy_prop.molar_membrane_acid_group_concentration,
    'e_0': op_con.open_circuit_voltage,
    'v_tn': phy_prop.v_thermo_neutral,
    'header_p_in_cat': op_con.p_manifold_cathode_out,
    'header_p_in_ano': op_con.p_manifold_anode_out,
    'cp_liq': phy_prop.heat_capacity_coolant,
    'vap_m_temp_coeff': phy_prop.vapour_mass_transport_coefficient
    }
