import input.operating_conditions as op_con
import input.simulation as sim
import input.physical_property as phy_prop
import input.geometry as geo
import numpy as np


dx = geo.channel_length / float(sim.elements)
k_cool_ch = np.pi * dx * phy_prop.convection_coefficient_coolant_channel \
            * geo.bipolar_plate_thickness * .5


temp_sys_const_dict = {'cell_num': op_con.cell_number,
                       'nodes': sim.elements + 1,
                       'cool_ch_bc': op_con.cooling_bc,
                       'temp_gas_in': [op_con.temp_air_in,
                                       op_con.temp_anode_gas_in],
                       'g_cool': phy_prop.heat_capacity_coolant
                                 * op_con.mass_flow_coolant,
                       'heat_pow': op_con.endplates_heat_power
                                   / float(sim.elements),
                       'temp_layer_init': op_con.temp_initial,
                       'k_cool_ch': k_cool_ch,
                       'temp_cool_in': op_con.temp_coolant_in}


def temp_sys_dyn_dict(k_alpha_ch, gamma, omega, v_los,
                   m_reac_flow_delta, g_gas, cp_h2, i):
    temp_dict = {'k_gas_ch': k_alpha_ch, 'gamma': gamma, 'omega': omega,
                 'v_los': v_los, 'm_reac_flow_delta': m_reac_flow_delta,
                 'g_gas': g_gas, 'cp_h2': cp_h2, 'i': i}
    return temp_dict



