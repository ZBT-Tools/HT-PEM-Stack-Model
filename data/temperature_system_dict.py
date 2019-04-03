import input.operating_conditions as op_con
import input.simulation as sim
import input.physical_properties as phy_prop
import input.geometry as geom

dict_temp_sys = {
    'cell_numb': op_con.cell_number,
    'nodes': sim.elements + 1,
    'channel_length': geom.channel_length,
    'channel_width': geom.coolant_channel_widht,
    'channel_height': geom.coolant_channel_height,
    'gas_ch_numb': geom.gas_channel_number,
    'cool_ch_numb': geom.coolant_channel_number,
    'cool_ch_bc': op_con.cooling_bc,
    'temp_gas_in': [op_con.temp_air_in, op_con.temp_anode_gas_in],
    'cool_cp': phy_prop.heat_capacity_coolant,
    'cool_m_flow': op_con.mass_flow_coolant,
    'cool_density': phy_prop.density_coolant,
    'cool_visc': phy_prop.dynamic_viscosity_coolant,
    'cool_th': geom.bipolar_plate_thickness * .5,
    'heat_pow': op_con.endplates_heat_power / float(sim.elements),
    'temp_layer_init': op_con.temp_initial,
    'cool_lambda':phy_prop.thermal_conductivity_coolant,
    'cool_temp_in': op_con.temp_coolant_in
    }
