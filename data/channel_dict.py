import input.geometry as geom
import input.physical_properties as phy_prop
import input.operating_conditions as op_con


dict_cathode_channel = {
    'channel_length': geom.channel_length,
    'p_out': op_con.p_manifold_cathode_out,
    'temp_in': op_con.temp_air_in,
    'hum_in': op_con.inlet_humidity_cathode,
    'flow_dir': geom.cathode_flow_direction,
    'channel_width': geom.channel_width,
    'channel_height': geom.channel_height,
    'bend_number': geom.channel_bends,
    'bend_fri_fac': geom.bend_pressure_loss_coefficient,
    'rib_width': geom.rib_width
    }

dict_anode_channel = {
    'channel_length': geom.channel_length,
    'p_out': op_con.p_manifold_anode_out,
    'temp_in': op_con.temp_anode_gas_in,
    'hum_in': op_con.inlet_humidity_anode,
    'flow_dir': geom.anode_flow_direction,
    'channel_width': geom.channel_width,
    'channel_height': geom.channel_height,
    'bend_number': geom.channel_bends,
    'bend_fri_fac': geom.bend_pressure_loss_coefficient,
    'rib_width': geom.rib_width
    }
