import input.geometry as geom
import input.physical_property as phy_prop
import input.operating_conditions as oper_con


dict_cathode_channel = {
    'channel_length': geom.channel_length,
    'p_in': oper_con.p_manifold_cathode_out,
    'temp_in': oper_con.temp_air_in,
    'hum_in': phy_prop.inlet_humidity_cathode,
    'flow_dir': phy_prop.cathode_channel_flow_direction,
    'channel_width': geom.channel_width,
    'channel_height': geom.channel_height,
    'bend_numb': phy_prop.channel_bends,
    'bend_fri_fac': phy_prop.bend_pressure_loss_coefficient,
    'rack_width': geom.rack_width
                        }

dict_anode_channel = {
    'channel_length': geom.channel_length,
    'p_in': oper_con.p_manifold_anode_out,
    'temp_in': oper_con.temp_anode_gas_in,
    'hum_in': phy_prop.inlet_humidity_anode,
    'flow_dir': phy_prop.anode_channel_flow_direction,
    'channel_width': geom.channel_width,
    'channel_height': geom.channel_height,
    'bend_numb': phy_prop.channel_bends,
    'bend_fri_fac': phy_prop.bend_pressure_loss_coefficient,
    'rack_width': geom.rack_width
                      }
