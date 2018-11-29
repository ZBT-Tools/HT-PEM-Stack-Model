import input.geometry as geo
import input.physical_property as phy_prop
import input.operating_conditions as oper_con


cathode_channel = {'length': geo.channel_length, 'p_in': oper_con.p_manifold_cathode_out,
                   't_in': oper_con.temp_air_in,
                   'hum_in': phy_prop.inlet_humidity_cathode,
                   'flow_dir': phy_prop.cathode_channel_flow_direction,
                   'width': geo.channel_width, 'height': geo.channel_height,
                   'num_bends': phy_prop.channel_bends,
                   'bend_fri_fac': phy_prop.bend_pressure_loss_coefficient,
                   'rack_width': geo.rack_width}
anode_channel = {'length': geo.channel_length, 'p_in': oper_con.p_manifold_anode_out,
                 't_in': oper_con.temp_anode_gas_in,
                 'hum_in': phy_prop.inlet_humidity_anode,
                 'flow_dir': phy_prop.anode_channel_flow_direction,
                 'width': geo.channel_width, 'height': geo.channel_height,
                 'num_bends': phy_prop.channel_bends,
                 'bend_fri_fac': phy_prop.bend_pressure_loss_coefficient,
                 'rack_width': geo.rack_width}
