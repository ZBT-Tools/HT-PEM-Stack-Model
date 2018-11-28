import input.geometry as geo
import input.physical_property as phy_prop
import input.operating_conditions as oper_con


cathode_channel = {'length': geo.channel_length, 'p_in': oper_con.p_cathode_basic,
                   't_in': oper_con.t_cathode_in,
                   'hum_in': phy_prop.cathode_inlet_humidity,
                   'flow_dir': phy_prop.cathode_channel_flow_direction,
                   'width': geo.channel_width, 'height': geo.channel_height,
                   'num_bends': phy_prop.channel_bends,
                   'bend_fri_fac': phy_prop.bend_pressure_loss_coefficient,
                   'rack_width': geo.rack_width}
anode_channel = {'length': geo.channel_length, 'p_in': oper_con.p_anode_basic,
                 't_in': oper_con.t_anode_in,
                 'hum_in': phy_prop.anode_inlet_humidity,
                 'flow_dir': phy_prop.anode_channel_flow_direction,
                 'width': geo.channel_width, 'height': geo.channel_height,
                 'num_bends': phy_prop.channel_bends,
                 'bend_fri_fac': phy_prop.bend_pressure_loss_coefficient,
                 'rack_width': geo.rack_width}
