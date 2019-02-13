import input.operating_conditions as op_con
import input.geometry as geom
import input.physical_properties as phy_prop
import input.simulation as sim

dict_stack = {
    'cell_numb': op_con.cell_number,
    'heat_power': op_con.endplates_heat_power,
    'header_height': geom.manifold_height,
    'header_width': geom.manifold_width,
    'dis_dis_fac': phy_prop.manifold_pressure_loss_coefficient,
    'stoi_cat': op_con.stoichiometry_cathode,
    'stoi_ano': op_con.stoichiometry_anode,
    'cool_ch_bc': op_con.cooling_bc,
    'alpha_env': phy_prop.convection_coefficient_stack_environment,
    'calc_temperature': sim.calc_temperature,
    'calc_current_density': sim.calc_current_density,
    'calc_flow_distribution': sim.calc_flow_distribution
    }
