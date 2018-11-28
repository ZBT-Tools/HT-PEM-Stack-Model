import input.operating_conditions as op_con
import input.geometry as geo
import input.simulation as sim

electrical_coupling_dict_const = {'cell_num': op_con.cell_number,
                                  'd_x': geo.channel_length / (sim.nodes - 1),
                                  'th_plate': geo.plate_thickness,
                                  'w_ch': geo.channel_width}


def electrical_coupling_dict_dyn(v_los, r_cell):
    return {'v_los': v_los, 'r_cell': r_cell}
