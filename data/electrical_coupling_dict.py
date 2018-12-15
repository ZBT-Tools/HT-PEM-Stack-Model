import input.operating_conditions as op_con
import input.geometry as geo
import input.simulation as sim


dict_electrical_coupling =\
    {
        'cell_num': op_con.cell_number,
        'dx': geo.channel_length / float(sim.elements),
        'th_bpp': geo.bipolar_plate_thickness,
        'channel_width': geo.channel_width
    }


def electrical_coupling(v_loss, r_cell):
    return {'v_loss': v_loss, 'r_cell': r_cell}
