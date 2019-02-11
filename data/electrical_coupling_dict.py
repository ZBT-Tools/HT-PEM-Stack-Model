import input.operating_conditions as op_con
import input.geometry as geo
import input.simulation as sim


dict_electrical_coupling =\
    {
        'cell_numb': op_con.cell_number,
        'dx': geo.channel_length / float(sim.elements),
        'th_bpp': geo.bipolar_plate_thickness,
        'width_channels': geo.channel_width * op_con.gas_channel_number
        + geo.rack_width * (op_con.gas_channel_number + 1)
    }


def electrical_coupling(v_loss, r_cell):
    return {'v_loss': v_loss, 'r_cell': r_cell}
