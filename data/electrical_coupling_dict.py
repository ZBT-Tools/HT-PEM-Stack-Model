import input.operating_conditions as op_con
import input.geometry as geom
import input.simulation as sim


dict_electrical_coupling =\
    {
        'cell_numb': op_con.cell_number,
        'dx': geom.channel_length / float(sim.elements),
        'th_bpp': geom.bipolar_plate_thickness,
        'width_channels': geom.channel_width * geom.gas_channel_number
                          + geom.rack_width * (geom.gas_channel_number + 1)
    }


def electrical_coupling(v_loss, r_cell):
    return {'v_loss': v_loss, 'r_cell': r_cell}
