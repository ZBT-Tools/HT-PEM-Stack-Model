import input.operating_conditions as op_con
import input.geometry as geom
import input.simulation as sim


dict_electrical_coupling =\
    {
        'cell_number': geom.cell_number,
        'dx': geom.channel_length / float(sim.elements),
        'th_bpp': geom.bipolar_plate_thickness,
        'width_channels': geom.channel_width * geom.gas_channel_number
                          + geom.rib_width * (geom.gas_channel_number + 1)
    }
