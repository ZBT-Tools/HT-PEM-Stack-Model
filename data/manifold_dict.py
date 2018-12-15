import copy as copy
import numpy as np
import input.operating_conditions as op_con
import input.geometry as geom
import input.physical_property as phy_prop


dict_mfold_cat = {
    'cell_num': op_con.cell_number,
    'header_width': geom.manifold_width,
    'header_height': geom.manifold_height,
    'kf': phy_prop.manifold_pressure_loss_coefficient,
    'cell_height': np.full(op_con.cell_number,
                           2. * (geom.bipolar_plate_thickness
                                 + geom.gas_diffusion_layer_thickness
                                 + geom.catalyst_layer_thickness)
                           + geom.membrane_thickness),
    'cell_channel_length': np.full(op_con.cell_number, geom.channel_length),
    'cell_channel_cross_area': np.full(op_con.cell_number,
                                       geom.channel_width
                                       * geom.channel_height),
    'p_out': op_con.p_manifold_cathode_out
                    }

dict_mfold_ano = copy.copy(dict_mfold_cat)
dict_mfold_ano['p_out'] = op_con.p_manifold_anode_out


def manifold(mol_flow, cell_temp, cell_cp, cell_visc,
             cell_p, cell_r, mass_flow):
    return {'mol_flow': mol_flow, 'cell_temp': cell_temp, 'cell_cp': cell_cp,
            'cell_visc': cell_visc, 'cell_p': cell_p, 'cell_r': cell_r,
            'mass_flow': mass_flow}
