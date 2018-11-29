import copy as copy
import numpy as np
import input.operating_conditions as op_con
import input.geometry as geo
import input.physical_property as phy_prop


m_fold_const_cathode = {'cell_num': op_con.cell_number,
                        'header_width': geo.manifold_width,
                        'header_height': geo.manifold_height,
                        'kf': phy_prop.manifold_pressure_loss_coefficient,
                        'cell_height': np.full(op_con.cell_number,
                                               2. * (geo.bipolar_plate_thickness
                                                     + geo.gas_diffusion_layer_thickness
                                                     + geo.catalyst_layer_thickness)
                                               + geo.membrane_thickness),
                        'cell_ch_length': np.full(op_con.cell_number,
                                                  geo.channel_length),
                        'cell_ch_ca': np.full(op_con.cell_number,
                                              geo.channel_width
                                              * geo.channel_height),
                        'p_out': op_con.p_manifold_cathode_out}
m_fold_const_anode = copy.copy(m_fold_const_cathode)
m_fold_const_anode['p_out'] = op_con.p_manifold_anode_out


def dict_manifold_dyn(mol_flow, cell_t, cell_cp,
                      cell_visc, cell_p, cell_r, mass_flow):
    return {'mol_flow': mol_flow, 'cell_t': cell_t, 'cell_cp': cell_cp,
            'cell_visc': cell_visc, 'cell_p': cell_p, 'cell_r': cell_r,
            'mass_flow': mass_flow}
