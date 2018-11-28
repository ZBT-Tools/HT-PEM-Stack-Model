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
                                               2. * (geo.plate_thickness
                                                     + geo.gdl_thickness
                                                     + geo.electrode_thickness)
                                               + geo.membrane_thickness),
                        'cell_ch_length': np.full(op_con.cell_number,
                                                  geo.channel_length),
                        'cell_ch_ca': np.full(op_con.cell_number,
                                              geo.channel_width
                                              * geo.channel_height),
                        'p_out': op_con.p_cathode_basic}
m_fold_const_anode = copy.copy(m_fold_const_cathode)
m_fold_const_anode['p_out'] = op_con.p_anode_basic


def dict_manifold_dyn(mol_flow, cell_t, cell_cp,
                      cell_visc, cell_p, cell_r, mass_flow):
    return {'mol_flow': mol_flow, 'cell_t': cell_t, 'cell_cp': cell_cp,
            'cell_visc': cell_visc, 'cell_p': cell_p, 'cell_r': cell_r,
            'mass_flow': mass_flow}
