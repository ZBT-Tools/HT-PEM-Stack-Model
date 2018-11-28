import numpy as np


def fwd_mat_reac_flow(nodes):
    return np.tril(np.full((nodes, nodes), 1.))


def bwd_mat_reac_flow(nodes):
    return np.triu(np.full((nodes, nodes), 1.))
