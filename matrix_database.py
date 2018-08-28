from scipy.linalg import block_diag
import numpy as np
from scipy import sparse


def b(n, M, d_x):
    m = np.full((n, n), 0.)
    for j in range(n):
        for l in range(n):
            if l is 0 and j is 0:
                m[l, j] = -1.
                m[l, j + 1] = 1.
            elif j is n - 1 and l is n - 1:
                m[l, j] = 1.
                m[l, j - 1] = -1.
            elif 0 < l < n - 1 and j is l:
                m[l, j - 1] = 1
                m[l, j] = -2
                m[l, j + 1] = 1
    b = [m] * M
    return block_diag(*b) / d_x ** 2.


def c(n, M, d_x):
    m = np.full((n * M, n * M), 0.)
    for j in range(n * M):
        for l in range(n * M):
            if j is l and (n <= j < M * n - n):
                m[l, j] = -2.
                m[l, j - n] = 1.
                m[l, j + n] = 1.
            elif j is l and j <= n:
                m[l, j] = -1
                m[l, j + n] = 1
            elif j is l and j >= M * n - n:
                m[l, j] = 1
                m[l, j - n] = -1
    return m / d_x ** 2.


def forward_matrix_reac_flow(n):
    return np.tril(np.full((n, n), 1.))


def backward_matrix_reac_flow(n):
    return np.triu(np.full((n,n), 1.))


def col_mat_m1(nodes):
    mat = np.full((nodes, nodes), 0.)
    vec = np.hstack(([0], np.tile([1], nodes-1)))
    vec2 = np.full(nodes, 1.)
    return np.asarray(mat + sparse.spdiags(vec, 0, nodes, nodes) + sparse.spdiags(vec2, -1, nodes, nodes))


def col_mat_m2(nodes,var):
    mat = np.full((nodes, nodes), 0.)
    vec = np.hstack(([0], np.tile([1], nodes-1)))
    vec2 = np.full(nodes, -var)
    return np.asarray(mat + sparse.spdiags(vec, 0, nodes, nodes) + sparse.spdiags(vec2, -1, nodes, nodes))
