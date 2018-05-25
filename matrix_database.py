from scipy.linalg import block_diag
import numpy as np


def b(n, M, d_x):
    m = np.full((n, n), 0.)
    for j in range(n):
        for l in range(n):
            if l is 0 and j is 0:  # ok
                m[l, j] = -1.
                m[l, j + 1] = 1.
            elif j is n - 1 and l is n - 1:  # ok
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


def forward_matrix_reac_flow_new(n):  # element forward
    m = np.full((n, n), 0.)
    for j in range(n):
        for i in range(n):
            if 1 <= i <= j - 1:
                m[j, i] = -1.
            elif i is 0 and j >= 1:
                m[j, i] = -0.5
            elif i is j and j >= 1:
                m[j, i] = -0.5
    return m


def backward_matrix_reac_flow_new(n):  # element backward
    return np.flip(np.flip(forward_matrix_reac_flow_new(n), 0), 1)


def temperature_matrix_conv(n, M, mu_p, mu_g, mu_m, a, alpha_1, alpha_2, alpha_3):
    m = np.full((5, 5), 0.)
    m[0, 0] = - mu_g
    m[0, 1] = mu_p + mu_g + alpha_1
    m[0, 2] = - mu_p
    m[1, 0] = mu_g + mu_m + alpha_2
    m[1, 1] = - mu_g
    m[1, 3] = - mu_m
    m[2, 0] = - mu_m
    m[2, 3] = mu_g + mu_m + alpha_2
    m[2, 4] = - mu_g
    m[3, 3] = - mu_g
    m[3, 4] = mu_p + mu_g + alpha_1
    m[4, 1] = mu_p
    m[4, 2] = -2. * mu_p - a - alpha_3
    return np.linalg.inv(np.kron(np.eye(n * M), m))
