# Global functions:
from numpy import reshape, exp


def dw(t):
    return 2.1 * 10. ** -7. * exp(-2436. / t)


def d_pos(t):
    return 1.6 * 10. ** -8. * exp(-1683. / t)


def toarray(var, m, n):
    return reshape(var.flatten(order='C'), (m, n))
