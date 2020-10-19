from mpmath import mpf, sqrt, power, mp, fabs
def one_to_one_binding(p, l, kdpl):
    p = mpf(p)
    l = mpf(l)
    kdpl = mpf(kdpl)
    return ((p + kdpl + l - sqrt(-4 * p * l + power(p + kdpl + l, 2))) / 2.0).real
    