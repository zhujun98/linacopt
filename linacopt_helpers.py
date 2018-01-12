#!/usr/bin/python


def quad_k2g(k, p):
    """Convert the K value of a quadrupole to gradient.

    Parameters
    ----------
    k: float
        Quadrupole strength (1/m^2)
    p: float
        Normalized momentum

    Returns
    -------
    Quadrupole gradient (T/m)
    """
    me = 9.10938291e-31
    v_light = 299792458
    qe = 1.60217657e-19

    return -1.0*k*p*me*v_light/qe