#cython

cimport cython

import numpy as np
cimport numpy as np

from libc.math cimport exp

def c_calculate_solvent_probability(double solvent):
    cdef list coeffs

    coeffs = [-14.105436736742137,
              -0.47015366358636385,
              -2.9151681976244639,
              -0.49308859741473005,
              0.90132625209729045,
              0.033529051311488103,
              0.088901407582105796,
              0.10749856607909694,
              0.055000918494099861,
              -0.052424473641668454,
              -0.045698882840119227,
              0.076048484096718036,
              -0.097645159906868589,
              0.03904454313991608,
              -0.072186667173865071]

    chebyshev_poly = np.polynomial.Chebyshev(coeffs, domain=[0, 1])
    return exp(chebyshev_poly(solvent))

def c_get_max_score(list scores):
    cdef double max = 0.0
    cdef double n_copies = 1
    cdef double solvent_fraction = 0.5

    for i in range(len(scores)):
        prob = scores[i][-1]
        if prob > max:
            max = prob
            n_copies = scores[i][0]
            solvent_fraction = scores[i][1]
    return solvent_fraction, n_copies