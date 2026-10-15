"""Dielectric constant of water, Born functions and Debye-Hueckel parameters.

The dielectric constant is the Johnson and Norton (1991, Am. J. Sci. 291, 541)
polynomial in density and temperature, the form SUPCRT92 uses and against which
the HKF parameters of aqueous species were regressed.  Density comes from
IAPWS-95 (`water.py`).  Derivatives are central finite differences of
w = -1/epsilon, which is smooth away from the critical point.

Units: T in K, P in bar, density g/cm3.
"""
from __future__ import annotations

import math

from . import water

T_REF = 298.15

# Johnson and Norton (1991), eqn. 9 coefficients
_K = (0.1470333593e2, 0.2128462733e3, -0.1154445173e3, 0.1955210915e2,
      -0.8330347980e2, 0.3213240048e2, -0.6694098645e1,
      -0.3786202045e2, 0.6887359646e2, -0.2729401652e2)


def epsilon_from_density(rho_g_cm3: float, t_k: float) -> float:
    k = _K
    tt = t_k / T_REF
    r = rho_g_cm3
    return (1.0 + k[0] * r / tt
            + (k[1] / tt + k[2] + k[3] * tt) * r ** 2
            + (k[4] / tt + k[5] * tt + k[6] * tt ** 2) * r ** 3
            + (k[7] / tt ** 2 + k[8] / tt + k[9]) * r ** 4)


def epsilon(t_k: float, p_bar: float, liquid: bool = True) -> float:
    """Static dielectric constant of water at T, P."""
    rho = water.density(t_k, p_bar, liquid) / 1000.0
    return epsilon_from_density(rho, t_k)


def born_z(t_k: float, p_bar: float, liquid: bool = True) -> float:
    """Born function Z = -1/epsilon - 1... returned as -1/epsilon (SUPCRT's ZBorn + 1)."""
    return -1.0 / epsilon(t_k, p_bar, liquid)


def born_functions(t_k: float, p_bar: float, liquid: bool = True,
                   h_t: float = 0.05, h_p: float = 0.5):
    """(Y, Q, X) = (dw/dT, dw/dP, d2w/dT2) of w = -1/epsilon, at constant P resp. T.

    Y in 1/K, Q in 1/bar, X in 1/K2.
    """
    w0 = born_z(t_k, p_bar, liquid)
    wp = born_z(t_k + h_t, p_bar, liquid)
    wm = born_z(t_k - h_t, p_bar, liquid)
    y = (wp - wm) / (2.0 * h_t)
    x = (wp - 2.0 * w0 + wm) / h_t ** 2
    q = (born_z(t_k, p_bar + h_p, liquid) - born_z(t_k, max(p_bar - h_p, 1e-3), liquid)) / (
        (p_bar + h_p) - max(p_bar - h_p, 1e-3))
    return y, q, x


def debye_huckel_a_b(t_k: float, p_bar: float, liquid: bool = True) -> tuple[float, float]:
    """A (kg^0.5 mol^-0.5, base-10 form) and B (kg^0.5 mol^-0.5 angstrom^-1).

    A = 1.8246e6 sqrt(rho) / (eps T)^1.5 and B = 50.29 sqrt(rho) / (eps T)^0.5
    with rho in g/cm3 (Helgeson and Kirkham, 1974; SUPCRT92 constants).
    """
    rho = water.density(t_k, p_bar, liquid) / 1000.0
    eps = epsilon_from_density(rho, t_k)
    a = 1.8246e6 * math.sqrt(rho) / (eps * t_k) ** 1.5
    b = 50.29 * math.sqrt(rho) / (eps * t_k) ** 0.5
    return a, b
