"""How far to trust the aqueous thermodynamics at a given temperature and pressure.

The HKF equations that give logK of aqueous species at any T and P are calibrated to about 5 kbar and
1000 C but lose accuracy where water becomes a poor solvent (low density, near the critical point).
There is no single reference that covers every species, but water itself is measured everywhere: the
ionization constant of water is given from 0 to 1000 C and to 1000 MPa by the IAPWS Revised Release
R11-24 (2024), the Bandura and Lvov (2006, J. Phys. Chem. Ref. Data 35, 15) equation with parameters
re-fitted in 2024 to newer near-critical and supercritical data.  It replaces, here, the older Marshall
and Franck (1981) equation, which is kept for comparison.
The engine's own logK for H2O = H+ + OH- depends on the same dielectric constant, g function and
Born coefficients as every charged species, so its distance from the measured value is a cheap,
independent indicator of how stretched the model is at that point.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from . import water

# IAPWS R11-24 (2024), Revised Release on the Ionization Constant of H2O, eqs. 2-4 and Table 1
# (T in K, rho in g/cm3); checked against all six verification values of its Table 2 in tests/test_pressure.py
_IAPWS_N = 6
_IAPWS_ALPHA = (-0.702132, 8681.05, -24145.1)
_IAPWS_BETA = (0.813876, -51.4471, -0.469920)
_IAPWS_MW = 18.015268
# Marshall and Franck (1981), J. Phys. Chem. Ref. Data 10, 295, eq. 8 (T in K, rho in g/cm3)
_MF = (-4.098, -3245.2, 2.2362e5, -3.984e7, 13.957, -1262.3, 8.5641e5)

# thresholds on |engine - measured| in log units
GOOD, CAUTION = 0.15, 0.35
LOW_DENSITY = 0.35          # g/cm3: below this the HKF equations are extrapolated
# one standard deviation of a tabulated log K where water ionization is reproduced: the usual spread
# between compilations for minerals and complexes at moderate temperature
BASE_SIGMA = 0.1


@dataclass(frozen=True)
class Reliability:
    t_c: float
    p_bar: float
    density: float           # g/cm3
    deviation: float         # engine log Kw minus IAPWS R11-24 log Kw
    level: str               # "good", "caution" or "poor"


def iapws_pkw(t_k: float, rho_g_cm3: float) -> float:
    """pKw = -log10 Kw (molal standard state) from IAPWS R11-24 eq. 2 at T (K) and density (g/cm3)."""
    a0, a1, a2 = _IAPWS_ALPHA
    b0, b1, b2 = _IAPWS_BETA
    rho = rho_g_cm3
    pkw_ideal = 0.61415 + 48251.33 / t_k - 67707.93 / t_k ** 2 + 10102100.0 / t_k ** 3     # eq. 4
    z = rho * math.exp(a0 + a1 / t_k + a2 / t_k ** 2 * rho ** (2.0 / 3.0)) if rho > 0 else 0.0  # eq. 3
    return (-2.0 * _IAPWS_N * (math.log10(1.0 + z) - z / (z + 1.0) * rho * (b0 + b1 / t_k + b2 * rho))
            + pkw_ideal + 2.0 * math.log10(_IAPWS_MW / 1000.0))


def iapws_log_kw(t_c: float, p_bar: float) -> tuple[float, float]:
    """(log10 Kw, water density in g/cm3) from IAPWS R11-24 with the IAPWS-95 density of liquid water."""
    t_k = t_c + 273.15
    rho = water.density(t_k, p_bar, liquid=True) / 1000.0
    return -iapws_pkw(t_k, rho), rho


def marshall_franck_log_kw(t_c: float, p_bar: float) -> tuple[float, float]:
    """(log10 Kw, water density in g/cm3) from the Marshall and Franck (1981) equation."""
    t_k = t_c + 273.15
    rho = water.density(t_k, p_bar, liquid=True) / 1000.0
    a, b, c, d, e, f, g = _MF
    log_kw = (a + b / t_k + c / t_k ** 2 + d / t_k ** 3
              + (e + f / t_k + g / t_k ** 2) * math.log10(rho))
    return log_kw, rho


def water_ionization_check(db, t_c: float, p_bar: float) -> Reliability:
    """Compare the database's logK of OH- (H2O = H+ + OH-) with IAPWS R11-24."""
    measured, rho = iapws_log_kw(t_c, p_bar)
    engine = db.species["OH-"].logk(float(t_c), float(p_bar))      # OH- + H+ = H2O, so -log Kw
    dev = engine + measured
    size = abs(dev)
    level = "good" if size <= GOOD and rho >= LOW_DENSITY else \
        "caution" if size <= CAUTION else "poor"
    return Reliability(t_c, p_bar, rho, dev, level)


def logk_sigma(db, t_c: float, p_bar: float) -> float:
    """One standard deviation (log units) assumed for every log K at T and P, for uncertainty bands:
    BASE_SIGMA with the water-ionization deviation added in quadrature, and at least the caution
    threshold at low water density.  Outside the liquid field the 'poor' half-width is used."""
    try:
        r = water_ionization_check(db, t_c, p_bar)
    except Exception:
        return 2.0 * CAUTION
    if not math.isfinite(r.deviation):
        return 2.0 * CAUTION
    sigma = math.hypot(BASE_SIGMA, r.deviation)
    return max(sigma, CAUTION) if r.density < LOW_DENSITY else sigma


LEVEL_TEXT = {
    "good": "Water ionization is reproduced to within {tol} log units: charged-species data are reliable here.",
    "caution": "Water ionization is off by up to {worst:.2f} log units here: treat charged-species "
               "results as approximate.",
    "poor": "Water ionization is off by up to {worst:.2f} log units here, near the critical point "
            "or at low water density: treat charged-species results as qualitative only.",
}
