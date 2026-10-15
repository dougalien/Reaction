"""How far to trust the aqueous thermodynamics at a given temperature and pressure.

The HKF equations that give logK of aqueous species at any T and P are calibrated to about 5 kbar and
1000 C but lose accuracy where water becomes a poor solvent (low density, near the critical point).
There is no single reference that covers every species, but water itself is measured everywhere: the
ionization constant of water is known from 0 to 1000 C and to 10 kbar (Marshall and Franck, 1981).
The engine's own logK for H2O = H+ + OH- depends on the same dielectric constant, g function and
Born coefficients as every charged species, so its distance from the measured value is a cheap,
independent indicator of how stretched the model is at that point.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from . import water

# Marshall and Franck (1981), J. Phys. Chem. Ref. Data 10, 295, eq. 8 (T in K, rho in g/cm3)
_MF = (-4.098, -3245.2, 2.2362e5, -3.984e7, 13.957, -1262.3, 8.5641e5)

# thresholds on |engine - measured| in log units
GOOD, CAUTION = 0.15, 0.35
LOW_DENSITY = 0.35          # g/cm3: below this the HKF equations are extrapolated


@dataclass(frozen=True)
class Reliability:
    t_c: float
    p_bar: float
    density: float           # g/cm3
    deviation: float         # engine log Kw minus Marshall-Franck log Kw
    level: str               # "good", "caution" or "poor"


def marshall_franck_log_kw(t_c: float, p_bar: float) -> tuple[float, float]:
    """(log10 Kw, water density in g/cm3) from the Marshall and Franck (1981) equation."""
    t_k = t_c + 273.15
    rho = water.density(t_k, p_bar, liquid=True) / 1000.0
    a, b, c, d, e, f, g = _MF
    log_kw = (a + b / t_k + c / t_k ** 2 + d / t_k ** 3
              + (e + f / t_k + g / t_k ** 2) * math.log10(rho))
    return log_kw, rho


def water_ionization_check(db, t_c: float, p_bar: float) -> Reliability:
    """Compare the database's logK of OH- (H2O = H+ + OH-) with the measured value."""
    measured, rho = marshall_franck_log_kw(t_c, p_bar)
    engine = db.species["OH-"].logk(float(t_c), float(p_bar))      # OH- + H+ = H2O, so -log Kw
    dev = engine + measured
    size = abs(dev)
    level = "good" if size <= GOOD and rho >= LOW_DENSITY else \
        "caution" if size <= CAUTION else "poor"
    return Reliability(t_c, p_bar, rho, dev, level)


LEVEL_TEXT = {
    "good": "Water ionization is reproduced to within {tol} log units: charged-species data are reliable here.",
    "caution": "Water ionization is off by up to {worst:.2f} log units here: treat charged-species "
               "results as approximate.",
    "poor": "Water ionization is off by up to {worst:.2f} log units here, near the critical point "
            "or at low water density: treat charged-species results as qualitative only.",
}
