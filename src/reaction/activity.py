"""Aqueous activity models (extended Debye-Hueckel family).

All ``log_*`` functions return base-10 logarithms.  Formulas follow the
EQ3/6 theory manual (Wolery, 1992) and the sources it cites; each was
checked numerically against EQ3NR output (tests/test_eq36_seawater.py).
"""
from __future__ import annotations

import math

from .constants import KELVIN, LN10, WATER_MOLES_PER_KG

HARD_CORE_WATER_ANGSTROM = 4.0   # fixed solute diameter in the water-activity equation
DEFAULT_ION_SIZE = 4.0           # angstrom, used when the database gives none


def log_gamma_bdot(z: float, a0: float, ionic_strength: float,
                   A: float, B: float, bdot: float) -> float:
    """Helgeson (1969) B-dot equation for a charged species."""
    s = math.sqrt(ionic_strength)
    return -A * z * z * s / (1.0 + a0 * B * s) + bdot * ionic_strength


def log_gamma_davies(z: float, ionic_strength: float, A: float) -> float:
    s = math.sqrt(ionic_strength)
    return -A * z * z * (s / (1.0 + s) - 0.3 * ionic_strength)


def log_gamma_nonpolar(ionic_strength: float, t_c: float, cco2) -> float:
    """Drummond (1981) salting-out polynomial for nonpolar neutrals (CO2, O2, H2...)."""
    t = t_c + KELVIN
    c1, c2, c3, c4, c5 = cco2
    ln_g = (c1 + c2 * t + c3 / t) * ionic_strength \
        - (c4 + c5 * t) * ionic_strength / (ionic_strength + 1.0)
    return ln_g / LN10


def log_gamma_chloride_bates_guggenheim(ionic_strength: float, A: float) -> float:
    """Bates-Guggenheim convention fixing the NBS pH scale."""
    s = math.sqrt(ionic_strength)
    return -A * s / (1.0 + 1.5 * s)


def _sigma(x: float) -> float:
    if x < 1e-3:
        return 1.0 - 1.5 * x + 1.8 * x * x
    return 3.0 / x ** 3 * (1.0 + x - 1.0 / (1.0 + x) - 2.0 * math.log(1.0 + x))


def log_activity_water_bdot(sum_molality: float, ionic_strength: float,
                            A: float, B: float, bdot: float) -> float:
    """log10 a_w, quasi-consistent with the B-dot equation (EQ3/6 eq. 94)."""
    x = HARD_CORE_WATER_ANGSTROM * B * math.sqrt(ionic_strength)
    return -(1.0 / WATER_MOLES_PER_KG) * (
        sum_molality / LN10
        - (2.0 / 3.0) * A * ionic_strength ** 1.5 * _sigma(x)
        + bdot * ionic_strength ** 2
    )
