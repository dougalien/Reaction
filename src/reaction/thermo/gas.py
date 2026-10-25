"""Fugacity coefficients of a gas mixture from the Peng-Robinson (1976) equation of state.

The gas phase is a real mixture: van der Waals one-fluid mixing rules with binary interaction
parameters k_ij for the pairs that matter most in natural fluids (water with CO2, CH4, N2 and H2S;
CO2 with CH4 and N2).  Critical constants are from Poling, Prausnitz and O'Connell (2001); the k_ij are
typical vapour-phase values from the Peng-Robinson literature (water pairs fitted to the gas side, not to
gas solubility in the liquid, which is handled by the aqueous HKF data and the Drummond salting-out
coefficient).  Gases without critical constants here are treated as ideal (phi = 1) and said so.

Accuracy: within a few percent for CO2- and CH4-rich gases to a few hundred bar; poorer for water-rich
gases near the critical point of water, where the aqueous model is itself stretched.
"""
from __future__ import annotations

import math

import numpy as np

R_BAR = 83.14462618          # cm3 bar / (mol K)

# Tc (K), Pc (bar), acentric factor
CRITICAL = {
    "H2O": (647.10, 220.64, 0.3443), "CO2": (304.13, 73.77, 0.2250), "CH4": (190.56, 45.99, 0.0114),
    "N2": (126.20, 33.98, 0.0377), "O2": (154.58, 50.43, 0.0222), "H2": (33.19, 13.13, -0.2160),
    "H2S": (373.10, 89.63, 0.0900), "CO": (132.86, 34.94, 0.0450), "NH3": (405.40, 113.53, 0.2526),
    "SO2": (430.80, 78.84, 0.2454), "HCl": (324.70, 83.10, 0.1316), "HF": (461.0, 64.8, 0.3830),
    "N2O": (309.57, 72.45, 0.1409), "NO": (180.15, 64.80, 0.5820), "S2": (1313.0, 182.0, 0.2620),
}
_K = {frozenset(p): k for p, k in {
    ("H2O", "CO2"): 0.19, ("H2O", "CH4"): 0.48, ("H2O", "N2"): 0.48, ("H2O", "H2S"): 0.08,
    ("H2O", "H2"): 0.0, ("CO2", "CH4"): 0.10, ("CO2", "N2"): -0.02, ("CH4", "N2"): 0.03,
    ("CO2", "H2S"): 0.10, ("CH4", "H2S"): 0.08,
}.items()}


def formula_of(gas_name: str) -> str:
    """'CO2(g)' -> 'CO2'."""
    return gas_name[:-3] if gas_name.endswith("(g)") else gas_name


def _cubic_largest_root(c2: float, c1: float, c0: float) -> float:
    roots = np.roots([1.0, c2, c1, c0])
    real = roots[np.abs(roots.imag) < 1e-9].real
    return float(real.max())


def log_phi(names: list[str], y: np.ndarray, t_k: float, p_bar: float) -> np.ndarray:
    """log10 fugacity coefficients of gases ``names`` (formulas) at mole fractions ``y``."""
    y = np.asarray(y, dtype=float)
    out = np.zeros(len(names))
    known = [i for i, n in enumerate(names) if n in CRITICAL]
    if not known or p_bar <= 0:
        return out
    ys = y[known]
    total = ys.sum()
    if total <= 0:
        ys = np.full(len(known), 1.0 / len(known))
    else:
        ys = ys / total
    a, b = np.zeros(len(known)), np.zeros(len(known))
    for k, i in enumerate(known):
        tc, pc, w = CRITICAL[names[i]]
        kappa = 0.37464 + 1.54226 * w - 0.26992 * w * w
        alpha = (1.0 + kappa * (1.0 - math.sqrt(t_k / tc))) ** 2
        a[k] = 0.45724 * R_BAR ** 2 * tc ** 2 / pc * alpha
        b[k] = 0.07780 * R_BAR * tc / pc
    kij = np.array([[_K.get(frozenset((names[i], names[j])), 0.0) if i != j else 0.0
                     for j in known] for i in known])
    aij = np.sqrt(np.outer(a, a)) * (1.0 - kij)
    am = float(ys @ aij @ ys)
    bm = float(ys @ b)
    A = am * p_bar / (R_BAR * t_k) ** 2
    B = bm * p_bar / (R_BAR * t_k)
    z = _cubic_largest_root(-(1.0 - B), A - 3 * B * B - 2 * B, -(A * B - B * B - B ** 3))
    z = max(z, B + 1e-12)
    s2 = math.sqrt(2.0)
    term = math.log((z + (1 + s2) * B) / (z + (1 - s2) * B))
    ln_phi = (b / bm * (z - 1.0) - math.log(z - B)
              - A / (2 * s2 * B) * (2.0 * (aij @ ys) / am - b / bm) * term)
    out[known] = ln_phi / math.log(10.0)
    return out
