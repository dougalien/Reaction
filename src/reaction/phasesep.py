"""Phase separation of a saline fluid: boiling and condensation in the NaCl-H2O system.

Above the critical point of water, and below the boiling curve, a saline fluid can split into a brine and a
vapour of different salinity.  Where that happens, and the compositions of the two phases, come from the
correlation formulae of Driesner and Heinrich (2007) in ``nacl_h2o.py`` (0 to 1000 C, to 5000 bar).  They put
the two-phase boundary of seawater within 3 bar of Bischoff and Rosenbauer's (1988) measurements from 196 to
450 C, and its critical point at 408 C and 301 bar (measured: 407 C, 298.5 bar).  The phases' densities come
from the lookup table in ``data/nacl_h2o_phases.tsv`` (see its header for the source; 1 to 996 C and 10 to
490 bar on a 5 C by 20 bar grid, interpolated bilinearly).  The fluid's salinity is taken as NaCl-equivalent
from its chloride, as is usual for vent fluids and fluid inclusions.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from . import nacl_h2o as dh

DATA = Path(__file__).with_name("data") / "nacl_h2o_phases.tsv"
NACL_G_PER_MOL = 58.443
T_RANGE = (0.0, 1000.0)        # the correlations' range
P_RANGE = (1.0, 5000.0)
COLUMNS = ("x_liq", "x_vap", "rho_liq", "rho_vap", "h_liq", "h_vap", "x_vap_halite", "x_liq_halite")


@lru_cache(maxsize=1)
def grid() -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """(temperatures, pressures, {column: array [P, T]}) with NaN where a field does not exist."""
    rows = [line.split("\t") for line in DATA.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#") and not line.startswith("t_c")]
    t = np.array(sorted({float(r[0]) for r in rows}))
    p = np.array(sorted({float(r[1]) for r in rows}))
    it = {v: i for i, v in enumerate(t)}
    ip = {v: i for i, v in enumerate(p)}
    cols = {c: np.full((len(p), len(t)), np.nan) for c in COLUMNS}
    for r in rows:
        i, j = ip[float(r[1])], it[float(r[0])]
        for c, v in zip(COLUMNS, r[2:]):
            if v:
                cols[c][i, j] = float(v)
    return t, p, cols


def _at(name: str, t_c: float, p_bar: float) -> float:
    """Bilinear interpolation; NaN unless the field exists at all four corners."""
    t, p, cols = grid()
    a = cols[name]
    j = int(np.clip(np.searchsorted(t, t_c) - 1, 0, len(t) - 2))
    i = int(np.clip(np.searchsorted(p, p_bar) - 1, 0, len(p) - 2))
    ft = (t_c - t[j]) / (t[j + 1] - t[j])
    fp = (p_bar - p[i]) / (p[i + 1] - p[i])
    c = a[i:i + 2, j:j + 2]
    if np.isnan(c).any():
        # the edge of the field: the nearest corner, if it is in the field
        near = a[i + int(round(fp)), j + int(round(ft))]
        return float(near)
    return float((c[0, 0] * (1 - ft) + c[0, 1] * ft) * (1 - fp) + (c[1, 0] * (1 - ft) + c[1, 1] * ft) * fp)


def salinity_from_chloride(cl_molal: float) -> float:
    """NaCl-equivalent mass fraction of a fluid with ``cl_molal`` mol of Cl per kg of water."""
    g = max(cl_molal, 0.0) * NACL_G_PER_MOL
    return g / (1000.0 + g)


@dataclass(frozen=True)
class PhaseState:
    kind: str                     # "one phase", "liquid + vapour", "vapour + halite", "liquid + halite", "outside"
    x_bulk: float                 # NaCl mass fraction of the whole fluid
    x_liq: float = math.nan       # of the brine
    x_vap: float = math.nan       # of the vapour
    vapour_mass: float = 0.0      # share of the fluid's mass in the vapour (lever rule)
    rho_liq: float = math.nan
    rho_vap: float = math.nan

    @property
    def separated(self) -> bool:
        return self.kind == "liquid + vapour"


def _grid_value(name: str, t_c: float, p_bar: float) -> float:
    t, p, _ = grid()
    if t[0] <= t_c <= t[-1] and p[0] <= p_bar <= p[-1]:
        return _at(name, t_c, p_bar)
    return math.nan


def _coexisting(t_c: float, p_bar: float) -> tuple[float, float]:
    """NaCl mass fractions of the coexisting liquid and vapour at T and P, or NaN outside the field."""
    if not dh.in_two_phase_field(t_c, p_bar):
        return math.nan, math.nan
    try:
        xl, xv = dh.x_liquid(t_c, p_bar), dh.x_vapour(t_c, p_bar)
    except (ValueError, ZeroDivisionError, OverflowError):
        return math.nan, math.nan
    if not (0.0 <= xv < xl <= 1.0):
        return math.nan, math.nan
    return dh.mass_fraction(xl), dh.mass_fraction(xv)


def phase_state(t_c: float, p_bar: float, x_bulk: float) -> PhaseState:
    """Whether a fluid of NaCl mass fraction ``x_bulk`` is one phase or splits at T and P."""
    if not (T_RANGE[0] <= t_c <= T_RANGE[1] and P_RANGE[0] <= p_bar <= P_RANGE[1]):
        return PhaseState("outside", x_bulk)
    xl, xv = _coexisting(t_c, p_bar)
    if math.isfinite(xl) and xv < x_bulk < xl:
        return PhaseState("liquid + vapour", x_bulk, xl, xv, (xl - x_bulk) / (xl - xv),
                          _grid_value("rho_liq", t_c, p_bar), _grid_value("rho_vap", t_c, p_bar))
    if t_c < dh.T_TRIPLE_NACL and p_bar <= dh.p_vlh(t_c):
        xvh = dh.mass_fraction(dh.x_vapour_halite(t_c, p_bar))
        if x_bulk > xvh:
            return PhaseState("vapour + halite", x_bulk, math.nan, xvh, 1.0)
        return PhaseState("one phase", x_bulk)
    xlh = dh.mass_fraction(dh.x_halite_liquid(t_c, p_bar)) if t_c < dh.t_halite_melting(p_bar) else 1.0
    if x_bulk > xlh:
        return PhaseState("liquid + halite", x_bulk, xlh)
    return PhaseState("one phase", x_bulk)
