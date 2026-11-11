"""The extended-term parameter b_gamma of the Helgeson, Kirkham and Flowers (1981) activity model.

HKF81 (Am. J. Sci. 281, 1249-1516) write the activity coefficient of an ion in an NaCl-dominated
solution as an extended Debye-Hueckel equation with one distance of closest approach for every ion
(a = 3.72 angstrom for NaCl, their Table 2), the conversion from mole fraction to molality, and an
extended term b_gamma * I that depends on temperature and pressure.

The b_gamma values and the way they are interpolated are those of CHNOSZ (Dick, 2019, Front. Earth
Sci. 7, 180; function ``bgamma`` in ``R/nonideal.R``, GPL-3.0), transcribed here number for number:

* along the water saturation curve (P = Psat): B-dot of Helgeson (1969, Am. J. Sci. 267, 729), read
  by CHNOSZ from his Figure (not his Table 2, whose spline turns sharply upward above 300 C);
* 0.5, 1, 2, 3, 4 and 5 kbar, 0 to 400 or 500 C: HKF81 Table 27;
* 550 and 600 C at 2 and 5 kbar, and 10, 20 and 30 kbar: read from Figure 11 of Manning et al. (2013,
  Rev. Mineral. Geochem. 75, chapter 5), so less precise than the tabulated values;
* 40, 50 and 60 kbar: CHNOSZ's own extrapolation of the 10-30 kbar values.
The 10-60 kbar rows are kept below for reference but are not used here (see below).

The numbers have not been checked here against the printed HKF81 table.  Each isobar is a cubic spline
in T (Forsythe, Malcolm and Moler end conditions, as R's ``splinefun``), and at a given T the
saturation row and the 0.5-5 kbar isobars are joined by a second such spline in P.

Departures from CHNOSZ, chosen for Reaction's range (to 600 C and 5 kbar): the 10-60 kbar isobars are
not used; the saturation node sits at the saturation pressure (1.013 bar to 100 C) instead of at 1, 16
or 86 bar; an isobar is not extended past its last tabulated temperature with its end cubic (CHNOSZ's
saturation row reaches -0.67 at 350 C that way) but held at its last value; and outside the pressure
nodes b_gamma is held at the end value.  These keep b_gamma continuous in T and P.  Every point that
relies on a held value is reported by :func:`extrapolated`.

Known limit: below 500 bar the pressure spline runs between the Psat row, which is Helgeson's (1969)
B-dot (fitted with individual ion sizes, without the mole-fraction term), and the 500 bar row of HKF81
(a = 3.72 angstrom, with that term).  The two are not the same parameter, so b_gamma at low pressure
inherits B-dot's values and the jump between 1 and 500 bar is partly an artefact of the mixed sources.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np

A_CIRC_NACL = 3.72          # angstrom, HKF81 Table 2 (via CHNOSZ)


class FmmSpline:
    """Cubic spline with the end conditions of Forsythe, Malcolm and Moler (1977), as R's ``splinefun``
    (method "fmm"): the third derivative at each end matches that of the cubic through the four points
    nearest that end.  Outside the knots the end cubic is extended."""

    def __init__(self, x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        n = len(x)
        self.x, self.y = x, y
        b, c, d = np.zeros(n), np.zeros(n), np.zeros(n)
        if n < 2:
            raise ValueError("need at least two points")
        if n < 3:
            b[:] = (y[1] - y[0]) / (x[1] - x[0])
            self.b, self.c, self.d = b, c, d
            return
        nm1 = n - 1
        d[0] = x[1] - x[0]
        c[1] = (y[1] - y[0]) / d[0]
        for i in range(1, nm1):
            d[i] = x[i + 1] - x[i]
            b[i] = 2.0 * (d[i - 1] + d[i])
            c[i + 1] = (y[i + 1] - y[i]) / d[i]
            c[i] = c[i + 1] - c[i]
        b[0] = -d[0]
        b[nm1] = -d[n - 2]
        c[0] = c[nm1] = 0.0
        if n > 3:
            c[0] = c[2] / (x[3] - x[1]) - c[1] / (x[2] - x[0])
            c[nm1] = c[n - 2] / (x[nm1] - x[n - 3]) - c[n - 3] / (x[n - 2] - x[n - 4])
            c[0] = c[0] * d[0] * d[0] / (x[3] - x[0])
            c[nm1] = -c[nm1] * d[n - 2] * d[n - 2] / (x[nm1] - x[n - 4])
        for i in range(1, n):
            t = d[i - 1] / b[i - 1]
            b[i] = b[i] - t * d[i - 1]
            c[i] = c[i] - t * c[i - 1]
        c[nm1] = c[nm1] / b[nm1]
        for i in range(n - 2, -1, -1):
            c[i] = (c[i] - d[i] * c[i + 1]) / b[i]
        b[nm1] = (y[nm1] - y[n - 2]) / d[n - 2] + d[n - 2] * (c[n - 2] + 2.0 * c[nm1])
        for i in range(nm1):
            b[i] = (y[i + 1] - y[i]) / d[i] - d[i] * (c[i + 1] + 2.0 * c[i])
            d[i] = (c[i + 1] - c[i]) / d[i]
            c[i] = 3.0 * c[i]
        c[nm1] = 3.0 * c[nm1]
        d[nm1] = d[n - 2]
        self.b, self.c, self.d = b, c, d

    def __call__(self, u: float) -> float:
        x = self.x
        i = int(np.searchsorted(x, u, side="right")) - 1
        i = min(max(i, 0), len(x) - 1)
        dx = u - x[i]
        return float(self.y[i] + dx * (self.b[i] + dx * (self.c[i] + dx * self.d[i])))


_25_500 = list(range(0, 525, 25))
# (pressure in bar, temperatures in C, b_gamma * 100), as in CHNOSZ R/nonideal.R
_ISOBARS = {
    "Psat": ([23.8, 49.4, 98.9, 147.6, 172.6, 197.1, 222.7, 248.1, 268.7],
             [4.07, 4.27, 4.30, 4.62, 4.86, 4.73, 4.09, 3.61, 1.56]),
    500: (list(range(0, 425, 25)),
          [5.6, 7.1, 7.8, 8.0, 7.8, 7.5, 7.0, 6.4, 5.7, 4.8, 3.8, 2.6, 1.0, -1.2, -4.1, -8.4, -15.2]),
    1000: (_25_500,
           [6.6, 7.7, 8.7, 8.3, 8.2, 7.9, 7.5, 7.0, 6.5, 5.9, 5.2, 4.4, 3.5, 2.5, 1.1, -0.6, -2.8, -5.7, -9.3,
            -13.7, -19.2]),
    2000: (_25_500 + [550, 600],
           [7.4, 8.3, 8.8, 8.9, 8.9, 8.7, 8.5, 8.1, 7.8, 7.4, 7.0, 6.6, 6.2, 5.8, 5.2, 4.6, 3.8, 2.9, 1.8, 0.5,
            -1.0, -3.93, -4.87]),
    3000: (_25_500,
           [6.5, 8.3, 9.2, 9.6, 9.7, 9.6, 9.4, 9.3, 9.2, 9.0, 8.8, 8.6, 8.3, 8.1, 7.8, 7.5, 7.1, 6.6, 6.0, 5.4,
            4.8]),
    4000: (_25_500,
           [4.0, 7.7, 9.5, 10.3, 10.7, 10.8, 10.8, 10.8, 10.7, 10.6, 10.5, 10.4, 10.3, 10.2, 10.0, 9.8, 9.6,
            9.3, 8.9, 8.5, 8.2]),
    5000: (_25_500 + [550, 600],
           [0.1, 6.7, 9.6, 11.1, 11.8, 12.2, 12.4, 12.4, 12.4, 12.4, 12.4, 12.3, 12.3, 12.2, 12.1, 11.9, 11.8,
            11.5, 11.3, 11.0, 10.8, 11.2, 12.52]),
    10000: ([25] + list(range(300, 1050, 50)),
            [12, 17.6, 17.8, 18, 18.2, 18.9, 21, 23.3, 26.5, 28.8, 31.4, 34.1, 36.5, 39.2, 41.6, 44.1]),
    20000: ([25] + list(range(300, 1050, 50)),
            [16, 21.2, 21.4, 22, 22.4, 23.5, 26.5, 29.2, 32.6, 35.2, 38.2, 41.4, 44.7, 47.7, 50.5, 53.7]),
    30000: ([25] + list(range(300, 1050, 50)),
            [19, 23.9, 24.1, 24.6, 25.2, 26.7, 30.3, 32.9, 36.5, 39.9, 43, 46.4, 49.8, 53.2, 56.8, 60]),
    40000: (list(range(300, 1050, 50)),
            [25.8, 26, 26.4, 27.2, 28.9, 33, 35.5, 39.2, 43.2, 46.4, 49.9, 53.4, 57.1, 61.2, 64.4]),
    50000: (list(range(300, 1050, 50)),
            [27.1, 27.3, 27.7, 28.5, 30.5, 34.8, 37.3, 41.1, 45.5, 48.7, 52.4, 55.9, 59.8, 64.3, 67.5]),
    60000: (list(range(300, 1050, 50)),
            [28, 28.2, 28.6, 29.5, 31.6, 36.1, 38.6, 42.5, 47.1, 50.4, 54.1, 57.6, 61.6, 66.5, 69.7]),
}


@lru_cache(maxsize=None)
def _isobar(key) -> FmmSpline:
    t, b = _ISOBARS[key]
    return FmmSpline(t, np.asarray(b) / 100.0)


T_PSAT_LAST = 268.7         # last temperature of the saturation row
_IN_RANGE = ("Psat", 500, 1000, 2000, 3000, 4000, 5000)


def _psat_node(t_c: float) -> float:
    """Pressure of the saturation-row node: 1.013 bar to 100 C, then the saturation pressure (capped
    just below the critical point), the curve Helgeson's (1969) B-dot refers to."""
    from .eosdb import reference_pressure_bar
    return reference_pressure_bar(min(t_c, 373.9))


def _isobar_value(key, t_c: float) -> float:
    """An isobar's spline inside its data; beyond its last temperature, its last value."""
    t = _ISOBARS[key][0]
    return _isobar(key)(min(t_c, t[-1]))


def _nodes(t_c: float):
    """(pressures, values) of the saturation row and the 0.5-5 kbar isobars at t_c."""
    pressures = [_psat_node(t_c)] + [float(k) for k in _IN_RANGE[1:]]
    values = [_isobar_value(k, t_c) for k in _IN_RANGE]
    return tuple(pressures), tuple(values)


class _Isotherm:
    def __init__(self, t_c: float):
        self.pressures, values = _nodes(t_c)
        self.lowest, self.highest = values[0], values[-1]
        self.spline = FmmSpline(self.pressures, values)

    def __call__(self, p_bar: float) -> float:
        # held constant outside the pressure nodes rather than extended with the end cubic
        if p_bar <= self.pressures[0]:
            return self.lowest
        if p_bar >= self.pressures[-1]:
            return self.highest
        return self.spline(p_bar)


@lru_cache(maxsize=4096)
def _isotherm(t_c: float) -> _Isotherm:
    return _Isotherm(t_c)


def bgamma(t_c: float, p_bar: float) -> float:
    """b_gamma (kg/mol) for NaCl-dominated solutions at t_c (degC) and p_bar (bar)."""
    if t_c < 0:
        raise ValueError("b_gamma is defined from 0 C")
    return _isotherm(round(float(t_c), 6))(float(p_bar))


def extrapolated(t_c: float, p_bar: float) -> str | None:
    """Why b_gamma at (T, P) rests on a held value rather than on tabulated values on both sides,
    or None when it does not."""
    if p_bar > 5000.0:
        return "above 5 kbar b_gamma is held at its 5 kbar value"
    pressures, _ = _nodes(t_c)
    below = [k for k, p in zip(_IN_RANGE, pressures) if p <= p_bar]
    above = [k for k, p in zip(_IN_RANGE, pressures) if p >= p_bar]
    near = ([below[-1]] if below else []) + ([above[0]] if above else [])
    held = [k for k in dict.fromkeys(near) if t_c > _ISOBARS[k][0][-1]]
    notes = []
    if held:
        names = " and ".join("saturation-curve" if k == "Psat" else f"{k / 1000:g} kbar" for k in held)
        last = ", ".join(f"{_ISOBARS[k][0][-1]:g} C" for k in held)
        notes.append(f"the {names} b_gamma data end at {last} and are held at their last value")
    if t_c > 500.0 and any(k in (2000, 5000) for k in near):
        notes.append("above 500 C the 2 and 5 kbar b_gamma values are read from a figure (Manning et al., 2013)")
    return "; ".join(notes) or None
