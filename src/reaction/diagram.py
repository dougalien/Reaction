"""Predominance and activity diagrams at any temperature and pressure.

A diagram is a map over two variables (pH, Eh, log fO2, the log activity of a basis species, or the log
of an activity ratio such as a(K+)/a(H+)) of the species or mineral of one element that is most stable,
the others held fixed.  It is computed from the same logK as every run, so it holds at any (T, P) the
data cover, not only at 25 C or along the saturation curve.

Each candidate is written over the basis (H2O, H+, O2(g), the element's master species and the masters of
any other elements asked for).  Per atom of the element, the candidate that forms with the largest
affinity wins (the "balanced on the element" rule of Helgeson and of CHNOSZ's ``diagram``):

    g_i = (sum_j nu_ij log a_j - logK_i - log a_i) / n_i

with log a_i the activity given to dissolved species of the element (1e-6 by default) and 0 for minerals.
The master's own activity drops out, because every candidate carries n_i of it.

Reliability is shown two ways.  A field is flagged when the data for its species are extrapolated at
this (T, P).  And a band along every boundary marks where the winner leads the runner-up by less than the
logK uncertainty at this point, taken from how well water ionization is reproduced here (see
``thermo/validity.py``): 0.15 log units where it is good, 0.35 where it is fair, 0.7 where it is poor.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .constants import FARADAY, LN10, R
from .database.model import Kind
from .network import OXYGEN, PROTON, WATER, Network
from .thermo import validity
from .thermo.eosdb import DEFAULT_MASTERS

AXIS_KINDS = ("pH", "Eh", "log_fO2", "log_a", "log_ratio")
BAND = {"good": 0.15, "caution": 0.35, "poor": 0.7}
DEFAULT_LOG_A = -6.0          # dissolved species of the element
DEFAULT_OTHER_LOG_A = -3.0    # masters of the other elements, when not an axis


class DiagramError(ValueError):
    pass


@dataclass
class Axis:
    kind: str                         # one of AXIS_KINDS
    lo: float
    hi: float
    species: str | None = None        # the basis species of log_a and log_ratio

    def label(self) -> str:
        if self.kind == "pH":
            return "pH"
        if self.kind == "Eh":
            return "Eh (V)"
        if self.kind == "log_fO2":
            return "log fO2"
        if self.kind == "log_a":
            return f"log a({self.species})"
        return f"log a({self.species})/a(H+)"


@dataclass
class Diagram:
    element: str
    t_c: float
    p_bar: float
    x: Axis
    y: Axis
    xs: np.ndarray
    ys: np.ndarray
    names: list[str]                  # candidates, the index of the grid
    kinds: list[str]                  # "aqueous" or "mineral"
    grid: np.ndarray                  # (ny, nx) index of the winner
    margin: np.ndarray                # (ny, nx) lead of the winner over the runner-up, log units
    band: float
    reliability: validity.Reliability | None
    extrapolated: set[str]
    water: dict[str, list[tuple[float, float]]] = field(default_factory=dict)
    limits: list[dict] = field(default_factory=list)      # saturation of a mineral of an axis species
    fixed: dict[str, float] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def present(self) -> list[str]:
        return [self.names[i] for i in sorted(set(self.grid.ravel().tolist()))]

    def at(self, x: float, y: float) -> str:
        i = int(np.argmin(abs(self.ys - y)))
        j = int(np.argmin(abs(self.xs - x)))
        return self.names[self.grid[i, j]]


def _eh_factor(t_c: float) -> float:
    """Volts per unit of pe."""
    return R * (t_c + 273.15) * LN10 / FARADAY


def predominance(db, element: str, x: Axis, y: Axis, *, t_c: float = 25.0, p_bar: float | None = None,
                 ph: float = 7.0, log_a: float = DEFAULT_LOG_A, log_fo2: float | None = None,
                 others: dict[str, float] | None = None, include: str = "all",
                 exclude=(), nx: int = 160, ny: int = 140) -> Diagram:
    """The most stable species of ``element`` over the plane of ``x`` and ``y``.

    ``others`` gives log activities of further basis species (or elements, by their master), which
    brings their compounds in (S for pyrite, C for siderite).  ``include`` is "all", "minerals" or
    "aqueous".  ``exclude`` lists minerals to leave out.
    """
    master = DEFAULT_MASTERS.get(element, element)
    if master not in db.species:
        raise DiagramError(f"no master species for {element}")
    others = dict(others or {})
    fixed: dict[str, float] = {}
    for key, value in others.items():
        name = DEFAULT_MASTERS.get(key, key)
        if name not in db.species:
            raise DiagramError(f"unknown species {key}")
        fixed[name] = float(value)
    for ax in (x, y):
        if ax.kind not in AXIS_KINDS:
            raise DiagramError(f"unknown axis {ax.kind}")
        if ax.kind in ("log_a", "log_ratio"):
            ax.species = DEFAULT_MASTERS.get(ax.species, ax.species)
            if ax.species not in db.species:
                raise DiagramError(f"unknown species {ax.species}")
            fixed.setdefault(ax.species, DEFAULT_OTHER_LOG_A)
    redox_axes = [ax for ax in (x, y) if ax.kind in ("Eh", "log_fO2")]
    if len(redox_axes) > 1:
        raise DiagramError("Eh and log fO2 cannot both be axes")
    if sum(ax.kind == "pH" for ax in (x, y)) > 1:
        raise DiagramError("pH cannot be both axes")
    basis = [WATER, PROTON, OXYGEN, master] + [n for n in fixed if n not in (WATER, PROTON, OXYGEN,
                                                                               master)]
    p = db.reference_pressure(t_c) if p_bar is None else float(p_bar)
    net = Network(db, basis)
    ev = net.evaluate(t_c, p)
    idx = {n: k for k, n in enumerate(net.basis)}
    im = idx[master]
    excluded = {e.lower() for e in exclude}

    rows, logk, assigned, names, kinds = [], [], [], [master], ["aqueous"]
    extrapolated = {n for n in ev.extrapolated}
    if include != "minerals":
        rows.append(np.eye(len(basis))[im]); logk.append(0.0); assigned.append(log_a)
    else:
        names, kinds = [], []
    if include != "minerals":
        for k, lk in zip(ev.dep_index, ev.logk_dep):
            nu = net.nu_dep[k]
            if nu[im] > 0:
                names.append(net.dep_names[k]); kinds.append("aqueous")
                rows.append(nu); logk.append(lk); assigned.append(log_a)
    if include != "aqueous":
        for k, lk in zip(ev.phase_index, ev.logk_phase):
            nu, name = net.nu_phase[k], net.phase_names[k]
            if net.phase_kind[k] == Kind.MINERAL and nu[im] > 0 and name.lower() not in excluded:
                names.append(name); kinds.append("mineral")
                rows.append(nu); logk.append(lk); assigned.append(0.0)
    if not names:
        raise DiagramError(f"nothing of {element} to draw")
    nu = np.array(rows)
    n_el = nu[:, im]
    logk = np.array(logk)
    assigned = np.array(assigned)

    xs = np.linspace(x.lo, x.hi, nx)
    ys = np.linspace(y.lo, y.hi, ny)
    X, Y = np.meshgrid(xs, ys)
    la = {name: np.full(X.shape, v) for name, v in fixed.items()}
    la[WATER] = np.zeros(X.shape)
    la[master] = np.zeros(X.shape)
    ph_grid = np.full(X.shape, float(ph))
    for ax, V in ((x, X), (y, Y)):
        if ax.kind == "pH":
            ph_grid = V
    la[PROTON] = -ph_grid
    k_eh = db.logk_eh(t_c, p)
    if log_fo2 is None and not redox_axes:
        log_fo2 = -0.7 if t_c < 100 else -30.0
    la[OXYGEN] = np.full(X.shape, float(log_fo2 if log_fo2 is not None else 0.0))
    for ax, V in ((x, X), (y, Y)):
        if ax.kind == "Eh":
            la[OXYGEN] = 4.0 * V / _eh_factor(t_c) + k_eh + 4.0 * ph_grid
        elif ax.kind == "log_fO2":
            la[OXYGEN] = V
        elif ax.kind == "log_a":
            la[ax.species] = V
        elif ax.kind == "log_ratio":
            la[ax.species] = V - db.species[ax.species].charge * ph_grid
    stack = np.stack([la[n] for n in net.basis], axis=-1)                 # (ny, nx, nb)
    g = (stack @ nu.T - logk - assigned) / n_el                            # (ny, nx, ncand)
    best = np.argmax(g, axis=-1)
    present = sorted(set(best.ravel().tolist()))
    shown = g[..., present]
    if len(present) > 1:              # the lead over the next field drawn: how sure each boundary is
        top = -np.sort(-shown, axis=-1)
        margin = top[..., 0] - top[..., 1]
    else:
        margin = np.full(best.shape, np.inf)
    winner = np.take_along_axis(g, best[..., None], axis=-1)[..., 0]

    try:
        rel = validity.water_ionization_check(db, t_c, p)
        band = BAND[rel.level]
    except Exception:                                  # outside the liquid field of water
        rel, band = None, BAND["poor"]
    notes = []
    if rel is not None and rel.level != "good":
        notes.append(validity.LEVEL_TEXT[rel.level].format(tol=validity.GOOD, worst=abs(rel.deviation)))
    for k in present:                 # candidates that never win but come within the uncertainty
        where = best == k
        close = [names[c] for c in range(len(names)) if c not in present
                 and np.max(g[..., c][where] - winner[where]) > -band]
        if close:
            notes.append(f"{', '.join(close)} {'comes' if len(close) == 1 else 'come'} within {band:g} "
                         f"log units of {names[k]} in its field: with slightly different data "
                         f"{'it' if len(close) == 1 else 'they'} could take part of it.")
    flagged = sorted(n for n in set(names) & extrapolated)
    if flagged:
        notes.append("Data extrapolated at this temperature for: " + ", ".join(flagged))
    if ev.dropped:
        gone = sorted({n for n in ev.dropped if element in getattr(db.species.get(n), "elements", {})})
        if gone:
            notes.append("No data at this temperature for: " + ", ".join(gone))

    dia = Diagram(element, t_c, p, x, y, xs, ys, names, kinds, best, margin, band, rel,
                  extrapolated & set(names), fixed=fixed, notes=notes)
    dia.water = _water_lines(net, ev, x, y, xs, t_c, k_eh, ph, log_fo2)
    dia.limits = _limits(net, ev, x, y)
    if not redox_axes:
        dia.fixed[OXYGEN] = float(log_fo2)
    if not any(ax.kind == "pH" for ax in (x, y)):
        dia.fixed[PROTON] = -float(ph)
    return dia


def _water_lines(net, ev, x: Axis, y: Axis, xs, t_c, k_eh, ph, log_fo2):
    """Where water is stable: fO2 = 1 above and fH2 = 1 below, on an Eh or log fO2 against pH plane."""
    if x.kind != "pH" or y.kind not in ("Eh", "log_fO2"):
        return {}
    lower_fo2 = None
    if "H2(g)" in net.phase_names:
        k = net.phase_names.index("H2(g)")
        pos = list(ev.phase_index).index(k) if k in ev.phase_index else None
        if pos is not None:
            nu, lk = net.nu_phase[k], ev.logk_phase[pos]
            # H2(g) = H2O - 0.5 O2(g): log fH2 = 0 at log fO2 = -2 logK
            lower_fo2 = -lk / -nu[net.basis.index(OXYGEN)] if nu[net.basis.index(OXYGEN)] else None
    out = {}
    for key, fo2 in (("upper", 0.0), ("lower", lower_fo2)):
        if fo2 is None:
            continue
        if y.kind == "log_fO2":
            out[key] = [(float(v), float(fo2)) for v in (xs[0], xs[-1])]
        else:
            out[key] = [(float(v), float((fo2 - k_eh - 4.0 * v) / 4.0 * _eh_factor(t_c)))
                        for v in (xs[0], xs[-1])]
    return out


def _limits(net, ev, x: Axis, y: Axis) -> list[dict]:
    """Where an axis species saturates with its own mineral (quartz for SiO2): beyond that line the
    fluid is supersaturated, so that part of the diagram cannot be reached at equilibrium."""
    out = []
    for name, ax in (("x", x), ("y", y)):
        if ax.kind != "log_a":
            continue
        j = net.basis.index(ax.species)
        best = None
        for k, lk in zip(ev.phase_index, ev.logk_phase):
            nu = net.nu_phase[k]
            others = [i for i, v in enumerate(nu) if v and i != j and net.basis[i] != WATER]
            if net.phase_kind[k] != Kind.MINERAL or others or nu[j] <= 0:
                continue
            value = lk / nu[j]
            if best is None or value < best[0]:
                best = (float(value), net.phase_names[k])
        if best is not None and ax.lo < best[0] < ax.hi:
            out.append({"axis": name, "value": best[0], "mineral": best[1]})
    return out


def label_points(dia: Diagram) -> dict[str, tuple[float, float]]:
    """A place inside each field for its label: the cell nearest the field's centroid."""
    out = {}
    X, Y = np.meshgrid(np.arange(len(dia.xs)), np.arange(len(dia.ys)))
    for k in sorted(set(dia.grid.ravel().tolist())):
        mask = dia.grid == k
        cx, cy = X[mask].mean(), Y[mask].mean()
        d = (X - cx) ** 2 + (Y - cy) ** 2
        d[~mask] = np.inf
        i, j = np.unravel_index(int(np.argmin(d)), d.shape)
        out[dia.names[k]] = (float(dia.xs[j]), float(dia.ys[i]))
    return out
