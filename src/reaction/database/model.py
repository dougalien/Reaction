"""Normalized thermodynamic-database model.

Every datum keeps its provenance, so a number in a result can be traced back
to the file and literature source it came from.  Reactions follow the EQ3/6
"dissociation" convention: the species itself has coefficient -1 and the
other participants are written in terms of other species, so that

    log a(species) = sum_k c_k log a_k  -  log K
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

import numpy as np
from scipy.interpolate import BarycentricInterpolator


class Kind(str, Enum):
    BASIS = "basis"
    AUX_BASIS = "aux_basis"
    AQUEOUS = "aqueous"
    MINERAL = "mineral"
    GAS = "gas"
    LIQUID = "liquid"


class ThermoRangeError(ValueError):
    """Requested temperature/pressure lies outside the range of the data."""


@dataclass(eq=False)
class GridLogK:
    """log10 K tabulated on a temperature grid (NaN where the datum is absent).

    Interpolation follows EQPT (EQ3/6): one exact polynomial through the valid
    points of the 0-100 C part of the grid and another through the valid
    points of the 100-300 C part, so continuity at 100 C holds when that grid
    point is valid.  A part with a single valid point is a constant; a part
    with none suppresses the species there (ThermoRangeError).
    """

    temps_c: np.ndarray
    values: np.ndarray
    t_min: float = 0.0
    t_max: float = 300.0
    _parts: Any = field(default=None, repr=False)

    SPLIT_C = 100.0

    def __post_init__(self) -> None:
        self.temps_c = np.asarray(self.temps_c, dtype=float)
        self.values = np.asarray(self.values, dtype=float)
        ok = np.isfinite(self.values)
        low = ok & (self.temps_c <= self.SPLIT_C)
        high = ok & (self.temps_c >= self.SPLIT_C)
        self._parts = []
        self._valid = []                 # (lowest, highest) valid grid T per part
        for mask in (low, high):
            t, v = self.temps_c[mask], self.values[mask]
            self._valid.append((float(t.min()), float(t.max())) if len(t) else None)
            if len(t) >= 2:
                self._parts.append(BarycentricInterpolator(t, v))
            elif len(t) == 1:
                self._parts.append(float(v[0]))
            else:
                self._parts.append(None)

    @property
    def available(self) -> bool:
        return any(p is not None for p in self._parts)

    def extrapolated(self, t_c: float) -> bool:
        """True if t_c lies outside the span of valid grid points in its part.

        EQ3/6 still evaluates its polynomial there (a lone 25 C value becomes
        a constant), so Reaction does too, but callers can flag the result.
        """
        span = self._valid[0 if t_c <= self.SPLIT_C else 1]
        if span is None:
            return True
        return t_c < span[0] - 1e-9 or t_c > span[1] + 1e-9

    def __call__(self, t_c: float, p_bar: float | None = None) -> float:
        if t_c < self.t_min - 1e-9 or t_c > self.t_max + 1e-9:
            raise ThermoRangeError(f"T={t_c} degC outside the data range "
                                   f"[{self.t_min}, {self.t_max}]")
        part = self._parts[0 if t_c <= self.SPLIT_C else 1]
        if part is None:
            raise ThermoRangeError(f"no logK data at T={t_c} degC")
        if isinstance(part, float):
            return part
        return float(part(t_c))


@dataclass
class Provenance:
    database: str
    file: str
    section: str
    revised: str | None = None
    sources: list[str] = field(default_factory=list)
    data_sets: list[str] = field(default_factory=list)


@dataclass
class Species:
    name: str
    kind: Kind
    formula: str | None
    charge: float
    elements: dict[str, float]
    mol_wt: float | None
    a0: float | None                 # ion-size parameter, angstrom
    activity_type: int               # 0 ion / polar neutral, -1 nonpolar neutral
    reaction: dict[str, float]       # empty for strict basis species
    logk: GridLogK | None
    molar_volume: float | None       # cm3/mol
    ref_state: dict[str, float]      # dG_f, dH_f (J/mol), S (J/mol/K), Cp coefs (J/mol/K)
    provenance: Provenance
    subtype: str = ""


@dataclass
class ThermoDatabase:
    title: str
    source: str
    temps_c: np.ndarray
    pressures_bar: np.ndarray
    adh: np.ndarray                  # Debye-Hueckel A (log10 form)
    bdh: np.ndarray                  # Debye-Hueckel B, (kg/mol)^0.5 / angstrom
    bdot: np.ndarray
    cco2: np.ndarray                 # Drummond (1981) polynomial coefficients
    logk_eh: GridLogK                # 2H2O = O2(g) + 4H+ + 4e-
    elements: dict[str, float]
    species: dict[str, Species] = field(default_factory=dict)
    skipped: list[str] = field(default_factory=list)   # import decisions worth auditing
    t_limits: tuple[float, float] = (0.0, 300.0)       # degC range the file claims

    def __getitem__(self, name: str) -> Species:
        return self.species[name]

    def __contains__(self, name: str) -> bool:
        return name in self.species

    def __len__(self) -> int:
        return len(self.species)

    def of_kind(self, *kinds: Kind) -> list[Species]:
        return [s for s in self.species.values() if s.kind in kinds]

    def _grid(self, values: np.ndarray, t_c: float) -> float:
        return GridLogK(self.temps_c, values, *self.t_limits)(t_c)

    def dh_params(self, t_c: float) -> tuple[float, float, float]:
        """Return (A, B, bdot) at t_c."""
        return (self._grid(self.adh, t_c), self._grid(self.bdh, t_c),
                self._grid(self.bdot, t_c))

    def molecular_weight(self, name: str) -> float:
        """Molecular weight from element weights (g/mol)."""
        sp = self.species[name]
        return sum(n * self.elements[el] for el, n in sp.elements.items())
