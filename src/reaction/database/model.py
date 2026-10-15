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
from scipy.interpolate import PchipInterpolator


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

    Evaluation passes exactly through the tabulated points (monotone cubic
    between them) and refuses to extrapolate.
    """

    temps_c: np.ndarray
    values: np.ndarray
    _interp: Any = field(default=None, repr=False)

    def __post_init__(self) -> None:
        self.temps_c = np.asarray(self.temps_c, dtype=float)
        self.values = np.asarray(self.values, dtype=float)
        ok = np.isfinite(self.values)
        self._t = self.temps_c[ok]
        self._v = self.values[ok]
        if len(self._t) >= 3:
            self._interp = PchipInterpolator(self._t, self._v, extrapolate=False)

    @property
    def available(self) -> bool:
        return len(self._t) > 0

    @property
    def t_range(self) -> tuple[float, float]:
        return float(self._t[0]), float(self._t[-1])

    def __call__(self, t_c: float, p_bar: float | None = None) -> float:
        if not self.available:
            raise ThermoRangeError("no logK data")
        lo, hi = self.t_range
        if t_c < lo - 1e-9 or t_c > hi + 1e-9:
            raise ThermoRangeError(f"T={t_c} degC outside logK data range [{lo}, {hi}]")
        for t, v in zip(self._t, self._v):
            if abs(t - t_c) < 1e-9:
                return float(v)
        if len(self._t) == 1:
            raise ThermoRangeError("single-point logK cannot be interpolated")
        if len(self._t) == 2:
            return float(np.interp(t_c, self._t, self._v))
        return float(self._interp(t_c))


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

    def __getitem__(self, name: str) -> Species:
        return self.species[name]

    def __contains__(self, name: str) -> bool:
        return name in self.species

    def __len__(self) -> int:
        return len(self.species)

    def of_kind(self, *kinds: Kind) -> list[Species]:
        return [s for s in self.species.values() if s.kind in kinds]

    def _grid(self, values: np.ndarray, t_c: float) -> float:
        return GridLogK(self.temps_c, values)(t_c)

    def dh_params(self, t_c: float) -> tuple[float, float, float]:
        """Return (A, B, bdot) at t_c."""
        return (self._grid(self.adh, t_c), self._grid(self.bdh, t_c),
                self._grid(self.bdot, t_c))

    def molecular_weight(self, name: str) -> float:
        """Molecular weight from element weights (g/mol)."""
        sp = self.species[name]
        return sum(n * self.elements[el] for el, n in sp.elements.items())
