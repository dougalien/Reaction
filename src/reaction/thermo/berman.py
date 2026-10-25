"""Minerals from the internally consistent dataset of Berman (1988) and its later additions.

Equations as implemented in CHNOSZ (Berman.R): heat capacity with k0..k6 terms, a volume
that is quadratic in T and P, lambda transitions (Eqs. 8-12) and order-disorder (Eqs. 15-20).
The parameter files are CHNOSZ's extdata/Berman/*.csv (GPL-3); when a mineral appears in
several files the youngest file wins, as in CHNOSZ.

Gibbs energies are returned on the Benson-Helgeson scale used by OBIGT and HKF
(apparent G of formation from the elements), so they can be mixed with aqueous species:

    G = Ga + Tr * S(elements),   Ga = H(T,P) - T S(T,P)   (Berman-Brown convention)

As in CHNOSZ, the tabulated G at 298 K is not used: G comes from H and S with the element
entropies, so it is consistent with the enthalpies and entropies of the dataset.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

T_REF, P_REF = 298.15, 1.0
CAL = 4.184

# Entropies of the elements at 298.15 K and 1 bar (Cox, Wagman and Medvedev, 1989; Fe from
# Wagman et al., 1982), cal/K per atom, as in CHNOSZ thermo/element.csv
_S_ELEMENT_CAL = {
    "O": 49.0325047801147 / 2, "H": 31.2332695984704 / 2, "F": 48.4682122370937 / 2,
    "Cl": 53.317638623327 / 2, "C": 1.37189292543021, "Si": 4.49569789674952,
    "Al": 6.76386233269599, "Ti": 7.34225621414914, "P": 9.82074569789675,
    "Mg": 7.80831739961759, "Ca": 9.94024856596558, "Na": 12.2609942638623,
    "K": 15.4588910133843, "Fe": 6.52007648183556, "Mn": 7.65057361376673,
    "Li": 6.95984703632887,
}

_V_SCALE = {"v1": 1e-5, "v2": 1e-5, "v3": 1e-5, "v4": 1e-8}


def element_entropy(elements: dict[str, float]) -> float:
    """Entropy of the elements in a formula at Tr, J/K."""
    try:
        return sum(n * _S_ELEMENT_CAL[el] for el, n in elements.items()) * CAL
    except KeyError as exc:
        raise ValueError(f"no entropy for element {exc.args[0]}") from None


@dataclass(frozen=True)
class BermanParams:
    name: str
    g0: float | None        # J/mol (Benson-Helgeson), None when not tabulated
    h0: float               # J/mol
    s0: float               # J/K/mol
    v0: float               # J/bar
    k: tuple[float, ...]    # k0..k6
    v: tuple[float, ...]    # v1..v4, unscaled
    lam: tuple | None       # (Tlambda, Tref, dTdP, l1, l2)
    dis: tuple | None       # (Tmin, Tmax, d0, d1, d2, d3, d4, Vad)
    source: str


def _f(row, key):
    t = (row.get(key) or "").strip()
    return None if t in ("", "NA") else float(t)


def load_parameters(folder: str | Path) -> dict[str, BermanParams]:
    """Read every csv in ``folder``; for a mineral in several files the youngest file wins."""
    folder = Path(folder)
    files = sorted(folder.glob("*.csv"), key=lambda p: (p.stem.rsplit("_", 1)[-1], p.stem))
    params: dict[str, BermanParams] = {}
    for path in files:                                   # oldest first, so later rows override
        with path.open(encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh):
                lam = None
                if _f(row, "Tlambda") is not None and _f(row, "Tref") is not None:
                    lam = tuple(_f(row, k) or 0.0 for k in ("Tlambda", "Tref", "dTdP", "l1", "l2"))
                dis = None
                if _f(row, "Tmin") is not None and _f(row, "Tmax") is not None:
                    dis = tuple(_f(row, k) or 0.0 for k in ("Tmin", "Tmax", "d0", "d1", "d2",
                                                            "d3", "d4", "Vad"))
                params[row["name"]] = BermanParams(
                    name=row["name"], g0=_f(row, "GfPrTr"), h0=_f(row, "HfPrTr"),
                    s0=_f(row, "SPrTr"), v0=_f(row, "VPrTr"),
                    k=tuple(_f(row, f"k{i}") or 0.0 for i in range(7)),
                    v=tuple((_f(row, k) or 0.0) * _V_SCALE[k] for k in ("v1", "v2", "v3", "v4")),
                    lam=lam, dis=dis, source=path.name)
    return params


@dataclass(frozen=True)
class BermanSpecies:
    """A mineral evaluated with the Berman equations; ``gibbs`` matches CrystalSpecies."""
    name: str
    state: str
    params: BermanParams
    s_elements: float        # J/K

    @property
    def v0(self) -> float:   # cm3/mol, like CrystalSpecies
        return self.params.v0 * 10.0

    @property
    def g_ref(self) -> float:
        return self.params.h0 - T_REF * (self.params.s0 - self.s_elements)

    def extrapolated(self, t_k: float) -> bool:
        return False

    def properties(self, t_k: float, p_bar: float) -> tuple[float, float, float]:
        """(G, H, S) at T and P: G on the Benson-Helgeson scale, J and J/K."""
        p = self.params
        k0, k1, k2, k3, k4, k5, k6 = p.k
        v1, v2, v3, v4 = p.v
        v0 = p.v0
        T, Tr, P, Pr = float(t_k), T_REF, float(p_bar), P_REF
        int_cp = (k0 * (T - Tr) + 2.0 * k1 * (T ** 0.5 - Tr ** 0.5) - k2 * (1 / T - 1 / Tr)
                  - k3 / 2.0 * (T ** -2 - Tr ** -2) + k4 * math.log(T / Tr)
                  + k5 / 2.0 * (T ** 2 - Tr ** 2) + k6 / 3.0 * (T ** 3 - Tr ** 3))
        int_cp_t = (k0 * math.log(T / Tr) - 2.0 * k1 * (T ** -0.5 - Tr ** -0.5)
                    - k2 / 2.0 * (T ** -2 - Tr ** -2) - k3 / 3.0 * (T ** -3 - Tr ** -3)
                    - k4 * (1 / T - 1 / Tr) + k5 * (T - Tr) + k6 / 2.0 * (T ** 2 - Tr ** 2))
        dt, dp = T - Tr, P - Pr
        # int V dP from Pr to P, and its temperature derivative
        int_vdp = v0 * ((1 + v1 * dt + v2 * dt ** 2) * dp + v3 / 2 * dp ** 2 + v4 / 3 * dp ** 3)
        d_int_vdp_dt = v0 * (v1 + 2 * v2 * dt) * dp
        h = p.h0 + int_cp + int_vdp - T * d_int_vdp_dt
        s = p.s0 + int_cp_t - d_int_vdp_dt
        if p.lam is not None and T > p.lam[1]:
            htr, str_ = self._transition(T, P)
            h, s = h + htr, s + str_
        if p.dis is not None and T > p.dis[0]:
            hds, sds = self._disorder(T, P)
            h, s = h + hds, s + sds
        ga = h - T * s
        return ga + Tr * self.s_elements, h, s

    def gibbs(self, t_k: float, p_bar: float) -> float:
        return self.properties(t_k, p_bar)[0]

    def _transition(self, T, P):
        t_lambda, t_ref, dtdp, l1, l2 = self.params.lam
        t_lambda_p = t_lambda + dtdp * (P - 1.0)
        td = t_lambda - t_lambda_p
        ttr = min(T, t_lambda_p)
        tref = t_ref - td
        x1 = l1 ** 2 * td + 2 * l1 * l2 * td ** 2 + l2 ** 2 * td ** 3
        x2 = l1 ** 2 + 4 * l1 * l2 * td + 3 * l2 ** 2 * td ** 2
        x3 = 2 * l1 * l2 + 3 * l2 ** 2 * td
        x4 = l2 ** 2
        htr = (x1 * (ttr - tref) + x2 / 2 * (ttr ** 2 - tref ** 2) + x3 / 3 * (ttr ** 3 - tref ** 3)
               + x4 / 4 * (ttr ** 4 - tref ** 4))
        str_ = (x1 * math.log(ttr / tref) + x2 * (ttr - tref) + x3 / 2 * (ttr ** 2 - tref ** 2)
                + x4 / 3 * (ttr ** 3 - tref ** 3))
        return htr, str_

    def _disorder(self, T, P):
        t_min, t_max, d0, d1, d2, d3, d4, vad = self.params.dis
        tds = min(T, t_max)
        hds = (d0 * (tds - t_min) + 2 * d1 * (tds ** 0.5 - t_min ** 0.5) - d2 * (1 / tds - 1 / t_min)
               + d3 / 2 * (tds ** 2 - t_min ** 2) + d4 / 3 * (tds ** 3 - t_min ** 3))
        sds = (d0 * math.log(tds / t_min) - 2 * d1 * (tds ** -0.5 - t_min ** -0.5)
               - d2 / 2 * (tds ** -2 - t_min ** -2) + d3 * (tds - t_min) + d4 / 2 * (tds ** 2 - t_min ** 2))
        if vad != 0:
            hds += hds / vad * (P - P_REF)
        if T > t_max:
            hds -= (T - t_max) * sds
        return hds, sds
