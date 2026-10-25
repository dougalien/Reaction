"""Pitzer's ion-interaction model for activity coefficients and the activity of water in brines.

The equations are those of Harvie, Moller and Weare (1984, GCA 48, 723), with unsymmetrical mixing
(E-theta) from Pitzer's (1975) approximation of the J integral.  Parameters are PHREEQC's pitzer.dat
(public domain; see ``data/pitzer_phreeqc.tsv`` for its sources), with their temperature functions.
They are fitted at the saturation pressure of water and are used here at the path's pressure without
a pressure term: the pressure effect on the equilibrium constants themselves comes from HKF.

In this model ions are mostly free: NaCl(aq), CaSO4(aq), NaSO4- and the like are not separate species,
their effect is in the interaction parameters.  ``keep_species`` says which dependent species of the
Pitzer system stay; species of other elements (Al, Fe(III) ...) keep their complexes and get only the
Debye-Hueckel part (z^2 F) of the Pitzer equations.
"""
from __future__ import annotations

import math
from functools import lru_cache
from pathlib import Path

import numpy as np

from .constants import WATER_MOLES_PER_KG

DATA = Path(__file__).with_name("data") / "pitzer_phreeqc.tsv"
T_REF = 298.15
B_DH = 1.2                          # the universal b of the Pitzer equations, kg^1/2 mol^-1/2
T_FITTED_C = 200.0                  # most temperature terms are fitted to about this temperature

# basis species that make up the Pitzer system, and the dependent species it keeps
COMPONENTS = {"H2O", "H+", "Na+", "K+", "Li+", "Mg+2", "Ca+2", "Sr+2", "Ba+2", "Mn+2", "Fe+2", "Cl-",
              "Br-", "SO4-2", "HCO3-", "B(OH)3", "SiO2"}
KEEP = {"OH-", "CO3-2", "CO2", "HSO4-", "BO2-", "MgOH+", "MgCO3", "HSiO3-"}


@lru_cache(maxsize=1)
def table() -> list[tuple[str, tuple[str, ...], tuple[float, ...]]]:
    rows = []
    for line in DATA.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("type\t"):
            continue
        cols = line.split("\t")
        rows.append((cols[0], tuple(cols[1].split(",")), tuple(float(v) for v in cols[2:8])))
    return rows


def at_temperature(a: tuple[float, ...], t_k: float) -> float:
    a0, a1, a2, a3, a4, a5 = a
    return (a0 + a1 * (1.0 / t_k - 1.0 / T_REF) + a2 * math.log(t_k / T_REF) + a3 * (t_k - T_REF)
            + a4 * (t_k * t_k - T_REF * T_REF) + a5 * (1.0 / (t_k * t_k) - 1.0 / (T_REF * T_REF)))


def _j(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Pitzer's (1975) approximation of J(x) and its derivative dJ/dx."""
    c1, c2, c3, c4 = 4.581, 0.7237, 0.0120, 0.528
    x = np.maximum(x, 1e-30)
    e = np.exp(-c3 * x ** c4)
    d = 4.0 + c1 * x ** (-c2) * e
    j = x / d
    dd = c1 * e * (-c2 * x ** (-c2 - 1.0) - c3 * c4 * x ** (c4 - c2 - 1.0))
    return j, (d - x * dd) / (d * d)


def _g(x):
    return 2.0 * (1.0 - (1.0 + x) * np.exp(-x)) / (x * x)


def _gp(x):
    return -2.0 * (1.0 - (1.0 + x + 0.5 * x * x) * np.exp(-x)) / (x * x)


class Pitzer:
    """Activity coefficients of ``names`` (charges ``z``) at one temperature."""

    def __init__(self, names: list[str], z: np.ndarray, temperature_c: float, a_phi: float):
        self.names, self.z = list(names), np.asarray(z, dtype=float)
        self.a_phi = a_phi
        t_k = temperature_c + 273.15
        idx = {n: i for i, n in enumerate(self.names)}
        n = len(self.names)
        self.cat = np.where(self.z > 0)[0]
        self.an = np.where(self.z < 0)[0]
        self.neu = np.where(self.z == 0)[0]
        b0, b1, b2, cphi = (np.zeros((n, n)) for _ in range(4))
        theta, lam = np.zeros((n, n)), np.zeros((n, n))
        psi, zeta = np.zeros((n, n, n)), np.zeros((n, n, n))
        self.with_params: set[str] = set()
        for kind, sp, a in table():
            if not all(s in idx for s in sp):
                continue
            v = at_temperature(a, t_k)
            ii = [idx[s] for s in sp]
            self.with_params.update(sp)
            if kind in ("B0", "B1", "B2", "C0"):
                i, j = ii
                target = {"B0": b0, "B1": b1, "B2": b2, "C0": cphi}[kind]
                target[i, j] = target[j, i] = v
            elif kind == "THETA":
                i, j = ii
                theta[i, j] = theta[j, i] = v
            elif kind == "LAMBDA":
                i, j = ii
                lam[i, j] = lam[j, i] = v
            elif kind == "PSI":
                for p in ((0, 1, 2), (1, 0, 2), (0, 2, 1), (2, 0, 1), (1, 2, 0), (2, 1, 0)):
                    psi[ii[p[0]], ii[p[1]], ii[p[2]]] = v
            elif kind == "ZETA":
                for p in ((0, 1, 2), (0, 2, 1)):           # neutral first, the ion pair either way
                    zeta[ii[p[0]], ii[p[1]], ii[p[2]]] = v
        az = np.abs(self.z)
        both2 = np.outer(az >= 2, az >= 2)
        self.alpha1 = np.where(both2, 1.4, 2.0)
        self.alpha2 = np.where(both2 & (np.outer(az, az) == 4), 12.0, 50.0)
        self.b0, self.b1, self.b2 = b0, b1, b2
        zz = np.sqrt(np.maximum(np.outer(az, az), 1e-30))
        self.c = cphi / (2.0 * zz)
        self.theta, self.lam, self.psi, self.zeta = theta, lam, psi, zeta
        same = np.outer(self.z, self.z) > 0
        self.unsym = same & (np.outer(az, np.ones(n)) != np.outer(np.ones(n), az))

    def _etheta(self, ionic: float):
        """E-theta and its derivative for pairs of ions of the same sign and different charge."""
        n = len(self.names)
        et, etp = np.zeros((n, n)), np.zeros((n, n))
        if ionic <= 0 or not self.unsym.any():
            return et, etp
        zz = np.outer(self.z, self.z)
        x = 6.0 * zz * self.a_phi * math.sqrt(ionic)
        jx, jpx = _j(x)
        dj = np.diag(jx)
        xd = np.diag(x)
        full = zz / (4.0 * ionic) * (jx - 0.5 * dj[:, None] - 0.5 * dj[None, :])
        fullp = -full / ionic + zz / (8.0 * ionic ** 2) * (
            x * jpx - 0.5 * (xd * np.diag(jpx))[:, None] - 0.5 * (xd * np.diag(jpx))[None, :])
        et = np.where(self.unsym, full, 0.0)
        etp = np.where(self.unsym, fullp, 0.0)
        return et, etp

    def compute(self, m: np.ndarray) -> tuple[np.ndarray, float]:
        """(ln gamma of every species, ln a(H2O)) at molalities ``m``."""
        m = np.maximum(np.asarray(m, dtype=float), 0.0)
        z, az = self.z, np.abs(self.z)
        ionic = 0.5 * float(np.sum(m * z * z))
        big_z = float(np.sum(m * az))
        s = math.sqrt(ionic)
        cat, an, neu = self.cat, self.an, self.neu
        mc, ma, mn = m[cat], m[an], m[neu]
        if s > 0:
            x1, x2 = self.alpha1 * s, self.alpha2 * s
            bmx = self.b0 + self.b1 * _g(x1) + self.b2 * _g(x2)
            bp = (self.b1 * _gp(x1) + self.b2 * _gp(x2)) / ionic
            bphi = self.b0 + self.b1 * np.exp(-x1) + self.b2 * np.exp(-x2)
        else:
            bmx, bp, bphi = self.b0.copy(), np.zeros_like(self.b0), self.b0 + self.b1 + self.b2
        et, etp = self._etheta(ionic)
        phi = self.theta + et
        phip = etp
        phiphi = self.theta + et + ionic * etp

        B_ca, Bp_ca, C_ca = bmx[np.ix_(cat, an)], bp[np.ix_(cat, an)], self.c[np.ix_(cat, an)]
        f = -self.a_phi * (s / (1.0 + B_DH * s) + 2.0 / B_DH * math.log(1.0 + B_DH * s))
        f += mc @ Bp_ca @ ma
        f += 0.5 * mc @ phip[np.ix_(cat, cat)] @ mc + 0.5 * ma @ phip[np.ix_(an, an)] @ ma
        sum_mcmac = mc @ C_ca @ ma
        ln_g = z * z * f + az * sum_mcmac
        if len(cat) and len(an):
            psi_cca = self.psi[np.ix_(cat, cat, an)]
            psi_aac = self.psi[np.ix_(an, an, cat)]
            ln_g[cat] += (bmx[np.ix_(cat, an)] * 2.0 + big_z * C_ca) @ ma
            ln_g[cat] += 2.0 * phi[np.ix_(cat, cat)] @ mc + np.einsum("ijk,j,k->i", psi_cca, mc, ma)
            ln_g[cat] += 0.5 * np.einsum("jki,j,k->i", psi_aac, ma, ma)
            ln_g[an] += (bmx[np.ix_(cat, an)] * 2.0 + big_z * C_ca).T @ mc
            ln_g[an] += 2.0 * phi[np.ix_(an, an)] @ ma + np.einsum("ijk,j,k->i", psi_aac, ma, mc)
            ln_g[an] += 0.5 * np.einsum("jki,j,k->i", psi_cca, mc, mc)
        lam = self.lam
        ln_g += 2.0 * lam[:, neu] @ mn                         # every species with the neutrals
        ln_g[neu] += 2.0 * lam[np.ix_(neu, cat)] @ mc + 2.0 * lam[np.ix_(neu, an)] @ ma
        zeta = self.zeta
        if len(neu):
            ln_g[neu] += np.einsum("ijk,j,k->i", zeta[np.ix_(neu, cat, an)], mc, ma)
            ln_g[cat] += np.einsum("nik,n,k->i", zeta[np.ix_(neu, cat, an)], mn, ma)
            ln_g[an] += np.einsum("nki,n,k->i", zeta[np.ix_(neu, cat, an)], mn, mc)

        # osmotic coefficient (Gibbs-Duhem consistent with the above)
        osm = -self.a_phi * ionic * s / (1.0 + B_DH * s)
        osm += mc @ (bphi[np.ix_(cat, an)] + big_z * C_ca) @ ma
        if len(cat) and len(an):
            osm += 0.5 * mc @ phiphi[np.ix_(cat, cat)] @ mc + 0.5 * np.einsum("i,j,k,ijk->", mc, mc, ma, psi_cca)
            osm += 0.5 * ma @ phiphi[np.ix_(an, an)] @ ma + 0.5 * np.einsum("i,j,k,ijk->", ma, ma, mc, psi_aac)
        osm += mn @ lam[np.ix_(neu, cat)] @ mc + mn @ lam[np.ix_(neu, an)] @ ma
        osm += 0.5 * mn @ lam[np.ix_(neu, neu)] @ mn
        if len(neu):
            osm += np.einsum("n,c,a,nca->", mn, mc, ma, zeta[np.ix_(neu, cat, an)])
        total = float(np.sum(m))
        ln_aw = -2.0 * osm / WATER_MOLES_PER_KG - total / WATER_MOLES_PER_KG
        return ln_g, ln_aw

    def osmotic(self, m: np.ndarray) -> float:
        total = float(np.sum(m))
        if total <= 0:
            return 1.0
        return -self.compute(m)[1] * WATER_MOLES_PER_KG / total
