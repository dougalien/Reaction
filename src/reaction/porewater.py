"""Steady pore water of a carbonate sediment: alkalinity and DIC with depth below the seafloor.

The model of Hales and Emerson (1997, GCA 61, 501-514; after Hales and Emerson 1996 and Emerson and Bender
1981, J. Mar. Res. 39, 139): one-dimensional, steady diffusion and reaction, in which calcite dissolves at

    R_cal = k (1 - phi) [C] (1 - Omega)^n          (their Eqn. 4; n = 4.5 there, after Keir 1980)

k a fraction of the calcite per day, [C] the calcite of the solids, and oxygen consumed by oxic respiration
makes beta CO2 per O2 (their Eqn. 3; beta 0.77 Redfield, 0.64 their middle case, 0.55 alkane-like).  Here the
O2 consumption is one exponential, scale L, that integrates to the O2 flux into the sediment (theirs has two
terms, about 70% of it with a 3 mm scale: their Eqn. 2).  Per kg of pore water, z down from the seafloor:

    D_s d2(DIC)/dz2 + R_cal + R_CO2 = 0
    D_s d2(TA)/dz2 + 2 R_cal - (17/106) R_CO2 = 0

with the bottom water's TA and DIC at z = 0 and no flux at the base.  The 17/106 is the alkalinity lost to the
nitrate and phosphate of Redfield organic matter oxidized (Wolf-Gladrow et al. 2007, Mar. Chem. 106, 287-300:
from memory, not checked).  Omega is the marine carbonate system of marine.py (Lueker et al. 2000 K1 and K2,
Mucci 1983 calcite solubility, Millero 1995 pressure terms); Hales and Emerson used the UNESCO (1987) solubility
and found that their Station G data need one at least 5% below it.

Diffusion: D_s = D0 / theta^2, theta^2 = 1 - ln(phi^2) (Boudreau 1996, GCA 60, 3139), D0 of HCO3- =
(5.06 + 0.275 T) 1e-6 cm2/s, T in C (Boudreau 1997, Diagenetic Models and Their Implementation, Table 4.8),
for both solutes.  Both of these from memory, not checked against the sources.  Hales and Emerson instead
measured the formation factor F (2.5-3 by 4-5 cm at the Ceara Rise) with phi = F^(-1/b), b about 2.3 (their
Eqn. 1): theta^2 = phi F, about 1.8 there, close to Boudreau's 1.86 at phi 0.65.

Left out: bioturbation and bioirrigation, suboxic respiration (O2 is not tracked: the respiration does not
stop when it runs out; Hales and Emerson find the Ceara Rise pore waters oxic for several cm), the calcite lost
from the column, aragonite in the sediment, CO2 hydration kinetics, and precipitation.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

from . import marine

CALCITE_G_MOL = 100.087
REDFIELD_ALK_PER_C = 17.0 / 106.0
SMOOTH_OMEGA = 1e-7


def d0_hco3_m2_s(t_c: float) -> float:
    """Tracer diffusion coefficient of HCO3- in water, m2/s (Boudreau 1997, Table 4.8: from memory)."""
    return (5.06 + 0.275 * t_c) * 1e-6 * 1e-4


def d_sediment_m2_s(t_c: float, porosity: float) -> float:
    """Bulk-sediment diffusion coefficient, D0 / (1 - ln phi^2) (Boudreau 1996)."""
    return d0_hco3_m2_s(t_c) / (1.0 - math.log(porosity ** 2))


@dataclass
class Sediment:
    """A carbonate sediment in the form of Hales and Emerson (1997, GCA 61, 501-514, Eqns. 2-4): calcite
    dissolves at R = k (1 - phi) [C] (1 - Omega)^n, k a fraction of the calcite per day; oxygen is consumed at
    a rate that integrates to the O2 flux into the sediment, falling off exponentially with depth, and makes
    beta CO2 per O2 (0.77 Redfield, 0.64 their middle case, 0.55 alkane-like)."""
    porosity: float                 # volume fraction of pore water
    caco3_wt_pct: float             # calcite, weight % of the dry solids
    k_pct_day: float                # dissolution rate constant, % of the calcite per day
    n: float                        # reaction order in (1 - Omega)
    o2_flux_umol_cm2_yr: float      # O2 consumed in the sediment, umol per cm2 of seafloor per year
    beta: float                     # CO2 made per O2 consumed
    resp_depth_cm: float            # e-folding depth of the O2 consumption, cm
    solid_density: float = 2.70     # g/cm3 of the solids (calcite 2.71)
    water_density: float = 1.025    # g/cm3 of the pore water

    def calcite_mol_per_kg(self) -> float:
        """Moles of sediment calcite per kg of pore water."""
        solids_g_per_cm3_water = (1.0 - self.porosity) / self.porosity * self.solid_density
        return solids_g_per_cm3_water * self.caco3_wt_pct / 100.0 / self.water_density * 1000.0 / CALCITE_G_MOL

    def k_per_s(self) -> float:
        return self.k_pct_day / 100.0 / 86400.0

    def resp_at_top(self) -> float:
        """CO2 made at the seafloor, mol per kg of pore water per s: beta x the O2 flux spread over the pore
        water of the exponential's depth scale (phi rho_w L)."""
        per_cm3_bulk_yr = self.beta * self.o2_flux_umol_cm2_yr * 1e-6 / self.resp_depth_cm   # mol cm-3 yr-1
        return per_cm3_bulk_yr / (self.porosity * self.water_density * 1e-3) / (365.25 * 86400.0)

    def problems(self) -> list[str]:
        out = []
        if not 0.05 <= self.porosity <= 0.98:
            out.append("porosity must be between 0.05 and 0.98")
        if not 0.0 <= self.caco3_wt_pct <= 100.0:
            out.append("CaCO3 must be 0 to 100 weight %")
        if not self.k_pct_day > 0.0:
            out.append("the sediment's calcite rate constant must be above 0")
        if not 0.0 < self.n <= 10.0:
            out.append("the reaction order must be above 0 and at most 10")
        if self.o2_flux_umol_cm2_yr < 0.0 or not 0.0 <= self.beta <= 1.5 or not self.resp_depth_cm > 0.0:
            out.append("the O2 flux must be 0 or more, CO2 per O2 0 to 1.5, the respiration depth above 0")
        return out


@dataclass
class Column:
    z_cm: np.ndarray
    ta: np.ndarray                  # mol/kg
    dic: np.ndarray
    omega: np.ndarray               # calcite
    r_cal: np.ndarray               # mol calcite per kg pore water per s
    flux_umol_cm2_yr: float         # calcite dissolved in the column, umol per cm2 of seafloor per year
    iterations: int

    def rows(self) -> list[list[float]]:
        """[z cm, TA umol/kg, DIC umol/kg] rows, as a card's porewater_profile."""
        return [[float(z), float(a * 1e6), float(d * 1e6)] for z, a, d in zip(self.z_cm, self.ta, self.dic)]

    def at(self, z_cm: float) -> tuple[float, float]:
        """(TA, DIC) mol/kg at this depth, linear between the nodes and held below the base."""
        return float(np.interp(z_cm, self.z_cm, self.ta)), float(np.interp(z_cm, self.z_cm, self.dic))


def _grid(depth_cm: float, first_cm: float, nodes: int) -> np.ndarray:
    """Nodes from the seafloor down, fine at the top (where the bottom water meets the calcite)."""
    ratio = 1.0
    if nodes > 2:                    # geometric spacing: first gap first_cm, all gaps summing to depth_cm
        lo, hi = 1.0, 2.0
        while first_cm * (hi ** (nodes - 1) - 1) / (hi - 1) < depth_cm:
            hi *= 2.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            total = first_cm * (nodes - 1) if mid == 1.0 else first_cm * (mid ** (nodes - 1) - 1) / (mid - 1)
            lo, hi = (mid, hi) if total < depth_cm else (lo, mid)
        ratio = 0.5 * (lo + hi)
    gaps = first_cm * ratio ** np.arange(nodes - 1)
    gaps *= depth_cm / gaps.sum()
    return np.concatenate([[0.0], np.cumsum(gaps)])


def solve_column(sed: Sediment, salinity: float, t_c: float, p_bar: float, ta_bw: float, dic_bw: float,
                 depth_cm: float = 50.0, nodes: int = 160, first_cm: float = 1e-3, tol: float = 1e-10,
                 max_iter: int = 200) -> Column:
    """The steady pore-water profile (see the module docstring).  ``ta_bw`` and ``dic_bw`` are the bottom
    water's, mol/kg.  Newton's method on the finite-volume equations, from the bottom water everywhere."""
    c = marine.constants(salinity, t_c, p_bar)

    def omega(ta, dic):
        # the alkalinity equation of marine.solve, vectorised: bisection in log10 H+
        lo = np.full_like(ta, -11.0)
        hi = np.full_like(ta, -5.0)
        for _ in range(64):
            mid = 0.5 * (lo + hi)
            h = 10.0 ** mid
            den = h * h + c.k1 * h + c.k1 * c.k2
            ex = dic * (c.k1 * h + 2.0 * c.k1 * c.k2) / den + c.total_boron * c.kb / (c.kb + h) + c.kw / h - h - ta
            up = ex > 0.0
            lo = np.where(up, mid, lo)
            hi = np.where(up, hi, mid)
        h = 10.0 ** (0.5 * (lo + hi))
        co3 = dic * c.k1 * c.k2 / (h * h + c.k1 * h + c.k1 * c.k2)
        return c.calcium * co3 / c.ksp_calcite

    z = _grid(depth_cm, first_cm, nodes) / 100.0                       # m
    d = d_sediment_m2_s(t_c, sed.porosity)
    k_area = sed.k_per_s() * sed.calcite_mol_per_kg()                  # mol per kg pore water per s at Omega 0
    r0 = sed.resp_at_top()
    r_org = r0 * np.exp(-z / (sed.resp_depth_cm / 100.0))

    # finite-volume diffusion operator on the interior and base nodes (node 0 is the bottom water)
    m = len(z)
    gap = np.diff(z)
    vol = np.empty(m)
    vol[1:-1] = 0.5 * (gap[:-1] + gap[1:])
    vol[-1] = 0.5 * gap[-1]
    up = np.zeros(m)
    dn = np.zeros(m)
    up[1:] = d / gap / vol[1:]                                          # to the node above
    dn[1:-1] = d / gap[1:] / vol[1:-1]                                  # to the node below
    lap = diags([-(up + dn)[1:], up[2:], dn[1:-1]], [0, -1, 1], shape=(m - 1, m - 1)).tocsr()

    def rate(ta, dic):
        # (1 - Omega) for Omega < 1 and 0 above, its corner rounded over 1e-7 in Omega (a softplus) so that
        # Newton's method does not cycle between nodes either side of saturation; the rate it changes is
        # below what the profile resolves
        w = omega(ta, dic)
        x = (1.0 - w) / SMOOTH_OMEGA
        under = SMOOTH_OMEGA * np.where(x > 30.0, x, np.log1p(np.exp(np.minimum(x, 30.0))))
        return k_area * under ** sed.n, w

    ta = np.full(m, ta_bw)
    dic = np.full(m, dic_bw)
    bc = np.zeros(m - 1)
    bc[0] = up[1]                                                      # the bottom water's term at node 1
    it = 0
    for it in range(1, max_iter + 1):
        r, w = rate(ta[1:], dic[1:])
        f_dic = lap @ dic[1:] + bc * dic_bw + r + r_org[1:]
        f_ta = lap @ ta[1:] + bc * ta_bw + 2.0 * r - REDFIELD_ALK_PER_C * r_org[1:]
        # the reaction's derivatives, by differences at each node
        h_ta = 1e-9 * np.maximum(ta[1:], 1e-6)
        h_dic = 1e-9 * np.maximum(dic[1:], 1e-6)
        dr_ta = (rate(ta[1:] + h_ta, dic[1:])[0] - r) / h_ta
        dr_dic = (rate(ta[1:], dic[1:] + h_dic)[0] - r) / h_dic
        jac = _block([[lap + diags(dr_dic), diags(dr_ta)], [diags(2.0 * dr_dic), lap + diags(2.0 * dr_ta)]])
        step = spsolve(jac, -np.concatenate([f_dic, f_ta]))
        s_dic, s_ta = step[:m - 1], step[m - 1:]
        # damp a step that would take a node far from its value (Omega's response is steep near saturation)
        big = max(np.max(np.abs(s_dic) / dic_bw), np.max(np.abs(s_ta) / ta_bw))
        lam = min(1.0, 0.3 / big) if big > 0 else 1.0
        dic[1:] += lam * s_dic
        ta[1:] += lam * s_ta
        if lam == 1.0 and big < tol:
            break
    else:
        raise RuntimeError(f"pore-water column did not converge in {max_iter} iterations")
    r, w = rate(ta[1:], dic[1:])
    r = np.concatenate([[rate(ta[:1], dic[:1])[0][0]], r])
    w = np.concatenate([[omega(ta[:1], dic[:1])[0]], w])
    # calcite dissolved in the column per cm2 of seafloor: the rate per kg of pore water over the pore water
    # in each cell (phi rho_w dz)
    cell = np.empty(m)
    cell[0] = 0.5 * gap[0]
    cell[1:-1] = vol[1:-1]
    cell[-1] = vol[-1]
    flux = float(np.sum(r * cell * 100.0) * sed.porosity * sed.water_density * 1e-3) * 1e6 * 365.25 * 86400.0
    return Column(z * 100.0, ta, dic, w, r, flux, it)


def _block(rows):
    from scipy.sparse import bmat
    return bmat(rows).tocsr()


def mixed_mass_left(z_cm: np.ndarray, lam_per_yr: np.ndarray, w_cm_yr: float, db_cm2_yr: float,
                    mix_depth_cm: float) -> np.ndarray:
    """Mean fraction of its mass left in a test now at each depth, for tests that rain onto the seafloor and
    are buried at w while biodiffusion D_b mixes the top ``mix_depth_cm`` (the standard biodiffusion model:
    Guinasso and Schink 1975, JGR 80, 3032; Boudreau 1997, Ch. 3; not checked against them here).  Each test
    loses mass at lam(z) per year in proportion to its mass (first order in the calcite, as R = k [C] (1 -
    Omega)^n), so the mean mass M of the tests at depth z (per test, their number being uniform) satisfies

        d/dz (w M - D_b dM/dz) = -lam M,    w M - D_b dM/dz = w at z = 0 (each test arrives whole),

    with D_b dropping to 0 below the mixed layer and M carried down by burial alone below the grid.
    Finite volumes on the nodes ``z_cm`` (cm), cell faces halfway between them."""
    z = np.asarray(z_cm, float)
    lam = np.asarray(lam_per_yr, float)
    m = len(z)
    face = 0.5 * (z[1:] + z[:-1])                           # interior faces
    edges = np.concatenate([[0.0], face, [z[-1]]])
    width = np.diff(edges)
    db = np.where(face < mix_depth_cm, db_cm2_yr, 0.0)       # D_b at the interior faces
    diag = lam * width
    lower = np.zeros(m)                                     # coefficient of M[i-1] in row i
    upper = np.zeros(m)                                     # coefficient of M[i+1] in row i
    b = np.zeros(m)
    b[0] = w_cm_yr                                          # the rain: each test arrives whole
    g = db / np.diff(z)                                     # diffusive conductance of each interior face
    # flux out through each cell's lower face minus in through its upper face = -lam M width; advection upwind
    diag[1:] += g
    lower[1:] = -(w_cm_yr + g)
    diag[:-1] += w_cm_yr + g
    upper[:-1] = -g
    diag[-1] += w_cm_yr                                     # out of the base by burial
    from scipy.linalg import solve_banded
    ab = np.vstack([np.concatenate([[0.0], upper[:-1]]), diag, np.concatenate([lower[1:], [0.0]])])
    return solve_banded((1, 1), ab, b)
