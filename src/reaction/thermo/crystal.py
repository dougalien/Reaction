"""Standard-state Gibbs energy of crystalline phases, gases and non-aqueous liquids.

Maier-Kelley style heat capacity Cp = a + bT + c/T^2 + d/sqrt(T) + eT^2 (+ f T^lambda,
not supported), constant molar volume, as in SUPCRT92 (Helgeson et al., 1978):

    G(T,P) = dGf - S (T - Tr) + int Cp dT - T int Cp dlnT + V (P - Pr)

Gases carry V = 0, so their G is the 1-bar ideal-gas standard state; fugacity
coefficients belong to the activity model.  Quartz and coesite use the SUPCRT92
volume functions that include the alpha-beta transition.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

CAL = 4.184
T_REF, P_REF = 298.15, 1.0
CM3BAR_TO_J_SUPCRT = 0.023901488 * CAL      # SUPCRT92's cm3 bar -> J (REAC92D.F)


@dataclass(frozen=True)
class CrystalSpecies:
    """Parameters in J, bar and K."""
    name: str
    state: str                 # cr, cr2, ..., gas, liq
    g0: float                  # dGf at Tr, Pr (J/mol)
    s0: float                  # J/mol/K
    v0: float                  # cm3/mol
    a: float
    b: float
    c: float
    d: float
    e: float
    t_max: float               # upper limit of validity of the Cp equation / transition (K)
    cp0: float | None = None   # constant Cp used when no coefficients are given

    def _has_cp_terms(self) -> bool:
        return any((self.a, self.b, self.c, self.d, self.e))

    def integrals(self, t_k: float) -> tuple[float, float]:
        """(int Cp dT, int Cp dlnT) from Tr to T."""
        if self._has_cp_terms():
            int_cp_dt = (self.a * (t_k - T_REF) + self.b * (t_k ** 2 - T_REF ** 2) / 2.0
                         - self.c * (1.0 / t_k - 1.0 / T_REF)
                         + 2.0 * self.d * (math.sqrt(t_k) - math.sqrt(T_REF))
                         + self.e * (t_k ** 3 - T_REF ** 3) / 3.0)
            int_cp_dlnt = (self.a * math.log(t_k / T_REF) + self.b * (t_k - T_REF)
                           - self.c / 2.0 * (t_k ** -2 - T_REF ** -2)
                           - 2.0 * self.d * (t_k ** -0.5 - T_REF ** -0.5)
                           + self.e * (t_k ** 2 - T_REF ** 2) / 2.0)
            return int_cp_dt, int_cp_dlnt
        cp = self.cp0 or 0.0
        return cp * (t_k - T_REF), cp * math.log(t_k / T_REF)

    def heat_capacity(self, t_k: float) -> float:
        if self._has_cp_terms():
            return (self.a + self.b * t_k + self.c * t_k ** -2 + self.d * t_k ** -0.5
                    + self.e * t_k ** 2)
        return self.cp0 or 0.0

    def gibbs(self, t_k: float, p_bar: float) -> float:
        int_cp_dt, int_cp_dlnt = self.integrals(t_k)
        if self.name in ("quartz", "coesite"):
            int_vdp = _quartz_coesite_int_vdp(self, t_k, p_bar)
        else:
            int_vdp = self.v0 * (p_bar - P_REF) * CM3BAR_TO_J_SUPCRT
        return self.g0 - self.s0 * (t_k - T_REF) + int_cp_dt - t_k * int_cp_dlnt + int_vdp

    def extrapolated(self, t_k: float) -> bool:
        return self.t_max > 0.0 and t_k > self.t_max


def _quartz_coesite_int_vdp(sp: CrystalSpecies, t_k: float, p_bar: float) -> float:
    """SUPCRT92 volume integral for quartz (cr = alpha, cr2 = beta) and coesite."""
    pr, tr = P_REF, T_REF
    t_tr_pr = 848.0
    aa, ba, ca = 549.824, 0.65995, -0.4973e-4
    v_pt_tta, v_pr_ttb = 23.348, 23.72
    s_tran = 0.342
    v_pr_tra = 22.688
    v_diff = 2.047
    k = 38.5
    if sp.state == "cr2":
        p_star = p_bar
        s_star = 0.0
    else:
        p_star = pr + k * (t_k - t_tr_pr)
        s_star = s_tran
    if t_k < t_tr_pr:
        p_star, s_star = pr, 0.0
    if sp.name == "coesite":
        v_pr_tra -= v_diff
        v_pr_ttb -= v_diff
    gv_term = 0.023901488 * (
        v_pr_tra * (p_bar - p_star) + v_pr_ttb * (p_star - pr)
        - 0.5 * ca * (2.0 * pr * (p_bar - p_star) - (p_bar ** 2 - p_star ** 2))
        - ca * k * (t_k - tr) * (p_bar - p_star)
        + k * (ba + aa * ca * k) * (t_k - tr) * math.log((aa + p_bar / k) / (aa + p_star / k)))
    return gv_term * CAL
