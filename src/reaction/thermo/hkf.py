"""Standard-state Gibbs energy of aqueous species: the revised Helgeson-Kirkham-Flowers
(HKF) equations of state (Helgeson et al., 1981; Tanger and Helgeson, 1988; Shock et
al., 1992 for the g function that makes the Born coefficient depend on T and P).

    G(T,P) = dGf - S (T - Tr)                                  (origin)
           + non-solvation part (c1, c2, a1..a4)               (Cp and V functions)
           - w(T,P)(Z + 1) + w (Zr + 1) + w Yr (T - Tr)       (solvation, Z = -1/eps)

with w(T,P) from the effective electrostatic radius r_e = r_ref + |z| g(T,P).

Parameters are kept in J and bar.  Water properties come from IAPWS-95 and the
Johnson-Norton dielectric constant (`water.py`, `electrostatics.py`).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

from . import electrostatics, water

CAL = 4.184
T_REF, P_REF = 298.15, 1.0
THETA = 228.0            # K
PSI = 2600.0             # bar
ETA = 1.66027e5 * CAL    # J angstrom / mol  (1.66027e5 cal A/mol)
CM3BAR_TO_J = 0.1


class WaterStateError(ValueError):
    """The requested T, P is not in the liquid/dense-fluid field of water."""


@dataclass(frozen=True)
class WaterState:
    t_k: float
    p_bar: float
    rho: float        # g/cm3
    eps: float


@lru_cache(maxsize=8192)
def water_state(t_k: float, p_bar: float) -> WaterState:
    """Liquid or supercritical water at (T, P); steam is rejected.

    A pressure up to 1% below the saturation pressure is accepted as metastable
    liquid, which the EQ3/6 grid needs at 100 C (1.0132 bar against 1.0142 bar).
    """
    if t_k < water.T_CRIT:
        psat = water.saturation_pressure_mpa(t_k) * water.BAR_PER_MPA
        if p_bar < 0.99 * psat:
            raise WaterStateError(f"P={p_bar} bar is below Psat={psat:.3f} bar at {t_k} K: "
                                  f"water is steam and aqueous species are undefined")
    rho = water.density(t_k, p_bar, liquid=True) / 1000.0
    return WaterState(t_k, p_bar, rho, electrostatics.epsilon_from_density(rho, t_k))


@lru_cache(maxsize=1)
def _reference_born() -> tuple[float, float]:
    """(Zr + 1 style constants) -> (-1/eps_r, Y_r) at 25 C, 1 bar."""
    ws = water_state(T_REF, P_REF)
    y, _, _ = electrostatics.born_functions(T_REF, P_REF)
    return -1.0 / ws.eps, y


# --- g function (Shock et al., 1992) -------------------------------------------------------
_AG = (-2.037662, 5.747000e-3, -6.557892e-6)
_BG = (6.107361, -1.074377e-2, 1.268348e-5)
_AF = (0.3666666e2, -0.1504956e-9, 0.5017997e-13)


def g_function(t_k: float, p_bar: float, rho_g_cm3: float) -> float:
    """Solvent function g (angstrom); zero at and above 1 g/cm3 and above 6 kbar."""
    if rho_g_cm3 >= 1.0 or p_bar > 6000.0:
        return 0.0
    tc = t_k - 273.15
    a_g = _AG[0] + _AG[1] * tc + _AG[2] * tc ** 2
    b_g = _BG[0] + _BG[1] * tc + _BG[2] * tc ** 2
    g = a_g * (1.0 - rho_g_cm3) ** b_g
    if tc > 155.0 and p_bar < 1000.0 and tc < 355.0:
        x = (tc - 155.0) / 300.0
        f = (x ** 4.8 + _AF[0] * x ** 16) * (_AF[1] * (1000.0 - p_bar) ** 3
                                             + _AF[2] * (1000.0 - p_bar) ** 4)
        g -= f
    return g


# --- species --------------------------------------------------------------------------------
@dataclass(frozen=True)
class HKFSpecies:
    """HKF parameters in J, bar and K (already unscaled)."""
    name: str
    charge: float
    g0: float          # dGf at Tr, Pr (J/mol)
    s0: float          # third-law entropy (J/mol/K)
    a1: float          # J/mol/bar
    a2: float          # J/mol
    a3: float          # J K/mol/bar
    a4: float          # J K/mol
    c1: float          # J/mol/K
    c2: float          # J K/mol
    omega: float       # J/mol
    cp0: float | None = None
    v0: float | None = None      # cm3/mol

    def omega_pt(self, t_k: float, p_bar: float, rho: float) -> float:
        """Born coefficient at (T, P) from the g function; H+ and neutral species keep omega."""
        z = self.charge
        if z == 0.0 or self.name == "H+":
            return self.omega
        g = g_function(t_k, p_bar, rho)
        re_ref = z * z / (self.omega / ETA + z / 3.082)
        re = re_ref + abs(z) * g
        return ETA * (z * z / re - z / (3.082 + g))

    def gibbs(self, t_k: float, p_bar: float) -> float:
        ws = water_state(t_k, p_bar)
        z_r, y_r = _reference_born()
        dt, dp = t_k - T_REF, p_bar - P_REF
        ln_p = math.log((PSI + p_bar) / (PSI + P_REF))
        non_solv = (-self.c1 * (t_k * math.log(t_k / T_REF) - t_k + T_REF)
                    - self.c2 * ((1.0 / (t_k - THETA) - 1.0 / (T_REF - THETA)) * ((THETA - t_k) / THETA)
                                 - (t_k / THETA ** 2) * math.log((T_REF * (t_k - THETA))
                                                                 / (t_k * (T_REF - THETA))))
                    + self.a1 * dp + self.a2 * ln_p
                    + (self.a3 * dp + self.a4 * ln_p) / (t_k - THETA))
        w_pt = self.omega_pt(t_k, p_bar, ws.rho)
        solv = -w_pt * (-1.0 / ws.eps + 1.0) + self.omega * (z_r + 1.0) + self.omega * y_r * dt
        return self.g0 - self.s0 * dt + non_solv + solv

    # --- diagnostics used to validate the parameters against tabulated V and Cp -----------
    def volume_298(self) -> float:
        """V(Tr, Pr) in cm3/mol from the HKF volume function."""
        _, q, _ = electrostatics.born_functions(T_REF, P_REF)
        v_j = (self.a1 + self.a2 / (PSI + P_REF)
               + (self.a3 + self.a4 / (PSI + P_REF)) / (T_REF - THETA) - self.omega * q)
        return v_j / CM3BAR_TO_J

    def heat_capacity_298(self) -> float:
        """Cp(Tr, Pr) in J/mol/K from the HKF heat-capacity function."""
        _, _, x = electrostatics.born_functions(T_REF, P_REF)
        return self.c1 + self.c2 / (T_REF - THETA) ** 2 + self.omega * T_REF * x
