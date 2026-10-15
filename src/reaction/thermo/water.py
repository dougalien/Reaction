"""Properties of ordinary water: the IAPWS-95 formulation (Wagner and Pruss, 2002).

The dimensionless Helmholtz energy phi(delta, tau) = phi0 + phir, with
delta = rho / rho_c and tau = T_c / T, is the single source of the equation of
state; pressure, density, Gibbs energy and their derivatives follow from it.
The 56 residual coefficients and the ideal-gas constants are transcribed from
Tables 6.1 and 6.2 of the IAPWS-95 release; the derivative formulas are derived
here and checked against finite differences and published verification values
(tests/test_water.py).

Units: T in K, P in bar (MPa where noted), rho in kg/m3, G/H/S per mole of H2O.
"""
from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq

T_CRIT = 647.096              # K
P_CRIT = 22.064               # MPa
RHO_CRIT = 322.0              # kg/m3
R_SPECIFIC = 0.46151805       # kJ/(kg K)
MOLAR_MASS = 18.015268        # g/mol
BAR_PER_MPA = 10.0

# --- ideal-gas part (Table 6.1) ----------------------------------------------
_N0 = (-8.32044648201, 6.6832105268, 3.00632)
_N0_EXP = (0.012436, 0.97315, 1.27950, 0.96956, 0.24873)
_G0_EXP = (1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105)

# --- residual part (Table 6.2) -------------------------------------------------
_N = np.array([
    0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1,
    0.31802509345418, -0.26145533859358, -0.78199751687981e-2,
    0.88089493102134e-2, -0.66856572307965, 0.20433810950965,
    -0.66212605039687e-4, -0.19232721156002, -0.25709043003438,
    0.16074868486251, -0.40092828925807e-1, 0.39343422603254e-6,
    -0.75941377088144e-5, 0.56250979351888e-3, -0.15608652257135e-4,
    0.11537996422951e-8, 0.36582165144204e-6, -0.13251180074668e-11,
    -0.62639586912454e-9, -0.10793600908932, 0.17611491008752e-1,
    0.22132295167546, -0.40247669763528, 0.58083399985759,
    0.49969146990806e-2, -0.31358700712549e-1, -0.74315929710341,
    0.47807329915480, 0.20527940895948e-1, -0.13636435110343,
    0.14180634400617e-1, 0.83326504880713e-2, -0.29052336009585e-1,
    0.38615085574206e-1, -0.20393486513704e-1, -0.16554050063734e-2,
    0.19955571979541e-2, 0.15870308324157e-3, -0.16388568342530e-4,
    0.43613615723811e-1, 0.34994005463765e-1, -0.76788197844621e-1,
    0.22446277332006e-1, -0.62689710414685e-4, -0.55711118565645e-9,
    -0.19905718354408, 0.31777497330738, -0.11841182425981,
    -0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4,
    -0.14874640856724, 0.31806110878444])
_D = np.array([1, 1, 1, 2, 2, 3, 4, 1, 1, 1, 2, 2, 3, 4,
               4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3, 4,
               4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12,
               3, 4, 4, 5, 14, 3, 6, 6, 6, 3, 3, 3], dtype=float)  # terms 1-54
_T = np.array([-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1, 4, 6, 12, 1, 5, 4, 2,
               13, 9, 3, 4, 11, 4, 13, 1, 7, 1, 9, 10, 10, 3,
               7, 10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8,
               16, 22, 23, 23, 10, 50, 44, 46, 50, 0, 1, 4], dtype=float)
_C = np.array([1] * 15 + [2] * 20 + [3] * 4 + [4] + [6] * 4, dtype=float)   # terms 8-51
_ALPHA, _BETA3, _GAMMA, _EPS = 20.0, np.array([150.0, 150.0, 250.0]), \
    np.array([1.21, 1.21, 1.25]), 1.0                                     # terms 52-54
# terms 55-56 (non-analytical, critical region)
_N4 = _N[54:56]
_A4 = np.array([3.5, 3.5])
_B4 = np.array([0.85, 0.95])
_BB4 = np.array([0.2, 0.2])
_C4 = np.array([28.0, 32.0])
_D4 = np.array([700.0, 800.0])
_AA4 = np.array([0.32, 0.32])
_BETA4 = np.array([0.3, 0.3])

_I1 = slice(0, 7)
_I2 = slice(7, 51)
_I3 = slice(51, 54)


def phi_ideal(delta: float, tau: float):
    """Ideal-gas part and its derivatives: (phi, phi_d, phi_dd, phi_t, phi_tt, phi_dt)."""
    n1, n2, n3 = _N0
    e = np.array(_N0_EXP)
    g = np.array(_G0_EXP)
    ex = np.exp(-g * tau)
    phi = math.log(delta) + n1 + n2 * tau + n3 * math.log(tau) + float(np.sum(e * np.log1p(-ex)))
    phi_t = n2 + n3 / tau + float(np.sum(e * g * (1.0 / (1.0 - ex) - 1.0)))
    phi_tt = -n3 / tau ** 2 - float(np.sum(e * g ** 2 * ex / (1.0 - ex) ** 2))
    return phi, 1.0 / delta, -1.0 / delta ** 2, phi_t, phi_tt, 0.0


def phi_residual(delta: float, tau: float):
    """Residual part and its derivatives: (phi, phi_d, phi_dd, phi_t, phi_tt, phi_dt)."""
    # --- terms 1-7: polynomial -------------------------------------------------
    d, t, n = _D[_I1], _T[_I1], _N[_I1]
    base = n * delta ** d * tau ** t
    phi = base.sum()
    phi_d = (base * d / delta).sum()
    phi_dd = (base * d * (d - 1) / delta ** 2).sum()
    phi_t = (base * t / tau).sum()
    phi_tt = (base * t * (t - 1) / tau ** 2).sum()
    phi_dt = (base * d * t / (delta * tau)).sum()

    # --- terms 8-51: exponential in delta^c ----------------------------------
    d, t, c, n = _D[_I2], _T[_I2], _C, _N[_I2]       # _C holds terms 8-51 only
    dc = delta ** c
    ex = np.exp(-dc)
    base = n * delta ** d * tau ** t * ex
    phi += base.sum()
    phi_d += (base * (d - c * dc) / delta).sum()
    phi_dd += (base * ((d - c * dc) * (d - 1 - c * dc) - c ** 2 * dc) / delta ** 2).sum()
    phi_t += (base * t / tau).sum()
    phi_tt += (base * t * (t - 1) / tau ** 2).sum()
    phi_dt += (base * t * (d - c * dc) / (delta * tau)).sum()

    # --- terms 52-54: Gaussian bell-shaped --------------------------------------
    d, t, n = _D[_I3], _T[_I3], _N[_I3]
    dd, dt = delta - _EPS, tau - _GAMMA
    base = n * delta ** d * tau ** t * np.exp(-_ALPHA * dd ** 2 - _BETA3 * dt ** 2)
    a_d = d / delta - 2.0 * _ALPHA * dd                     # d ln(term) / d delta
    a_t = t / tau - 2.0 * _BETA3 * dt                        # d ln(term) / d tau
    phi += base.sum()
    phi_d += (base * a_d).sum()
    phi_dd += (base * (a_d ** 2 - d / delta ** 2 - 2.0 * _ALPHA)).sum()
    phi_t += (base * a_t).sum()
    phi_tt += (base * (a_t ** 2 - t / tau ** 2 - 2.0 * _BETA3)).sum()
    phi_dt += (base * a_d * a_t).sum()

    # --- terms 55-56: critical-region non-analytical terms ------------------------
    a, b, B, C, D, A, beta, n = _A4, _B4, _BB4, _C4, _D4, _AA4, _BETA4, _N4
    dm1 = delta - 1.0
    u = dm1 ** 2
    u_p = u ** (1.0 / (2.0 * beta) - 1.0) if u > 0.0 else np.zeros(2)     # u^(1/(2 beta) - 1)
    theta = (1.0 - tau) + A * u ** (1.0 / (2.0 * beta))
    Delta = theta ** 2 + B * u ** a
    # derivatives of Delta
    Dd = dm1 * (A * theta * (2.0 / beta) * u_p + 2.0 * B * a * u ** (a - 1.0))
    if u > 0.0:
        Ddd = Dd / dm1 + u * (4.0 * B * a * (a - 1.0) * u ** (a - 2.0)
                              + 2.0 * A ** 2 * (1.0 / beta) ** 2 * u_p ** 2
                              + A * theta * (4.0 / beta) * (1.0 / (2.0 * beta) - 1.0)
                              * u ** (1.0 / (2.0 * beta) - 2.0))
    else:
        Ddd = np.zeros(2)
    # Delta^b and its derivatives
    Db = Delta ** b
    Db_d = b * Delta ** (b - 1.0) * Dd
    Db_dd = b * (Delta ** (b - 1.0) * Ddd + (b - 1.0) * Delta ** (b - 2.0) * Dd ** 2)
    Db_t = -2.0 * theta * b * Delta ** (b - 1.0)
    Db_tt = 2.0 * b * Delta ** (b - 1.0) + 4.0 * theta ** 2 * b * (b - 1.0) * Delta ** (b - 2.0)
    Db_dt = (-A * b * (2.0 / beta) * Delta ** (b - 1.0) * dm1 * u_p
             - 2.0 * theta * b * (b - 1.0) * Delta ** (b - 2.0) * Dd)
    # psi and its derivatives
    psi = np.exp(-C * u - D * (tau - 1.0) ** 2)
    ps_d = -2.0 * C * dm1 * psi
    ps_dd = (2.0 * C * u - 1.0) * 2.0 * C * psi
    ps_t = -2.0 * D * (tau - 1.0) * psi
    ps_tt = (2.0 * D * (tau - 1.0) ** 2 - 1.0) * 2.0 * D * psi
    ps_dt = 4.0 * C * D * dm1 * (tau - 1.0) * psi
    phi += (n * Db * delta * psi).sum()
    phi_d += (n * (Db * (psi + delta * ps_d) + Db_d * delta * psi)).sum()
    phi_dd += (n * (Db * (2.0 * ps_d + delta * ps_dd) + 2.0 * Db_d * (psi + delta * ps_d)
                    + Db_dd * delta * psi)).sum()
    phi_t += (n * delta * (Db_t * psi + Db * ps_t)).sum()
    phi_tt += (n * delta * (Db_tt * psi + 2.0 * Db_t * ps_t + Db * ps_tt)).sum()
    phi_dt += (n * (Db * (ps_t + delta * ps_dt) + delta * Db_d * ps_t
                    + Db_t * (psi + delta * ps_d) + Db_dt * delta * psi)).sum()
    return float(phi), float(phi_d), float(phi_dd), float(phi_t), float(phi_tt), float(phi_dt)


# --- thermodynamic functions of (rho, T) -----------------------------------------
def pressure_mpa(rho: float, t_k: float) -> float:
    delta, tau = rho / RHO_CRIT, T_CRIT / t_k
    return rho * R_SPECIFIC * t_k * (1.0 + delta * phi_residual(delta, tau)[1]) / 1000.0


def _dp_drho(rho: float, t_k: float) -> float:
    """(dP/drho)_T in MPa per kg/m3."""
    delta, tau = rho / RHO_CRIT, T_CRIT / t_k
    res = phi_residual(delta, tau)
    return R_SPECIFIC * t_k * (1.0 + 2.0 * delta * res[1] + delta ** 2 * res[2]) / 1000.0


# --- saturation (auxiliary equations of Wagner and Pruss, 2002) ---------------------
def saturation_pressure_mpa(t_k: float) -> float:
    v = 1.0 - t_k / T_CRIT
    a = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)
    p = (1.0, 1.5, 3.0, 3.5, 4.0, 7.5)
    s = sum(ai * v ** pi for ai, pi in zip(a, p))
    return P_CRIT * math.exp(T_CRIT / t_k * s)


def saturated_liquid_density(t_k: float) -> float:
    v = 1.0 - t_k / T_CRIT
    b = (1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -6.74694450e5)
    e = (1 / 3, 2 / 3, 5 / 3, 16 / 3, 43 / 3, 110 / 3)
    return RHO_CRIT * (1.0 + sum(bi * v ** ei for bi, ei in zip(b, e)))


def saturated_vapor_density(t_k: float) -> float:
    v = 1.0 - t_k / T_CRIT
    c = (-2.03150240, -2.68302940, -5.38626492, -17.2991605, -44.7586581, -63.9201063)
    e = (2 / 6, 4 / 6, 8 / 6, 18 / 6, 37 / 6, 71 / 6)
    return RHO_CRIT * math.exp(sum(ci * v ** ei for ci, ei in zip(c, e)))


class WaterError(ValueError):
    pass


def density(t_k: float, p_bar: float, liquid: bool = False) -> float:
    """Density (kg/m3) of the stable fluid at T (K) and P (bar).

    On the saturation curve the liquid is returned.  With ``liquid=True`` the
    liquid branch is used below the critical temperature even when P is a little
    under the saturation pressure (metastable liquid; the EQ3/6 databases use
    1.0132 bar at 100 C, where Psat is 1.0142 bar).
    """
    if not (273.15 - 1e-9 <= t_k <= 1273.15):
        raise WaterError(f"T={t_k} K outside 273.15-1273.15 K")
    if p_bar <= 0:
        raise WaterError("P must be positive")
    p_mpa = p_bar / BAR_PER_MPA
    if p_mpa > 1000.0:
        raise WaterError("P above 10 kbar is outside the IAPWS-95 range")

    def f(rho):
        return pressure_mpa(rho, t_k) - p_mpa

    if t_k < T_CRIT:
        psat = saturation_pressure_mpa(t_k)
        rho_l = saturated_liquid_density(t_k)
        rho_v = saturated_vapor_density(t_k)
        if liquid or p_mpa >= psat * (1.0 - 1e-9):
            lo, hi = 0.93 * rho_l, 1500.0              # liquid branch (monotonic above the spinodal)
            while f(lo) > 0.0 and lo > 0.5 * rho_l:
                lo *= 0.97
        else:
            lo, hi = 1e-9, min(rho_v * 1.05, RHO_CRIT)  # vapour branch
    else:
        lo, hi = 1e-10, 1500.0
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0.0:
        raise WaterError(f"no density bracket at T={t_k} K, P={p_bar} bar")
    rho = brentq(f, lo, hi, xtol=1e-12, rtol=1e-14, maxiter=200)
    for _ in range(3):                                  # Newton polish to machine precision
        rho -= f(rho) / _dp_drho(rho, t_k)
    return rho


# --- Gibbs energy and derivatives ---------------------------------------------------------
def _molar(x_kj_per_kg: float) -> float:
    """kJ/kg -> J/mol."""
    return x_kj_per_kg * MOLAR_MASS


def _state(t_k: float, p_bar: float, liquid: bool):
    rho = density(t_k, p_bar, liquid)
    delta, tau = rho / RHO_CRIT, T_CRIT / t_k
    return rho, delta, tau, phi_ideal(delta, tau), phi_residual(delta, tau)


def gibbs_iapws(t_k: float, p_bar: float, liquid: bool = False) -> float:
    """Molar Gibbs energy (J/mol) in the IAPWS convention (u = s = 0 for the
    saturated liquid at the triple point)."""
    _, delta, tau, i, r = _state(t_k, p_bar, liquid)
    return _molar((1.0 + i[0] + r[0] + delta * r[1]) * R_SPECIFIC * t_k)


def enthalpy_iapws(t_k: float, p_bar: float, liquid: bool = False) -> float:
    _, delta, tau, i, r = _state(t_k, p_bar, liquid)
    return _molar((1.0 + tau * (i[3] + r[3]) + delta * r[1]) * R_SPECIFIC * t_k)


def entropy_iapws(t_k: float, p_bar: float, liquid: bool = False) -> float:
    _, delta, tau, i, r = _state(t_k, p_bar, liquid)
    return _molar((tau * (i[3] + r[3]) - i[0] - r[0]) * R_SPECIFIC)


def heat_capacity_p(t_k: float, p_bar: float, liquid: bool = False) -> float:
    """Isobaric heat capacity, J/(mol K)."""
    _, delta, tau, i, r = _state(t_k, p_bar, liquid)
    x = -tau ** 2 * (i[4] + r[4]) + (1.0 + delta * r[1] - delta * tau * r[5]) ** 2 / (
        1.0 + 2.0 * delta * r[1] + delta ** 2 * r[2])
    return _molar(x * R_SPECIFIC)


def expansivity_and_compressibility(t_k: float, p_bar: float, liquid: bool = False):
    """(alpha, beta): isobaric expansivity (1/K) and isothermal compressibility (1/bar)."""
    rho, delta, tau, _, r = _state(t_k, p_bar, liquid)
    dp_drho = _dp_drho(rho, t_k) * BAR_PER_MPA                 # bar per kg/m3
    dp_dt = rho * R_SPECIFIC * (1.0 + delta * r[1] - delta * tau * r[5]) / 100.0   # bar/K
    beta = 1.0 / (rho * dp_drho)
    alpha = dp_dt / (rho * dp_drho)
    return alpha, beta


# --- apparent Gibbs energy: the SUPCRT92 / HKF convention --------------------------------------
# G(T,P) = dGf(Tr,Pr) + [G(T,P) - G(Tr,Pr)] on an absolute-entropy scale, so that
# dG/dT = -S(Tr,Pr) at the reference state.  IAPWS-95 puts s = 0 at the triple point,
# so its entropy is shifted by a constant to the third-law value used by SUPCRT.
_CAL = 4.184
_T_REF, _P_REF = 298.15, 1.0
_GF_REF = -56688.0 * _CAL          # J/mol, formation from the elements
_S_REF = 16.712 * _CAL             # J/(mol K), third-law entropy


@lru_cache(maxsize=1)
def _entropy_shift() -> float:
    """Third-law entropy minus IAPWS-95 entropy at 25 C, 1 bar (J/(mol K))."""
    return _S_REF - entropy_iapws(_T_REF, _P_REF)


def gibbs(t_k: float, p_bar: float, liquid: bool = False) -> float:
    """Apparent molar Gibbs energy of formation of liquid/fluid water, J/mol."""
    return (_GF_REF + gibbs_iapws(t_k, p_bar, liquid) - gibbs_iapws(_T_REF, _P_REF)
            - _entropy_shift() * (t_k - _T_REF))
