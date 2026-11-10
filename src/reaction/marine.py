"""The carbonate system of seawater from the oceanographic constants, as CO2SYS and PyCO2SYS compute it
with their default settings, for marine waters (salinity about 20 to 43, 0 to 40 C, to about 1000 bar).

These are stoichiometric constants on the total pH scale, in mol per kg of seawater, measured in seawater
itself, so no activity coefficients enter: they are the constants that the seawater dissolution rates of
calcite and aragonite (Naviaux et al. 2019, Dong et al. 2019) were measured and fitted with.  Inside that
range they are better ground truth than an HKF database with B-dot or Pitzer activities, which puts
calcite's saturation in seawater 0.06 to 0.16 log units low and aragonite's solubility only 1.2 times
calcite's rather than 1.5 (scripts/seawater_omega_benchmark.py).

  K1, K2      Lueker, Dickson and Keeling (2000), Mar. Chem. 70, 105-119
  KB          Dickson (1990), Deep-Sea Res. 37, 755-766
  Kw          Millero (1995), Geochim. Cosmochim. Acta 59, 661-677
  Ksp         Mucci (1983), Am. J. Sci. 283, 780-799 (calcite and aragonite)
  pressure    Millero (1995), Ksp after Ingle (1975) as in CO2SYS
  boron       Uppstrom (1974), Deep-Sea Res. 21, 161-162: 0.0004157 S / 35 mol/kg
  calcium     Riley and Tongudai (1967): 0.02128 / 40.087 S / 1.80655 mol/kg

Alkalinity counts carbonate, borate, hydroxide and H+ (phosphate, silicate and the small sulfate and
fluoride terms are left out: a few umol/kg at most in deep water).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

R_BAR = 83.14462618                 # cm3 bar / (mol K)

# Millero (1995) as in CO2SYS: delta V = a0 + a1 t + a2 t^2 (cm3/mol), delta kappa = (b0 + b1 t) 1e-3
PRESSURE_TERMS = {
    "K1": (-25.50, 0.1271, 0.0, -3.08, 0.0877),
    "K2": (-15.82, -0.0219, 0.0, 1.13, -0.1475),
    "KB": (-29.48, 0.1622, -0.002608, -2.84, 0.0),
    "Kw": (-20.02, 0.1119, -0.001409, -5.13, 0.0794),
    "calcite": (-48.76, 0.5304, 0.0, -11.76, 0.3692),
    "aragonite": (-45.96, 0.5304, 0.0, -11.76, 0.3692),
}


def _pressure_factor(name: str, t_c: float, p_bar: float) -> float:
    """ln(K at P / K at 1 atm); ``p_bar`` is the absolute pressure."""
    a0, a1, a2, b0, b1 = PRESSURE_TERMS[name]
    p = max(p_bar - 1.01325, 0.0)                       # gauge pressure, bar
    dv = a0 + a1 * t_c + a2 * t_c ** 2
    dk = (b0 + b1 * t_c) * 1e-3
    rt = R_BAR * (t_c + 273.15)
    return -dv / rt * p + 0.5 * dk / rt * p ** 2


@dataclass(frozen=True)
class Constants:
    k1: float
    k2: float
    kb: float
    kw: float
    ksp_calcite: float
    ksp_aragonite: float
    total_boron: float
    calcium: float


def constants(salinity: float, t_c: float, p_bar: float = 1.01325) -> Constants:
    s, tk = salinity, t_c + 273.15
    lnt = math.log(tk)
    pk1 = 3633.86 / tk - 61.2172 + 9.67770 * lnt - 0.011555 * s + 0.0001152 * s * s
    pk2 = 471.78 / tk + 25.9290 - 3.16967 * lnt - 0.01781 * s + 0.0001122 * s * s
    ln_kb = ((-8966.90 - 2890.53 * s ** 0.5 - 77.942 * s + 1.728 * s ** 1.5 - 0.0996 * s * s) / tk
             + 148.0248 + 137.1942 * s ** 0.5 + 1.62142 * s
             - (24.4344 + 25.085 * s ** 0.5 + 0.2474 * s) * lnt + 0.053105 * s ** 0.5 * tk)
    ln_kw = (148.9802 - 13847.26 / tk - 23.6521 * lnt
             + (-5.977 + 118.67 / tk + 1.0495 * lnt) * s ** 0.5 - 0.01615 * s)
    log_t = math.log10(tk)
    log_cal = (-171.9065 - 0.077993 * tk + 2839.319 / tk + 71.595 * log_t
               + (-0.77712 + 0.0028426 * tk + 178.34 / tk) * s ** 0.5 - 0.07711 * s + 0.0041249 * s ** 1.5)
    log_arag = (-171.945 - 0.077993 * tk + 2903.293 / tk + 71.595 * log_t
                + (-0.068393 + 0.0017276 * tk + 88.135 / tk) * s ** 0.5 - 0.10018 * s + 0.0059415 * s ** 1.5)
    f = lambda name: math.exp(_pressure_factor(name, t_c, p_bar))         # noqa: E731
    return Constants(10.0 ** -pk1 * f("K1"), 10.0 ** -pk2 * f("K2"), math.exp(ln_kb) * f("KB"),
                     math.exp(ln_kw) * f("Kw"), 10.0 ** log_cal * f("calcite"),
                     10.0 ** log_arag * f("aragonite"), 0.0004157 * s / 35.0, 0.02128 / 40.087 * s / 1.80655)


@dataclass(frozen=True)
class CarbonateState:
    ph_total: float
    co2: float                      # mol/kg of seawater
    hco3: float
    co3: float
    omega_calcite: float
    omega_aragonite: float


def solve(salinity: float, t_c: float, p_bar: float, alkalinity: float, dic: float) -> CarbonateState:
    """The carbonate system from total alkalinity and dissolved inorganic carbon (mol/kg of seawater)."""
    c = constants(salinity, t_c, p_bar)

    def excess(log_h: float) -> float:
        h = 10.0 ** log_h
        denom = h * h + c.k1 * h + c.k1 * c.k2
        carb = dic * (c.k1 * h + 2.0 * c.k1 * c.k2) / denom
        return carb + c.total_boron * c.kb / (c.kb + h) + c.kw / h - h - alkalinity

    lo, hi = -11.0, -5.0                    # excess falls as H+ rises
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if excess(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    h = 10.0 ** (0.5 * (lo + hi))
    denom = h * h + c.k1 * h + c.k1 * c.k2
    co3 = dic * c.k1 * c.k2 / denom
    return CarbonateState(-math.log10(h), dic * h * h / denom, dic * c.k1 * h / denom, co3,
                          c.calcium * co3 / c.ksp_calcite, c.calcium * co3 / c.ksp_aragonite)


def in_range(salinity: float, t_c: float, p_bar: float) -> bool:
    """Where the constants are measured (Lueker et al.: 2-35 C, S 19-43; Mucci 5-40 C, S 5-44, used down to
    0 C as CO2SYS does; the pressure terms to about 1000 bar)."""
    return 19.0 <= salinity <= 43.0 and -2.0 <= t_c <= 40.0 and p_bar <= 1100.0
