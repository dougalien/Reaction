"""Phase relations of NaCl-H2O from the correlation formulae of Driesner and Heinrich (2007, GCA 71, 4880-4901).

Temperatures in C, pressures in bar, compositions as NaCl mole fraction X (``mass_fraction`` and ``mole_fraction``
convert).  Equation and table numbers are the paper's.  Pure water's vapour pressure is IAPWS-IF97's (the paper
used IAPS-84; the difference does not matter here).

The liquid + vapour field is bounded below by vapour + liquid + halite (``p_vlh``) and above by the critical
curve (``p_crit``) or, below water's critical temperature, by water's boiling curve.  Inside it, ``x_liquid`` and
``x_vapour`` give the coexisting compositions.  Above the field the fluid is one phase; below ``p_vlh`` vapour
coexists with halite (``x_vapour_halite``), and above it a liquid may be saturated with halite
(``x_halite_liquid``, the halite liquidus).
"""
from __future__ import annotations

import math

T_CRIT_H2O = 373.976          # C (IAPS-84, as the paper)
P_CRIT_H2O = 220.54915        # bar
T_TRIPLE_NACL = 800.7         # C
P_TRIPLE_NACL = 5.0e-4        # bar
M_NACL, M_H2O = 58.443, 18.015

# Table 3
A_MELT, B_SUBL, B_BOIL = 2.47260e-2, 1.18061e4, 0.941812e4
# Table 4: critical pressure (c, exponents cA) and critical composition (d)
C = {1: -2.36, 2: 1.28534e-1, 3: -2.3707e-2, 4: 3.20089e-3, 5: -1.38917e-4, 6: 1.02789e-7, 7: -4.8376e-11,
     8: 2.36, 9: -1.31417e-2, 10: 2.98491e-3, 11: -1.30114e-4, 14: -4.88336e-4}
CA = {1: 1.0, 2: 1.5, 3: 2.0, 4: 2.5, 5: 3.0, 6: 4.0, 7: 5.0, 8: 1.0, 9: 2.0, 10: 2.5, 11: 3.0}
D = {1: 8.00000e-05, 2: 1.00000e-05, 3: -1.37125e-07, 4: 9.46822e-10, 5: -3.50549e-12, 6: 6.57369e-15,
     7: -4.89423e-18, 8: 7.77761e-2, 9: 2.7042e-4, 10: -4.244821e-7, 11: 2.580872e-10}
# Table 5: halite liquidus, e_i(P)
E = [(0.0989944, 3.30796e-6, -4.71759e-10), (0.00947257, -8.66460e-6, 1.69417e-9),
     (0.610863, -1.51716e-5, 1.19290e-8), (-1.64994, 2.03441e-4, -6.46015e-8),
     (3.36474, -1.54023e-4, 8.17048e-8)]
# Table 6: pressure of vapour + liquid + halite
F = [4.64e-3, 5.0e-7, 1.69078e1, -2.69148e2, 7.63204e3, -4.95636e4, 2.33119e5, -5.13556e5, 5.49708e5, -2.84628e5]
# Table 7: liquid branch of V + L
H = {1: 1.68486e-3, 2: 2.19379e-4, 3: 4.3858e2, 4: 1.84508e1, 5: -5.6765e-10, 6: 6.73704e-6, 7: 1.44951e-7,
     8: 3.84904e2, 9: 7.07477e0, 10: 6.06896e-5, 11: 7.62859e-3}
# Table 8: vapour branch of V + L (and halite-saturated vapour)
K = {0: -0.235694, 1: -0.188838, 2: 0.004, 3: 0.0552466, 4: 0.66918, 5: 396.848, 6: 45.0, 7: -3.2719e-7,
     8: 141.699, 9: -0.292631, 10: -0.00139991, 11: 1.95965e-6, 12: -7.3653e-10, 13: 0.904411, 14: 0.000769766,
     15: -1.18658e-6}


def mole_fraction(w: float) -> float:
    """NaCl mole fraction from mass fraction."""
    if w <= 0:
        return 0.0
    if w >= 1:
        return 1.0
    return (w / M_NACL) / (w / M_NACL + (1.0 - w) / M_H2O)


def mass_fraction(x: float) -> float:
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    return x * M_NACL / (x * M_NACL + (1.0 - x) * M_H2O)


def p_water_sat(t_c: float) -> float:
    """Vapour pressure of pure water (IAPWS-IF97 region 4), bar; valid to its critical point."""
    t = t_c + 273.15
    n = (0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5, -0.32325550322333e7,
         0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849, 0.65017534844798e3)
    th = t + n[8] / (t - n[9])
    a = th * th + n[0] * th + n[1]
    b = n[2] * th * th + n[3] * th + n[4]
    c = n[5] * th * th + n[6] * th + n[7]
    return (2.0 * c / (-b + math.sqrt(b * b - 4.0 * a * c))) ** 4 * 10.0     # MPa -> bar


def t_halite_melting(p_bar: float) -> float:                                  # eq. 1
    return T_TRIPLE_NACL + A_MELT * (p_bar - P_TRIPLE_NACL)


def p_nacl(t_c: float) -> float:
    """Vapour pressure of NaCl: sublimation below its triple point, boiling above (eqs. 2, 3)."""
    b = B_SUBL if t_c < T_TRIPLE_NACL else B_BOIL
    return 10.0 ** (math.log10(P_TRIPLE_NACL) + b * (1.0 / (T_TRIPLE_NACL + 273.15) - 1.0 / (t_c + 273.15)))


def _p_crit_5b(t_c: float) -> float:
    dt = t_c - T_CRIT_H2O
    return P_CRIT_H2O + sum(C[n] * dt ** CA[n] for n in range(8, 12))


def p_crit(t_c: float) -> float:
    """Critical pressure (eqs. 5a-c); below water's critical temperature the hypothetical extension 5a."""
    if t_c < T_CRIT_H2O:
        dt = T_CRIT_H2O - t_c
        return P_CRIT_H2O + sum(C[n] * dt ** CA[n] for n in range(1, 8))
    if t_c <= 500.0:
        return _p_crit_5b(t_c)
    dt = 500.0 - T_CRIT_H2O
    c12 = _p_crit_5b(500.0)
    c13 = sum(C[n] * CA[n] * dt ** (CA[n] - 1.0) for n in range(8, 12))
    return c12 + c13 * (t_c - 500.0) + C[14] * (t_c - 500.0) ** 2


def x_crit(t_c: float) -> float:
    """Critical composition (eqs. 7a, 7b), for T at or above water's critical temperature."""
    if t_c <= 600.0:
        dt = t_c - T_CRIT_H2O
        return sum(D[i] * dt ** i for i in range(1, 8))
    dt = t_c - 600.0
    return sum(D[i] * dt ** (i - 8) for i in range(8, 12))


def x_halite_liquid(t_c: float, p_bar: float) -> float:
    """The halite liquidus: liquid saturated with halite (eq. 8 with table 5)."""
    e = [a + b * p_bar + c * p_bar ** 2 for a, b, c in E]
    e.append(1.0 - sum(e))
    r = t_c / t_halite_melting(p_bar)
    return sum(ei * r ** i for i, ei in enumerate(e))


def p_vlh(t_c: float) -> float:
    """Pressure of vapour + liquid + halite (eq. 10), up to the triple point of NaCl."""
    f = F + [P_TRIPLE_NACL - sum(F)]
    r = t_c / T_TRIPLE_NACL
    return sum(fi * r ** i for i, fi in enumerate(f))


def _g12(t_c: float) -> tuple[float, float]:
    g1 = H[2] + (H[1] - H[2]) / (1.0 + math.exp((t_c - H[3]) / H[4])) + H[5] * t_c ** 2
    g2 = H[7] + (H[6] - H[7]) / (1.0 + math.exp((t_c - H[8]) / H[9])) + H[10] * math.exp(-H[11] * t_c)
    return g1, g2


def _liquid_branch(t_c: float) -> tuple[float, float, float, float, float]:
    """(Pcrit, Xcrit, g0, g1, g2) of eq. 11 at t_c, with g0 (and below water's critical point Xcrit) set by the
    branch ending on the halite liquidus at vapour + liquid + halite (or at NaCl's boiling curve, X = 1, above
    its triple point) and, below water's critical point, at X = 0 on water's boiling curve."""
    pc = p_crit(t_c)
    g1, g2 = _g12(t_c)
    if t_c < T_TRIPLE_NACL:
        p_end = p_vlh(t_c)
        x_end = x_halite_liquid(t_c, p_end)
    else:
        p_end, x_end = p_nacl(t_c), 1.0
    q_end = pc - p_end
    if t_c >= T_CRIT_H2O:
        xc = x_crit(t_c)
        g0 = (x_end - xc - g1 * q_end - g2 * q_end ** 2) / math.sqrt(q_end)
        return pc, xc, g0, g1, g2
    q_w = pc - p_water_sat(t_c)
    # X(q) = xc + g0 sqrt(q) + g1 q + g2 q^2 with X(q_w) = 0 and X(q_end) = x_end: linear in xc and g0
    r_w = -(g1 * q_w + g2 * q_w ** 2)
    r_e = x_end - (g1 * q_end + g2 * q_end ** 2)
    g0 = (r_e - r_w) / (math.sqrt(q_end) - math.sqrt(q_w))
    xc = r_w - g0 * math.sqrt(q_w)
    return pc, xc, g0, g1, g2


def p_top(t_c: float) -> float:
    """Upper pressure of the liquid + vapour field: the critical curve, or water's boiling curve below its
    critical temperature."""
    return p_crit(t_c) if t_c >= T_CRIT_H2O else p_water_sat(t_c)


def p_bottom(t_c: float) -> float:
    """Lower pressure of the liquid + vapour field: vapour + liquid + halite, or NaCl's boiling curve above its
    triple point."""
    return p_vlh(t_c) if t_c < T_TRIPLE_NACL else p_nacl(t_c)


def x_liquid(t_c: float, p_bar: float) -> float:
    """NaCl mole fraction of the liquid on the V + L surface (eq. 11)."""
    pc, xc, g0, g1, g2 = _liquid_branch(t_c)
    q = max(pc - p_bar, 0.0)
    return xc + g0 * math.sqrt(q) + g1 * q + g2 * q * q


def _log_k_bar(t_c: float, p_bar: float) -> float:                            # eqs. 16, 17 with table 8
    j0 = K[0] + K[1] * math.exp(-K[2] * t_c)
    j1 = K[4] + (K[3] - K[4]) / (1.0 + math.exp((t_c - K[5]) / K[6])) + K[7] * (t_c + K[8]) ** 2
    j2 = K[9] + K[10] * t_c + K[11] * t_c ** 2 + K[12] * t_c ** 3
    j3 = K[13] + K[14] * t_c + K[15] * t_c ** 2
    pn = p_nacl(t_c)
    pb = (p_bar - pn) / (p_crit(t_c) - pn)
    u = max(1.0 - pb, 0.0)
    return 1.0 + j0 * u ** j1 + j2 * u + j3 * u ** 2 - (1.0 + j0 + j2 + j3) * u ** 3


def _x_vapour_from(x_l: float, t_c: float, p_bar: float) -> float:
    """The vapour paired with liquid x_l through the distribution coefficient (eqs. 14, 15)."""
    pn = p_nacl(t_c)
    log_xls0 = math.log10(x_halite_liquid(t_c, pn)) if t_c < T_TRIPLE_NACL else 0.0
    log_k = _log_k_bar(t_c, p_bar) * (math.log10(pn / p_crit(t_c)) - log_xls0) + log_xls0
    return x_l / 10.0 ** log_k * pn / p_bar


def x_vapour(t_c: float, p_bar: float) -> float:
    """NaCl mole fraction of the vapour on the V + L surface (eqs. 13-17)."""
    return _x_vapour_from(x_liquid(t_c, p_bar), t_c, p_bar)


def x_vapour_halite(t_c: float, p_bar: float) -> float:
    """Vapour saturated with halite (eq. 9), below the V + L + H pressure."""
    return _x_vapour_from(x_halite_liquid(t_c, p_bar), t_c, p_bar)


def in_two_phase_field(t_c: float, p_bar: float) -> bool:
    """Whether liquid and vapour can coexist at T and P (for some composition)."""
    if t_c >= T_TRIPLE_NACL + 300.0:
        return False
    return p_bottom(t_c) < p_bar < p_top(t_c)
