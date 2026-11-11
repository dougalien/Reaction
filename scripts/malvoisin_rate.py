"""Derive the olivine serpentinization rate law of kinetics.py (MALVOISIN_LOG_K25, MALVOISIN_E_KJ).

Malvoisin et al. (2012), J. Geophys. Res. 117, B04102, section 3.4, para. 19: stage I rates of 38-50 um San
Carlos olivine (Fo91) at 500 bar, 4.91, 5.96, 6.86 and 0.15 % per 1000 h at 250, 270, 300 and 350 C.
Each is put on the surface of smooth 44 um spheres of 3.3 g/cm3 (assumptions, stated in kinetics.py),
divided by 1 - Omega of forsterite in a fluid buffered by chrysotile + brucite (this database, a_H2O = 1),
and 250-300 C are fitted as ln k = a - E/RT.  350 C is left out and predicted.

    python scripts/malvoisin_rate.py
"""
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from reaction.thermo.eosdb import load_eos_database          # noqa: E402

R = 8.314462618
P_BAR = 500.0
MEASURED = {250.0: 4.91, 270.0: 5.96, 300.0: 6.86, 350.0: 0.15}     # % per 1000 h
MOLAR_MASS = 0.91 * 140.69 + 0.09 * 203.77                          # g/mol, Fo91 from Mg2SiO4 and Fe2SiO4
RHO, DIAMETER = 3.3e6, 44e-6                                        # g/m3, m


def main():
    db = load_eos_database(ROOT / "src/reaction/data/reference/OBIGT", "current", eq36_database="table")

    def log_k_r(t_c):        # 2 forsterite + 3 H2O = chrysotile + brucite
        return 2 * db["forsterite"].logk(t_c, P_BAR) - db["chrysotile"].logk(t_c, P_BAR) - db["brucite"].logk(t_c, P_BAR)

    ssa = 6.0 / (RHO * DIAMETER)
    rows = []
    for t_c, pct in MEASURED.items():
        rate = pct / 100.0 / (1000.0 * 3600.0) / MOLAR_MASS / ssa
        drive = 1.0 - 10.0 ** (-log_k_r(t_c) / 2.0)
        rows.append((t_c, rate, drive))
        print(f"{t_c:5.0f} C  rate {rate:.3e} mol/m2/s  log K_r {log_k_r(t_c):6.3f}  1 - Omega {drive:.4f}")
    fit = [(1.0 / (t + 273.15), math.log(r / d)) for t, r, d in rows if t <= 300.0]
    n = len(fit)
    sx, sy = sum(x for x, _ in fit), sum(y for _, y in fit)
    sxx, sxy = sum(x * x for x, _ in fit), sum(x * y for x, y in fit)
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    a = (sy - slope * sx) / n
    print(f"ln k = {a:.3f} {slope:+.1f}/T: E = {-slope * R / 1000:.2f} kJ/mol, "
          f"log k(25 C) = {(a + slope / 298.15) / math.log(10):.3f}")
    t, r, d = rows[-1]
    pred = math.exp(a + slope / (t + 273.15)) * d
    print(f"350 C predicted {pred:.3e}, measured {r:.3e}, ratio {pred / r:.2f}")


if __name__ == "__main__":
    main()
