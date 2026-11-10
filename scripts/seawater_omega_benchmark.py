"""Calcite and aragonite saturation of S = 35 seawater (TA 2300, DIC 2000 umol/kg) in Reaction (B-dot and
Pitzer) against PyCO2SYS (pip install PyCO2SYS; its defaults: Lueker et al. 2000 K1/K2, Mucci 1983 Ksp,
Millero 1995 pressure terms).  Reaction's pH is found by matching the total alkalinity of its speciation.
Major ions: the reference composition of Millero et al. (2008), mol per kg of seawater.
Run: python scripts/seawater_omega_benchmark.py"""
import PyCO2SYS as pyco2
from reaction.app.engine import Engine
from reaction.speciation import speciate
from reaction.network import PROTON
from reaction.speciation import Molal, PH
db = Engine().db
W = 0.9648                                    # kg water per kg seawater, S = 35
sw = {"Na+": 0.4690, "Mg+2": 0.05282, "Ca+2": 0.01028, "K+": 0.01021, "Sr+2": 9.07e-5, "Cl-": 0.5457,
      "SO4-2": 0.02824, "Br-": 8.42e-4, "F-": 6.83e-5, "B(OH)3": 4.16e-4}
TA, DIC = 2300e-6, 2000e-6
def alk(sp):
    a = 0.0
    for n, st in sp.species.items():
        s = db.species[n]
        if n == "HCO3-": a += st.molality
        elif n == PROTON: a -= st.molality
        elif s.reaction:
            a += st.molality * (s.reaction.get("HCO3-", 0.0) - s.reaction.get(PROTON, 0.0))
    return a
def reaction_omega(t, p, model):
    cons = {k: Molal(v / W) for k, v in sw.items()}
    cons["HCO3-"] = Molal(DIC / W)
    target = TA / W
    lo, hi = 7.0, 9.0
    for _ in range(40):
        ph = 0.5 * (lo + hi); cons[PROTON] = PH(ph)
        sp = speciate(db, t, cons, pressure_bar=p, activity_model=model)
        if alk(sp) > target: hi = ph
        else: lo = ph
    co3 = sp.species["CO3-2"].molality
    return ph, 10 ** sp.saturation["calcite"], 10 ** sp.saturation["aragonite"], co3 * W * 1e6, sp.species["CO3-2"].log_gamma
for t, p in ((25, 1.013), (2, 1.013), (2, 400.0)):
    r = pyco2.sys(par1=TA * 1e6, par2=DIC * 1e6, par1_type=1, par2_type=2, salinity=35, temperature=t,
                  pressure=(p - 1.013) * 10)
    print(f"T {t} P {p:g} bar | PyCO2SYS: pH_T {r['pH_total']:.3f} pH_free {r['pH_free']:.3f} CO3 {r['carbonate']:.1f} umol/kg "
          f"Om_cal {r['saturation_calcite']:.3f} Om_arag {r['saturation_aragonite']:.3f}")
    for model in ("bdot", "pitzer"):
        try:
            ph, oc, oa, co3, lg = reaction_omega(t, p, model)
            print(f"   Reaction {model:6s}: pH {ph:.3f} free CO3 {co3:.1f} umol/kg (log gamma {lg:.2f}) Om_cal {oc:.3f} Om_arag {oa:.3f}")
        except Exception as e:
            print("   Reaction", model, "failed:", str(e)[:200])
