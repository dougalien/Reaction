"""What would move the model onto the Cullen and Curry (1997) Ceara core tops?  Each core top on the regional O2
trend at 1.85 cm/kyr (ceara_surface.py), with one change at a time:
  ssa      the foram's surface area 2 -> 8 m2/g (in the sediment column the test dissolves at the sediment's rate
           constant per mass, so its own surface area should cancel; this checks that)
  factor   the buried test's rate constant per mass as a fraction of the sediment's calcite: 0.5, 0.25, 0.1
           (Hales and Emerson find the dissolving calcite mostly in the 8-64 um fraction, finer than a test)
  k        the sediment's rate constant at the ends of Hales and Emerson's range, 25 and 300 %/day (best 125)
    python scripts/ceara_sensitivity.py out.csv      (reruns skip cases already in out.csv)
"""
from __future__ import annotations

import csv
import dataclasses
import sys
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from ceara_bracket import boudreau_db, kept_pct, pressure_bar  # noqa: E402
from ceara_surface import CORETOPS, trend_site  # noqa: E402

W = 1.85
CASES = ([("ssa", 8.0)] + [("factor", f) for f in (0.5, 0.25, 0.1)] + [("k", k) for k in (25.0, 300.0)])


def run(case: tuple) -> tuple:
    from reaction.app.engine import Engine
    from reaction.app.scenario import example_card
    (kind, value), z = case
    _, _, caco3, o2 = trend_site(z)
    p = pressure_bar(z)
    card = example_card("ceara")
    extra = {}
    if kind == "ssa":
        extra["reactants"] = [dict(r, ssa_m2_g=value) if r["name"] == "forams" else r for r in card.reactants]
    if kind == "factor":
        extra["shell_k_factor"] = value
    card = dataclasses.replace(
        card, sinking=False, depth_axis=False, burial=True, p_start_bar=p, p_end_bar=p, steps=21,
        sedimentation_cm_kyr=W, burial_depth_cm=20.0, sediment_column=True, porosity=0.65, caco3_wt_pct=caco3,
        sediment_k_pct_day=value if kind == "k" else 125.0, sediment_n=4.5, o2_flux_umol_cm2_yr=o2, co2_per_o2=0.64,
        respiration_depth_cm=0.5, mixed_layer_cm=10.0, biodiffusion_cm2_yr=boudreau_db(W), **extra)
    return kind, value, z, o2, kept_pct(Engine().run(card), "calcite")


def main(out: str) -> int:
    done = set()
    try:
        with open(out, newline="") as f:
            done = {(r[0], float(r[1]), float(r[2])) for r in list(csv.reader(f))[1:]}
    except FileNotFoundError:
        with open(out, "w", newline="") as f:
            csv.writer(f).writerow(["change", "value", "depth_m", "o2_flux_umol_cm2_yr", "foram_kept_pct"])
    cases = [(c, z) for c in CASES for z, _ in CORETOPS if (c[0], c[1], z) not in done]
    with ProcessPoolExecutor(4) as pool, open(out, "a", newline="") as f:
        wr = csv.writer(f)
        for row in pool.map(run, cases):
            wr.writerow([row[0], row[1], row[2], round(row[3], 3), round(row[4], 3)])
            f.flush()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "ceara_sensitivity.csv"))
