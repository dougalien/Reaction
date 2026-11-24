"""Bracket the conditions a Ceara Rise foram and pteropod met: three measured water columns and three
sediment sites, each sediment site at three sedimentation rates.  Prints the mass each shell keeps.

Water column: a 20 ug foram (calcite) and a 20 ug pteropod (aragonite), 2 m2/g, sink from 1 to 4.6 km
through GLODAPv2.2023 stations 198 (1994), 484 (2015) and 39 (1972), the shallowest, middle and deepest
calcite saturation horizons among the eight deep Ceara stations (flush mode, seawater rates, the sinking speed
of Takahashi and Be 1984).

Sediment: the same shells buried 20 cm through the pore water of Hales and Emerson (1997, GCA 61, 501) at their
Stations B (3280 m), A (4000 m) and G (4685 m): porosity 0.65, k 125 %/day, n 4.5, CO2/O2 0.64, respiration
depth 0.5 cm, CaCO3 65/65/40 wt%, O2 flux 20/24/13 umol/cm2/yr (their Tables 1-2).  Bottom water from station
198 at each depth.  Bioturbation: a 10 cm mixed layer (Boudreau 1994, GCA 58, 1243, global mean 9.8 cm) with
Db = 15.7 w^0.69 cm2/yr (w in cm/yr, the same paper).  Sedimentation 0.5, 1, 2, 4 and 8 cm/kyr, the rates of
the earlier runs.

    python scripts/ceara_bracket.py            (about seven minutes)
"""
from __future__ import annotations

import dataclasses
import sys

from reaction.app.engine import Engine
from reaction.app.scenario import (WATER_PROFILE_CEARA, WATER_PROFILE_CEARA_1972, WATER_PROFILE_CEARA_2015,
                                   example_card)

COLUMNS = (("Station 198, 1994 (horizon about 4275 m)", WATER_PROFILE_CEARA),
           ("Station 484, 2015 (about 4380 m)", WATER_PROFILE_CEARA_2015),
           ("Station 39, 1972 (about 4470 m)", WATER_PROFILE_CEARA_1972))
# Hales and Emerson (1997) sites: depth m, CaCO3 wt%, O2 flux umol/cm2/yr
SITES = (("Station B, 3280 m", 3280.0, 65.0, 20.0), ("Station A, 4000 m", 4000.0, 65.0, 24.0),
         ("Station G, 4685 m", 4685.0, 40.0, 13.0))
RATES_CM_KYR = (0.5, 1.0, 2.0, 4.0, 8.0)
MIXED_LAYER_CM, BURIAL_CM = 10.0, 20.0
UG = {"calcite": 100.087e6, "aragonite": 100.087e6}


def pressure_bar(z_m: float) -> float:
    return 1.01325 + z_m * 1033.0 * 9.81 / 1e5


def boudreau_db(w_cm_kyr: float) -> float:
    """Biodiffusion coefficient from the sedimentation rate, Boudreau (1994): Db = 15.7 w^0.69, cm2/yr."""
    return 15.7 * (w_cm_kyr / 1000.0) ** 0.69


def kept_pct(res: dict, mineral: str) -> float:
    run = res["modes"]["flush"]
    mixed = run.get("kinetic", {}).get(f"mixed:{mineral}")
    if mixed:
        return float(mixed[-1])
    return 100.0 * run["minerals"][mineral][-1] * UG[mineral] / 20.0


def water_cases():
    for name, prof in COLUMNS:
        yield name, dataclasses.replace(example_card("ceara"), water_profile=prof)


def sediment_cases():
    for name, z, caco3, o2 in SITES:
        for w in RATES_CM_KYR:
            p = pressure_bar(z)
            yield f"{name}, {w:g} cm/kyr", dataclasses.replace(
                example_card("ceara"), sinking=False, depth_axis=False, burial=True, p_start_bar=p, p_end_bar=p,
                steps=21, sedimentation_cm_kyr=w, burial_depth_cm=BURIAL_CM, sediment_column=True, porosity=0.65,
                caco3_wt_pct=caco3, sediment_k_pct_day=125.0, sediment_n=4.5, o2_flux_umol_cm2_yr=o2,
                co2_per_o2=0.64, respiration_depth_cm=0.5, mixed_layer_cm=MIXED_LAYER_CM,
                biodiffusion_cm2_yr=boudreau_db(w))


def main() -> int:
    engine = Engine()
    print("| Case | Foram kept % | Pteropod kept % |\n|---|---|---|")
    for name, card in [*water_cases(), *sediment_cases()]:
        assert not card.problems(), (name, card.problems())
        res = engine.run(card)
        print(f"| {name} | {kept_pct(res, 'calcite'):.2f} | {kept_pct(res, 'aragonite'):.2f} |", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
