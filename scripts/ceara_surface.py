"""A surface of foram mass remaining over water depth and sedimentation rate at the Ceara Rise.

Sites with in-situ O2 fluxes: every station of Hales and Emerson (1997, GCA 61, 501-514), Tables 1-2: B1 3280 m, B2 3290 m,
B3 3270 m (Station B sampled three times), A 4000 m, C 4120 m, G 4685 m, each with the mean of its own O2-flux
electrode readings (B1 10, 20; B2 15, 38; B3 13, 24; A 24; C 12; G 11, 13, 14 umol/cm2/yr) and CaCO3 of 65 wt%
(B, A; C takes A's, as they did) or 40 wt% (G).  The four further sites of Martin and Sayles (1996, GCA 60,
243-263, Table 1 and Fig. 5C), E 3772 m, F 4164 m, H 4267 m and J 4342 m, have no O2 flux; they get their organic
carbon oxidation rate (E 26, F 18.5, H 10, J 15 umol/cm2/yr, read from Fig. 5C) times 1.1, the mean ratio of Hales
and Emerson's O2 flux to Martin and Sayles's oxidation rate at the three sites both studied (B 20/22, A 24/16,
G 12.7/14), and 65 wt% CaCO3 (60-67% at all their sites but G).  Between the sites, a 100 m grid from 3300 to 4600 m with the O2 flux
and CaCO3 interpolated linearly in depth between the neighbouring stations (an assumption: there are no
measurements there).  All else as in ceara_bracket.py: k 125 %/day, n 4.5, beta 0.64, porosity 0.65,
respiration depth 0.5 cm, bottom water from GLODAP station 198, 20 cm burial, a 10 cm mixed layer with the
Db of Boudreau (1994).  Sedimentation 0.5, 1, 2, 4 and 8 cm/kyr.  Writes a CSV.

    python scripts/ceara_surface.py out.csv [--trend] [--likely] [--family] [--coretops]     (about 40 minutes, the trend grids 70 more on 4 cores; reruns skip cases already in out.csv)
"""
from __future__ import annotations

import csv
import dataclasses
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from ceara_bracket import RATES_CM_KYR, boudreau_db, kept_pct, pressure_bar  # noqa: E402

# name, depth m, CaCO3 wt%, O2 flux umol/cm2/yr (mean of the electrodes)
STATIONS = (("B3", 3270.0, 65.0, 18.5), ("B1", 3280.0, 65.0, 15.0), ("B2", 3290.0, 65.0, 26.5),
            ("E", 3772.0, 65.0, 28.6), ("A", 4000.0, 65.0, 24.0), ("C", 4120.0, 65.0, 12.0),
            ("F", 4164.0, 65.0, 20.35), ("H", 4267.0, 65.0, 11.0), ("J", 4342.0, 65.0, 16.5),
            ("G", 4685.0, 40.0, 12.67))
GRID_M = tuple(float(z) for z in range(3300, 4601, 100))


def interpolated(z: float) -> tuple[str, float, float, float]:
    """A grid site: CaCO3 and O2 flux linear in depth between the sites (Station B as the mean of B1-B3)."""
    sites = [("B", 3280.0, 65.0, 20.0)] + [s for s in STATIONS if not s[0].startswith("B")]
    zs = [s[1] for s in sites]
    return (f"grid {z:.0f} m", z, float(np.interp(z, zs, [s[2] for s in sites])),
            float(np.interp(z, zs, [s[3] for s in sites])))


# The regional trend: the O2 flux fitted linearly in depth through the eight site depths (B as the mean of B1-B3),
# O2 = 52.63 - 0.008456 z (r2 0.32; residual standard deviation 5.6 umol/cm2/yr); organic carbon oxidation falls
# with depth on the Rise as Martin and Sayles found.  TREND_SHIFT moves it by that deviation for low and high rain.
TREND = (52.6305, -0.00845614)
TREND_SD = 5.6
TREND_M = tuple(float(z) for z in range(3300, 4701, 100))


def trend_site(z: float, shift: float = 0.0) -> tuple[str, float, float, float]:
    """A site on the regional trend: O2 flux from the fit (plus ``shift`` deviations), CaCO3 as interpolated."""
    o2 = TREND[0] + TREND[1] * z + shift * TREND_SD
    caco3 = interpolated(z)[2] if z <= 4685.0 else 40.0
    return f"trend{shift:+g} {z:.0f} m", z, caco3, round(o2, 3)


def run(case: tuple) -> tuple:
    from reaction.app.engine import Engine
    from reaction.app.scenario import example_card
    (name, z, caco3, o2), w = case
    p = pressure_bar(z)
    card = dataclasses.replace(
        example_card("ceara"), sinking=False, depth_axis=False, burial=True, p_start_bar=p, p_end_bar=p, steps=21,
        sedimentation_cm_kyr=w, burial_depth_cm=20.0, sediment_column=True, porosity=0.65, caco3_wt_pct=caco3,
        sediment_k_pct_day=125.0, sediment_n=4.5, o2_flux_umol_cm2_yr=o2, co2_per_o2=0.64, respiration_depth_cm=0.5,
        mixed_layer_cm=10.0, biodiffusion_cm2_yr=boudreau_db(w))
    res = Engine().run(card)
    return name, z, caco3, o2, w, kept_pct(res, "calcite"), kept_pct(res, "aragonite")


HEADER = ["site", "depth_m", "caco3_wt_pct", "o2_flux_umol_cm2_yr", "sedimentation_cm_kyr", "foram_kept_pct",
          "pteropod_kept_pct"]


# the likely Holocene sedimentation on the Rise: 1.2-2.5 cm/kyr (Martin and Sayles 1996, Table 4; Curry and
# Lohmann 1990) and its midpoint
LIKELY_CM_KYR = (1.2, 1.85, 2.5)


# the family of curves: a fixed O2 flux at every depth (0 = no organic rain, the bottom-water-only master curve)
FAMILY_O2 = (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0)
FAMILY_M = tuple(float(z) for z in range(3300, 4701, 25))     # every 25 m, to resolve the measured profile


def family_site(z: float, o2: float) -> tuple[str, float, float, float]:
    return f"o2={o2:g} {z:.0f} m", z, (interpolated(z)[2] if z <= 4685.0 else 40.0), o2


# the 14 Ceara Rise core tops of Cullen and Curry (1997, Proc. ODP Sci. Results 154, ch. 14, Fig. 5): water depth (m)
# and % whole foraminifers, read by eye from the figure (about +-30 m, +-2 %).  Each is run on the regional O2 trend
# (extrapolated above 3270 m, the shallowest measured site) and one deviation either side.
CORETOPS = ((2807.0, 92.0), (3057.0, 95.0), (3159.0, 96.0), (3295.0, 88.0), (3364.0, 85.0), (3534.0, 92.0),
            (3636.0, 88.0), (3727.0, 85.0), (3795.0, 91.0), (4000.0, 89.0), (4284.0, 81.0), (4330.0, 88.0),
            (4557.0, 29.0), (4580.0, 20.0))


def coretop_cases() -> list[tuple]:
    """The core tops on the central trend at the likely rates, and at 1.85 cm/kyr one O2 deviation either side."""
    cases = [(trend_site(z), w) for z, _ in CORETOPS for w in LIKELY_CM_KYR]
    return cases + [(trend_site(z, k), 1.85) for z, _ in CORETOPS for k in (-1.0, 1.0)]


def main(out: str, trend: bool = False, likely: bool = False, family: bool = False, coretops: bool = False) -> int:
    """Runs the cases not already in ``out`` (same site, depth, CaCO3, O2 flux and rate) and appends them.
    ``trend``: the regional-trend grid at the fitted O2 flux and one deviation below and above it.
    ``likely``: only the sites and the central trend, at the likely sedimentation rates."""
    sites = list(STATIONS) + [interpolated(z) for z in GRID_M]
    if trend:
        sites += [trend_site(z, k) for k in (0.0, -1.0, 1.0) for z in TREND_M]
    done = set()
    try:
        with open(out, newline="") as f:
            done = {(r[0], *(round(float(x), 3) for x in r[1:5])) for r in list(csv.reader(f))[1:]}
    except FileNotFoundError:
        with open(out, "w", newline="") as f:
            csv.writer(f).writerow(HEADER)
    if likely:
        sites = list(STATIONS) + [trend_site(z) for z in TREND_M]
    if family:
        sites = [family_site(z, o2) for o2 in FAMILY_O2 for z in FAMILY_M]
    rates = (1.85,) if family else LIKELY_CM_KYR if likely else RATES_CM_KYR
    cases = coretop_cases() if coretops else [(s, w) for s in sites for w in rates]
    cases = [(s, w) for s, w in cases if (s[0], *(round(x, 3) for x in (*s[1:], w))) not in done]
    with ProcessPoolExecutor(4) as pool, open(out, "a", newline="") as f:
        wr = csv.writer(f)
        for row in pool.map(run, cases):
            wr.writerow([row[0], *(round(x, 3) for x in row[1:])])
            f.flush()
    return 0


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    sys.exit(main(args[0] if args else "ceara_surface.csv", "--trend" in sys.argv, "--likely" in sys.argv,
                  "--family" in sys.argv, "--coretops" in sys.argv))
