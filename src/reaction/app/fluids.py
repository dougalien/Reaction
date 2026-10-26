"""Named starting fluids.  Each is an analysis (as a user would type it) plus a pH at a measurement
temperature, read by ``analysis.py`` exactly as a custom analysis is."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FluidPreset:
    key: str
    label: str
    description: str
    measured_t_c: float
    analysis: str
    ph: float
    source: str
    units: str = "mg/kg"
    charge_balance: str = "auto"            # "auto", "none" or the ion to adjust
    caveats: tuple[str, ...] = ()
    redox: str = "auto"                     # see REDOX_CHOICES
    redox_value: float | None = None        # log fO2, Eh (V) or pe when ``redox`` asks for one
    measured_p_bar: float | None = None     # pressure of the measurement; None = 1 bar or on the boiling curve


CUSTOM = "custom"
CUSTOM_LABEL = "Custom fluid"

# How the redox state of a fluid is set.  A species name (O2, H2, Fe+3, H2S, NO3-) uses its
# concentration in the analysis as the redox couple, with its partner in the basis.
REDOX_CHOICES = {
    "auto": "from a redox species in the analysis (dissolved O2 or H2, Fe(III), sulfide), else not tracked",
    "none": "not tracked: only the basis oxidation states (Fe(II), sulfate, ...) are present",
    "air": "in equilibrium with air (O2 fugacity 0.21 bar)",
    "log_fo2": "a given log oxygen fugacity",
    "eh": "a given Eh, volts",
    "pe": "a given pe",
}
AIR_LOG_FO2 = -0.678         # log10(0.2095)

PRESETS: dict[str, FluidPreset] = {
    "seawater": FluidPreset(
        "seawater", "Seawater",
        "Standard seawater, major ions only (no silica, iron or trace metals).",
        25.0,
        "Na+ 10768\nK+ 399.1\nCa+2 412.3\nMg+2 1291.8\nCl- 19353\nSO4-2 2712\nHCO3- 2.022 mmolal",
        8.22,
        "composition of the EQ3/6 test case swmaj (pH and alkalinity measured at 25 C); "
        "oxygen from equilibrium with air",
        charge_balance="none", redox="air",
        caveats=(
            "The magnesium hydroxide sulfate hydrate that forms when seawater is heated in "
            "experiments is not in the database, so magnesium removal and the acid it releases "
            "are underestimated here.",
            "The fluid contains no silica, iron or aluminium, so only carbonate, sulfate and "
            "hydroxide minerals can form from it alone.",
        )),
    "rainwater": FluidPreset(
        "rainwater", "Rainwater",
        "Dilute rain in equilibrium with atmospheric CO2, with a little sea salt and sulfate.",
        25.0,
        "Na+ 0.4\nK+ 0.1\nCa+2 0.46\nMg+2 0.1\nCl- 0.7\nSO4-2 1.5\nHCO3- 0.016 mmolal",
        5.6,
        "illustrative values of typical magnitude for continental rain, not a specific sample; "
        "carbonate from CO2 at 10^-3.4 bar, oxygen from equilibrium with air",
        redox="air"),
    "vent_fluid": FluidPreset(
        "vent_fluid", "Vent fluid",
        "Black-smoker end-member fluid: acid, reduced, magnesium- and sulfate-free, rich in silica, "
        "iron, manganese and hydrogen sulfide.",
        25.0,
        "Na+ 432\nK+ 23.2\nLi+ 0.89\nCa+2 15.6\nSr+2 0.08\nMn+2 0.96\nFe+2 1.66\nAl+3 0.005\n"
        "SiO2 17.3\nCl- 489\nHCO3- 5.7\nH2S 6.5\nH2 0.5",
        3.4,
        "rounded values similar to East Pacific Rise 21 N end-member fluids (Von Damm and others, "
        "1985); pH measured at 25 C after quenching. Check the paper before quoting numbers. "
        "Dissolved H2 (which sets the redox state) is a typical magnitude for such fluids, not a "
        "measured value",
        units="mmol/kg",
        caveats=("The redox state comes from an illustrative 0.5 mmol/kg of dissolved hydrogen; "
                 "change it on the card if you have a measurement.",)),
    "formation_brine": FluidPreset(
        "formation_brine", "Formation brine",
        "A sodium-calcium chloride basin brine, about 160 g/kg of dissolved solids.",
        25.0,
        "Na+ 45000\nK+ 1000\nCa+2 12000\nMg+2 1500\nSr+2 600\nCl- 98000\nBr- 600\nSO4-2 300\n"
        "HCO3- 150",
        6.0,
        "illustrative composition of a sedimentary basin brine, not a specific sample",
        caveats=("The ionic strength is about 3 molal, beyond what the B-dot activity model is "
                 "calibrated for; treat activities and saturation states as rough until Pitzer "
                 "equations are added.",
                 "The redox state is not tracked (no redox species in the analysis), so only "
                 "ferrous iron and sulfate are present.")),
    "freshwater": FluidPreset(
        "freshwater", "Fresh water",
        "World-average river water in equilibrium with the atmosphere (air oxygen, 420 ppm CO2).",
        25.0,
        "Ca+2 15\nMg+2 4.1\nNa+ 6.3\nK+ 2.3\nHCO3- 58.4\nSO4-2 11.2\nCl- 7.8\nSiO2 13.1",
        8.14,
        "mean composition of the world's rivers (Livingstone 1963, USGS Professional Paper 440-G; "
        "mg/L); nitrate and iron left out. The pH is the one that puts the fluid at a CO2 "
        "fugacity of 10^-3.38 bar (420 ppm) at 25 C with this bicarbonate; oxygen from equilibrium with air",
        units="mg/l", charge_balance="none", redox="air",
        caveats=("World-average river water is a mean of many rivers, not a sample; local waters differ "
                 "widely, so adjust the analysis to the water you have in mind.",
                 "The silica is above quartz saturation at 25 C, as in most rivers; block quartz on the "
                 "card if it should stay dissolved.")),
    "gsl_brine": FluidPreset(
        "gsl_brine", "Great Salt Lake brine",
        "A sodium chloride-sulfate lake brine near halite saturation, about 250 g/kg of dissolved "
        "solids, in equilibrium with the atmosphere.",
        25.0,
        "Na+ 85600\nK+ 4070\nMg+2 7200\nCa+2 241\nCl- 140000\nSO4-2 16400\nHCO3- 251",
        8.2,
        "Great Salt Lake brine as tabulated by Drever (1997, The Geochemistry of Natural Waters, 3rd "
        "edition; mg/kg), a brine of about 25% salinity like the lake's concentrated stages and its "
        "north arm (Jones and others 2009, Aquatic Geochemistry 15, 123-195); silica (48 mg/kg) left "
        "out. The pH, not part of that analysis, puts the brine close to the atmosphere's CO2 "
        "(10^-3.36 bar against 10^-3.38) at 25 C with Pitzer activities; oxygen from equilibrium with air",
        charge_balance="auto", redox="air",
        caveats=("At an ionic strength near 6.5 molal only the Pitzer activity model is meaningful; the "
                 "brine card uses it.",
                 "Silica is left out: at 25 C it would make magnesium silicates (talc, sepiolite) that do "
                 "not form in the lake on these time scales. Add \"SiO2 48\" to the analysis to put it "
                 "back.")),
    "hydrothermal": FluidPreset(
        "hydrothermal", "Hydrothermal fluid",
        "Seawater after the default heating run: heated from 2 to 350 C at 500 bar, with the anhydrite "
        "and brucite it formed left behind. Acid (pH 5.1 at 350 C), calcium- and sulfate-poor.",
        350.0,
        "Na+ 485.431\nK+ 10.5791\nCa+2 0.38103\nMg+2 52.3797\nCl- 565.752\nSO4-2 18.978\n"
        "HCO3- 2.02214\nO2 0.228445",
        5.0755,
        "the fluid at the last point of this app's seawater heating card (seawater preset, 2 to 350 C "
        "at 500 bar, 36 steps, static, B-dot), recomputed by the test suite; element totals per kg of "
        "water, pH and dissolved O2 at 350 C and 500 bar",
        units="mmolal", charge_balance="none", measured_p_bar=500.0,
        caveats=("This is heated seawater, not a vent fluid: it has not reacted with rock, so it carries "
                 "no silica, iron or metals (the Vent fluid preset is a natural black-smoker fluid).",
                 "Its pH and dissolved oxygen are at 350 C and 500 bar; the run must stay at a pressure "
                 "where water is liquid.")),
}


def preset(key: str) -> FluidPreset:
    try:
        return PRESETS[key]
    except KeyError:
        raise KeyError(f"unknown fluid {key!r}; known fluids: {', '.join(PRESETS)}") from None


def label(key: str) -> str:
    return PRESETS[key].label if key in PRESETS else CUSTOM_LABEL if key == CUSTOM else key
