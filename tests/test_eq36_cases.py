"""Regression: every EQ3NR speciation case in the EQ3/6 test library (3tlib_cmp).

Each case is run through Reaction from its ``.3i`` input and compared with the
stored ``.3o`` output: ionic strength, water activity, osmotic coefficient, pH,
log molality and log gamma of every printed species, log Q/K of every solid and
log fugacity of every gas.  EQ3NR prints five significant figures, so
agreement within 1e-4 means agreement to the precision of the reference.

The cases exercise: charge balancing (on H+, Na+, Cl-), mg/L and mg/kg.sol
inputs, equilibrium with minerals and gases, a special basis switch
(NH3(aq) -> NO3-), redox from Eh, log fO2 and the NO3-/NH3 couple, decoupled
components (O2(aq), Fe+++, NO2-), suppressed species and phases, the Davies
equation, temperatures of 9.5, 10, 25.5, 70 and 250 C, and multi-problem files.
"""
import os
from pathlib import Path

import pytest

from reaction.app import config
from reaction.validation.eq3nr_run import run_case

# the test library: REACTION_EQ36_TESTS, else the app's data folder (see conftest.py)
_ROOT = os.environ.get("REACTION_EQ36_TESTS") or str(config.external_dir() / "EQ36_tests")
_LIB = Path(_ROOT) / "3tlib_cmp"
CASES = sorted(p.stem for p in _LIB.glob("*.3i")) if _LIB.is_dir() else []

TOLERANCE = 1e-4


@pytest.mark.parametrize("case", CASES or [pytest.param(None, marks=pytest.mark.skip(
    reason="set REACTION_EQ36_TESTS to the EQ3/6 test library folder"))])
def test_case_matches_eq3nr(data0_cmp, case):
    for rep in run_case(data0_cmp, _LIB / f"{case}.3i"):
        assert not rep.missing_species, f"{rep.name}: species absent from Reaction"
        assert not rep.missing_solids, f"{rep.name}: solids absent from Reaction"
        worst = rep.worst(3)
        assert rep.max_delta < TOLERANCE, (
            f"{rep.name}: " + "; ".join(f"{d.what}: {d.ours:.6g} vs {d.reference:.6g}" for d in worst))


def test_whole_library_is_covered(eq36_tests):
    assert len(CASES) == 41
