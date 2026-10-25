"""The phase-rule bookkeeping of the active set, on stoichiometry alone (no database needed)."""
import numpy as np

from reaction.equilibrium import _dependence, _independent, _phase_rule_swap

# columns: H2O, H+, Ca++, Mg++, HCO3-, SO4--
NAMES = ["calcite", "dolomite", "dolomite-ord", "gypsum", "anhydrite", "magnesite"]
NU = np.array([
    [0, -1, 1, 0, 1, 0],         # calcite      CaCO3 + H+ = Ca++ + HCO3-
    [0, -2, 1, 1, 2, 0],         # dolomite
    [0, -2, 1, 1, 2, 0],         # dolomite-ord, the same formula
    [2, 0, 1, 0, 0, 1],          # gypsum       CaSO4.2H2O = Ca++ + SO4-- + 2 H2O
    [0, 0, 1, 0, 0, 1],          # anhydrite
    [0, -1, 0, 1, 1, 0],         # magnesite
], dtype=float)
UNK = [1, 2, 3, 4, 5]            # water activity is not a Newton unknown


def test_same_formula_is_dependent():
    lam = _dependence(NU, [1], 2, UNK)
    assert lam is not None and np.allclose(lam, [1.0])


def test_hydrate_and_anhydrous_phase_are_dependent_over_the_unknowns():
    assert _dependence(NU, [4], 3, UNK) is not None


def test_independent_minerals_are_not():
    assert _dependence(NU, [0], 5, UNK) is None
    assert _dependence(NU, [], 0, UNK) is None


def test_combination_of_two_is_dependent():
    lam = _dependence(NU, [0, 5], 1, UNK)            # dolomite = calcite + magnesite
    assert lam is not None and np.allclose(lam, [1.0, 1.0])


def test_swap_hands_the_amount_to_the_newcomer():
    leave, theta, lam = _phase_rule_swap(NU, [0, 1], 2, UNK, {"calcite": 0.3, "dolomite": 0.1}, NAMES)
    assert leave == 1 and abs(theta - 0.1) < 1e-12


def test_ratio_test_picks_the_mineral_that_runs_out_first():
    # dolomite joining calcite + magnesite: each moves by one mole per mole of dolomite
    leave, theta, _ = _phase_rule_swap(NU, [0, 5], 1, UNK, {"calcite": 0.3, "magnesite": 0.05}, NAMES)
    assert leave == 5 and abs(theta - 0.05) < 1e-12


def test_independent_keeps_the_first_of_each_dependent_set():
    assert _independent(NU, [2, 1, 0, 3, 4], UNK) == [2, 0, 3]


def test_a_phase_on_fixed_components_only_is_never_swapped_in():
    # its row over the unknowns is empty: nothing to pivot on, the active set's blocking handles it
    nu = np.vstack([NU, [1, 0, 0, 0, 0, 0]])
    assert _phase_rule_swap(nu, [0], 6, UNK, {"calcite": 0.3}, NAMES + ["ice"])[0] is None
