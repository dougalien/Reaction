"""Reading water analyses, presets and custom fluids on the card (no thermodynamic data needed)."""
import pytest

from reaction.app import analysis
from reaction.app.assistant import interpret_rules
from reaction.app.fluids import PRESETS
from reaction.app.scenario import Scenario


def test_names_in_any_common_notation():
    for written in ("Ca", "Ca2+", "Ca++", "Ca+2", "calcium", "CALCIUM"):
        assert analysis.lookup(written)[0] == "Ca+2", written
    for written in ("SO4", "SO4--", "SO42-", "SO4-2", "SO4=", "sulphate", "SO₄²⁻"):
        assert analysis.lookup(written)[0] == "SO4-2", written
    assert analysis.lookup("Si") == ("SiO2", "Si")
    assert analysis.lookup("nonsense") is None


def test_units_per_line_and_heading_units():
    got = analysis.parse("Na: 10,768 mg/L\nCa2+ (mmol/kg) 10.3\nSiO2 12\nMg <0.1")
    assert [(e.species, e.value, e.unit) for e in got.entries] == [
        ("Na+", 10768.0, "mg/l"), ("Ca+2", 10.3, "mmol/kg"), ("SiO2", 12.0, "mg/kg")]
    assert any("below detection" in n for n in got.notes)


def test_the_card_reports_what_it_cannot_read():
    got = analysis.parse("Na 10\nunobtainium 3\nNa 4\nalkalinity 120")
    assert any("unobtainium" in p for p in got.problems)
    assert any("given twice" in p for p in got.problems)
    assert any("alkalinity" in p for p in got.problems)


def test_seawater_molality_uses_dissolved_solids():
    """10768 mg/kg of solution is 0.4855 mol/kg of water once the 35 g/kg of salt is taken out."""
    got = analysis.parse(PRESETS["seawater"].analysis)
    molal, tds = analysis.molalities(got.entries)
    assert tds == pytest.approx(35055, rel=1e-3)
    assert molal["Na+"] == pytest.approx(10768 / 1000 / 22.98977 / (1 - tds / 1e6), rel=1e-9)
    assert molal["HCO3-"] == pytest.approx(2.022e-3, rel=1e-6)          # molal entries stay as given


def test_mass_of_an_element_is_weighed_as_that_element():
    si, = analysis.molalities(analysis.parse("Si 28.0855 mg/kg").entries)[0].values()
    sio2, = analysis.molalities(analysis.parse("SiO2 60.0843 mg/kg").entries)[0].values()
    assert si == pytest.approx(sio2, rel=1e-4) == pytest.approx(1e-3, rel=1e-3)


def test_milliequivalents_divide_by_the_charge():
    molal, _ = analysis.molalities(analysis.parse("Ca 2 meq/l").entries)
    assert molal["Ca+2"] == pytest.approx(1e-3, rel=1e-3)


def test_balance_ion_adds_to_the_short_side():
    assert analysis.balance_ion({"Na+": 0.5, "Cl-": 0.4, "SO4-2": 0.01}) == "Cl-"
    assert analysis.balance_ion({"Na+": 0.3, "Ca+2": 0.01, "Cl-": 0.4}) == "Na+"


def test_every_preset_reads_cleanly():
    for key, p in PRESETS.items():
        card = Scenario(fluid=key)
        assert not card.problems(), (key, card.problems())
        got = analysis.parse(p.analysis, p.units)
        assert got.entries and not got.problems, key


def test_editing_a_preset_starts_a_custom_fluid_from_it():
    card, problems = Scenario().changed({"ph": 7.0})
    assert not problems and card.fluid == "custom" and card.ph == 7.0
    assert card.analysis == PRESETS["seawater"].analysis and card.charge_balance == "none"
    card, _ = Scenario().changed({"fluid": "vent_fluid", "ph": 4.0})   # preset named in the same change
    assert card.fluid == "custom" and "Fe+2" in card.analysis and card.units == "mmol/kg"


def test_custom_fluid_problems():
    card = Scenario(fluid="custom", analysis="", ph=None)
    text = " ".join(card.problems())
    assert "empty" in text and "pH" in text
    card = Scenario(fluid="custom", analysis="Na 10\nCl 15", ph=7.0, charge_balance="SO4-2")
    assert any("SO4-2" in p for p in card.problems())
    card, problems = card.changed({"charge_balance": "chloride"})
    assert not problems and card.charge_balance == "Cl-" and not card.problems()


def test_a_pasted_analysis_in_the_chat_fills_a_custom_fluid():
    got = interpret_rules("Na 120 mg/L, K 4, Ca 40, Mg 10, Cl 180, SO4 50, HCO3 120, pH 7.6", Scenario())
    ch = got.scenario_changes
    assert ch["fluid"] == "custom" and ch["ph"] == 7.6 and not got.run
    assert analysis.parse(ch["analysis"]).entries[1].unit == "mg/l"      # the one unit applies to all
    card, _ = Scenario().changed(ch)
    assert not card.problems()


def test_a_pasted_analysis_without_ph_asks_for_it():
    got = interpret_rules("Na 10 mmol/kg, Cl 12, Ca 1", Scenario())
    assert got.question and got.scenario_changes["ph"] is None
    card, _ = Scenario().changed(got.scenario_changes)
    assert any("pH" in p for p in card.problems())


def test_ordinary_requests_are_not_read_as_analyses():
    for text in ("heat seawater to 300 C at 500 bar, both", "cool the vent fluid from 350 to 25 C",
                 "plot Ca++ and Mg++ against temperature"):
        assert "analysis" not in interpret_rules(text, Scenario()).scenario_changes, text
    assert interpret_rules("use the vent fluid", Scenario()).scenario_changes["fluid"] == "vent_fluid"
    assert interpret_rules("a formation brine", Scenario()).scenario_changes["fluid"] == "formation_brine"


def test_reactants_from_text_and_from_the_chat():
    from reaction.app.scenario import parse_reactants
    assert parse_reactants("10 g albite, quartz 0.1 mol; 5 mmol HCl") == [
        {"name": "albite", "amount": 10.0, "unit": "g"}, {"name": "quartz", "amount": 0.1, "unit": "mol"},
        {"name": "HCl", "amount": 0.005, "unit": "mol"}]
    assert parse_reactants("2 g dolomite,ordered")[0]["name"] == "dolomite,ordered"
    with pytest.raises(ValueError):
        parse_reactants("lots of albite")
    got = interpret_rules("react 10 g of albite with seawater at 25 C", Scenario())
    card, problems = Scenario().changed(got.scenario_changes)
    assert not problems and card.reactants == [{"name": "albite", "amount": 10.0, "unit": "g"}]
    assert card.t_start_c == card.t_end_c == 25.0 and got.run
    card.p_start_bar = card.p_end_bar = 1.0
    assert not card.problems()                        # nothing else varies, but the reactant does
    assert "adding 10 g albite" in card.describe()


def test_mixing_on_the_card_and_in_the_chat():
    card = Scenario(fluid="vent_fluid", t_start_c=350.0, mix_with="seawater", mix_ratio=10.0, mix_t_c=2.0,
                    p_start_bar=250.0, p_end_bar=250.0, steps=11)
    assert not card.problems() and card.varies == ["temperature_c"]
    assert card.points()[-1][0] == pytest.approx(370.0 / 11)
    assert "mixed with up to 10 kg of seawater at 2 °C" in card.describe()
    assert Scenario(mix_with="custom").problems()
    got = interpret_rules("mix vent fluid at 350 C with up to 10 kg of seawater at 2 C, 250 bar", Scenario())
    assert got.scenario_changes == {"pressure_bar": 250.0, "fluid": "vent_fluid", "mix_with": "seawater",
                                    "mix_ratio": 10.0, "mix_t_c": 2.0, "t_start_c": 350.0}
    assert interpret_rules("no mixing", Scenario()).scenario_changes == {"mix_with": ""}
