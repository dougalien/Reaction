"""Importer test on a small synthetic file written in the data0 layout."""
import math

import pytest

from reaction.database import Kind, ThermoRangeError, load_data0

SYNTHETIC = """\
data0.syn.R0
SYNTHETIC TEST DATABASE (not real data)
+--------------------------------------------------------------------
Miscellaneous parameters
+--------------------------------------------------------------------
Temperature limits (degC)
         0.0000  300.0000
temperatures
         0.0100   25.0000   60.0000  100.0000
       150.0000  200.0000  250.0000  300.0000
pressures
         1.0132    1.0132    1.0132    1.0132
         4.7572   15.5365   39.7365   85.8378
debye huckel a (adh)
         0.4939    0.5114    0.5465    0.5995
         0.6855    0.7994    0.9593    1.2180
debye huckel b (bdh)
         0.3253    0.3288    0.3346    0.3421
         0.3525    0.3639    0.3766    0.3925
bdot
         0.0394    0.0410    0.0438    0.0460
         0.0470    0.0470    0.0340    0.0000
cco2   (coefficients for the Drummond (1981) polynomial)
        -1.0312              0.0012806
          255.9                 0.4445
      -0.001606
log k for eh reaction
       -91.0448  -83.1049  -74.0534  -65.8641
       -57.8929  -51.6848  -46.7256  -42.6828
+--------------------------------------------------------------------
bdot parameters
+--------------------------------------------------------------------
*  species name                 azer0  neutral ion type
CO2(aq)                         3.0000    -1
+--------------------------------------------------------------------
elements
+--------------------------------------------------------------------
O         15.99940
H          1.00794
C         12.01100
Ca        40.07800
+--------------------------------------------------------------------
basis species
+--------------------------------------------------------------------
H2O
     sp.type =  basis
*    EQ3/6   =  com
     revised =  08-feb-1994
*    mol.wt. =  18.015 g/mol
*    DHazero =   3.0
     charge  =   0.0
****
     2 element(s):
      2.0000 H              1.0000 O
****
*    gflag = 1 [reported delG0f used]
*    ref-state data  [source: syn/one ]
*         delG0f =    -56.688 kcal/mol     [reported]
*         delH0f =    -68.317 kcal/mol     [reported]
*         S0PrTr =     16.712 cal/(mol*K)  [reported]
+--------------------------------------------------------------------
O2(g)
     sp.type =  basis
     revised =  30-jun-1993
*    mol.wt. =  31.999 g/mol
*    DHazero =   3.0
     charge  =   0.0
****
     1 element(s):
      2.0000 O
****
+--------------------------------------------------------------------
auxiliary basis species
+--------------------------------------------------------------------
+--------------------------------------------------------------------
aqueous species
+--------------------------------------------------------------------
CO2(aq)
     sp.type =  aqueous
     revised =  01-jan-2000
*    DHazero =   3.0
     charge  =   0.0
****
     3 element(s):
      1.0000 C              2.0000 O
****
     4 species in aqueous dissociation reaction:
    -1.0000  CO2(aq)                     -1.0000  H2O
     1.0000  H+                           1.0000  HCO3-
*
**** logK grid [0-25-60-100C @1.0132bar; 150-200-250-300C @Psat-H2O]:
        -6.5804   -6.3447   -6.3000   -6.4000
        -6.5000   -6.7000   500.0000  500.0000
*
+--------------------------------------------------------------------
solids
+--------------------------------------------------------------------
Calcite                 CaCO3
     sp.type =  solid   polymorph
     revised =  02-aug-1989
*    mol.wt. =   100.087 g/mol
     V0PrTr  =    36.934 cm**3/mol [source: syn/two]
****
     3 element(s):
      1.0000 C              1.0000 Ca             3.0000 O
****
     3 species in aqueous dissociation reaction:
    -1.0000  Calcite                     -1.0000  H+
     1.0000  Ca++                         1.0000  HCO3-
*
**** logK grid [0-25-60-100C @1.0132bar; 150-200-250-300C @Psat-H2O]:
         2.2257    1.8487    1.3330    0.7743
         0.0999   -0.5838   -1.3262   -2.2154
*
*    ref-state data  [source: syn/two, syn/three ]
*         delG0f =   -269.880 kcal/mol     [reported]
*         delG0f =   -270.000 kcal/mol     [calculated]
*    Cp coefficients [source: syn/two               ]
*         T**0   =  0.24980000E+02
*         T**1   =  0.52400000E-02
*         T**-2  = -0.62000000E+06
*         Tlimit =  926.85C
+--------------------------------------------------------------------
gases
+--------------------------------------------------------------------
O2(g)                   O2
     sp.type =  gas     refstate
     revised =  18-may-1990
*    mol.wt. =    31.999 g/mol
****
     1 element(s):
      2.0000 O
****
     2 species in aqueous dissociation reaction:
    -1.0000  O2(g)                        1.0000  O2(g)
*
**** logK grid [0-25-60-100C @1.0132bar; 150-200-250-300C @Psat-H2O]:
         0.0000    0.0000    0.0000    0.0000
         0.0000    0.0000    0.0000    0.0000
*
+--------------------------------------------------------------------
solid solutions
+--------------------------------------------------------------------
+--------------------------------------------------------------------
references
+--------------------------------------------------------------------
"""


@pytest.fixture(scope="module")
def db(tmp_path_factory):
    path = tmp_path_factory.mktemp("data0") / "data0.syn"
    path.write_text(SYNTHETIC, encoding="latin-1")
    return load_data0(path)


def test_header(db):
    assert db.title == "data0.syn.R0"
    assert db.dh_params(25.0) == pytest.approx((0.5114, 0.3288, 0.0410))
    assert db.logk_eh(25.0) == pytest.approx(-83.1049)
    assert db.elements["Ca"] == pytest.approx(40.078)
    assert db.molecular_weight("H2O") == pytest.approx(18.0153, abs=1e-3)


def test_species_kinds_and_reaction(db):
    assert db["H2O"].kind is Kind.BASIS and db["H2O"].reaction == {}
    assert db["CO2(aq)"].kind is Kind.AQUEOUS
    assert db["CO2(aq)"].reaction == {"CO2(aq)": -1.0, "H2O": -1.0, "H+": 1.0, "HCO3-": 1.0}
    assert db["Calcite"].kind is Kind.MINERAL
    assert db["Calcite"].molar_volume == pytest.approx(36.934)


def test_nonpolar_neutral_type_comes_from_the_bdot_table(db):
    assert db["CO2(aq)"].activity_type == -1


def test_logk_grid_values_and_missing_data_flag(db):
    co2 = db["CO2(aq)"].logk
    assert co2(25.0) == pytest.approx(-6.3447)
    assert math.isnan(co2.values[-1])                 # 500.0 means "no data"
    assert co2(150.0) == pytest.approx(-6.5)          # exact at grid points
    # EQ3/6 fits one polynomial per part of the grid and evaluates it everywhere,
    # so values beyond the last valid point are extrapolated; Reaction flags them
    assert not co2.extrapolated(200.0) and co2.extrapolated(250.0)
    assert math.isfinite(co2(250.0))
    with pytest.raises(ThermoRangeError):
        co2(350.0)                                    # outside the file's 0-300 C limits
    # interpolation stays between neighbours in the well-sampled low-T part
    assert -6.4 < co2(80.0) < -6.3


def test_reference_state_units_are_converted_to_si(db):
    ref = db["Calcite"].ref_state
    assert ref["dG_f"] == pytest.approx(-269.880 * 4184.0)      # reported value preferred
    assert ref["cp_T0"] == pytest.approx(24.98 * 4.184)
    assert ref["cp_Tlimit_C"] == pytest.approx(926.85)


def test_provenance_is_kept(db):
    prov = db["Calcite"].provenance
    assert prov.section == "solids" and prov.revised == "02-aug-1989"
    assert prov.sources == ["syn/three", "syn/two"] or set(prov.sources) >= {"syn/two"}


def test_identity_duplicate_is_skipped_and_recorded(db):
    assert db["O2(g)"].kind is Kind.BASIS
    assert any("O2(g)" in note for note in db.skipped)
