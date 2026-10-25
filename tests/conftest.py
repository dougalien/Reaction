"""Shared fixtures.

External reference material is never stored in this repository.  By default the tests look in the
app's data folder (``%LOCALAPPDATA%\\Reaction\\external``, or ``REACTION_DATA_DIR``), which holds

    EQ36_DB/data0.cmp     the EQ3/6 v8.0a database (and the other data0.* files)
    EQ36_tests/           the EQ3/6 test library (contains 3tlib_cmp/ and 6tlib_cmp/)
    OBIGT/                the CHNOSZ OBIGT csv files

Each can be overridden to point at another copy:

    REACTION_EQ36_DIR    an Eq3_6v8.0a folder (contains DB/data0.cmp)
    REACTION_EQ36_TESTS  the EQ3/6 test library
    REACTION_OBIGT_DIR   the OBIGT csv files
"""
import os
from pathlib import Path

import pytest

from reaction.app import config
from reaction.database import load_data0


def _dir(var: str) -> Path | None:
    value = os.environ.get(var)
    return Path(value) if value and Path(value).is_dir() else None


@pytest.fixture(scope="session")
def data0_cmp():
    root = _dir("REACTION_EQ36_DIR")
    path = root / "DB" / "data0.cmp" if root else config.external_dir() / "EQ36_DB" / "data0.cmp"
    if not path.exists():
        pytest.skip("set REACTION_EQ36_DIR to an EQ3/6 v8.0a folder, or fill the data folder")
    return load_data0(path)


@pytest.fixture(scope="session")
def eq36_tests() -> Path:
    root = _dir("REACTION_EQ36_TESTS") or config.external_dir() / "EQ36_tests"
    if not (root / "3tlib_cmp").is_dir():
        pytest.skip("set REACTION_EQ36_TESTS to the EQ3/6 test library folder, or fill the data folder")
    return root


@pytest.fixture(scope="session")
def obigt_dir() -> Path:
    root = _dir("REACTION_OBIGT_DIR") or config.external_dir() / "OBIGT"
    if not (root / "inorganic_aq.csv").exists():
        pytest.skip("set REACTION_OBIGT_DIR to a folder with the CHNOSZ OBIGT csv files")
    return root


@pytest.fixture(scope="session")
def eos_db_supcrt92(obigt_dir, data0_cmp):
    from reaction.thermo.eosdb import load_eos_database
    return load_eos_database(obigt_dir, "supcrt92", eq36_database=data0_cmp)
