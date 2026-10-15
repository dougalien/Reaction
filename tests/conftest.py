"""Shared fixtures.

External reference material is never stored in this repository.  Point the
environment at local copies to enable the validation tests:

    REACTION_EQ36_DIR    the Eq3_6v8.0a folder (contains DB/data0.cmp)
    REACTION_EQ36_TESTS  the EQ3/6 test library (contains 3tlib_cmp/)
"""
import os
from pathlib import Path

import pytest

from reaction.database import load_data0


def _dir(var: str) -> Path | None:
    value = os.environ.get(var)
    return Path(value) if value and Path(value).is_dir() else None


@pytest.fixture(scope="session")
def data0_cmp():
    root = _dir("REACTION_EQ36_DIR")
    if root is None or not (root / "DB" / "data0.cmp").exists():
        pytest.skip("set REACTION_EQ36_DIR to an EQ3/6 v8.0a folder")
    return load_data0(root / "DB" / "data0.cmp")


@pytest.fixture(scope="session")
def eq36_tests() -> Path:
    root = _dir("REACTION_EQ36_TESTS")
    if root is None:
        pytest.skip("set REACTION_EQ36_TESTS to the EQ3/6 test library folder")
    return root
