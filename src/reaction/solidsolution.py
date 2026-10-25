"""Solid solutions: minerals whose composition varies between endmembers.

Each one is an ideal mixture on ``sites`` equivalent sites, so the activity of endmember i is x_i^sites
(EQ3/6's type 1).  At equilibrium with the fluid x_i = (Q_i/K_i)^(1/sites), and the solid solution is
stable when these add up to one or more: a phase whose amount is an unknown and whose condition is
sum x_i = 1, solved the same way as the gas phase (see ``equilibrium.py``).

Ideal mixing ignores the miscibility gaps of real carbonates and feldspars at low temperature, where it
overstates how much of a minor endmember goes in.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

DATA = Path(__file__).with_name("data") / "solid_solutions.tsv"


@dataclass(frozen=True)
class SolidSolution:
    name: str
    sites: float
    endmembers: tuple[str, ...]
    source: str


@lru_cache(maxsize=1)
def table() -> tuple[SolidSolution, ...]:
    rows = []
    for line in DATA.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("name\t"):
            continue
        name, sites, members, source = line.split("\t")
        rows.append(SolidSolution(name, float(sites), tuple(members.split(";")), source))
    return tuple(rows)


def by_name() -> dict[str, SolidSolution]:
    return {s.name: s for s in table()}
