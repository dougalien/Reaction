"""Reader for the CSV tables EQ6 writes at every print point of a reaction path.

Rows are keyed by table letter (B1, C1, D1, E1, J, P, Q ...); each table has the
columns Xi, t(days) and then its own.  Values EQ6 could not compute are -99999.
Short rows (the product-mineral table P lists only minerals that are present) are padded
with zeros.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

MISSING = -99999.0


@dataclass
class Table:
    key: str
    title: str
    columns: list[str] = field(default_factory=list)        # after Xi, t(days)
    xi: list[float] = field(default_factory=list)
    time: list[float] = field(default_factory=list)
    rows: list[list[float | None]] = field(default_factory=list)

    def column(self, name: str) -> list[float | None]:
        k = self.columns.index(name)
        return [r[k] for r in self.rows]

    def row_at(self, i: int) -> dict[str, float | None]:
        return dict(zip(self.columns, self.rows[i]))


@dataclass
class Eq6Csv:
    tables: dict[str, Table]

    def __getitem__(self, key: str) -> Table:
        return self.tables[key]

    def __contains__(self, key: str) -> bool:
        return key in self.tables

    @property
    def xi(self) -> list[float]:
        return self.tables["B1"].xi


# tables whose blank entries mean "not reported" (phases far from saturation are omitted)
BLANK_IS_MISSING = {"Q", "K", "T"}


def _number(text: str, blank: float | None = 0.0) -> float | None:
    text = text.strip()
    if text == "":
        return blank
    if "*" in text:                                  # Fortran field overflow
        return None
    m = re.fullmatch(r"([+-]?\d+\.?\d*)([+-]\d{3})", text)   # 1.4494-128: no 'E' for 3-digit exponents
    value = float(f"{m.group(1)}E{m.group(2)}") if m else float(text)
    return None if value == MISSING else value


def read_eq6_csv(path: str | Path) -> Eq6Csv:
    tables: dict[str, Table] = {}
    current: Table | None = None
    with Path(path).open(encoding="latin-1", newline="") as fh:
        for row in csv.reader(fh):
            if not row or not row[0]:
                continue
            key = row[0]
            if key in ("A", "Z"):
                continue
            if len(row) > 2 and row[1].startswith("Table "):
                current = Table(key, row[2])
                tables[key] = current
                continue
            if current is None or current.key != key:
                continue
            if len(row) > 1 and row[1] in ("EndTable:", "Endtable:"):
                current = None
                continue
            if len(row) > 1 and row[1] == "Xi":
                current.columns = [c for c in row[3:] if c != ""]
                continue
            if len(row) < 3:
                continue
            try:
                xi, t = float(row[1]), float(row[2])
            except ValueError:
                continue
            blank = None if key in BLANK_IS_MISSING else 0.0
            values = [_number(x, blank) for x in row[3:3 + len(current.columns)]]
            values += [blank] * (len(current.columns) - len(values))
            current.xi.append(xi)
            current.time.append(t)
            current.rows.append(values)
    return Eq6Csv(tables)
