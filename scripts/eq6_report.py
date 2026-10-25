"""Run the EQ6 reaction-path test cases through Reaction and list what disagrees with EQ6's own output.

    python scripts/eq6_report.py                 # every case
    python scripts/eq6_report.py methane pyrsw   # selected cases

Reads the EQ3/6 data and test library from the app's data folder (``%LOCALAPPDATA%\\Reaction\\external``,
or ``REACTION_DATA_DIR``): ``EQ36_DB/data0.cmp`` and ``EQ36_tests/6tlib_cmp``.  For each case it prints the
number of compared values and how many fall outside the tolerances in ``validation/eq6_run.py``, then the
worst few.  EQ6 is the oracle for the numerics only; it does not shape the design (docs/PRODUCT.md).
"""
from __future__ import annotations

import glob
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from reaction.app import config                                      # noqa: E402
from reaction.database import load_data0                              # noqa: E402
from reaction.validation.eq3nr_run import Unsupported                 # noqa: E402
from reaction.validation.eq6_run import exceeds, run_case             # noqa: E402


def main(only: list[str]) -> None:
    ext = config.external_dir()
    db = load_data0(ext / "EQ36_DB" / "data0.cmp")
    root = ext / "EQ36_tests" / "6tlib_cmp"
    for f in sorted(glob.glob(str(root / "*.6i"))):
        name = Path(f).stem
        if only and name not in only:
            continue
        t0 = time.time()
        try:
            report, _ = run_case(db, f)
        except Unsupported as exc:
            print(f"{name:10} UNSUPPORTED: {exc}")
            continue
        except Exception as exc:                                     # keep going through the other cases
            tail = traceback.format_exc().strip().splitlines()
            print(f"{name:10} ERROR {type(exc).__name__}: {str(exc)[:100]}  [{tail[-3].strip()[:90]}]")
            continue
        bad = [d for d in report.deviations if exceeds(d.what, d)]
        kinds: dict[str, int] = {}
        for d in bad:
            kinds[d.what.split(" ")[0]] = kinds.get(d.what.split(" ")[0], 0) + 1
        print(f"{name:10} {report.n_points:3d} pts, {len(report.deviations):5d} numbers, {len(bad):4d} outside "
              f"tolerance  ({time.time() - t0:.1f}s)  " + "; ".join(f"{k}:{v}" for k, v in kinds.items()))
        for d in sorted(bad, key=lambda d: -d.delta / (abs(d.reference) + 1e-9))[:6]:
            print(f"      {d.what} at xi={d.xi:.4g}: ours {d.ours:.6g} ref {d.reference:.6g}")


if __name__ == "__main__":
    main(sys.argv[1:])
