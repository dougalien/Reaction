"""Diagrams for the results window: what the chat asks for, and the drawing data for the plot.

A diagram request is a small dict kept by the session (``spec``).  ``compute`` turns it into the
fields, boundaries' uncertainty, water and saturation lines and, when a run is open, that run's path
drawn on top.  Unless the request says otherwise, the diagram is drawn at the run's last temperature and
pressure, and anything held fixed (pH, log fO2, the activities of other species) comes from the run's
fluid there, so the fields are those the path actually crosses.
"""
from __future__ import annotations

import math
import re
from typing import Any

from .. import diagram as dg
from ..network import OXYGEN, PROTON
from ..thermo.eosdb import DEFAULT_MASTERS
from .names import pretty_mineral, pretty_species
from .scenario import ALWAYS_EXCLUDED

KINDS = {"eh_ph": "Eh-pH", "fo2_ph": "log fO2-pH", "activity": "activity"}
CATIONS = {"na": "Na+", "sodium": "Na+", "k": "K+", "potassium": "K+", "ca": "Ca+2", "calcium": "Ca+2",
           "mg": "Mg+2", "magnesium": "Mg+2"}
_ELEMENT_WORDS = {"iron": "Fe", "manganese": "Mn", "copper": "Cu", "sulfur": "S", "sulphur": "S",
                  "carbon": "C", "aluminum": "Al", "aluminium": "Al", "zinc": "Zn", "lead": "Pb",
                  "nitrogen": "N", "uranium": "U", "silver": "Ag", "nickel": "Ni", "cobalt": "Co",
                  "mercury": "Hg", "arsenic": "As", "silica": "Si", "silicon": "Si"}
_GRID = (160, 140)


def default_spec(kind: str = "eh_ph", element: str = "Fe") -> dict:
    if kind == "activity":
        return {"kind": "activity", "element": "Al", "cation": "K+", "x": [-6.0, -1.0], "y": [0.0, 12.0],
                "t_c": None, "p_bar": None, "log_a": dg.DEFAULT_LOG_A, "with": [], "include": "minerals"}
    y = [-1.0, 1.2] if kind == "eh_ph" else [-80.0, 0.0]
    return {"kind": kind, "element": element, "x": [0.0, 14.0], "y": y, "t_c": None, "p_bar": None,
            "log_a": dg.DEFAULT_LOG_A, "with": [], "include": "all"}


def _element(word: str) -> str | None:
    w = word.strip()
    if w.lower() in _ELEMENT_WORDS:
        return _ELEMENT_WORDS[w.lower()]
    cap = w[:1].upper() + w[1:].lower()
    return cap if cap in DEFAULT_MASTERS else None


def parse_request(text: str) -> dict | None:
    """A diagram request in plain words, or None if the text does not ask for one.

    "Eh-pH diagram for Fe at 300 C", "fO2-pH diagram of copper with sulfur", "activity diagram for K at
    200 C", "close the diagram" (returns {"close": True})."""
    low = text.lower()
    if re.search(r"\b(close|hide|remove)\b.*\bdiagram\b", low):
        return {"close": True}
    if not re.search(r"\bdiagrams?\b|\bpourbaix\b", low):
        return None
    if re.search(r"activity|stability field|\bsilicate|feldspar|\bclay", low):
        kind = "activity"
    elif re.search(r"f\s*o2|fugacity|oxygen", low):
        kind = "fo2_ph"
    else:
        kind = "eh_ph"
    spec: dict[str, Any] = {"kind": kind}
    words = re.findall(r"[A-Za-z]+", text)
    m = re.search(r"\b(?:for|of)\s+([A-Za-z]+)", text)
    if kind == "activity":
        cat = None
        for w in words:
            if w.lower() in CATIONS and (len(w) > 2 or w[0].isupper()):
                cat = CATIONS[w.lower()]
                break
        if cat:
            spec["cation"] = cat
    elif m and _element(m.group(1)):
        spec["element"] = _element(m.group(1))
    else:
        for k, w in enumerate(words):
            if w.lower() in _ELEMENT_WORDS or (w[:1].isupper() and len(w) <= 2 and w != "I"):
                if w == "C" and re.search(r"\d\s*C\b", text):
                    continue                          # degrees Celsius
                e = _element(w)
                if e:
                    spec["element"] = e
                    break
    extra = []
    for mm in re.finditer(r"\b(?:with|and)\s+([A-Za-z]+)", text):
        e = _element(mm.group(1))
        if e and e != spec.get("element") and e not in extra:
            extra.append(e)
    if extra:
        spec["with"] = extra
    t = re.search(r"(-?\d+(?:\.\d+)?)\s*(?:°\s*)?(?:c\b|deg|celsius)", low)
    if t:
        spec["t_c"] = float(t.group(1))
    p = re.search(r"(\d+(?:\.\d+)?)\s*(kbar|bar|mpa)\b", low)
    if p:
        v = float(p.group(1))
        spec["p_bar"] = v * 1000 if p.group(2) == "kbar" else v * 10 if p.group(2) == "mpa" else v
    a = re.search(r"(?:activity|activities)\s+(?:of\s+)?(?:10\^?|1e)(-\d+)", low)
    if a:
        spec["log_a"] = float(a.group(1))
    return spec


def merge(old: dict | None, request: dict) -> dict:
    """A request applied over the diagram already open: a new kind starts from that kind's defaults,
    keeping the temperature and pressure unless the request names new ones."""
    kind = request.get("kind") if request.get("kind") in KINDS else (old or {}).get("kind", "eh_ph")
    if old is None or old.get("kind") != kind:
        base = default_spec(kind, request.get("element") or (old or {}).get("element", "Fe")
                            if kind != "activity" else "Al")
        if old is not None:
            base["t_c"], base["p_bar"] = old.get("t_c"), old.get("p_bar")
    else:
        base = dict(old)
    for k, v in request.items():
        if k in ("element", "cation", "t_c", "p_bar", "log_a", "with", "include", "x", "y"):
            base[k] = _clean(k, v, base.get(k))
    return base


def _clean(key: str, value, old):
    """A request value of the right type, or the old value when it is not one."""
    try:
        if key in ("t_c", "p_bar"):
            return None if value in (None, "") else float(value)
        if key == "log_a":
            return float(value)
        if key in ("x", "y"):
            lo, hi = (float(v) for v in value)
            return [lo, hi] if hi > lo else old
        if key == "element":
            return _element(str(value)) or old
        if key == "cation":
            return value if value in CATIONS.values() else old
        if key == "include":
            return value if value in ("all", "minerals", "aqueous") else old
        if key == "with":
            items = value.replace(",", " ").split() if isinstance(value, str) else list(value)
            return [e for e in (_element(str(v)) for v in items) if e]
    except (TypeError, ValueError):
        pass
    return old


def _axes(spec: dict) -> tuple[dg.Axis, dg.Axis]:
    x0, x1 = spec["x"]
    y0, y1 = spec["y"]
    if spec["kind"] == "activity":
        return dg.Axis("log_a", x0, x1, "SiO2"), dg.Axis("log_ratio", y0, y1, spec.get("cation", "K+"))
    return dg.Axis("pH", x0, x1), dg.Axis("Eh" if spec["kind"] == "eh_ph" else "log_fO2", y0, y1)


def _last(values):
    for v in reversed(values or []):
        if v is not None:
            return v
    return None


def _point_values(result: dict | None, t_c: float) -> tuple[dict, int | None, str | None]:
    """The run's fluid at the point nearest ``t_c``: pH, log fO2 and log activities of its basis."""
    if not result:
        return {}, None, None
    temps = result["axes"]["temperature_c"]["values"]
    for mode in ("static", "flow_through"):
        run = result["modes"].get(mode)
        if not run:
            continue
        good = [i for i, v in enumerate(run["properties"]["pH"]) if v is not None]
        if not good:
            continue
        i = min(good, key=lambda k: (abs(temps[k] - t_c), -k))
        out = {"pH": run["properties"]["pH"][i]}
        fo2 = run["properties"].get("log_fO2")
        if fo2 and fo2[i] is not None:
            out[OXYGEN] = fo2[i]
        for name, vals in run.get("log_a", {}).items():
            if vals[i] is not None:
                out[name] = vals[i]
        return out, i, mode
    return {}, None, None


def compute(db, spec: dict, result: dict | None = None) -> dict:
    """Everything the results window needs to draw ``spec``."""
    spec = dict(spec)
    x, y = _axes(spec)
    notes: list[str] = []
    points = None
    if result:
        temps = result["axes"]["temperature_c"]["values"]
        press = result["axes"]["pressure_bar"]["values"]
    t_c = spec.get("t_c")
    p_bar = spec.get("p_bar")
    if t_c is None:
        t_c = _last(temps) if result else 25.0
        if p_bar is None and result:
            p_bar = _last(press)
    t_c = float(t_c)
    fluid, at, mode = _point_values(result, t_c)
    others: dict[str, float] = {}
    for el in spec.get("with", []):
        master = DEFAULT_MASTERS.get(el, el)
        if master in fluid:
            others[master] = fluid[master]
        else:
            others[master] = dg.DEFAULT_OTHER_LOG_A
            notes.append(f"log a({pretty_species(master, db.species[master].charge if master in db.species else None)}) "
                         f"held at {dg.DEFAULT_OTHER_LOG_A:g}.")
    ph = fluid.get("pH", 7.0)
    log_fo2 = fluid.get(OXYGEN)
    exclude = list(ALWAYS_EXCLUDED)
    if result:
        exclude += result.get("scenario", {}).get("blocked_minerals", [])
    dia = dg.predominance(db, spec["element"], x, y, t_c=t_c, p_bar=p_bar, ph=ph, log_a=spec["log_a"],
                          log_fo2=log_fo2, others=others, include=spec.get("include", "all"),
                          exclude=exclude, nx=_GRID[0], ny=_GRID[1])
    if fluid and at is not None:
        held = []
        if x.kind != "pH" and y.kind != "pH":
            held.append(f"pH {ph:.2f}")
        if spec["kind"] == "activity" and log_fo2 is not None:
            held.append(f"log fO2 {log_fo2:.1f}")
        if others:
            held.append(", ".join(f"log a({pretty_species(k, db.species[k].charge)}) {v:.2f}"
                                  for k, v in others.items() if k in fluid))
        if held:
            notes.append(f"Held at the run's values at {temps[at]:g} C: " + "; ".join(h for h in held if h) + ".")
    else:
        if x.kind != "pH" and y.kind != "pH":
            notes.append(f"pH held at {ph:g}.")
        if spec["kind"] == "activity" and spec["element"] == "Al":
            pass
        elif dia.fixed.get(OXYGEN) is not None and y.kind not in ("Eh", "log_fO2"):
            notes.append(f"log fO2 held at {dia.fixed[OXYGEN]:g}.")
    labels = dg.label_points(dia)
    present = set(dia.grid.ravel().tolist())
    fields = []
    for k, (name, kind) in enumerate(zip(dia.names, dia.kinds)):
        if k not in present:
            continue
        sp = db.species.get(name)
        label = pretty_mineral(name) if kind == "mineral" else pretty_species(name, sp.charge if sp else None)
        cells = int((dia.grid == k).sum())
        fields.append({"index": k, "name": name, "label": label, "kind": kind,
                       "extrapolated": name in dia.extrapolated, "at": labels[name],
                       "small": cells < 0.004 * dia.grid.size})
    small = [f["label"] for f in fields if f["small"]]
    if small:
        notes.append("Too small to label (hover to see them): " + ", ".join(small) + ".")
    uncertain = ["".join("1" if m < dia.band else "0" for m in row) for row in dia.margin]
    if result:
        points = _path(result, x, y)
        if points:
            ts = [p[2] for pts in points.values() for p in pts]
            if max(abs(t - t_c) for t in ts) > 10:
                notes.append(f"The run's path runs from {min(ts):g} to {max(ts):g} C; the fields are "
                             f"drawn for {t_c:g} C, so read the path against them only near that "
                             f"temperature (filled dots).")
    title = (f"{KINDS[spec['kind']]} diagram for {spec['element']}" if spec["kind"] != "activity"
             else f"Activity diagram, {spec['element']}-{pretty_species(x.species, 0)}-"
                  f"{pretty_species(y.species, db.species[y.species].charge)}")
    rel = dia.reliability
    return {
        "spec": spec, "title": f"{title} at {t_c:g} C, {dia.p_bar:.0f} bar",
        "t_c": t_c, "p_bar": dia.p_bar,
        "x": {"label": _axis_label(x, db), "lo": x.lo, "hi": x.hi},
        "y": {"label": _axis_label(y, db), "lo": y.lo, "hi": y.hi},
        "nx": len(dia.xs), "ny": len(dia.ys),
        "grid": dia.grid.tolist(), "uncertain": uncertain, "fields": fields,
        "water": dia.water, "limits": [{**lim, "label": pretty_mineral(lim["mineral"]) + " saturation"}
                                       for lim in dia.limits],
        "reliability": {"level": rel.level if rel else "poor",
                        "deviation": None if rel is None or not math.isfinite(rel.deviation)
                        else round(rel.deviation, 3), "band": dia.band},
        "path": points or {}, "notes": dia.notes + notes,
        "log_a": spec["log_a"], "include": spec.get("include", "all"),
    }


def _axis_label(ax: dg.Axis, db) -> str:
    if ax.kind == "log_a":
        return f"log a({pretty_species(ax.species, db.species[ax.species].charge)})"
    if ax.kind == "log_ratio":
        z = db.species[ax.species].charge
        hp = "a(H+)" if z == 1 else f"a(H+)^{z:g}"
        return f"log a({pretty_species(ax.species, z)})/{hp}"
    return ax.label()


def _path(result: dict, x: dg.Axis, y: dg.Axis) -> dict[str, list]:
    """The run's points in the diagram's coordinates, with their temperatures."""
    temps = result["axes"]["temperature_c"]["values"]
    out = {}
    for mode, run in result["modes"].items():
        props, la = run["properties"], run.get("log_a", {})

        def coord(ax, i):
            if ax.kind == "pH":
                return props["pH"][i]
            if ax.kind == "Eh":
                return (props.get("Eh") or [None] * len(temps))[i]
            if ax.kind == "log_fO2":
                return (props.get("log_fO2") or [None] * len(temps))[i]
            vals = la.get(ax.species)
            if vals is None or vals[i] is None:
                return None
            if ax.kind == "log_a":
                return vals[i]
            z = (result.get("catalog", {}).get("species", {}).get(ax.species) or {}).get("charge") or 1.0
            return vals[i] + z * props["pH"][i] if props["pH"][i] is not None else None
        pts = []
        for i in range(len(temps)):
            a, b = coord(x, i), coord(y, i)
            if a is not None and b is not None:
                pts.append([a, b, temps[i]])
        if pts:
            out[mode] = pts
    return out
