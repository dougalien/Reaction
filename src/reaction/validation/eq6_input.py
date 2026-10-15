"""Reader for EQ6 menu-style input files (``.6i``).

Only what a reaction-path run needs is extracted: the initial composition of the system
(`Mass Balance Totals`), reactants and their rate laws, temperature and pressure options,
the run limits, fixed fugacities, mineral suppression and the system-model switches.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Eq6Reactant:
    name: str
    kind: str                       # Pure mineral, Special reactant, Aqueous species, Gas species ...
    status: str
    remaining: float
    destroyed: float
    composition: dict[str, float] = field(default_factory=dict)    # element -> stoichiometric number
    reaction: dict[str, float] = field(default_factory=dict)       # basis species -> coefficient
    forward_law: str = "Relative rate equation"
    rkb: tuple[float, float, float] = (1.0, 0.0, 0.0)
    backward_law: str = "Partial equilibrium"


@dataclass
class Eq6Input:
    path: str
    description: str
    temperature: tuple[int, float, float, dict]      # (option, base C, derivative, extras)
    pressure: tuple[int, float, float]
    reactants: list[Eq6Reactant]
    xi_start: float
    xi_max: float
    iopt: dict[int, int]
    fixed_fugacity: list[tuple[str, float, float]]   # (gas, moles to add, log fugacity)
    suppress: list[str]                              # nxmod "Suppress" entries
    other_alter: list[tuple[str, str, float]]
    subset_suppression: list[tuple[str, str]]        # nxopt (option, defining species)
    subset_exceptions: list[str]
    totals: dict[str, tuple[float, float]]           # component -> (equilibrium system, aqueous) moles
    log_moles: dict[str, float]                      # starting values of the basis variables
    activity_model: str
    ph_scale: str
    switches: list[tuple[str, str]]
    es_solids: dict[str, float] = field(default_factory=dict)
    # species listed among the mass-balance rows but flagged "Make non-basis" (O2(aq), H2(aq)):
    # their totals are reported separately from the basis species whose reactions contain them
    non_basis: list[str] = field(default_factory=list)


def _num(text: str) -> float:
    return float(text.strip().replace("D", "E"))


def _checked(lines: list[str], start: int, span: int) -> int | None:
    for ln in lines[start:start + span]:
        m = re.search(r"\[x\]\s*\(\s*(-?\d+)\)", ln)
        if m:
            return int(m.group(1))
    return None


def _find(lines: list[str], pattern: str, start: int = 0) -> int:
    for i in range(start, len(lines)):
        if re.search(pattern, lines[i]):
            return i
    raise ValueError(f"pattern {pattern!r} not found")


def read_eq6_input(path: str | Path) -> Eq6Input:
    path = Path(path)
    lines = path.read_text(encoding="latin-1").splitlines()
    text = "\n".join(lines)
    m = re.search(r'Description=\s*"([^"]*)"', text)
    description = m.group(1) if m else ""

    # --- temperature and pressure -----------------------------------------------------------------
    it = _find(lines, r"Temperature option \(jtemp\)")
    ip = _find(lines, r"Pressure option \(jpress\)")
    block = lines[it:ip]
    jt = _checked(lines, it, 14)
    labelled = [(m.group(1), _num(m.group(2))) for m in re.finditer(
        r"\|\s*(Value \(C\)|Base Value \(C\)|Derivative|T of fluid 1 \(C\)|T of fluid 2 \(C\)"
        r"|Mass ratio factor)\s*\|\s*(\S+)\|", "\n".join(lines[it:ip]))]
    base = deriv = 0.0
    extras: dict = {}
    if jt == 0:
        base = next(v for k, v in labelled if k == "Value (C)")
    elif jt in (1, 2):
        bases = [v for k, v in labelled if k == "Base Value (C)"]
        derivs = [v for k, v in labelled if k == "Derivative"]
        base, deriv = bases[jt - 1], derivs[jt - 1]
    elif jt == 3:
        got = dict(labelled)
        extras = {"t1": got["T of fluid 1 (C)"], "t2": got["T of fluid 2 (C)"],
                  "ratio": got["Mass ratio factor"]}
    temperature = (jt, base, deriv, extras)

    ir = _find(lines, r"Reactants \(Irreversible Reactions\)")
    jp = _checked(lines, ip, 14)
    plabelled = [(m.group(1), _num(m.group(2))) for m in re.finditer(
        r"\|\s*(Value \(bars\)|Base Value \(bars\)|Derivative)\s*\|\s*(\S+)\|",
        "\n".join(lines[ip:ir]))]
    pbase = pderiv = 0.0
    if jp == 2:
        pbase = next(v for k, v in plabelled if k == "Value (bars)")
    elif jp in (3, 4):
        pbases = [v for k, v in plabelled if k == "Base Value (bars)"]
        pderivs = [v for k, v in plabelled if k == "Derivative"]
        pbase, pderiv = pbases[jp - 3], pderivs[jp - 3]
    pressure = (jp, pbase, pderiv)

    # --- reactants --------------------------------------------------------------------------------------
    ie = _find(lines, r"Starting, minimum, and maximum values", ir)
    reactants: list[Eq6Reactant] = []
    current: Eq6Reactant | None = None
    section = None
    for ln in lines[ir + 1:ie]:
        m = re.match(r"\|Reactant\s*\|(.*?)\s*\|\s*\(ureac", ln)
        if m:
            current = Eq6Reactant(m.group(1).strip(), "", "", 0.0, 0.0)
            reactants.append(current)
            section = None
            continue
        if current is None:
            continue
        m = re.match(r"\|->\|\s*(Type|Status)\s*\|(.*?)\s*\|", ln)
        if m:
            setattr(current, "kind" if m.group(1) == "Type" else "status", m.group(2).strip())
            continue
        m = re.match(r"\|->\|\s*Amount (remaining|destroyed) \(moles\)\s*\|\s*(\S+)\|", ln)
        if m:
            setattr(current, "remaining" if m.group(1) == "remaining" else "destroyed", _num(m.group(2)))
            continue
        if re.match(r"\|->\|\s*Composition", ln):
            section = "composition"
            continue
        if re.match(r"\|->\|\s*Reaction\s*\|", ln):
            section = "reaction"
            continue
        m = re.match(r"\|--->\|\s*(\S.*?)\s*\|\s*(\S+)\|\s*\(", ln)
        if m and section in ("composition", "reaction"):
            if section == "composition":
                current.composition[m.group(1)] = _num(m.group(2))
            else:
                current.reaction[m.group(1)] = _num(m.group(2))
            continue
        m = re.match(r"\|->\|\s*Forward rate law\s*\|(.*?)\s*\|", ln)
        if m:
            current.forward_law = m.group(1).strip()
            section = "forward"
            rkb = []
            continue
        m = re.match(r"\|--->\|\s*d(\d?)Xi\(n\)/dXi\d?.*?\|\s*(\S+)\|", ln)
        if m and section == "forward":
            rkb.append(_num(m.group(2)))
            if len(rkb) == 3:
                current.rkb = tuple(rkb)
            continue
        m = re.match(r"\|->\|\s*Backward rate law\s*\|(.*?)\s*\|", ln)
        if m:
            current.backward_law = m.group(1).strip()
            section = None

    reactants = [r for r in reactants if r.name != "None"]

    def scalar(label: str) -> float:
        return _num(re.search(rf"\|{label}\s*\|\s*(\S+)\|", text).group(1))

    xi_start, xi_max = scalar("Starting Xi value"), scalar("Maximum Xi value")

    # --- iopt switches -------------------------------------------------------------------------------
    iopt: dict[int, int] = {}
    for i, ln in enumerate(lines):
        m = re.match(r"\|iopt\((\d+)\)\s*-", ln)
        if m:
            sel = _checked(lines, i + 1, 8)
            if sel is not None:
                iopt[int(m.group(1))] = sel

    # --- fixed fugacity, suppression ------------------------------------------------------------------------
    iff = _find(lines, r"Fixed Fugacity Options")
    inp = _find(lines, r"Numerical Parameters", iff)
    fixed = []
    for ln in lines[iff:inp]:
        m = re.match(r"\|(\S[^|]*?)\s*\|\s*(\S+)\|\s*(\S+)\|\s*--", ln)
        if m and m.group(1) not in ("None", "Gas") and not m.group(1).startswith("("):
            fixed.append((m.group(1), _num(m.group(2)), _num(m.group(3))))

    ix = _find(lines, r"Mineral Sub-Set Selection Suppression Options")
    ixe = _find(lines, r"Exceptions to the Mineral Sub-Set", ix)
    subset = []
    for ln in lines[ix:ixe]:
        m = re.match(r"\|(None|All|Alwith|Allwith)\s*\|(.*?)\s*\|", ln)
        if m and m.group(1) != "None":
            subset.append((m.group(1), m.group(2).strip()))
    exceptions = []
    ixe_end = _find(lines, r"Fixed Fugacity Options", ixe)
    for ln in lines[ixe:ixe_end]:
        m = re.match(r"\|(\S[^|]*?)\s*\|\s*\(uxopex", ln)
        if m and m.group(1) != "None":
            exceptions.append(m.group(1))

    suppress, other = [], []
    ia = _find(lines, r"Alter/Suppress Options")
    iae = _find(lines, r"^\|Iopg Activity", ia)
    for ln in lines[ia:iae]:
        m = re.match(r"\|(\S[^|]*?)\s*\|(Suppress|Replace|AugmentLogK|AugmentG)\s*\|\s*(\S+)\|", ln)
        if m:
            if m.group(2) == "Suppress":
                suppress.append(m.group(1))
            else:
                other.append((m.group(1), m.group(2), _num(m.group(3))))

    act = {-1: "davies", 0: "bdot", 1: "pitzer", 2: "hc"}[_checked(lines, _find(lines, r"iopg\(1\)") + 1, 6)]
    scale = {-1: "internal", 0: "nbs", 1: "mesmer"}[_checked(lines, _find(lines, r"iopg\(2\)") + 1, 5)]

    # --- initial composition --------------------------------------------------------------------------------------
    totals: dict[str, tuple[float, float]] = {}
    it0 = _find(lines, r"Mass Balance Totals \(moles\)")
    it1 = _find(lines, r"Electrical imbalance", it0)
    for ln in lines[it0:it1]:
        m = re.match(r"\|(\S+)\s+Aqueous\s*\|\s*(\S+)\|\s*(\S+)\|", ln)
        if m:
            totals[m.group(1)] = (_num(m.group(2)), _num(m.group(3)))
    log_moles: dict[str, float] = {}
    ic = _find(lines, r"Matrix Column Variables and Values")
    ic1 = _find(lines, r"Phases and Species in the PRS", ic)
    for ln in lines[ic:ic1]:
        m = re.match(r"\|(\S+)\s+Aqueous solution\s*\|\s*(\S+)\|", ln)
        if m:
            log_moles[m.group(1)] = _num(m.group(2))

    non_basis = []
    im = _find(lines, r"Mass Balance Species \(Matrix Row Variables\)")
    for ln in lines[im:it0]:
        m = re.match(r"\|(\S+)\s+Aqueous solution\s*\|\s*Make non-basis", ln)
        if m:
            non_basis.append(m.group(1))

    switches = []
    i0 = _find(lines, r"Special Basis Switches")
    for a, b in zip(re.findall(r"\|Replace\s*\|(.*?)\s*\|", "\n".join(lines[i0:i0 + 6])),
                    re.findall(r"\|\s*with\s*\|(.*?)\s*\|", "\n".join(lines[i0:i0 + 6]))):
        if a != "None":
            switches.append((a, b))

    return Eq6Input(
        path=str(path), description=description, temperature=temperature, pressure=pressure,
        reactants=reactants, xi_start=xi_start, xi_max=xi_max, iopt=iopt, fixed_fugacity=fixed,
        suppress=suppress, other_alter=other, subset_suppression=subset,
        subset_exceptions=exceptions, totals=totals, log_moles=log_moles,
        activity_model=act, ph_scale=scale, switches=switches, non_basis=non_basis)
