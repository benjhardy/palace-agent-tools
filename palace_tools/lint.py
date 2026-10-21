"""Lint rules for Palace configs. Each rule encodes a failure seen in real projects or a documented Palace rule.

``lint(config, geometry=None, context=None)`` returns findings with an id, severity (error | warning | info), a
message, a fix, the lesson behind it and where that lesson came from. Rules that need geometry (distances,
patch shapes, which faces are outer) use the ``*.tags.json`` written by ``palace_tools.mesh``; without it they are
skipped and listed as such.

Severity: ``error`` = Palace will abort, silently ignore your intent, or give a wrong answer; ``warning`` = very
likely wrong or very costly; ``info`` = advice.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

from . import rf
from .config import count_frequencies, excitations, frequency_range

DOCS = "https://awslabs.github.io/palace/v0.16.0/"  # versioned docs matching the target release
RULES: dict[str, dict] = {}


def rule(rid: str, title: str, lesson: str, source: str, docs: str = "", needs: str = ""):
    def deco(fn):
        RULES[rid] = {"id": rid, "title": title, "lesson": lesson, "source": source,
                      "docs": DOCS + docs if docs else "", "needs": needs, "fn": fn}
        return fn
    return deco


@dataclass
class Ctx:
    cfg: dict
    geometry: dict | None = None
    version: tuple = (0, 16)
    unknowns: int | None = None
    ram_gb: float | None = None
    config_dir: Path | None = None
    spec: dict | None = None
    findings: list = field(default_factory=list)

    @property
    def ptype(self):
        return self.cfg.get("Problem", {}).get("Type")

    @property
    def bnd(self):
        return self.cfg.get("Boundaries", {})

    @property
    def lin(self):
        return self.cfg.get("Solver", {}).get("Linear", {})

    def attr_geom(self) -> dict[int, dict]:
        """attribute -> geometry entry (surfaces and volumes) from tags.json."""
        if not self.geometry:
            return {}
        out = {}
        for grp in ("surfaces", "volumes"):
            for name, e in self.geometry.get(grp, {}).items():
                if isinstance(e, dict) and "tag" in e:
                    out[int(e["tag"])] = dict(e, name=name, kind=grp[:-1])
        return out


def _f(ctx: Ctx, rid: str, severity: str, message: str, fix: str = "", **data):
    r = RULES[rid]
    ctx.findings.append({"id": rid, "severity": severity, "title": r["title"], "message": message, "fix": fix,
                         "lesson": r["lesson"], "source": r["source"], "docs": r["docs"], **data})


def _lumped(ctx):
    return ctx.bnd.get("LumpedPort", [])


def _wave(ctx):
    return ctx.bnd.get("WavePort", [])


def _attrs_of(entry) -> list[int]:
    if "Attributes" in entry:
        return list(entry["Attributes"])
    return [a for el in entry.get("Elements", []) for a in el.get("Attributes", [])]


def _is_direct(ctx) -> bool:
    return ctx.lin.get("MGMaxLevels", 100) == 1


def _closed(ctx) -> bool:
    b = ctx.bnd
    return not (b.get("Absorbing") or b.get("WavePort"))


def _excited_lumped(ctx):
    return [p for p in _lumped(ctx) if p.get("Excitation") not in (None, False, 0)]


# ============================================================================================ solver

@rule("SOL001", "Sparse direct solver: use MUMPS, not SuperLU / 'Default'",
      "SuperLU (what 'Default' resolves to when Palace is built with it, as the 0.16 container is) gave "
      "partition-dependent wrong points (S21 ~ -45 dB at 2 of 21 frequencies on 6 ranks) and 'singular matrix' "
      "aborts in closed boxes. MUMPS gave identical results on any rank count.",
      "gds-bench (G01, validation V4)", "config/solver/")
def _sol001(ctx: Ctx):
    if ctx.ptype not in ("Driven", "Eigenmode"):
        return
    t = ctx.lin.get("Type", "Default")
    if t in ("SuperLU", "Default"):
        _f(ctx, "SOL001", "warning",
           f"Solver.Linear.Type is '{t}'" + (" (resolves to SuperLU in builds with SuperLU_DIST, e.g. the 0.16 "
                                               "container)" if t == "Default" else ""),
           'Set "Solver": {"Linear": {"Type": "MUMPS"}}.')


@rule("SOL002", "Direct factorisation vs p-multigrid",
      "Eigenmode with lumped elements: p-multigrid took ~1 min per shift-and-invert solve at order 2; a direct "
      "MUMPS factorisation (MGMaxLevels 1) did the whole 484k-unknown run in 6.6 min. Driven, 102k unknowns: "
      "full complex MUMPS gave the same S to 7e-9, 5.6x faster than p-multigrid, at 2x the memory. Above roughly "
      "0.6-1 M unknowns the factorisation memory (~15-25 kB per unknown) decides.",
      "birdcage, port-tuning", "config/solver/")
def _sol002(ctx: Ctx):
    if ctx.ptype not in ("Driven", "Eigenmode"):
        return
    direct = _is_direct(ctx)
    has_lumped_lc = any(p.get("L") or p.get("C") or p.get("Ls") or p.get("Cs") for p in _lumped(ctx))
    if ctx.ptype == "Eigenmode" and not direct:
        sev = "warning" if has_lumped_lc else "info"
        _f(ctx, "SOL002", sev, "Eigenmode solve with p-multigrid preconditioning (MGMaxLevels not 1)"
           + (" and lumped L/C elements" if has_lumped_lc else ""),
           'Use "Linear": {"Type": "MUMPS", "MGMaxLevels": 1} unless the factorisation does not fit in RAM.')
    if ctx.unknowns and ctx.ptype == "Driven":
        if not direct and ctx.unknowns < 300_000:
            _f(ctx, "SOL002", "info", f"{ctx.unknowns:,} unknowns with p-multigrid: a direct factorisation is "
               "usually several times faster at this size", 'Try "MGMaxLevels": 1 with MUMPS.')
    if ctx.unknowns and direct:
        per = 22e3 if ctx.ptype == "Driven" else 15e3
        need = ctx.unknowns * per / 1e9
        if ctx.ram_gb and need > 0.8 * ctx.ram_gb:
            _f(ctx, "SOL002", "warning", f"Direct factorisation of {ctx.unknowns:,} unknowns needs ~{need:.0f} GB "
               f"(budget {ctx.ram_gb:g} GB)", "Drop MGMaxLevels (p-multigrid with MUMPS coarse solve) or use "
               "symmetry planes / a coarser mesh.")


@rule("SOL003", "Error-estimator cost when AMR is off",
      "Palace always runs the flux-projection error estimate. With the defaults (Jacobi PCG, EstimatorTol 1e-6) it "
      "took 1850 of 3660 s on an inductor; EstimatorMG true + EstimatorTol 1e-4 gave bit-identical S-parameters "
      "5.8x faster overall.",
      "gds-bench (validation V1/V3)", "config/solver/")
def _sol003(ctx: Ctx):
    if ctx.ptype not in ("Driven", "Eigenmode"):
        return
    amr = ctx.cfg.get("Model", {}).get("Refinement", {}).get("MaxIts", 0)
    if amr:
        return
    if not ctx.lin.get("EstimatorMG", False) or ctx.lin.get("EstimatorTol", 1e-6) < 1e-4:
        _f(ctx, "SOL003", "info", "Default error-estimator settings without AMR (slow, and not needed for S or f)",
           'Add "EstimatorMG": true, "EstimatorTol": 1e-4 to Solver.Linear.')


@rule("SOL004", "Complex coarse solve for driven problems",
      "With MGMaxLevels 1 the 'coarse' solve is the whole problem. ComplexCoarseSolve true factorises the true "
      "complex matrix, so GMRES converges in 1-2 iterations instead of tens.",
      "port-tuning", "config/solver/")
def _sol004(ctx: Ctx):
    if ctx.ptype == "Driven" and _is_direct(ctx) and not ctx.lin.get("ComplexCoarseSolve", False):
        _f(ctx, "SOL004", "info", "Direct preconditioner on the real approximation of a complex system",
           'Set "ComplexCoarseSolve": true.')


@rule("SOL005", "Krylov iteration limit",
      "An iterative solve that stops at MaxIts gives a wrong answer at that frequency; Palace only warns in the log.",
      "Palace docs", "config/solver/")
def _sol005(ctx: Ctx):
    mi = ctx.lin.get("MaxIts", 100)
    if not _is_direct(ctx) and mi < 100:
        _f(ctx, "SOL005", "info", f"Solver.Linear.MaxIts = {mi} with an iterative preconditioner",
           "Use 200-500 and check the log for 'did NOT converge'.")


@rule("SOL006", "Element order",
      "Order 1 shifted a birdcage k=1 mode by 1.25 % against order 2 on the same mesh; order 3 costs ~3x the "
      "unknowns per element of order 2 (21 vs 7 per tet).",
      "birdcage, microstrip-amr", "config/solver/")
def _sol006(ctx: Ctx):
    o = ctx.cfg.get("Solver", {}).get("Order", 1)
    if o == 1 and ctx.ptype in ("Driven", "Eigenmode"):
        _f(ctx, "SOL006", "warning", "Solver.Order is 1 (Palace's default)", "Use Order 2 for production runs.")
    if o >= 4:
        _f(ctx, "SOL006", "warning", f"Solver.Order = {o}: memory grows fast", "Check the cost estimate first.")


# ============================================================================================ ports

@rule("PRT001", "Mixed wave and lumped ports",
      "Palace 0.16 writes no port-S.csv when lumped and wave ports are both present (PostOperator::MeasureSParameter "
      "returns early). Fixed in 0.17. Workarounds: model a resistor as an Impedance sheet (Rs = R w / l) and keep "
      "wave ports, or make every port lumped.",
      "wilkinson, port-tuning", "guide/boundaries/")
def _prt001(ctx: Ctx):
    if ctx.ptype != "Driven" or not (_lumped(ctx) and _wave(ctx)):
        return
    if ctx.version < (0, 17):
        _f(ctx, "PRT001", "error", "LumpedPort and WavePort both present: Palace 0.16 computes no S-parameters",
           "Use an Impedance boundary for resistors/loads, or make all ports lumped (or run Palace >= 0.17).")
    else:
        _f(ctx, "PRT001", "info", "Mixed lumped and wave ports (S-parameters supported from Palace 0.17)")


@rule("PRT002", "Excitation index must equal port index",
      "For one-port-per-excitation runs Palace assumes excitation index == port index; anything else gives no or "
      "mislabelled S columns.",
      "Palace docs (guide/boundaries), source", "guide/boundaries/")
def _prt002(ctx: Ctx):
    if ctx.ptype != "Driven":
        return
    for ex, ports in excitations(ctx.cfg).items():
        if len(ports) == 1 and ports[0][1] != ex:
            _f(ctx, "PRT002", "error", f"{ports[0][0]} {ports[0][1]} is excitation {ex}",
               f'Set "Excitation": {ports[0][1]} (or true).')


@rule("PRT003", "Several ports in one excitation",
      "Ports sharing an excitation index are driven together: Palace writes port-V/port-I but no port-S. Useful for "
      "grouped/combiner drive; compute the metric from port-V.csv and port-I.csv.",
      "gds-bench (G02 combiner)", "guide/boundaries/")
def _prt003(ctx: Ctx):
    if ctx.ptype != "Driven":
        return
    for ex, ports in excitations(ctx.cfg).items():
        if len(ports) > 1:
            _f(ctx, "PRT003", "info", f"Excitation {ex} drives {len(ports)} ports together: no port-S.csv",
               "Intended? Otherwise give each port its own excitation index.")


@rule("PRT004", "Surface current source disables S-parameters",
      "Surface currents are added to every excitation, so no excitation is 'simple' and S is not computed.",
      "Palace source (portexcitations)", "config/boundaries/")
def _prt004(ctx: Ctx):
    if ctx.ptype == "Driven" and ctx.bnd.get("SurfaceCurrent") and (_lumped(ctx) or _wave(ctx)):
        _f(ctx, "PRT004", "warning", "SurfaceCurrent together with ports: no port-S.csv",
           "Excite through ports only if you need S-parameters.")


@rule("PRT005", "Port indices must be unique across port kinds",
      "Index values must be unique across LumpedPort, WavePort, SurfaceCurrent and Terminal (Palace >= 0.14).",
      "Palace docs", "guide/boundaries/")
def _prt005(ctx: Ctx):
    seen = {}
    for kind in ("LumpedPort", "WavePort", "SurfaceCurrent", "Terminal"):
        for p in ctx.bnd.get(kind, []):
            i = p.get("Index")
            if i in seen:
                _f(ctx, "PRT005", "error", f"Index {i} used by {seen[i]} and {kind}", "Renumber the ports.")
            seen[i] = kind


@rule("PRT006", "Driven problem without an excitation",
      "A driven run needs at least one excited port, a SurfaceCurrent or a CurrentDipole.",
      "Palace docs", "guide/problem/")
def _prt006(ctx: Ctx):
    if ctx.ptype != "Driven":
        return
    if not excitations(ctx.cfg) and not ctx.bnd.get("SurfaceCurrent") and not ctx.cfg.get("Domains", {}).get("CurrentDipole"):
        _f(ctx, "PRT006", "error", "No port has Excitation set", 'Set "Excitation": <port index> on the driven port(s).')


@rule("PRT007", "Lumped port values: circuit or per-square, not both, not neither",
      "Palace aborts if a LumpedPort has both R/L/C and Rs/Ls/Cs, or none of them.",
      "Palace source (lumpedportoperator)", "config/boundaries/")
def _prt007(ctx: Ctx):
    for p in _lumped(ctx):
        circ = any(p.get(k) for k in ("R", "L", "C"))
        surf = any(p.get(k) for k in ("Rs", "Ls", "Cs"))
        if circ and surf:
            _f(ctx, "PRT007", "error", f"LumpedPort {p.get('Index')} has both circuit and per-square values",
               "Keep R/L/C (circuit values; Palace converts with the patch l/w).")
        if not circ and not surf:
            _f(ctx, "PRT007", "error", f"LumpedPort {p.get('Index')} has no R/L/C or Rs/Ls/Cs", "Set R, L or C.")


@rule("PRT008", "Excited lumped port must be resistive",
      "Palace <= 0.17 aborts on an excited lumped port with R = 0 or with L/C. Put reactive parts on separate, "
      "non-excited lumped elements (0.18 allows reactive excited ports).",
      "Palace source / changelog 0.18", "config/boundaries/")
def _prt008(ctx: Ctx):
    for p in _excited_lumped(ctx):
        r = p.get("R") or p.get("Rs") or 0
        reac = any(p.get(k) for k in ("L", "C", "Ls", "Cs"))
        if ctx.version < (0, 18) and (r <= 0 or reac):
            _f(ctx, "PRT008", "error", f"Excited LumpedPort {p.get('Index')}: R = {r}, reactive = {reac}",
               "Give the excited port R > 0 only.")


def _deviation_deg(d, axes):
    d = [float(x) for x in d]
    n = math.sqrt(sum(x * x for x in d)) or 1
    d = [x / n for x in d]
    out = []
    for ax in axes:
        c = abs(sum(a * b for a, b in zip(d, ax)))
        out.append(math.degrees(math.acos(min(1.0, c))))
    return out


DIRS = {"+X": [1, 0, 0], "-X": [-1, 0, 0], "+Y": [0, 1, 0], "-Y": [0, -1, 0], "+Z": [0, 0, 1], "-Z": [0, 0, -1],
        "X": [1, 0, 0], "Y": [0, 1, 0], "Z": [0, 0, 1]}


@rule("PRT009", "Lumped element direction vs patch shape",
      "Palace measures a lumped element along its Direction on the patch's bounding box: the direction must be "
      "within 1 degree of a box axis (abort) and the patch must be flat along it (projected-length check). Curved "
      "ring-gap patches failed until made planar; a direction along the patch normal gives zero length.",
      "birdcage", "config/boundaries/", needs="geometry")
def _prt009(ctx: Ctx):
    geo = ctx.attr_geom()
    if not geo:
        return
    for kind in ("LumpedPort", "SurfaceCurrent"):
        for p in ctx.bnd.get(kind, []):
            elems = p.get("Elements") or [{"Attributes": p.get("Attributes", []), "Direction": p.get("Direction")}]
            for el in elems:
                d = el.get("Direction")
                if d is None or (isinstance(d, str) and d.upper().endswith("R")):
                    continue
                vec = DIRS.get(d.upper()) if isinstance(d, str) else d
                if vec is None:
                    _f(ctx, "PRT009", "error", f"{kind} {p.get('Index')}: unknown direction {d!r}")
                    continue
                for a in el.get("Attributes", []):
                    g = geo.get(a)
                    if not g or "axes" not in g:
                        continue
                    if g.get("planar") is False:
                        _f(ctx, "PRT009", "error", f"{kind} {p.get('Index')} (attr {a}, '{g['name']}') is not planar",
                           "Make the lumped patch a flat rectangle.")
                    dev = _deviation_deg(vec, g["axes"])
                    k = min(range(3), key=lambda i: dev[i])
                    lengths = g.get("lengths", [1, 1, 1])
                    if lengths[k] <= 1e-9 * max(lengths):
                        _f(ctx, "PRT009", "error", f"{kind} {p.get('Index')} direction {d} is along the normal of "
                           f"patch '{g['name']}'", "Set Direction along the patch, from one conductor to the other.")
                    elif dev[k] > 1.0:
                        _f(ctx, "PRT009", "error", f"{kind} {p.get('Index')} direction {d} is {dev[k]:.2f} deg off the "
                           f"axes of patch '{g['name']}' (Palace aborts above 1 deg)", "Align Direction with an edge.")
                    elif dev[k] > 0.1:
                        _f(ctx, "PRT009", "warning", f"{kind} {p.get('Index')} direction {d} is {dev[k]:.2f} deg off "
                           f"the patch axes", "Align Direction with an edge.")


@rule("PRT010", "Lumped element units (H, F, ohm; circuit vs per-square)",
      "R/L/C are circuit values in ohm/H/F; Palace converts them to per-square values with the patch: Rs = R w/l, "
      "Ls = L w/l, Cs = C l/w (the log prints the per-square values). A value of 1.8 meant as pF is 1.8 F.",
      "birdcage, port-tuning", "config/boundaries/")
def _prt010(ctx: Ctx):
    for p in _lumped(ctx):
        i = p.get("Index")
        c, l, r = p.get("C"), p.get("L"), p.get("R")
        if c and (c > 1e-6 or c < 1e-19):
            _f(ctx, "PRT010", "error", f"LumpedPort {i}: C = {c:g} F", "C is in farads (1.8 pF = 1.8e-12).")
        if l and (l > 1e-3 or l < 1e-15):
            _f(ctx, "PRT010", "error", f"LumpedPort {i}: L = {l:g} H", "L is in henries (3 nH = 3e-9).")
        if r and (r > 1e6 or r < 1e-3):
            _f(ctx, "PRT010", "warning", f"LumpedPort {i}: R = {r:g} ohm", "Check the value.")
        if any(p.get(k) for k in ("Rs", "Ls", "Cs")):
            _f(ctx, "PRT010", "info", f"LumpedPort {i} uses per-square values (Rs/Ls/Cs)",
               "Prefer circuit values; per-square values do not follow geometry edits.")
    for z in ctx.bnd.get("Impedance", []):
        if z.get("Rs") and z["Rs"] > 1e5:
            _f(ctx, "PRT010", "warning", f"Impedance Rs = {z['Rs']:g} ohm/sq", "Impedance values are per square.")


@rule("PRT011", "Lumped ports with different R in Palace 0.16",
      "Palace 0.16 reports |S_ij|^2 between lumped ports of different R low by exactly R_j/R_i (measured: 50 -> 75 "
      "ohm gave |S21|^2 = 2/3 of the delivered power; a 400-ohm node port needed |S21|^2 x 8). Palace 0.17 fixed it "
      "('incorrect S-parameters when lumped port pairs have different R values', PR 743), so the x R_i/R_j "
      "correction must NOT be applied to 0.17+ results.",
      "lumped_port_line example, gds-bench (G02), Palace changelog 0.17", "config/boundaries/")
def _prt011(ctx: Ctx):
    if ctx.ptype != "Driven" or ctx.version >= (0, 17):
        return
    rs = {p.get("Index"): p.get("R") for p in _lumped(ctx) if p.get("R")}
    if len(set(rs.values())) > 1:
        _f(ctx, "PRT011", "warning", f"Lumped ports with different R: {rs}",
           "Use equal R (e.g. 50 ohm) and renormalise with scikit-rf, or correct |S_ij|^2 by R_i/R_j "
           "(0.16 only), or work from port-V/port-I.")


def _flat_sides(g: dict, dom: list, tol: float) -> list[str]:
    """Domain sides ('-x', '+z', ...) on which a group lies flat. Uses the per-entity 'sides' from tags.json when
    present (a group can span several walls), else the group's bounding box."""
    if g.get("sides") is not None:
        return list(g["sides"])
    bb = g.get("bbox")
    if not bb or not dom:
        return []
    out = []
    for k, ax in enumerate("xyz"):
        if abs(bb[k + 3] - bb[k]) < tol:
            if abs(bb[k] - dom[k]) < tol:
                out.append("-" + ax)
            if abs(bb[k] - dom[k + 3]) < tol:
                out.append("+" + ax)
    return out


@rule("PRT012", "Edge lumped port next to an absorbing wall",
      "With first-order absorbing end walls an edge lumped port lost ~20 % of its power into the wall beside it "
      "(eigenvalues of S^H S 0.45-0.54). A PEC end wall shorts the strip; a PMC (natural, unlisted) end wall is "
      "lossless and does not touch the strip.",
      "port-tuning", "config/boundaries/", needs="geometry")
def _prt012(ctx: Ctx):
    geo = ctx.attr_geom()
    dom = (ctx.geometry or {}).get("domain_bbox")
    ab = set(ctx.bnd.get("Absorbing", {}).get("Attributes", []))
    if not geo or not dom or not ab:
        return
    tol = 1e-6 * max(dom[3] - dom[0], dom[4] - dom[1], dom[5] - dom[2])
    ab_sides = {sd for a in ab if a in geo for sd in _flat_sides(geo[a], dom, tol)}
    for p in _lumped(ctx):
        for a in _attrs_of(p):
            if a in geo:
                hit = [sd for sd in _flat_sides(geo[a], dom, tol) if sd in ab_sides]
                if hit:
                    _f(ctx, "PRT012", "warning", f"LumpedPort {p.get('Index')} lies on the {hit[0]} wall, which has "
                       "absorbing faces", "Make that end wall PMC (leave it unlisted) or move the port inboard.")


@rule("PRT013", "Wave port that does not span the domain face",
      "A 5w x 6h wave port on an open face next to absorbing walls did not hold the whole microstrip field: 3-6 % "
      "loss on a bare line that fell with frequency. Ports spanning the full domain end face (in a shielded box) "
      "fixed it; Palace's 2-D port solve treats unlabelled port edges as PMC.",
      "gds-bench (G01)", "guide/boundaries/", needs="geometry")
def _prt013(ctx: Ctx):
    geo = ctx.attr_geom()
    dom = (ctx.geometry or {}).get("domain_bbox")
    if not geo or not dom or not ctx.bnd.get("Absorbing"):
        return
    for p in _wave(ctx):
        for a in _attrs_of(p):
            g = geo.get(a)
            if not g or not g.get("bbox"):
                continue
            bb = g["bbox"]
            ext = [bb[k + 3] - bb[k] for k in range(3)]
            dext = [dom[k + 3] - dom[k] for k in range(3)]
            n = min(range(3), key=lambda k: ext[k])
            frac = [ext[k] / dext[k] for k in range(3) if k != n]
            if min(frac) < 0.98:
                _f(ctx, "PRT013", "warning", f"WavePort {p.get('Index')} covers {100 * frac[0]:.0f} % x "
                   f"{100 * frac[1]:.0f} % of its domain face, with absorbing walls present",
                   "Make the port the whole end face, or check the loss with a bare-line run first.")


@rule("PRT014", "Wave ports with nonconformal AMR",
      "In Palace 0.16 the port eigen-solve switched to a spurious mode during nonconformal AMR (to the "
      "substrate-filled mode n_eff = sqrt(4.4) at 39k unknowns on a CPW; to n_eff = 1 at 295k on a waveguide tee). "
      "'Nonconformal': false ran to 470k unknowns with a steady kn. Palace 0.17 fixed 'nonconformal AMR erroneous "
      "results when waveports are present'.",
      "microstrip-amr", "config/model/")
def _prt014(ctx: Ctx):
    ref = ctx.cfg.get("Model", {}).get("Refinement", {})
    if _wave(ctx) and ref.get("MaxIts", 0) > 0 and ref.get("Nonconformal", True) and ctx.version < (0, 17):
        _f(ctx, "PRT014", "warning", "Wave ports with nonconformal AMR (Refinement.Nonconformal defaults to true)",
           'Set "Nonconformal": false, and check the kn per pass in the log (palace-tools check).')


@rule("PRT015", "Wave-port de-embedding ignores conductor loss in 0.16",
      "Palace 0.16's port-mode solve treats Conductivity boundaries as PEC (lossless kn). A de-embedded (Offset) "
      "length therefore carries no conductor loss: a 1-inch CPW read -0.001 dB instead of -0.30 dB, and beta is "
      "low by about alpha. 0.17 changed the wave-port formulation to support impedance boundaries.",
      "microstrip-amr", "config/boundaries/")
def _prt015(ctx: Ctx):
    if ctx.version >= (0, 17) or not ctx.bnd.get("Conductivity"):
        return
    off = [p for p in _wave(ctx) if p.get("Offset")]
    if off:
        _f(ctx, "PRT015", "info", f"{len(off)} wave port(s) de-embedded with Offset while conductors are lossy",
           "Model the physical length if you need its conductor loss.")


@rule("PRT016", "Field output with several excitations",
      "Palace 0.16 kept the ParaView fields of the last excitation only (fixed in 0.18). Run one excitation per field "
      "case. ParaView output is ~230 MB per excitation and frequency on a 1.3 M-unknown model.",
      "wilkinson", "config/solver/")
def _prt016(ctx: Ctx):
    if ctx.ptype != "Driven":
        return
    drv = ctx.cfg.get("Solver", {}).get("Driven", {})
    saving = drv.get("Save") or drv.get("SaveStep") or any(s.get("SaveStep") for s in drv.get("Samples", []))
    if saving and len(excitations(ctx.cfg)) > 1 and ctx.version < (0, 18):
        _f(ctx, "PRT016", "warning", "Fields are saved and there are several excitations",
           "Run a separate case per excitation for field plots.")


@rule("PRT017", "Wave ports must be on the outer boundary",
      "Wave ports are one-sided: mesh elements on one side only. An internal face needs a lumped port.",
      "Palace docs", "guide/boundaries/", needs="geometry")
def _prt017(ctx: Ctx):
    geo = ctx.attr_geom()
    for p in _wave(ctx):
        for a in _attrs_of(p):
            if a in geo and geo[a].get("outer") is False:
                _f(ctx, "PRT017", "error", f"WavePort {p.get('Index')} attribute {a} is an internal face",
                   "Put the wave port on the domain boundary, or use a lumped port.")


# ============================================================================================ boundaries

@rule("BND001", "Absorbing boundary too close (distance in wavelengths)",
      "First-order absorbing walls in the reactive near field act as a resistive load and fake a nearly "
      "frequency-flat loss: Wilkinson walls at 8/6 mm (0.21/0.16 lambda0 at 8 GHz) added 0.03 dB; 12/10 and 16/14 mm "
      "agreed within 0.003 dB. Walls 1 mm from a microstrip at 10 GHz gave 5 % flat loss and S11 != S22.",
      "wilkinson, gds-bench (G01)", "config/boundaries/", needs="geometry")
def _bnd001(ctx: Ctx):
    geo = ctx.attr_geom()
    g = ctx.geometry or {}
    sb = g.get("structure_bbox")
    ab = ctx.bnd.get("Absorbing", {}).get("Attributes", [])
    if not geo or not sb or not ab:
        return
    fmin, fmax = frequency_range(ctx.cfg)
    if not fmax:
        return
    l0 = ctx.cfg.get("Model", {}).get("L0", 1e-6)
    lam = rf.wavelength(fmax) / l0  # mesh units
    worst = None
    for a in ab:
        e = geo.get(a)
        if not e or not e.get("bbox"):
            continue
        bb = e["bbox"]
        ext = [bb[k + 3] - bb[k] for k in range(3)]
        n = min(range(3), key=lambda k: ext[k])
        plane = bb[n]
        d = min(abs(plane - sb[n]), abs(plane - sb[n + 3]))
        if worst is None or d < worst[0]:
            worst = (d, e["name"], "xyz"[n])
    if worst is None:
        return
    d, name, ax = worst
    ratio = d / lam
    msg = (f"Absorbing face '{name}' is {d:.4g} mesh units ({ratio:.2f} lambda0 at {fmax:g} GHz) from the structure "
           f"along {ax}")
    if ratio < 0.15:
        _f(ctx, "BND001", "warning", msg, "Move absorbing walls to >= 0.3-0.5 lambda0 (or use a PEC box) and compare "
           "two distances.", distance_lambda=ratio)
    elif ratio < 0.35:
        _f(ctx, "BND001", "info", msg, "Check by moving the walls out ~50 % and comparing loss/S.", distance_lambda=ratio)


@rule("BND002", "Adaptive sweep in a closed lossless box",
      "Palace 0.16's adaptive sweep aborted after the first sample in a closed PEC box ('SuperLUSolver: Found a "
      "singular matrix') on 4 and 6 ranks; direct point solves worked. With absorbing walls the adaptive sweep "
      "matched direct solves within 0.2 dB.",
      "gds-bench (G01)", "config/solver/")
def _bnd002(ctx: Ctx):
    drv = ctx.cfg.get("Solver", {}).get("Driven", {})
    if ctx.ptype != "Driven" or not drv.get("AdaptiveTol"):
        return
    if not ctx.bnd.get("Absorbing"):
        t = ctx.lin.get("Type", "Default")
        sev = "warning" if t in ("SuperLU", "Default") else "info"
        _f(ctx, "BND002", sev, "Adaptive frequency sweep without absorbing boundaries (closed box)",
           "Use MUMPS; if the sweep aborts, use direct point solves (AdaptiveTol 0) or validate a few points.")


@rule("BND003", "Finite-conductivity boundary needs thickness >> skin depth",
      "The Conductivity boundary is a surface-impedance model for conductors much thicker than the skin depth; "
      "below ~3 skin depths it misestimates the loss.",
      "Palace docs (guide/boundaries)", "guide/boundaries/")
def _bnd003(ctx: Ctx):
    fmin, _ = frequency_range(ctx.cfg)
    if not fmin:
        return
    l0 = ctx.cfg.get("Model", {}).get("L0", 1e-6)
    for c in ctx.bnd.get("Conductivity", []):
        t = c.get("Thickness")
        if t:
            dlt = rf.skin_depth(fmin, c["Conductivity"], c.get("Permeability", 1.0))
            if t * l0 < 3 * dlt:
                _f(ctx, "BND003", "warning", f"Conductor thickness {t * l0 * 1e6:.2f} um < 3 skin depths "
                   f"({3 * dlt * 1e6:.2f} um at {fmin:g} GHz)", "Mesh the conductor as a volume, or accept the error.")


@rule("BND004", "Symmetry planes: PEC/PMC choice and what they keep",
      "A PEC plane forces tangential E = 0 (B1 lies in it; currents odd), a PMC plane tangential H = 0 (B1 crosses "
      "it; currents even). PMC is Palace's natural BC, so an unlisted outer face is already PMC. A symmetry model "
      "keeps only modes of matching parity, and a lumped element cut by the plane must be rescaled: cut along its "
      "width -> C/2 (full C gave -7.6 %), cut across its length -> 2C. A driven port must respect the symmetry.",
      "birdcage", "guide/boundaries/")
def _bnd004(ctx: Ctx):
    sym = (ctx.spec or {}).get("symmetry") or []
    geo = ctx.attr_geom()
    names = [g["name"] for g in geo.values() if re.match(r"sym", g.get("name", ""), re.I)]
    if not sym and not names:
        return
    if ctx.ptype == "Driven":
        _f(ctx, "BND004", "warning", "Driven model with symmetry planes",
           "Only symmetric (PMC) or antisymmetric (PEC) excitation is modelled: check that the port drive has that "
           "parity, and that half/quarter ports carry 2x/4x the impedance accounting.")
    else:
        _f(ctx, "BND004", "info", "Eigenmode model with symmetry planes: the spectrum is filtered by parity",
           "Solve the other PEC/PMC combinations if you need the full spectrum; rescale cut lumped elements.")


@rule("BND005", "Attributes that are not in the mesh / volumes without material",
      "Palace only warns 'Unknown ... boundary attributes! Solver will just ignore them!' and runs on: a typo "
      "silently removes a boundary condition. Every volume needs a material.",
      "Palace source (iodata)", "config/boundaries/", needs="geometry")
def _bnd005(ctx: Ctx):
    geo = ctx.attr_geom()
    if not geo:
        return
    surf = {a for a, g in geo.items() if g["kind"] == "surface"}
    vol = {a for a, g in geo.items() if g["kind"] == "volume"}
    used = []
    for k, v in ctx.bnd.items():
        if k == "Postprocessing":
            continue
        if isinstance(v, dict):
            used += [(k, a) for a in v.get("Attributes", [])]
        elif isinstance(v, list):
            for e in v:
                used += [(k, a) for a in _attrs_of(e)]
    for k, a in used:
        if a not in surf:
            _f(ctx, "BND005", "error", f"Boundaries.{k} uses attribute {a}, which is not a surface group in the mesh",
               "Fix the tag (see the *.tags.json).")
    mats = {a for m in ctx.cfg.get("Domains", {}).get("Materials", []) for a in m.get("Attributes", [])}
    for a in mats - vol:
        _f(ctx, "BND005", "error", f"Material attribute {a} is not a volume group in the mesh")
    for a in vol - mats:
        _f(ctx, "BND005", "error", f"Volume '{geo[a]['name']}' (attribute {a}) has no material")
    unassigned = (ctx.geometry or {}).get("outer_unassigned") or []
    if unassigned:
        _f(ctx, "BND005", "info", f"Outer faces with no boundary condition (Palace treats them as PMC): {unassigned}",
           "Intended for PMC walls / symmetry planes; otherwise assign PEC or Absorbing.")


@rule("BND006", "Second-order absorbing boundary support",
      "Order-2 absorbing boundaries are only available for frequency-domain problems; transient supports order 1.",
      "Palace docs", "config/boundaries/")
def _bnd006(ctx: Ctx):
    ab = ctx.bnd.get("Absorbing", {})
    if ab.get("Order", 1) > 1 and ctx.ptype not in ("Driven", "Eigenmode"):
        _f(ctx, "BND006", "error", f"Absorbing Order {ab['Order']} with {ctx.ptype}", "Use Order 1.")


@rule("BND007", "One attribute, two boundary conditions",
      "An attribute listed under two conditions (e.g. PEC and Absorbing) gives an ill-defined model.",
      "Palace docs", "config/boundaries/")
def _bnd007(ctx: Ctx):
    owner: dict[int, list] = {}
    for k in ("PEC", "PMC", "Absorbing", "Impedance", "Conductivity", "LumpedPort", "WavePort"):
        v = ctx.bnd.get(k)
        if not v:
            continue
        for e in (v if isinstance(v, list) else [v]):
            for a in _attrs_of(e):
                owner.setdefault(a, []).append(k)
    for a, ks in owner.items():
        if len(set(ks)) > 1:
            _f(ctx, "BND007", "error", f"Attribute {a} appears in {sorted(set(ks))}", "Give each surface one condition.")


# ============================================================================================ model

@rule("MOD001", "Mesh units (Model.L0)",
      "Palace's default L0 is 1e-6: mesh coordinates are read as micrometres unless L0 is set. A mm mesh without L0 "
      "is 1000x too small and every resonance lands 1000x too high.",
      "Palace docs (config/model)", "config/model/")
def _mod001(ctx: Ctx):
    m = ctx.cfg.get("Model", {})
    if "L0" not in m:
        _f(ctx, "MOD001", "warning", "Model.L0 not set: mesh units are taken as micrometres",
           'Set "L0": 1e-3 for a mm mesh (1.0 for metres).')
    dom = (ctx.geometry or {}).get("domain_bbox")
    _, fmax = frequency_range(ctx.cfg)
    if dom and fmax:
        size = max(dom[3] - dom[0], dom[4] - dom[1], dom[5] - dom[2]) * m.get("L0", 1e-6)
        r = size / rf.wavelength(fmax)
        if r > 60:
            _f(ctx, "MOD001", "warning", f"Domain is {r:.0f} wavelengths across at {fmax:g} GHz: units or frequency "
               "wrong, or a very large model", "Check L0 and the frequency (GHz).")
        if r < 1e-4:
            _f(ctx, "MOD001", "warning", f"Domain is {r:.1e} wavelengths across: check L0 and the frequency unit (GHz)")


@rule("MOD002", "Mesh file exists",
      "A relative mesh path is resolved from Palace's working directory, so run Palace from the config folder.",
      "Palace docs (config/model)", "config/model/")
def _mod002(ctx: Ctx):
    if ctx.config_dir is None:
        return
    mp = ctx.cfg.get("Model", {}).get("Mesh")
    if mp and not Path(mp).is_absolute() and not (ctx.config_dir / mp).exists():
        _f(ctx, "MOD002", "error", f"Mesh '{mp}' not found next to the config", "Fix the path or build the mesh.")


@rule("MOD003", "Output folder",
      "Palace 0.14-0.16 abort with 'Invalid output directory' when Problem.Output is missing (0.17 defaults to "
      "'postpro').", "Palace changelog 0.17", "config/problem/")
def _mod003(ctx: Ctx):
    if not ctx.cfg.get("Problem", {}).get("Output") and ctx.version < (0, 17):
        _f(ctx, "MOD003", "error", "Problem.Output not set", 'Set "Output": "postpro".')


@rule("MOD004", "Field output size",
      "ParaView output is written for every saved step and is large (~230 MB per excitation and frequency at 1.3 M "
      "unknowns). It is not needed for S-parameters or eigenfrequencies.",
      "wilkinson", "config/problem/")
def _mod004(ctx: Ctx):
    sol = ctx.cfg.get("Solver", {})
    n = 0
    if ctx.ptype == "Eigenmode":
        n = sol.get("Eigenmode", {}).get("Save", 0)
    elif ctx.ptype == "Driven":
        drv = sol.get("Driven", {})
        n = len(drv.get("Save", []))
        nf = count_frequencies(ctx.cfg)
        for s in drv.get("Samples", []):
            if s.get("SaveStep"):
                n += max(1, nf // s["SaveStep"])
        n *= max(1, len(excitations(ctx.cfg)))
    if n > 10:
        _f(ctx, "MOD004", "info", f"About {n} field snapshots will be written", "Save only the steps you will plot.")


# ============================================================================================ AMR

@rule("AMR001", "AMR without a size cap",
      "Palace's Dorfler marking (UpdateFraction 0.7) grew meshes 1.3-2.8x per pass; MaxSize is checked after "
      "refining, so a 450k cap ended at 623k unknowns (14 GB). Without MaxSize the last pass can exhaust RAM.",
      "microstrip-amr", "config/model/")
def _amr001(ctx: Ctx):
    ref = ctx.cfg.get("Model", {}).get("Refinement", {})
    if ref.get("MaxIts", 0) > 0 and not ref.get("MaxSize"):
        _f(ctx, "AMR001", "warning", "Refinement.MaxIts > 0 with no MaxSize", "Set MaxSize to ~0.6x the unknowns your "
           "RAM allows (the cap can be overshot by ~1.4x).")


@rule("AMR002", "AMR over a whole frequency sweep",
      "Every AMR pass re-solves every frequency point. The cheaper pattern (HFSS-like): adapt at one or a few "
      "frequencies, then sweep on the adapted mesh (load the saved mesh, or reuse the refinement in a second run).",
      "microstrip-amr", "config/model/")
def _amr002(ctx: Ctx):
    ref = ctx.cfg.get("Model", {}).get("Refinement", {})
    if ctx.ptype == "Driven" and ref.get("MaxIts", 0) > 0:
        n = count_frequencies(ctx.cfg)
        drv = ctx.cfg.get("Solver", {}).get("Driven", {})
        if n > 5 and not drv.get("AdaptiveTol"):
            _f(ctx, "AMR002", "warning", f"AMR with {n} direct frequency points per pass",
               "Adapt at the design frequency, then sweep.")


@rule("AMR003", "AMR convergence is slow on resonant structures",
      "On an 8-rod interdigital filter the band-centre phase of S21 still moved 7 deg between the last two AMR passes "
      "at 623k (p2) and 673k (p3) unknowns; seeding the mesh along the resonators and order 3 got closest per "
      "unknown. Judge convergence on the quantity you need, between passes.",
      "microstrip-amr", "config/model/")
def _amr003(ctx: Ctx):
    ref = ctx.cfg.get("Model", {}).get("Refinement", {})
    if ref.get("MaxIts", 0) > 0 and ref.get("UpdateFraction", 0.7) >= 0.7 and ref.get("Tol", 1e-2) >= 1e-2:
        _f(ctx, "AMR003", "info", "AMR stops on the global indicator norm (Tol 1e-2 default), not on your figure of "
           "merit", "Track S (or f) per pass from postpro/iterationN (palace-tools check --convergence).")


# ============================================================================================ eigen

@rule("EIG001", "Eigenmode target",
      "Palace 0.16 requires Solver.Eigenmode.Target (GHz) and finds the N modes ABOVE it (shift-and-invert). Put it "
      "below the lowest mode you need (the birdcage used 0.55 x f0 for a crowded spectrum).",
      "birdcage, Palace changelog 0.16", "config/solver/")
def _eig001(ctx: Ctx):
    if ctx.ptype != "Eigenmode":
        return
    e = ctx.cfg.get("Solver", {}).get("Eigenmode", {})
    t = e.get("Target")
    if t is None:
        _f(ctx, "EIG001", "error", "Solver.Eigenmode.Target missing", "Set Target (GHz) below the lowest mode wanted.")
        return
    if t > 3000:
        _f(ctx, "EIG001", "error", f"Target = {t:g}: Palace frequencies are in GHz", "Convert to GHz.")
    exp = (ctx.spec or {}).get("eigen", {}).get("expected_ghz")
    if exp and t >= exp:
        _f(ctx, "EIG001", "error", f"Target {t:g} GHz is at/above the expected mode ({exp:g} GHz)",
           "Modes below Target are not found: lower Target.")


@rule("EIG002", "Lossy / nonlinear eigenproblems are slow",
      "A conducting phantom made the birdcage eigenproblem complex: GMRES ~19 iterations per solve and hours for a "
      "dozen modes. Wave ports and Conductivity boundaries make it nonlinear (Palace's Hybrid solver).",
      "birdcage", "config/solver/")
def _eig002(ctx: Ctx):
    if ctx.ptype != "Eigenmode":
        return
    mats = ctx.cfg.get("Domains", {}).get("Materials", [])
    lossy = [m for m in mats if m.get("Conductivity") or m.get("LossTan")]
    nonlin = ctx.bnd.get("WavePort") or ctx.bnd.get("Conductivity") or ctx.bnd.get("Absorbing", {}).get("Order", 1) > 1
    if lossy or nonlin:
        _f(ctx, "EIG002", "info", "Eigenproblem is complex" + (" and nonlinear" if nonlin else ""),
           "Expect several times the cost of the lossless run; estimate first, keep N small.")


@rule("EIG003", "Excitations are ignored in eigenmode runs",
      "Lumped ports in an eigenmode run are circuit elements (R adds port-Q, L/C enter the eigenproblem); Excitation "
      "has no effect. A lumped port with only C set is a lumped capacitor.",
      "birdcage", "guide/problem/")
def _eig003(ctx: Ctx):
    if ctx.ptype == "Eigenmode" and any(p.get("Excitation") for p in _lumped(ctx) + _wave(ctx)):
        _f(ctx, "EIG003", "info", "Excitation set on ports of an eigenmode run (ignored)")


# ============================================================================================ driven sweep

@rule("DRV001", "Frequency values and units",
      "Palace frequencies are in GHz. Zero or negative frequencies are invalid for driven runs.",
      "Palace docs", "config/solver/")
def _drv001(ctx: Ctx):
    if ctx.ptype != "Driven":
        return
    fmin, fmax = frequency_range(ctx.cfg)
    if fmin is None:
        _f(ctx, "DRV001", "error", "No frequencies in Solver.Driven", "Add Samples.")
        return
    if fmin <= 0:
        _f(ctx, "DRV001", "error", f"MinFreq = {fmin:g} GHz", "Start above 0 GHz.")
    if fmax > 3000:
        _f(ctx, "DRV001", "error", f"MaxFreq = {fmax:g}: frequencies are GHz", "Convert to GHz.")


@rule("DRV002", "Adaptive sweep settings",
      "AdaptiveTol 1e-3 kept S within ~1e-3 of direct solves (0.005 dB in a passband) but deep stopbands at "
      "-32...-65 dB differed by up to 0.33 dB. Palace needs more than a handful of points, or reverts to a uniform sweep.",
      "port-tuning, wilkinson", "config/solver/")
def _drv002(ctx: Ctx):
    if ctx.ptype != "Driven":
        return
    drv = ctx.cfg.get("Solver", {}).get("Driven", {})
    tol = drv.get("AdaptiveTol", 0.0)
    if not tol:
        return
    if tol > 1e-2:
        _f(ctx, "DRV002", "warning", f"AdaptiveTol = {tol:g}", "Use 1e-3 (1e-4 for deep stopbands/nulls).")
    if count_frequencies(ctx.cfg) < 6:
        _f(ctx, "DRV002", "warning", "Adaptive sweep with very few points", "Use direct points (AdaptiveTol 0).")


@rule("DRV003", "Adaptive sweep online cost",
      "After building the reduced model, Palace 0.16 re-assembles operators at every frequency of every excitation: "
      "~2.6 s each at 258k unknowns (97 of 128 min for 317 points x 7 excitations); with wave ports it also "
      "re-solves the port mode per point (a quarter of a 90 min run).",
      "port-tuning, wilkinson", "config/solver/")
def _drv003(ctx: Ctx):
    drv = ctx.cfg.get("Solver", {}).get("Driven", {})
    if ctx.ptype != "Driven" or not drv.get("AdaptiveTol"):
        return
    n = count_frequencies(ctx.cfg) * max(1, len(excitations(ctx.cfg)))
    if n > 600:
        _f(ctx, "DRV003", "warning", f"{n} frequency x excitation evaluations in the online phase",
           "Use fewer points (coarser step) or point solves at the frequencies you need.", evaluations=n)


@rule("DRV004", "Many direct frequency points",
      "Each direct point is a full factorisation/solve per excitation; the adaptive sweep needs ~5-12 full solves "
      "per excitation for a smooth band.",
      "wilkinson, microstrip-amr", "config/solver/")
def _drv004(ctx: Ctx):
    drv = ctx.cfg.get("Solver", {}).get("Driven", {})
    if ctx.ptype == "Driven" and not drv.get("AdaptiveTol"):
        n = count_frequencies(ctx.cfg)
        if n > 40:
            _f(ctx, "DRV004", "info", f"{n} direct frequency points x {max(1, len(excitations(ctx.cfg)))} excitations",
               'Consider "AdaptiveTol": 1e-3 (open models; see BND002 for closed boxes).')


# ============================================================================================ driver

def parse_version(v) -> tuple:
    if v is None:
        return (0, 16)
    if isinstance(v, tuple):
        return v
    m = re.search(r"(\d+)\.(\d+)", str(v))
    return (int(m.group(1)), int(m.group(2))) if m else (0, 16)


def lint(config: dict, geometry: dict | None = None, context: dict | None = None) -> dict:
    """Run every rule. context: {palace_version, unknowns, ram_gb, config_dir, spec}."""
    context = context or {}
    ctx = Ctx(cfg=config, geometry=geometry, version=parse_version(context.get("palace_version")),
              unknowns=context.get("unknowns") or ((geometry or {}).get("mesh", {}).get("unknowns", {}) or {}).get(
                  str(config.get("Solver", {}).get("Order", 1))),
              ram_gb=context.get("ram_gb"),
              config_dir=Path(context["config_dir"]) if context.get("config_dir") else None,
              spec=context.get("spec"))
    skipped = []
    for rid, r in RULES.items():
        if r["needs"] == "geometry" and not geometry:
            skipped.append(rid)
            continue
        r["fn"](ctx)
    sev = {"error": 0, "warning": 1, "info": 2}
    ctx.findings.sort(key=lambda f: (sev[f["severity"]], f["id"]))
    counts = {s: sum(1 for f in ctx.findings if f["severity"] == s) for s in sev}
    return {"ok": counts["error"] == 0, "counts": counts, "findings": ctx.findings,
            "skipped_rules_need_geometry": skipped, "palace_version_assumed": ".".join(map(str, ctx.version))}


def rules_table() -> list[dict]:
    return [{k: v for k, v in r.items() if k != "fn"} for r in RULES.values()]
