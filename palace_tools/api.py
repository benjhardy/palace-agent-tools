"""JSON-in / JSON-out operations shared by the CLI and the MCP server.

Every function takes one dict of arguments and returns a JSON-serialisable dict. ``TOOLS`` holds the name,
description and input schema of each, which the MCP server publishes as-is.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from . import __version__, checks, config as cfgmod, estimate as est, lint as lintmod, logparse, results, schema


def _load_json(x):
    if x is None or isinstance(x, dict):
        return x
    return json.loads(Path(x).read_text(encoding="utf-8"))


def _clean(o):
    """Make numpy / complex / NaN values JSON-safe."""
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return None if math.isnan(f) or math.isinf(f) else f
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, complex):
        return {"re": o.real, "im": o.imag}
    return o


# ------------------------------------------------------------------------------------------------ config

def build_config(a: dict) -> dict:
    spec = _load_json(a["spec"])
    base = Path(a["out"]).parent if a.get("out") else (Path(a["base_dir"]) if a.get("base_dir") else None)
    cfg = cfgmod.build_config(spec, base_dir=base)
    res = {"config": cfg}
    if a.get("out"):
        res["config_path"] = str(cfgmod.write_config(cfg, Path(a["out"])))
    geo = None
    tags = spec.get("tags")
    if isinstance(tags, str) and base is not None and (base / tags).exists():
        geo = _load_json(base / tags)
    res["lint"] = lintmod.lint(cfg, geo, {"palace_version": a.get("palace_version"), "spec": spec,
                                          "config_dir": str(base) if base else None, "ram_gb": a.get("ram_gb")})
    res["schema"] = schema.validate(cfg)
    return res


def lint_config(a: dict) -> dict:
    cfg = _load_json(a["config"])
    cdir = a.get("config_dir") or (str(Path(a["config"]).parent) if isinstance(a["config"], str) else None)
    geo = _load_json(a.get("geometry"))
    res = lintmod.lint(cfg, geo, {"palace_version": a.get("palace_version"), "unknowns": a.get("unknowns"),
                                  "ram_gb": a.get("ram_gb"), "config_dir": cdir, "spec": _load_json(a.get("spec"))})
    res["schema"] = schema.validate(cfg)
    res["ok"] = res["ok"] and res["schema"]["valid"] is not False
    return res


def validate_config(a: dict) -> dict:
    return schema.validate(_load_json(a["config"]))


def list_rules(a: dict) -> dict:
    return {"rules": lintmod.rules_table()}


def spec_doc(a: dict) -> dict:
    return {"spec_keys": cfgmod.SPEC_DOC, "example": EXAMPLE_SPEC}


EXAMPLE_SPEC = {
    "problem": "driven", "mesh": "line.msh", "units": "mm", "tags": "line.tags.json", "order": 2,
    "materials": [{"tags": ["substrate"], "eps_r": 3.55, "tan_d": 0.0027}, {"tags": ["air"], "eps_r": 1.0}],
    "pec": ["ground", "strip"],
    "lumped_ports": [{"index": 1, "tags": ["port1"], "direction": "+Z", "R": 50},
                     {"index": 2, "tags": ["port2"], "direction": "+Z", "R": 50}],
    "frequency": {"min_ghz": 1, "max_ghz": 10, "step_ghz": 0.1, "adaptive_tol": 1e-3},
}


# ------------------------------------------------------------------------------------------------ mesh

def mesh_template(a: dict) -> dict:
    from . import mesh

    t = a.pop("template")
    fn = {"rectangular_cavity": mesh.rectangular_cavity, "microstrip_line": mesh.microstrip_line,
          "parallel_plate_line": mesh.parallel_plate_line}[t]
    out = a.pop("out")
    info = fn(out, **a)
    return {"tags_file": str(Path(out).with_suffix(".tags.json")), "mesh_file": str(Path(out).with_suffix(".msh")),
            "groups": {g: {n: v["tag"] for n, v in info[g].items()} for g in ("volumes", "surfaces")},
            "mesh": info["mesh"], "domain_bbox": info["domain_bbox"], "structure_bbox": info["structure_bbox"]}


def mesh_report(a: dict) -> dict:
    from . import mesh

    g = mesh._gmsh()
    if not g.isInitialized():
        g.initialize()
    g.option.setNumber("General.Terminal", 0)
    try:
        g.open(str(a["mesh"]))
        rep = mesh.mesh_report(g)
        rep["physical_groups"] = [{"dim": d, "tag": t, "name": g.model.getPhysicalName(d, t)}
                                  for d, t in g.model.getPhysicalGroups()]
    finally:
        g.finalize()
    return rep


# ------------------------------------------------------------------------------------------------ cost

def estimate_cost(a: dict) -> dict:
    unknowns = a.get("unknowns")
    tags = _load_json(a.get("tags"))
    cfg = _load_json(a.get("config"))
    if unknowns is None and tags and cfg:
        unknowns = tags.get("mesh", {}).get("unknowns", {}).get(str(cfg.get("Solver", {}).get("Order", 1)))
    if cfg:
        r = est.estimate_from_config(cfg, unknowns=unknowns, tets=a.get("tets"), ranks=a.get("ranks", 4),
                                     ram_gb=a.get("ram_gb"))
    else:
        r = est.estimate(unknowns=unknowns, tets=a.get("tets"), problem=a.get("problem", "Driven"),
                         order=a.get("order", 2), solver_mode=a.get("solver_mode", "direct"),
                         n_freq=a.get("n_freq", 1), n_excitations=a.get("n_excitations", 1),
                         adaptive=a.get("adaptive", False), eigen_n=a.get("eigen_n", 10), ranks=a.get("ranks", 4),
                         ram_gb=a.get("ram_gb"))
    if a.get("calibration"):
        r["fits"] = est.calibration_summary()
    return r


# ------------------------------------------------------------------------------------------------ run

def run(a: dict) -> dict:
    from . import run as runmod

    est_mem = a.get("est_mem_gb")
    case = Path(a["case_dir"])
    cfgname = a.get("config", "config.json")
    if est_mem is None and (case / cfgname).exists():
        cfg = _load_json(case / cfgname)
        tags = next(iter(case.glob("*.tags.json")), None)
        unk = None
        if tags:
            unk = _load_json(tags).get("mesh", {}).get("unknowns", {}).get(str(cfg.get("Solver", {}).get("Order", 1)))
        if unk:
            e = est.estimate_from_config(cfg, unknowns=unk, ranks=a.get("ranks", 4))
            est_mem = e["memory_gb"]["typ"]
    return runmod.run_case(case, cfgname, a.get("ranks", 4), a.get("backend"), a.get("distro"), a.get("image"),
                           a.get("cmd"), a.get("force", False), est_mem, a.get("timeout_s"), a.get("dry_run", False))


def status(a: dict) -> dict:
    from . import run as runmod

    res = {"host": runmod.host_status(a.get("backend"), a.get("distro"))}
    if a.get("case_dir"):
        res["case"] = runmod.progress(Path(a["case_dir"]))
    if a.get("ranks"):
        res["guard"] = runmod.guard(int(a["ranks"]), a.get("est_mem_gb"), res["host"])
    return res


def doctor(a: dict) -> dict:
    from . import run as runmod

    mods = {}
    for m in ("gmsh", "skrf", "jsonschema"):
        try:
            mod = __import__(m)
            mods[m] = getattr(mod, "__version__", "ok")
        except ImportError:
            mods[m] = None
    return {"palace_tools": __version__, "python_modules": mods, "palace": runmod.which_palace()}


# ------------------------------------------------------------------------------------------------ results

def _out_dir(a: dict) -> tuple[Path, dict | None]:
    case = Path(a["case_dir"]) if a.get("case_dir") else None
    cfg = _load_json(a.get("config")) if a.get("config") else None
    if cfg is None and case and (case / "config.json").exists():
        cfg = _load_json(case / "config.json")
    if a.get("output_dir"):
        od = Path(a["output_dir"])
    elif case:
        od = case / ((cfg or {}).get("Problem", {}).get("Output") or "postpro")
    else:
        raise ValueError("give output_dir or case_dir")
    return od, cfg


def parse_results(a: dict) -> dict:
    od, cfg = _out_dir(a)
    res = results.summarize_output(od, cfg)
    case = Path(a["case_dir"]) if a.get("case_dir") else od.parent
    run = logparse.summarize_run(case)
    if "log" in run:
        lg = run["log"]
        res["log"] = {k: lg.get(k) for k in ("palace_version", "problem", "order", "unknowns", "ranks", "peak_mem_gb",
                                             "elapsed_s", "finished", "adaptive", "amr", "messages", "krylov")}
    if "time" in run:
        res["time"] = run["time"]
    if a.get("touchstone") and (od / "port-S.csv").exists():
        sp = results.read_port_s(od / "port-S.csv")
        z = a.get("z0")
        if z is None and cfg:
            zr = results.port_impedances(cfg)
            vals = [zr.get(p) for p in sp["ports"]]
            z = [v if v else 50.0 for v in vals]
            if any(v is None for v in vals):
                res["touchstone_note"] = ("wave-port S is normalised to each port's mode, not to an ohmic reference; "
                                          "written with z0 = 50 as a label only")
        res["touchstone"] = results.to_touchstone(sp, Path(a["touchstone"]), z)
    return _clean(res)


def check(a: dict) -> dict:
    od, cfg = _out_dir(a)
    out = []
    exp = a.get("expect", {})
    if (od / "port-S.csv").exists():
        sp = results.read_port_s(od / "port-S.csv")
        ref = None
        if a.get("reference_touchstone"):
            import skrf as rf

            ref = rf.Network(a["reference_touchstone"]).s
        out += checks.run_all_sparams(sp, a.get("symmetric_pairs"), ref, a.get("reference_tol", 0.02))
    elif cfg and cfg.get("Problem", {}).get("Type") == "Driven":
        out.append({"name": "port_s_present", "status": "fail", "value": None,
                    "detail": "no port-S.csv (mixed lumped/wave ports in 0.16, grouped excitation, or the run failed)"})
    if (od / "eig.csv").exists():
        modes = results.read_eig(od / "eig.csv")
        out.append(checks.eigen_errors(modes))
        if exp.get("eigen_ghz"):
            out.append(checks.eigen_expected(modes, exp["eigen_ghz"], exp.get("rel_tol", 0.01), exp.get("q_min")))
    case = Path(a["case_dir"]) if a.get("case_dir") else od.parent
    if (case / "palace.log").exists():
        lg = logparse.parse_log((case / "palace.log").read_text(encoding="utf-8", errors="replace"))
        errs = [m for m in lg["messages"] if m["severity"] in ("error", "warning")]
        out.append({"name": "log", "status": "fail" if lg["n_errors"] else ("warn" if errs else "pass"),
                    "value": {"errors": lg["n_errors"], "warnings": lg["n_warnings"]},
                    "detail": "; ".join(m["line"] for m in errs[:5]) or "no errors or warnings",
                    "finished": lg.get("finished")})
        if lg.get("wave_port_kn"):
            eps_max = max((m.get("Permittivity", 1.0) if not isinstance(m.get("Permittivity"), list)
                           else max(m["Permittivity"]) for m in (cfg or {}).get("Domains", {}).get("Materials", [])),
                          default=None)
            out.append(checks.wave_port_modes(lg["wave_port_kn"], eps_r_max=eps_max))
    if a.get("convergence"):
        c = a["convergence"]
        out.append(checks.convergence(c["values"], c.get("unknowns"), c.get("tol", 0.01), c.get("kind", "relative")))
    hist = results.amr_history(od)
    if len(hist) > 1 and any("port_s" in h for h in hist):
        vals = [h["port_s"]["s"] for h in hist if "port_s" in h]
        if all(v.shape == vals[0].shape for v in vals):
            r = checks.convergence(vals, [h.get("unknowns") for h in hist], a.get("amr_tol", 0.02), kind="absolute")
            r["name"] = "amr_convergence_max_dS"
            out.append(r)
    return _clean({"overall": checks.overall(out), "checks": out})


def gotchas(a: dict) -> dict:
    return {"gotchas": [{"id": r["id"], "title": r["title"], "lesson": r["lesson"], "source": r["source"],
                         "docs": r["docs"]} for r in lintmod.rules_table()],
            "also": OTHER_GOTCHAS}


OTHER_GOTCHAS = [
    {"title": "WSL kills detached jobs", "lesson": "When the last wsl.exe exits, WSL stops the VM and kills nohup'd or "
     "&-backgrounded Palace runs. Keep a foreground wsl.exe alive for the whole run (palace-tools run blocks; run "
     "that command itself in the background).", "source": "gds-bench, rfpa-open"},
    {"title": "wsl.exe re-parses bash -lc arguments", "lesson": "Loop variables in 'wsl bash -lc \"for ...\"' expand "
     "early; put loops in a script file.", "source": "wilkinson"},
    {"title": "Windows default encoding", "lesson": "Python on Windows writes cp1252 by default and breaks on Greek "
     "letters in Palace logs (omega, kn): always pass encoding='utf-8'.", "source": "wilkinson"},
    {"title": "Unequal lumped-port R in Palace 0.16", "lesson": "Between lumped ports of different R, 0.16 reports "
     "|S21|^2 low by R1/R2 (50 -> 400 ohm: power gain = |S21|^2 x 8). This is a 0.16 bug fixed in 0.17; after the "
     "upgrade the x8 correction would be wrong. Prefer equal R and renormalise in scikit-rf.",
     "source": "gds-bench (G02), lumped_port_line example"},
    {"title": "Unresolved vias create gain", "lesson": "A single 0.9 um via under a 70 um pad gave S22 = +0.8 dB; a "
     "solid post fixed it. Any |S| > 1 means a modelling error.", "source": "gds-bench (G01)"},
    {"title": "Floating ports and keep-outs", "lesson": "A via port inside a ground keep-out floats: |S21| ~ -50 dB and "
     "no error from the solver. Check port grounding geometrically before a 70-minute solve.", "source": "rfpa-open"},
    {"title": "Holes need booleans", "lesson": "Multi-loop addPlaneSurface gave an invalid solid (holes counted as "
     "material); cut holes with OCC booleans.", "source": "gds-bench"},
    {"title": "Wave-port log frequency", "lesson": "Palace 0.16 prints 'boundary modes at wave ports for w/2pi = X GHz' "
     "with X = 2 pi f (rad/ns); divide by 2 pi.", "source": "Palace 0.16 source"},
    {"title": "Coarse meshes shift resonances down", "lesson": "A coarse Wilkinson mesh put the minima 1-2 % low; "
     "converge on the figure of merit, not only on the AMR indicator.", "source": "wilkinson"},
    {"title": "Symmetry pays", "lesson": "Quarter / eighth birdcage models: 4.5x / 7.9x less memory, 5x faster, and "
     "the filtered spectrum makes mode identification easier.", "source": "birdcage"},
    {"title": "Internal lumped ports are not calibrated", "lesson": "A lumped sheet across a pad gap has its own "
     "self-inductance (~0.2 nH for 0.5 mm of 0.6 mm-wide line): port-tuned values carry that error.",
     "source": "port-tuning"},
]


# ------------------------------------------------------------------------------------------------ registry

def _obj(props: dict, required=()):
    return {"type": "object", "properties": props, "required": list(required), "additionalProperties": True}


S_STR = {"type": "string"}
S_INT = {"type": "integer"}
S_NUM = {"type": "number"}
S_BOOL = {"type": "boolean"}
S_ANY = {}

TOOLS = {
    "palace_spec": {"fn": spec_doc, "description": "Describe the compact config spec accepted by palace_build_config, "
                    "with an example.", "schema": _obj({})},
    "palace_build_config": {
        "fn": build_config,
        "description": "Build a Palace JSON config from a compact spec (physical-group names, units, ports, sweep) "
                       "with safe defaults (MUMPS, direct solve for eigen/small driven, cheap estimator, "
                       "Excitation = Index). Writes it if 'out' is given, and returns lint + schema results.",
        "schema": _obj({"spec": {"type": ["object", "string"], "description": "spec dict or path to a JSON file"},
                        "out": {"type": "string", "description": "config.json path to write"},
                        "palace_version": S_STR, "ram_gb": S_NUM}, ["spec"])},
    "palace_lint_config": {
        "fn": lint_config,
        "description": "Validate a Palace config against the 0.16 schema and run the lint rules (solver choice, ports, "
                       "boundaries, units, AMR, sweeps). Pass geometry (the *.tags.json) for distance/patch checks.",
        "schema": _obj({"config": {"type": ["object", "string"]}, "geometry": {"type": ["object", "string"]},
                        "palace_version": S_STR, "unknowns": S_INT, "ram_gb": S_NUM, "spec": S_ANY}, ["config"])},
    "palace_lint_rules": {"fn": list_rules, "description": "List every lint rule with its lesson and source.",
                          "schema": _obj({})},
    "palace_gotchas": {"fn": gotchas, "description": "The gotchas list (lint rules plus workflow lessons).",
                       "schema": _obj({})},
    "palace_mesh_template": {
        "fn": mesh_template,
        "description": "Build and mesh a template model with Gmsh: rectangular_cavity(a,b,d,h), "
                       "parallel_plate_line(w,h,length,eps_r,h_mesh), microstrip_line(w,h,length,eps_r,side,top,"
                       "h_mesh,h_edge,ports='lumped'|'wave',f_ghz). Writes .msh + .tags.json, returns groups and a "
                       "mesh report with exact unknowns per order.",
        "schema": _obj({"template": {"type": "string", "enum": ["rectangular_cavity", "parallel_plate_line",
                                                                 "microstrip_line"]},
                        "out": S_STR}, ["template", "out"])},
    "palace_mesh_report": {"fn": mesh_report,
                           "description": "Quality (minSICN, gamma), tets and exact Nedelec unknowns of a .msh file.",
                           "schema": _obj({"mesh": S_STR}, ["mesh"])},
    "palace_estimate_cost": {
        "fn": estimate_cost,
        "description": "Estimate peak RAM and wall time (quiet/typical/busy) from unknowns or a config + tags file, "
                       "calibrated on ~200 recorded Palace 0.16 runs.",
        "schema": _obj({"config": {"type": ["object", "string"]}, "tags": {"type": ["object", "string"]},
                        "unknowns": S_INT, "tets": S_INT, "problem": S_STR, "order": S_INT, "solver_mode": S_STR,
                        "n_freq": S_INT, "n_excitations": S_INT, "adaptive": S_BOOL, "eigen_n": S_INT,
                        "ranks": S_INT, "ram_gb": S_NUM, "calibration": S_BOOL})},
    "palace_status": {"fn": status, "description": "Running Palace jobs, CPUs and free RAM on the Palace host; "
                      "progress of a case folder; whether a run with 'ranks' may start now.",
                      "schema": _obj({"case_dir": S_STR, "ranks": S_INT, "backend": S_STR, "distro": S_STR})},
    "palace_run": {
        "fn": run,
        "description": "Run Palace on case_dir/config (blocks until done). Refuses if another Palace job would be "
                       "oversubscribed or RAM is short, unless force. Logs palace.log, time.log, run_record.json.",
        "schema": _obj({"case_dir": S_STR, "config": S_STR, "ranks": S_INT, "backend": S_STR, "distro": S_STR,
                        "image": S_STR, "cmd": S_STR, "force": S_BOOL, "dry_run": S_BOOL, "timeout_s": S_NUM},
                       ["case_dir"])},
    "palace_results": {
        "fn": parse_results,
        "description": "Summarise a Palace output folder (S-parameters, eigenmodes, EPR/Q tables, AMR, fields, log) "
                       "and optionally write Touchstone with the right per-port reference impedances.",
        "schema": _obj({"case_dir": S_STR, "output_dir": S_STR, "config": {"type": ["object", "string"]},
                        "touchstone": S_STR, "z0": S_ANY})},
    "palace_check": {
        "fn": check,
        "description": "Sanity checks: passivity, reciprocity, energy balance (fake wall loss), port symmetry, "
                       "eigenmode at the expected frequency, log errors, wave-port spurious modes, convergence.",
        "schema": _obj({"case_dir": S_STR, "output_dir": S_STR, "config": {"type": ["object", "string"]},
                        "expect": _obj({"eigen_ghz": S_NUM, "rel_tol": S_NUM, "q_min": S_NUM}),
                        "symmetric_pairs": {"type": "array"}, "reference_touchstone": S_STR,
                        "convergence": _obj({"values": {"type": "array"}, "unknowns": {"type": "array"},
                                             "tol": S_NUM})})},
    "palace_doctor": {"fn": doctor, "description": "Which optional modules and Palace backends are available.",
                      "schema": _obj({})},
}


def call(name: str, args: dict | None) -> dict:
    if name not in TOOLS:
        raise KeyError(f"unknown tool {name}; known: {sorted(TOOLS)}")
    return _clean(TOOLS[name]["fn"](dict(args or {})))
