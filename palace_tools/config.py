"""Build Palace JSON configs from a compact, agent-friendly spec.

The spec uses physical-group *names* from the mesh (resolved through the ``*.tags.json`` file written by
``palace_tools.mesh``), SI-ish units with explicit suffixes, and safe defaults learned the hard way:

* MUMPS, never the "Default" sparse solver (SuperLU in most builds gave partition-dependent wrong answers);
* a direct factorisation (MGMaxLevels 1) for eigenmode problems and for small/medium driven problems;
* the cheap error-estimator settings (EstimatorMG, EstimatorTol 1e-4) when AMR is off;
* each excited port gets its own excitation index equal to its port index (the only layout that gives S);
* lumped elements cut by a symmetry plane are rescaled (half-width -> C/2, L*2, R*2; half-length -> 2C, L/2, R/2);
* resistive sheets are given as R with the sheet's length (along the current) and width: Rs = R * w / l.

See ``SPEC_DOC`` (or ``palace-tools config spec``) for every key.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

from .rf import UNITS

SPEC_DOC = {
    "problem": "driven | eigenmode | electrostatic | magnetostatic | transient (builder fills driven/eigenmode fully)",
    "mesh": "mesh file (Gmsh .msh v2.2 recommended), relative to the config file",
    "units": "mesh length unit: m, cm, mm, um, mil, in (or a number = L0 in metres). Palace's default L0 is 1e-6 (um)!",
    "tags": "name -> attribute map, or path to the *.tags.json written by palace_tools.mesh",
    "order": "FE order (default 2)",
    "output": "output folder (default 'postpro')",
    "materials": "[{tags, eps_r, tan_d, mu_r, sigma}] volume materials",
    "pec / pmc": "[tags] boundaries (PMC is also Palace's natural BC: unlisted outer faces are PMC)",
    "absorbing": "{tags, order: 1|2} first/second-order absorbing (scattering) boundary",
    "conductivity": "[{tags, sigma, mu_r, thickness}] finite-conductivity surface (thick conductor, >> skin depth)",
    "impedance": "[{tags, Rs|Ls|Cs per square}] or [{tags, R, length, width}] -> Rs = R * width / length",
    "lumped_ports": "[{index, tags, direction ('+X'..'-Z', '+R', or [x,y,z]), R (ohm, circuit value), excite}]",
    "lumped_elements": "[{index, tags, direction, R|L|C (circuit values, H/F), cut: none|half-width|half-length}]",
    "wave_ports": "[{index, tags, mode, offset (mesh units, de-embedding), excite}]",
    "wave_port_pec": "[tags] extra PEC edges for the 2-D port-mode solve",
    "symmetry": "[{tags, type: pec|pmc}] symmetry planes (PEC plane: tangential E = 0; PMC: tangential H = 0)",
    "frequency": "{min_ghz, max_ghz, step_ghz | n, points_ghz: [...], adaptive_tol (0 = direct solves), "
                 "adaptive_max_samples, save_ghz: [...]}",
    "eigen": "{target_ghz (required; modes are searched ABOVE it), n, tol, save, expected_ghz}",
    "solver": "{preset: auto|direct|multigrid, linear_type (MUMPS), tol, max_its, unknowns (for 'auto'), "
              "estimator_fast: true}",
    "amr": "{max_its, tol, max_size, nonconformal, update_fraction}",
    "save_fields": "bool or int (eigen: modes to save; driven: SaveStep); default off",
    "probes": "[[x, y, z], ...] field probes in mesh units",
    "energy_domains": "[tags] domain-energy postprocessing (participation)",
    "extra": "raw Palace keys deep-merged last (escape hatch)",
}

DIRECT_UNKNOWNS_LIMIT = 600_000  # auto preset: above this, p-multigrid with a MUMPS coarse solve (memory)


def load_tags(tags) -> dict:
    """Return {name: int}. Accepts a dict (name -> int or -> {"tag": int}) or a *.tags.json path."""
    if tags is None:
        return {}
    if isinstance(tags, (str, Path)):
        tags = json.loads(Path(tags).read_text(encoding="utf-8"))
    flat = {}
    if "surfaces" in tags or "volumes" in tags:
        for grp in ("volumes", "surfaces"):
            for k, v in tags.get(grp, {}).items():
                flat[k] = int(v["tag"] if isinstance(v, dict) else v)
        return flat
    for k, v in tags.items():
        if isinstance(v, dict) and "tag" in v:
            flat[k] = int(v["tag"])
        elif isinstance(v, int):
            flat[k] = v
    return flat


def _attrs(x, tags: dict, what: str) -> list[int]:
    if x is None:
        return []
    if isinstance(x, (int, str)):
        x = [x]
    out = []
    for t in x:
        if isinstance(t, bool):
            raise ValueError(f"{what}: bad tag {t!r}")
        if isinstance(t, int):
            out.append(t)
        elif isinstance(t, str) and t.isdigit():
            out.append(int(t))
        elif t in tags:
            out.append(tags[t])
        else:
            raise KeyError(f"{what}: unknown physical group '{t}'. Known: {sorted(tags)}")
    return sorted(dict.fromkeys(out))


def _l0(units) -> float:
    if isinstance(units, (int, float)):
        return float(units)
    if units not in UNITS:
        raise ValueError(f"units must be one of {sorted(UNITS)} or a number (metres)")
    return UNITS[units]


def _deep_merge(a: dict, b: dict) -> dict:
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(a.get(k), dict):
            _deep_merge(a[k], v)
        else:
            a[k] = copy.deepcopy(v)
    return a


CUT_SCALE = {"none": {"R": 1.0, "L": 1.0, "C": 1.0},
             # a plane through the element along its current: half the width at the same voltage, half the current
             "half-width": {"R": 2.0, "L": 2.0, "C": 0.5},
             # a plane across the element (normal to the current): each half sees half the voltage, full current
             "half-length": {"R": 0.5, "L": 0.5, "C": 2.0}}


def _direction(d):
    if isinstance(d, str):
        return d.upper() if d[-1].lower() != "r" else d[:-1] + "R"
    return [float(v) for v in d]


def build_config(spec: dict, base_dir: Path | None = None) -> dict:
    """Compact spec -> Palace config dict. Raises ValueError/KeyError with a readable message on bad input."""
    sp = copy.deepcopy(spec)
    tags = load_tags(sp.get("tags") if not isinstance(sp.get("tags"), str) or base_dir is None
                     else (Path(base_dir) / sp["tags"]))
    ptype = {"driven": "Driven", "eigenmode": "Eigenmode", "eigen": "Eigenmode", "electrostatic": "Electrostatic",
             "magnetostatic": "Magnetostatic", "transient": "Transient"}[str(sp.get("problem", "driven")).lower()]
    if "mesh" not in sp:
        raise ValueError("spec needs 'mesh'")
    if "units" not in sp:
        raise ValueError("spec needs 'units' (Palace assumes micrometres when Model.L0 is missing)")
    order = int(sp.get("order", 2))
    cfg: dict = {
        "Problem": {"Type": ptype, "Verbose": int(sp.get("verbose", 2)), "Output": sp.get("output", "postpro")},
        "Model": {"Mesh": sp["mesh"], "L0": _l0(sp["units"])},
        "Domains": {"Materials": []},
        "Boundaries": {},
        "Solver": {"Order": order, "Device": "CPU"},
    }
    save = sp.get("save_fields", False)
    if not save:
        cfg["Problem"]["OutputFormats"] = {"Paraview": False, "GridFunction": False}

    # ---- materials
    for m in sp.get("materials", []):
        mat = {"Attributes": _attrs(m.get("tags"), tags, "materials"), "Permittivity": float(m.get("eps_r", 1.0)),
               "Permeability": float(m.get("mu_r", 1.0))}
        if m.get("tan_d"):
            mat["LossTan"] = float(m["tan_d"])
        if m.get("sigma"):
            mat["Conductivity"] = float(m["sigma"])
        cfg["Domains"]["Materials"].append(mat)
    if not cfg["Domains"]["Materials"]:
        raise ValueError("spec needs at least one entry in 'materials'")

    b = cfg["Boundaries"]
    pec = _attrs(sp.get("pec"), tags, "pec")
    pmc = _attrs(sp.get("pmc"), tags, "pmc")
    for sym in sp.get("symmetry", []):
        a = _attrs(sym.get("tags"), tags, "symmetry")
        if sym.get("type", "").lower() == "pec":
            pec += a
        elif sym.get("type", "").lower() == "pmc":
            pmc += a
        else:
            raise ValueError("symmetry entries need type 'pec' or 'pmc'")
    if pec:
        b["PEC"] = {"Attributes": sorted(set(pec))}
    if pmc:
        b["PMC"] = {"Attributes": sorted(set(pmc))}
    ab = sp.get("absorbing")
    if ab:
        if isinstance(ab, list):
            ab = {"tags": ab}
        b["Absorbing"] = {"Attributes": _attrs(ab.get("tags"), tags, "absorbing"), "Order": int(ab.get("order", 1))}
    if sp.get("conductivity"):
        b["Conductivity"] = []
        for c in sp["conductivity"]:
            e = {"Attributes": _attrs(c.get("tags"), tags, "conductivity"), "Conductivity": float(c["sigma"])}
            if c.get("mu_r"):
                e["Permeability"] = float(c["mu_r"])
            if c.get("thickness"):
                e["Thickness"] = float(c["thickness"])
            b["Conductivity"].append(e)
    if sp.get("impedance"):
        b["Impedance"] = []
        for z in sp["impedance"]:
            e = {"Attributes": _attrs(z.get("tags"), tags, "impedance")}
            if "R" in z:  # circuit value on a rectangular sheet: R = Rs * l / w
                if not (z.get("length") and z.get("width")):
                    raise ValueError("impedance with R needs 'length' (along the current) and 'width'")
                e["Rs"] = float(z["R"]) * float(z["width"]) / float(z["length"])
            for k in ("Rs", "Ls", "Cs"):
                if k in z:
                    e[k] = float(z[k])
            b["Impedance"].append(e)
    if sp.get("wave_port_pec"):
        b["WavePortPEC"] = {"Attributes": _attrs(sp["wave_port_pec"], tags, "wave_port_pec")}

    # ---- ports and lumped elements (Index unique across all port kinds; Excitation = Index)
    lps = []
    for p in sp.get("lumped_ports", []):
        e = {"Index": int(p["index"]), "Direction": _direction(p.get("direction", "+Z"))}
        if p.get("elements"):
            e.pop("Direction")
            e["Elements"] = [{"Attributes": _attrs(el["tags"], tags, "lumped port element"),
                              "Direction": _direction(el.get("direction", "+Z"))} for el in p["elements"]]
        else:
            e["Attributes"] = _attrs(p.get("tags"), tags, f"lumped port {p['index']}")
        if p.get("coordinate_system"):
            e["CoordinateSystem"] = p["coordinate_system"]
        for k in ("R", "L", "C", "Rs", "Ls", "Cs"):
            if p.get(k) is not None:
                e[k] = float(p[k])
        if "R" not in e and "Rs" not in e and ptype in ("Driven", "Transient"):
            e["R"] = 50.0
        if p.get("excite", ptype in ("Driven", "Transient")):
            e["Excitation"] = int(p["index"])
        if p.get("active") is False:
            e["Active"] = False
        lps.append(e)
    for el in sp.get("lumped_elements", []):
        cut = el.get("cut", "none")
        if cut not in CUT_SCALE:
            raise ValueError(f"lumped element cut must be one of {list(CUT_SCALE)}")
        e = {"Index": int(el["index"]), "Attributes": _attrs(el.get("tags"), tags, f"lumped element {el['index']}"),
             "Direction": _direction(el.get("direction", "+Z"))}
        for k in ("R", "L", "C"):
            if el.get(k) is not None:
                e[k] = float(el[k]) * CUT_SCALE[cut][k]
        for k in ("Rs", "Ls", "Cs"):
            if el.get(k) is not None:
                e[k] = float(el[k])  # per-square values do not change when the element is cut
        lps.append(e)
    if lps:
        b["LumpedPort"] = sorted(lps, key=lambda e: e["Index"])
    wps = []
    for p in sp.get("wave_ports", []):
        e = {"Index": int(p["index"]), "Attributes": _attrs(p.get("tags"), tags, f"wave port {p['index']}"),
             "Mode": int(p.get("mode", 1))}
        if p.get("offset"):
            e["Offset"] = float(p["offset"])
        if p.get("excite", ptype == "Driven"):
            e["Excitation"] = int(p["index"])
        wps.append(e)
    if wps:
        b["WavePort"] = sorted(wps, key=lambda e: e["Index"])

    # ---- postprocessing
    if sp.get("probes"):
        cfg["Domains"]["Postprocessing"] = {"Probe": [{"Index": i + 1, "Center": [float(c) for c in p]}
                                                      for i, p in enumerate(sp["probes"])]}
    if sp.get("energy_domains"):
        cfg["Domains"].setdefault("Postprocessing", {})["Energy"] = [
            {"Index": i + 1, "Attributes": _attrs(t, tags, "energy_domains")} for i, t in enumerate(sp["energy_domains"])]

    # ---- solver
    sv = sp.get("solver", {})
    amr = sp.get("amr", {})
    preset = sv.get("preset", "auto")
    unknowns = sv.get("unknowns")
    if preset == "auto":
        preset = "direct" if (ptype == "Eigenmode" or not unknowns or unknowns <= DIRECT_UNKNOWNS_LIMIT) else "multigrid"
    lin: dict = {"Type": sv.get("linear_type", "MUMPS"), "KSPType": "GMRES",
                 "Tol": float(sv.get("tol", 1e-10 if ptype == "Eigenmode" else 1e-8)),
                 "MaxIts": int(sv.get("max_its", 200 if preset == "direct" else 400))}
    if preset == "direct":
        lin["MGMaxLevels"] = 1
    if ptype == "Driven":
        lin["ComplexCoarseSolve"] = True
    if sv.get("estimator_fast", True) and not amr.get("max_its"):
        lin["EstimatorMG"] = True
        lin["EstimatorTol"] = 1e-4
    cfg["Solver"]["Linear"] = lin

    if ptype == "Driven":
        fr = sp.get("frequency")
        if not fr:
            raise ValueError("driven spec needs 'frequency'")
        samples = []
        if fr.get("points_ghz"):
            samples.append({"Type": "Point", "Freq": [float(x) for x in fr["points_ghz"]], "SaveStep": 0})
        if fr.get("min_ghz") is not None:
            s = {"Type": "Linear", "MinFreq": float(fr["min_ghz"]), "MaxFreq": float(fr["max_ghz"]), "SaveStep": 0}
            if fr.get("step_ghz"):
                s["FreqStep"] = float(fr["step_ghz"])
            elif fr.get("n"):
                s["NSample"] = int(fr["n"])
            else:
                raise ValueError("frequency range needs step_ghz or n")
            if isinstance(save, int) and not isinstance(save, bool) and save > 0:
                s["SaveStep"] = int(save)
            samples.append(s)
        if not samples:
            raise ValueError("frequency needs points_ghz or min_ghz/max_ghz")
        drv = {"Samples": samples, "AdaptiveTol": float(fr.get("adaptive_tol", 0.0))}
        if fr.get("adaptive_max_samples"):
            drv["AdaptiveMaxSamples"] = int(fr["adaptive_max_samples"])
        if fr.get("save_ghz"):
            drv["Save"] = [float(x) for x in fr["save_ghz"]]
        cfg["Solver"]["Driven"] = drv
    elif ptype == "Eigenmode":
        eg = sp.get("eigen", {})
        tgt = eg.get("target_ghz")
        if tgt is None and eg.get("expected_ghz"):
            tgt = 0.5 * float(eg["expected_ghz"])
        if tgt is None:
            raise ValueError("eigen spec needs target_ghz (modes are searched above it) or expected_ghz")
        e = {"Target": float(tgt), "N": int(eg.get("n", 5)), "Tol": float(eg.get("tol", 1e-8)),
             "Save": int(save) if isinstance(save, int) and not isinstance(save, bool) else (1 if save else 0)}
        if eg.get("type"):
            e["Type"] = eg["type"]
        cfg["Solver"]["Eigenmode"] = e

    if amr.get("max_its"):
        ref = {"MaxIts": int(amr["max_its"]), "Tol": float(amr.get("tol", 1e-2)),
               "MaxSize": int(amr.get("max_size", 0)), "Nonconformal": bool(amr.get("nonconformal", False)),
               "UpdateFraction": float(amr.get("update_fraction", 0.7))}
        cfg["Model"]["Refinement"] = ref

    if sp.get("extra"):
        _deep_merge(cfg, sp["extra"])
    return cfg


def write_config(cfg: dict, path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
    return path


def load_config(cfg) -> dict:
    if isinstance(cfg, dict):
        return cfg
    return json.loads(Path(cfg).read_text(encoding="utf-8"))


def count_frequencies(cfg: dict) -> int:
    drv = cfg.get("Solver", {}).get("Driven", {})
    n = 0
    for smp in drv.get("Samples", []):
        if smp.get("Freq"):
            n += len(smp["Freq"])
        elif smp.get("FreqStep"):
            n += int(round((smp["MaxFreq"] - smp["MinFreq"]) / smp["FreqStep"])) + 1
        elif smp.get("NSample"):
            n += int(smp["NSample"])
    if drv.get("FreqStep") and drv.get("MaxFreq") is not None:
        n += int(round((drv["MaxFreq"] - drv["MinFreq"]) / drv["FreqStep"])) + 1
    return n


def frequency_range(cfg: dict) -> tuple[float | None, float | None]:
    drv = cfg.get("Solver", {}).get("Driven", {})
    fs = []
    for smp in drv.get("Samples", []):
        fs += list(smp.get("Freq", []))
        for k in ("MinFreq", "MaxFreq"):
            if smp.get(k) is not None:
                fs.append(smp[k])
    for k in ("MinFreq", "MaxFreq"):
        if drv.get(k) is not None:
            fs.append(drv[k])
    if not fs:
        eg = cfg.get("Solver", {}).get("Eigenmode", {})
        if eg.get("Target"):
            return float(eg["Target"]), float(eg.get("TargetUpper", 3 * eg["Target"]))
        return None, None
    return float(min(fs)), float(max(fs))


def excitations(cfg: dict) -> dict[int, list[tuple[str, int]]]:
    """{excitation index: [(port kind, port index), ...]} following Palace's bool/int rules."""
    bnd = cfg.get("Boundaries", {})
    ports = [("LumpedPort", p) for p in bnd.get("LumpedPort", [])] + [("WavePort", p) for p in bnd.get("WavePort", [])]
    excited = [(k, p) for k, p in ports if p.get("Excitation") not in (None, False, 0)]
    res: dict[int, list] = {}
    bools = [x for x in excited if x[1]["Excitation"] is True]
    for k, p in excited:
        ex = p["Excitation"]
        if ex is True:
            ex = p["Index"] if len(bools) == 1 else 1
        res.setdefault(int(ex), []).append((k, int(p["Index"])))
    return res
