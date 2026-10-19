"""Harvest solve-cost records from finished Palace case folders into palace_tools/data/cost_records.json.

    python tools/harvest_costs.py wilkinson=<dir> birdcage=<dir> ... [--out palace_tools/data/cost_records.json]

Each case folder needs a ``palace.log`` (and ideally ``config.json``, ``time.log`` and the output ``palace.json``).
Only numbers and run names are stored, never paths.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from palace_tools.logparse import summarize_run  # noqa: E402


def find_config(d: Path) -> dict | None:
    cands = [d / "config.json"] + sorted(d.glob("*.json"))
    for c in cands:
        if c.exists():
            try:
                j = json.loads(c.read_text())
            except ValueError:
                continue
            if isinstance(j, dict) and "Solver" in j and "Problem" in j:
                return j
    return None


def record(label: str, d: Path) -> dict | None:
    s = summarize_run(d)
    lg = s.get("log", {})
    pj = s.get("palace_json", {})
    if not lg.get("finished"):
        return None
    cfg = find_config(d) or {}
    sol = cfg.get("Solver", {})
    lin = sol.get("Linear", {})
    bnd = cfg.get("Boundaries", {})
    ptype = cfg.get("Problem", {}).get("Type") or lg.get("problem")
    mats = cfg.get("Domains", {}).get("Materials", [])
    lossy_mat = any(m.get("LossTan") or m.get("Conductivity") for m in mats)
    lossy = bool(bnd.get("Absorbing") or bnd.get("Conductivity") or bnd.get("WavePort") or lossy_mat
                 or any(p.get("R") or p.get("Rs") for p in bnd.get("LumpedPort", [])))
    ports = bnd.get("LumpedPort", []) + bnd.get("WavePort", [])
    exc = {p.get("Excitation") for p in ports if p.get("Excitation") not in (None, False, 0)}
    drv = sol.get("Driven", {})
    nfreq = 0
    for smp in drv.get("Samples", []):
        if smp.get("Freq"):
            nfreq += len(smp["Freq"])
        elif smp.get("FreqStep"):
            nfreq += int(round((smp["MaxFreq"] - smp["MinFreq"]) / smp["FreqStep"])) + 1
        elif smp.get("NSample"):
            nfreq += int(smp["NSample"])
    if drv.get("FreqStep"):
        nfreq += int(round((drv["MaxFreq"] - drv["MinFreq"]) / drv["FreqStep"])) + 1
    mg = lin.get("MGMaxLevels", 100)
    lt = lin.get("Type", "Default")
    order = sol.get("Order") or lg.get("order")
    direct = (mg == 1 and lt in ("MUMPS", "SuperLU", "STRUMPACK", "Default")) or lg.get("multigrid_levels") == 1
    amr_its = cfg.get("Model", {}).get("Refinement", {}).get("MaxIts", 0) or 0
    rec = {
        "project": label, "run": d.name, "problem": ptype, "order": order,
        "linear_type": lt, "mg_max_levels": mg, "solver_mode": "direct" if direct else "multigrid",
        "complex": bool(ptype == "Driven" or lossy), "wave_ports": bool(bnd.get("WavePort")),
        "n_ports": len(ports), "n_excitations": max(1, len(exc)) if ptype == "Driven" else 0,
        "adaptive_tol": drv.get("AdaptiveTol", 0.0) if ptype == "Driven" else None,
        "n_freq": nfreq if ptype == "Driven" else None,
        "eigen_n": sol.get("Eigenmode", {}).get("N") if ptype == "Eigenmode" else None,
        "amr_its": amr_its,
        "estimator_mg": lin.get("EstimatorMG", False), "estimator_tol": lin.get("EstimatorTol", 1e-6),
        "unknowns": pj.get("unknowns") or lg.get("unknowns"), "elements": pj.get("elements") or lg.get("elements"),
        "ranks": pj.get("ranks") or lg.get("ranks"),
        "peak_mem_gb": pj.get("peak_mem_gb") or lg.get("peak_mem_gb"),
        "elapsed_s": pj.get("elapsed_s") or lg.get("elapsed_s"),
        "wall_s": s.get("time", {}).get("wall_s"),
        "linear_solves": pj.get("linear_solves") or lg.get("krylov", {}).get("solves"),
        "krylov_its_mean": lg.get("krylov", {}).get("mean_its"),
        "palace_version": pj.get("palace_version") or lg.get("palace_version"),
    }
    if not cfg:  # no config next to the log (e.g. gds2palace runs): fall back to what the log says
        rec["adaptive_tol"] = 1e-3 if lg.get("adaptive_sweep") else 0.0
        rec["n_freq"] = (lg.get("progress") or {}).get("of")
        rec["n_excitations"] = max(1, len(lg.get("adaptive") or []))
        rec["solver_mode"] = "direct" if lg.get("multigrid_levels") == 1 else "multigrid"
        rec["config_found"] = False
    if not rec["unknowns"] or not rec["elapsed_s"]:
        return None
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("roots", nargs="+", help="label=directory")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "palace_tools/data/cost_records.json"))
    a = ap.parse_args()
    recs = []
    for spec in a.roots:
        label, root = spec.split("=", 1)
        for log in sorted(Path(root).rglob("palace.log")):
            r = record(label, log.parent)
            if r:
                recs.append(r)
    note = ("Palace 0.16 solve costs recorded on one desktop (8 threads, 32 GB; WSL2 sees 25 GB). Some runs shared "
            "the CPU with another Palace job, so wall times are upper bounds. peak_mem_gb is Palace's own estimate "
            "summed over ranks; elapsed_s is Palace's total; wall_s includes container start-up.")
    Path(a.out).write_text(json.dumps({"note": note, "machine": {"threads": 8, "ram_gb": 32, "wsl_ram_gb": 25},
                                       "records": recs}, indent=1))
    print(f"{len(recs)} records -> {a.out}")


if __name__ == "__main__":
    main()
