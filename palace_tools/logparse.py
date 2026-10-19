"""Parse Palace stdout logs (``palace.log``), ``/usr/bin/time -v`` output and Palace's own ``palace.json``.

The patterns were taken from Palace 0.16 logs. Unknown lines are ignored, so newer versions degrade gracefully.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ANSI = re.compile(r"\x1b\[[0-9;]*m")
NUM = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"

# message fragment -> (severity, meaning / what to do). Order matters: first match wins.
KNOWN_MESSAGES = [
    (r"Unknown (\w+) boundary attributes!", "error",
     "Palace ignored boundary attributes that are not in the mesh: a tag typo or a physical group that was not "
     "written. The boundary condition you asked for is missing from the model."),
    (r"Found a singular matrix|singular matrix", "error",
     "Sparse direct solver hit a singular matrix. Seen with SuperLU (Palace's 'Default' when built with it) in "
     "closed boxes and adaptive sweeps: switch Solver.Linear.Type to MUMPS."),
    (r"ERROR RETURN \*\* FROM DMUMPS INFO\(1\)=\s*-9", "error",
     "MUMPS workspace too small (INFO(1) = -9). Usually tiny meshes split over too many ranks; use fewer ranks."),
    (r"ERROR RETURN \*\* FROM DMUMPS INFO\(1\)=\s*-1[03]", "error",
     "MUMPS failed (numerically singular matrix or memory allocation). Check RAM headroom and the boundary setup."),
    (r"ERROR RETURN \*\* FROM DMUMPS", "error", "MUMPS returned an error code; see the INFO(1)/INFO(2) lines."),
    (r"Specified direction does not align sufficiently with bounding box axes", "error",
     "Lumped port / surface current direction is more than 1 degree off every axis of the patch's bounding box. "
     "Make the patch planar and set Direction along one of its edges."),
    (r"should match projected length", "error",
     "Lumped element is not a flat rectangle along its Direction (curved or tilted patch). Make it planar with "
     "Direction along an edge."),
    (r"Excited lumped port must have nonzero resistance", "error", "Excited lumped ports need R > 0 (Palace <= 0.17)."),
    (r"Lumped port excitations do not support nonzero reactance", "error",
     "Palace <= 0.17: an excited lumped port must be purely resistive; put L/C on a separate, non-excited element."),
    (r"has both R/L/C and Rs/Ls/Cs defined", "error", "Use circuit values (R/L/C) or per-square values, not both."),
    (r"Verification failed|MFEM abort|MFEM_ABORT|terminate called|MPI_ABORT|Abort\(\d+\)|Segmentation fault",
     "error", "Palace aborted; read the lines around this message."),
    (r"User specified direction .* does not align with either bounding box", "warning",
     "Lumped element direction is 0.1-1 degree off the patch axes; Palace continues but the length may be off."),
    (r"Linear solver did not converge", "warning",
     "The Krylov solver stopped at MaxIts before Tol. Results at that step may be wrong: raise MaxIts or use a "
     "stronger preconditioner (MUMPS, MGMaxLevels 1 for small/medium problems)."),
    (r"GMRES solver did NOT converge|FGMRES solver did NOT converge", "warning",
     "GMRES did not reach its tolerance in MaxIts iterations."),
    (r"PCG solver did NOT converge", "info",
     "The error-estimator projection (PCG) stopped early; harmless for S-parameters, it only affects AMR marking."),
    (r"Adaptive frequency sweep requires > \d+ total frequency samples", "warning",
     "Too few frequency points for the adaptive sweep; Palace reverted to a uniform sweep."),
    (r"Minimal rational interpolation encountered rank-deficient matrix", "warning",
     "Adaptive sweep (PROM) basis is rank deficient; usually harmless, check the result against a few direct points."),
    (r"One or more external boundary attributes has no associated boundary", "info",
     "Some outer faces have no boundary condition: Palace treats them as PMC (natural). Fine for symmetry planes "
     "and PMC end walls, a bug if you meant PEC or absorbing."),
    (r"mesh faces with no associated boundary element exist on the domain boundary", "warning",
     "Outer faces have no boundary elements in the mesh (not in any physical group). They become PMC. Check that "
     "every outer face is in a physical group (Gmsh: Mesh.SaveAll = 0 drops untagged entities)."),
    (r"Probe \d+ at .* could not be found", "warning", "A probe point lies outside the mesh (or in a removed conductor)."),
    (r"ARPACK eigenvalue solver found only", "warning", "Fewer eigenmodes converged than requested."),
    (r"SLEPc EPS solve did not converge", "warning", "An internal SLEPc solve did not converge."),
    (r"Output folder is not empty", "info", "Output folder is reused; old CSVs are overwritten."),
    (r"is a nonconformal mesh, assumed from previous AMR", "info", "Restarted from an AMR mesh; preprocessing skipped."),
]


def strip_ansi(text: str) -> str:
    return ANSI.sub("", text)


def _last_int(pattern: str, text: str):
    m = re.findall(pattern, text)
    return int(m[-1]) if m else None


def _mem_to_gb(val: str, unit: str) -> float:
    return float(val) * {"K": 1 / 1024 ** 2, "M": 1 / 1024, "G": 1.0, "T": 1024.0}[unit]


def parse_log(text: str) -> dict:
    """Extract version, size, cost, warnings, errors, progress and solver diagnostics from a Palace log."""
    t = strip_ansi(text)
    lines = t.splitlines()
    out: dict = {}
    m = re.search(r"Git changeset ID: (\S+)", t)
    out["palace_version"] = m.group(1) if m else None
    m = re.search(r"Running with (\d+) MPI process", t)
    out["ranks"] = int(m.group(1)) if m else None
    m = re.search(r"(\d+) OpenMP thread", t)
    out["omp_threads"] = int(m.group(1)) if m else None
    # problem type
    if "Configuring SLEPc eigenvalue solver" in t or "Configuring ARPACK eigenvalue solver" in t:
        out["problem"] = "Eigenmode"
    elif re.search(r"ω/2π = ", t) or "Computing adaptive fast frequency response" in t:
        out["problem"] = "Driven"
    elif "terminal-C" in t or "Electrostatic" in t:
        out["problem"] = "Electrostatic"
    else:
        out["problem"] = None
    out["adaptive_sweep"] = "Computing adaptive fast frequency response" in t
    # sizes
    nd = re.findall(r"ND \(p = (\d+)\): (\d+)", t)
    out["order"] = int(nd[-1][0]) if nd else None
    out["unknowns"] = int(nd[-1][1]) if nd else _last_int(r"global unknowns = (\d+)", t)
    levels = re.findall(r"Level (\d+) \(p = (\d+)\): (\d+) unknowns", t)
    out["multigrid_levels"] = len({lv for lv, _, _ in levels}) if levels else None
    m = re.search(r"elements\s+\d+\s+[\d.]+\s+\d+\s+(\d+)", t)
    out["elements"] = int(m.group(1)) if m else None
    m = re.findall(r"Estimated peak per-node memory usage is: .*?Total ([\d.]+)([KMGT])", t)
    out["peak_mem_gb"] = round(_mem_to_gb(*m[-1]), 3) if m else None
    m = re.findall(r"Estimated peak per-rank memory usage is: .*?Max\. ([\d.]+)([KMGT])", t)
    out["peak_mem_rank_max_gb"] = round(_mem_to_gb(*m[-1]), 3) if m else None
    tot = re.findall(r"^Total\s+(" + NUM + r")\s+(" + NUM + r")", t, re.M)
    out["elapsed_s"] = float(tot[-1][1]) if tot else None
    out["timing_s"] = _timing_table(t)
    out["finished"] = bool(tot)
    # linear solver behaviour
    its = [int(x) for x in re.findall(r"(?:GMRES|FGMRES|CG) solver converged in (\d+) iteration", t)]
    out["krylov"] = {"solves": len(its), "max_its": max(its) if its else None,
                     "mean_its": round(sum(its) / len(its), 1) if its else None,
                     "not_converged": len(re.findall(r"(?:GMRES|FGMRES) solver did NOT converge", t))}
    m = re.search(r"Dropping \d+ small entries", t)
    out["direct_solver_used"] = bool(m) or "DMUMPS" in t or "SuperLU" in t or "STRUMPACK" in t
    # driven progress
    steps = re.findall(r"It (\d+)/(\d+): ω/2π = (" + NUM + r") GHz", t)
    if steps:
        out["progress"] = {"step": int(steps[-1][0]), "of": int(steps[-1][1]), "f_ghz": float(steps[-1][2])}
    ad = re.findall(r"Adaptive sampling (converged|did not converge|\w[\w ]*?) with (\d+) frequency samples", t)
    out["adaptive"] = [{"status": s, "samples": int(n)} for s, n in ad]
    out["adaptive_errors"] = [float(x) for x in re.findall(r"Greedy iteration \d+ \(n = \d+\): .*?error = (" + NUM + r")", t)]
    # wave ports: kn per port per frequency (spurious-mode detection)
    out["wave_port_kn"] = _wave_port_kn(lines)
    # eigen
    m = re.search(r"Shift-and-invert σ = (" + NUM + r") GHz", t)
    out["eigen_shift_ghz"] = float(m.group(1)) if m else None
    # AMR
    amr = re.findall(r"Indicator norm = (" + NUM + r"), global unknowns = (\d+)", t)
    out["amr"] = [{"indicator_norm": float(a), "unknowns": int(b)} for a, b in amr]
    m = re.search(r"Completed (\d+) iterations? of adaptive mesh refinement", t)
    out["amr_iterations"] = int(m.group(1)) if m else None
    # lumped elements as Palace understood them (per-square values; the check for unit mistakes)
    out["lumped_per_square"] = re.findall(r"^\s*(\d+): ((?:Rs|Ls|Cs) = .*?), n = \(([^)]*)\)", t, re.M)
    out["messages"] = _messages(lines)
    out["n_errors"] = sum(1 for x in out["messages"] if x["severity"] == "error")
    out["n_warnings"] = sum(1 for x in out["messages"] if x["severity"] == "warning")
    return out


def _timing_table(t: str) -> dict:
    m = re.search(r"Elapsed Time Report \(s\).*?\n=+\n(.*?)\n-+\nTotal", t, re.S)
    if not m:
        return {}
    res = {}
    for ln in m.group(1).splitlines():
        mm = re.match(r"^(\s*)(.+?)\s+(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*$", ln)
        if mm:
            key = ("  " if mm.group(1) else "") + mm.group(2).strip()
            res[key.strip() if not mm.group(1) else "." + mm.group(2).strip()] = float(mm.group(4))
    return res


def _wave_port_kn(lines: list[str]) -> list[dict]:
    res, f = [], None
    for ln in lines:
        m = re.search(r"boundary modes at wave ports for ω/2π = (" + NUM + r") GHz", ln)
        if m:
            # Palace 0.16 labels this value ω/2π but prints ω in rad/ns (2π x the frequency in GHz)
            f = round(float(m.group(1)) / (2 * 3.141592653589793), 6)
            continue
        m = re.search(r"Port (\d+), mode (\d+): kₙ = (" + NUM + r")([-+]" + r"\d*\.?\d+(?:[eE][-+]?\d+)?" + r")?i?", ln)
        if m:
            res.append({"f_ghz": f, "port": int(m.group(1)), "mode": int(m.group(2)), "kn_re": float(m.group(3)),
                        "kn_im": float(m.group(4)) if m.group(4) else 0.0})
    return res


def _messages(lines: list[str]) -> list[dict]:
    found: dict[str, dict] = {}
    for i, ln in enumerate(lines):
        text = ln.strip()
        if text == "--> Warning!" and i + 1 < len(lines):
            text = lines[i + 1].strip()
        for pat, sev, meaning in KNOWN_MESSAGES:
            if re.search(pat, text):
                key = pat
                if key in found:
                    found[key]["count"] += 1
                else:
                    found[key] = {"severity": sev, "line": text[:300], "line_no": i + 1, "meaning": meaning, "count": 1}
                break
    order = {"error": 0, "warning": 1, "info": 2}
    return sorted(found.values(), key=lambda d: (order[d["severity"]], d["line_no"]))


def parse_time_v(text: str) -> dict:
    """GNU ``/usr/bin/time -v`` output: wall time, max RSS of the largest child (one MPI rank, not the sum)."""
    out = {}
    m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", text)
    if m:
        parts = [float(p) for p in m.group(1).split(":")]
        sec = 0.0
        for p in parts:
            sec = sec * 60 + p
        out["wall_s"] = round(sec, 2)
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    if m:
        out["max_rss_one_process_gb"] = round(int(m.group(1)) / 1024 ** 2, 3)
    m = re.search(r"Exit status: (\d+)", text)
    if m:
        out["exit_status"] = int(m.group(1))
    m = re.search(r"Percent of CPU this job got: (\d+)%", text)
    if m:
        out["cpu_percent"] = int(m.group(1))
    return out


def parse_palace_json(path: Path) -> dict:
    """Palace's own run summary (written to the output folder): DOFs, ranks, peak memory, timings, solves."""
    d = json.loads(Path(path).read_text())
    pr = d.get("Problem", {})
    pk = d.get("PeakMemoryMegabytes", {})
    el = d.get("ElapsedTime", {}).get("Durations", {})
    ls = d.get("LinearSolver", {})
    return {"palace_version": d.get("GitTag"), "unknowns": pr.get("DegreesOfFreedom"), "ranks": pr.get("MPISize"),
            "elements": pr.get("MeshElements"), "multigrid_dofs": pr.get("MultigridDegreesOfFreedom"),
            "peak_mem_gb": round(pk.get("Total", 0) / 1024, 3) if pk else None,
            "elapsed_s": el.get("Total"), "linear_solves": ls.get("TotalSolves"), "krylov_its": ls.get("TotalIts"),
            "timing_s": el, "excitations": d.get("Excitations")}


def summarize_run(case_dir: Path, log_name: str = "palace.log") -> dict:
    """Everything known about a finished or running case folder (log + time + palace.json)."""
    case_dir = Path(case_dir)
    res: dict = {"case": str(case_dir)}
    lp = case_dir / log_name
    if lp.exists():
        res["log"] = parse_log(lp.read_text(errors="replace"))
    tp = case_dir / "time.log"
    if tp.exists():
        res["time"] = parse_time_v(tp.read_text(errors="replace"))
    pj = sorted(case_dir.rglob("palace.json"), key=lambda p: p.stat().st_mtime)
    if pj:
        try:
            res["palace_json"] = parse_palace_json(pj[-1])
        except (ValueError, KeyError):
            pass
    return res
