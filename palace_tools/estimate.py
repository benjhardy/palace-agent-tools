"""Cost estimate (unknowns, peak RAM, wall time) for a Palace run, calibrated on recorded Palace 0.16 runs.

The calibration set (``data/cost_records.json``) holds ~200 finished runs from six projects on one desktop
(8 threads, 32 GB, WSL2 with 25 GB): eigenmode and driven, direct (MUMPS, MGMaxLevels 1) and p-multigrid, orders
1-3, 1k-1.6M unknowns. Power laws ``y = a (N / 1e5)^b`` are fitted per class in log space. Wall time on a shared
machine varies a lot (the same 156k eigen run took 41 s idle and 239 s next to another job), so the estimate gives
a quiet-machine value (10th percentile of the residuals), a typical value and a busy-machine value (85th).

Rule-of-thumb anchors from the records: 1.12 M unknowns, order 2, eigenmode, MUMPS direct, 4 ranks: 16.4 GB and
11-19 min for 20 modes. Driven order 2 with p-multigrid: ~8-12 kB per unknown; full complex MUMPS: ~25-30 kB.
"""
from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent / "data" / "cost_records.json"

# unknowns per tetrahedron for Nedelec order p (interior ratio E/T ~ 1.2, F/T ~ 2.1 for Gmsh meshes);
# exact counts come from the mesh report (palace_tools.mesh.count_unknowns)
DOF_PER_TET = {1: 1.25, 2: 6.7, 3: 20.5, 4: 46.0}


def nedelec_unknowns(n_edges: int, n_faces: int, n_tets: int, order: int) -> int:
    """Exact number of Nedelec (first kind) DOFs of order p on a tetrahedral mesh: p E + p(p-1) F + p(p-1)(p-2)/2 T.
    Matches Palace's 'ND (p = 2)' line (e.g. 212,703 edges and 349,706 faces -> 1,124,818)."""
    p = order
    return p * n_edges + p * (p - 1) * n_faces + p * (p - 1) * (p - 2) // 2 * n_tets


@lru_cache(maxsize=1)
def records() -> list[dict]:
    return json.loads(DATA.read_text(encoding="utf-8"))["records"]


def _fit(rows, y, x=lambda r: r["unknowns"], extras=()):
    """Least squares in log space: log y = log a + b log(N/1e5) + sum_k c_k log(extra_k / ref_k).
    extras: tuples (name, function of record, reference value)."""
    N = np.array([x(r) for r in rows], float)
    Y = np.array([y(r) for r in rows], float)
    ok = np.isfinite(Y) & (Y > 0) & (N > 0)
    E = []
    for _, fn, _ref in extras:
        e = np.array([fn(r) or np.nan for r in rows], float)
        ok &= np.isfinite(e) & (e > 0)
        E.append(e)
    if ok.sum() < 3 + len(extras):
        return None
    cols = [np.ones(ok.sum()), np.log(N[ok] / 1e5)]
    for (_, _, ref), e in zip(extras, E):
        cols.append(np.log(e[ok] / ref))
    A = np.column_stack(cols)
    c, *_ = np.linalg.lstsq(A, np.log(Y[ok]), rcond=None)
    res = np.log(Y[ok]) - A @ c
    return {"log_a": float(c[0]), "b": float(c[1]),
            "extras": [{"name": n, "ref": ref, "c": float(ck)} for (n, _, ref), ck in zip(extras, c[2:])],
            "q20": float(np.quantile(res, 0.1)), "q85": float(np.quantile(res, 0.85)), "sd": float(res.std()),
            "n": int(ok.sum()), "N_range": [int(N[ok].min()), int(N[ok].max())]}


def _eval(fit, N, **extra):
    base = fit["log_a"] + fit["b"] * math.log(N / 1e5)
    for e in fit.get("extras", []):
        base += e["c"] * math.log(max(extra.get(e["name"], e["ref"]), 1e-9) / e["ref"])
    return {"low": math.exp(base + fit["q20"]), "typ": math.exp(base), "high": math.exp(base + fit["q85"])}


RANKS = ("ranks", lambda r: r.get("ranks"), 4.0)
MODES = ("modes", lambda r: r.get("eigen_n"), 10.0)


@lru_cache(maxsize=None)
def fits() -> dict:
    R = [r for r in records() if r.get("unknowns") and r.get("elapsed_s")]
    out = {"memory": {}, "time": {}}
    for prob in ("Driven", "Eigenmode"):
        for mode in ("direct", "multigrid"):
            for order in (1, 2, 3):
                rows = [r for r in R if r["problem"] == prob and r["solver_mode"] == mode and r["order"] == order
                        and r.get("peak_mem_gb")]
                f = _fit(rows, lambda r: r["peak_mem_gb"])
                if f:
                    out["memory"][f"{prob}/{mode}/p{order}"] = f
            rows = [r for r in R if r["problem"] == prob and r["solver_mode"] == mode and r.get("peak_mem_gb")]
            f = _fit(rows, lambda r: r["peak_mem_gb"])
            if f:
                out["memory"][f"{prob}/{mode}/any"] = f
    for mode in ("direct", "multigrid"):
        pts = [r for r in R if r["problem"] == "Driven" and r["solver_mode"] == mode and not r["amr_its"]
               and not r["adaptive_tol"] and r.get("linear_solves")]
        f = _fit(pts, lambda r: r["elapsed_s"] / r["linear_solves"], extras=())
        if f:
            out["time"][f"Driven/{mode}/per_solve"] = f
    # adaptive sweeps: HDM solves cost as point solves; the rest is the online phase per (frequency x excitation)
    for mode in ("direct", "multigrid"):
        t1 = out["time"].get(f"Driven/{mode}/per_solve") or out["time"].get("Driven/direct/per_solve")
        ad = [r for r in R if r["problem"] == "Driven" and r["solver_mode"] == mode and not r["amr_its"]
              and r["adaptive_tol"] and r.get("linear_solves") and r.get("n_freq")]
        rows = []
        for r in ad:
            hdm = r["linear_solves"] * _eval(t1, r["unknowns"], ranks=r.get("ranks") or 4)["typ"] if t1 else 0
            on = (r["elapsed_s"] - hdm) / (r["n_freq"] * max(1, r["n_excitations"]))
            rows.append(dict(r, _on=max(on, 0.02)))  # small sweeps: online cost below the HDM-fit noise
        f = _fit(rows, lambda r: r["_on"], extras=())
        if f:
            out["time"][f"Driven/{mode}/online_per_point"] = f
        f = _fit(ad, lambda r: r["linear_solves"] / max(1, r["n_excitations"]), x=lambda r: r["unknowns"])
        if f:
            out["time"][f"Driven/{mode}/hdm_per_excitation"] = f
    eig = [r for r in R if r["problem"] == "Eigenmode" and not r["amr_its"]]
    f = _fit(eig, lambda r: r["elapsed_s"], extras=(MODES,))
    if f:
        out["time"]["Eigenmode/total"] = f
    return out


def _pick(table: dict, *keys):
    for k in keys:
        if k in table:
            return k, table[k]
    return None, None


def estimate(unknowns: int | None = None, problem: str = "Driven", order: int = 2, solver_mode: str = "direct",
             n_freq: int = 1, n_excitations: int = 1, adaptive: bool = False, eigen_n: int = 10,
             ranks: int = 4, ram_gb: float | None = None, tets: int | None = None, amr_passes: int = 0,
             cpu_threads: int | None = None) -> dict:
    """Estimate peak RAM (GB, summed over ranks) and wall time (minutes)."""
    problem = {"driven": "Driven", "eigenmode": "Eigenmode", "eigen": "Eigenmode"}.get(problem.lower(), problem)
    notes = []
    if unknowns is None:
        if not tets:
            raise ValueError("give unknowns, or tets + order")
        unknowns = int(tets * DOF_PER_TET.get(order, DOF_PER_TET[2]))
        notes.append(f"unknowns estimated from {tets:,} tets at order {order} ({DOF_PER_TET.get(order)} per tet)")
    F = fits()
    mk, mf = _pick(F["memory"], f"{problem}/{solver_mode}/p{order}", f"{problem}/{solver_mode}/any",
                   f"Driven/{solver_mode}/p{order}", f"Driven/{solver_mode}/any")
    mem = _eval(mf, unknowns)
    mem = {k: round(v + 0.075 * max(0, ranks - 4), 2) for k, v in mem.items()}  # ~75 MB baseline per extra rank
    res: dict = {"unknowns": int(unknowns), "problem": problem, "order": order, "solver_mode": solver_mode,
                 "ranks": ranks, "memory_gb": mem, "memory_model": mk}
    if mf and not (mf["N_range"][0] <= unknowns <= mf["N_range"][1] * 1.5):
        notes.append(f"{unknowns:,} unknowns is outside the calibrated range {mf['N_range']} for {mk}: extrapolating")
    if problem == "Eigenmode":
        tk, tf = "Eigenmode/total", F["time"].get("Eigenmode/total")
        t = _eval(tf, unknowns, modes=max(1, eigen_n), ranks=ranks)
        res["solves"] = {"modes": eigen_n}
        if solver_mode != "direct":
            t = {k: v * 5 for k, v in t.items()}
            notes.append("eigenmode with p-multigrid: calibrated on direct runs, x5 applied (birdcage: ~1 min per "
                         "shift-and-invert solve with p-multigrid)")
    else:
        tk, t1f = _pick(F["time"], f"Driven/{solver_mode}/per_solve", "Driven/direct/per_solve")
        t1 = _eval(t1f, unknowns, ranks=ranks)
        if adaptive:
            _, hf = _pick(F["time"], f"Driven/{solver_mode}/hdm_per_excitation", "Driven/direct/hdm_per_excitation")
            hdm = _eval(hf, unknowns)["typ"] if hf else 8.0
            hdm = min(max(hdm, 4.0), 20.0)
            n_hdm = hdm * n_excitations
            _, of = _pick(F["time"], f"Driven/{solver_mode}/online_per_point", "Driven/direct/online_per_point")
            on = _eval(of, unknowns, ranks=ranks) if of else {"low": 0.1, "typ": 0.5, "high": 2.0}
            n_on = n_freq * n_excitations
            t = {k: n_hdm * t1[k] + n_on * on[k] for k in t1}
            res["solves"] = {"full_solves_expected": round(n_hdm, 1), "online_evaluations": n_on,
                             "online_s_per_point": round(on["typ"], 2)}
            tk = f"{tk} + online"
        else:
            n = n_freq * n_excitations
            t = {k: n * v for k, v in t1.items()}
            res["solves"] = {"full_solves": n, "s_per_solve": round(t1["typ"], 1)}
    if amr_passes:
        # each pass solves on a mesh ~1.3-2.8x larger than the last; the final pass dominates
        t = {k: v * 1.8 for k, v in t.items()}
        notes.append(f"{amr_passes} AMR passes: time scaled x1.8 for the earlier, smaller passes; give the final "
                     "(capped) size as unknowns")
    rank_factor = 1.0
    notes.append("wall time is calibrated on 2-8 ranks of one 8-thread desktop, where extra ranks bought little; "
                 "it is not scaled with the rank count")
    if unknowns < 50_000 and ranks > 4:
        rank_factor = 1.5
        notes.append("small problem on many ranks: 8 ranks were slower than 4 on an 8-thread desktop")
    res["wall_min"] = {k: round(v * rank_factor / 60, 2) for k, v in t.items()}
    res["time_model"] = tk
    advice = []
    if ram_gb:
        if mem["high"] > 0.85 * ram_gb:
            advice.append(f"peak RAM may reach {mem['high']:.1f} GB of {ram_gb:g} GB: "
                          + ("use p-multigrid (drop MGMaxLevels 1), " if solver_mode == "direct" else "")
                          + "symmetry planes or a coarser mesh")
        res["fits_in_ram"] = mem["typ"] < 0.85 * ram_gb
    if cpu_threads and ranks > cpu_threads:
        advice.append(f"{ranks} ranks on {cpu_threads} threads oversubscribes the CPU")
    if unknowns < 50_000 and ranks > 2:
        advice.append("tiny problem: 1-2 ranks are enough (MUMPS can also fail with INFO(1) = -9 on tiny meshes "
                      "split over many ranks)")
    if solver_mode == "multigrid" and unknowns < 300_000:
        advice.append("direct MUMPS (MGMaxLevels 1) is usually faster at this size if RAM allows")
    if res["wall_min"]["typ"] > 60:
        advice.append("long run: write a plan with this estimate first, check the mesh on a coarse version")
    res["advice"] = advice
    res["notes"] = notes
    res["calibration"] = {"records": len(records()), "machine": "8-thread desktop, 32 GB (WSL2 25 GB), Palace 0.16"}
    return res


def estimate_from_config(cfg: dict, unknowns: int | None = None, tets: int | None = None, ranks: int = 4,
                         ram_gb: float | None = None) -> dict:
    from .config import count_frequencies, excitations

    sol = cfg.get("Solver", {})
    lin = sol.get("Linear", {})
    prob = cfg.get("Problem", {}).get("Type", "Driven")
    order = sol.get("Order", 1)
    mode = "direct" if lin.get("MGMaxLevels", 100) == 1 else "multigrid"
    drv = sol.get("Driven", {})
    amr = cfg.get("Model", {}).get("Refinement", {}).get("MaxIts", 0)
    return estimate(unknowns=unknowns, tets=tets, problem=prob, order=order, solver_mode=mode,
                    n_freq=max(1, count_frequencies(cfg)), n_excitations=max(1, len(excitations(cfg))),
                    adaptive=bool(drv.get("AdaptiveTol")), eigen_n=sol.get("Eigenmode", {}).get("N", 1),
                    ranks=ranks, ram_gb=ram_gb, amr_passes=amr)


def calibration_summary() -> dict:
    F = fits()
    tidy = lambda d: {k: {"a": round(math.exp(v["log_a"]), 4), "b": round(v["b"], 3),
                          "extras": {e["name"]: round(e["c"], 3) for e in v.get("extras", [])},
                          "n": v["n"], "N_range": v["N_range"], "spread_x": round(math.exp(v["q85"] - v["q20"]), 2)}
                      for k, v in d.items()}
    return {"memory_gb": tidy(F["memory"]), "time_s": tidy(F["time"]), "records": len(records()),
            "form": "y = a * (N/1e5)^b * prod (extra/ref)^c  (ranks ref 4, modes ref 10)"}
