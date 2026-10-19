"""Rectangular cavity eigenmodes vs the analytic TE/TM resonances, and Q vs the exact dielectric-loss Q.

WR-90 cross-section cavity, a = 22.86 mm (x), b = 10.16 mm (y), d = 30 mm (z), PEC walls, filled with a
dielectric eps_r = 1, tan_d = 1e-3 (so Q = 1 / tan_d to first order, exactly 1 / (2 tan(atan(tan_d) / 2))).

    python examples/cavity_te101/run_example.py [--ranks 2] [--h 2.0] [--no-run]

Writes the mesh, config and Palace output to examples/cavity_te101/run/ and a summary to results.json.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from palace_tools import api, rf
from palace_tools.mesh import rectangular_cavity

HERE = Path(__file__).resolve().parent
A, B, D, TAN_D = 22.86, 10.16, 30.0, 1e-3


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ranks", type=int, default=2)
    ap.add_argument("--h", type=float, default=2.5, help="mesh size, mm")
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--no-run", action="store_true")
    a = ap.parse_args()
    run = HERE / "run"
    run.mkdir(exist_ok=True)

    tags = rectangular_cavity(run / "cavity", A, B, D, a.h, units="mm", verbosity=0)
    spec = {"problem": "eigenmode", "mesh": "cavity.msh", "units": "mm", "tags": "cavity.tags.json",
            "order": a.order, "materials": [{"tags": ["cavity"], "eps_r": 1.0, "tan_d": TAN_D}],
            "pec": ["walls"], "eigen": {"target_ghz": 5.0, "n": 4, "tol": 1e-8, "expected_ghz": 8.24}}
    built = api.call("palace_build_config", {"spec": spec, "out": str(run / "config.json")})
    lint = built["lint"]
    cost = api.call("palace_estimate_cost", {"config": str(run / "config.json"),
                                             "tags": str(run / "cavity.tags.json"), "ranks": a.ranks})
    print("lint:", lint["counts"], [f["id"] for f in lint["findings"]])
    print("estimate:", cost["unknowns"], "unknowns", cost["memory_gb"], "GB", cost["wall_min"], "min")
    if a.no_run:
        return
    rec = api.call("palace_run", {"case_dir": str(run), "ranks": a.ranks})
    if not rec.get("ok"):
        raise SystemExit(json.dumps(rec, indent=1))
    res = api.call("palace_results", {"case_dir": str(run)})
    modes = res["eig"]
    ana = rf.cavity_modes(A * 1e-3, B * 1e-3, D * 1e-3, n=6)
    q_exact = rf.cavity_q_dielectric(TAN_D)
    rows = []
    for m, ref in zip(modes, ana):
        rows.append({"mode": ref["mode"], "f_palace_ghz": m["f_ghz"], "f_analytic_ghz": ref["f_ghz"],
                     "rel_err": m["f_ghz"] / ref["f_ghz"] - 1, "q_palace": m["q"], "q_exact": q_exact})
    chk = api.call("palace_check", {"case_dir": str(run), "expect": {"eigen_ghz": ana[0]["f_ghz"], "rel_tol": 1e-3,
                                                                    "q_min": 0.99 * q_exact}})
    out = {"example": "cavity_te101", "geometry_mm": {"a": A, "b": B, "d": D}, "tan_d": TAN_D,
           "mesh": {"h_mm": a.h, "tets": tags["mesh"]["tets"], "order": a.order,
                    "unknowns": res["log"]["unknowns"]},
           "modes": rows, "check": chk["overall"],
           "cost": {"ranks": a.ranks, "wall_s": rec["wall_s"], "palace_elapsed_s": res["log"]["elapsed_s"],
                    "peak_mem_gb": res["log"]["peak_mem_gb"], "estimate": {"memory_gb": cost["memory_gb"],
                                                                           "wall_min": cost["wall_min"]}},
           "palace_version": res["log"]["palace_version"]}
    (HERE / "results.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    for r in rows:
        print(f"{r['mode']:6s} {r['f_palace_ghz']:.6f} GHz vs {r['f_analytic_ghz']:.6f} ({100 * r['rel_err']:+.4f} %)"
              f"  Q {r['q_palace']:.2f} vs {r['q_exact']:.2f}")
    print("check:", chk["overall"])


if __name__ == "__main__":
    main()
