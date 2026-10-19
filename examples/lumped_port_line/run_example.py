"""Lumped ports and a lumped load on an exact TEM line: Palace vs the closed-form transmission-line S-matrix.

A parallel-plate line (PEC plates at z = 0 and z = h, PMC side walls) carries an exact TEM mode with
Z0 = eta0 / sqrt(eps_r) * h / w, here 75 ohm in eps_r = 2.2 (h = 1 mm, w = 3.386 mm, 30 mm long). A lumped port
that covers the whole end face has the same uniform field as the TEM mode, so it is an exact port: no port
parasitics, unlike an edge port on a microstrip.

Two Palace runs on the same mesh:
  a) both ends are 50-ohm lumped ports (each its own excitation) -> full 2x2 S, compared with the analytic
     75-ohm line between 50-ohm ports (|dS| should be ~1e-3);
  b) port 2 is a passive lumped element R = 75 ohm (not excited) -> the line is terminated in Z0, so S11 at the
     50-ohm port 1 is the constant (75 - 50) / (75 + 50) = 0.2 in magnitude (-13.98 dB), with no ripple.
     It also shows a Palace 0.16 bug (fixed in 0.17): with ports of different R, the reported |S21|^2 is
     smaller than the delivered power by R1/R2, and palace-tools check flags the "missing" energy.

    python examples/lumped_port_line/run_example.py [--ranks 2]
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from palace_tools import api, rf
from palace_tools.mesh import parallel_plate_line
from palace_tools.results import read_port_s

HERE = Path(__file__).resolve().parent
ER, H, Z0, L = 2.2, 1.0, 75.0, 30.0
W = rf.ETA0 / math.sqrt(ER) * H / Z0


def run_case(name: str, spec: dict, ranks: int) -> tuple[dict, dict, dict]:
    run = HERE / "run"
    cfgname = f"config_{name}.json"
    spec = dict(spec, output=f"postpro_{name}")
    built = api.call("palace_build_config", {"spec": spec, "out": str(run / cfgname)})
    print(name, "lint:", built["lint"]["counts"], [f["id"] for f in built["lint"]["findings"]])
    rec = api.call("palace_run", {"case_dir": str(run), "config": cfgname, "ranks": ranks})
    if not rec.get("ok"):
        raise SystemExit(json.dumps(rec, indent=1))
    (run / f"palace_{name}.log").write_bytes((run / "palace.log").read_bytes())
    chk = api.call("palace_check", {"case_dir": str(run), "output_dir": str(run / f"postpro_{name}"),
                                    "config": str(run / cfgname)})
    return rec, chk, read_port_s(run / f"postpro_{name}" / "port-S.csv")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ranks", type=int, default=2)
    ap.add_argument("--h", type=float, default=0.5, help="mesh size, mm")
    a = ap.parse_args()
    run = HERE / "run"
    run.mkdir(exist_ok=True)
    tags = parallel_plate_line(run / "line", W, H, L, ER, a.h, units="mm", verbosity=0)
    base = {"problem": "driven", "mesh": "line.msh", "units": "mm", "tags": "line.tags.json", "order": 2,
            "materials": [{"tags": ["dielectric"], "eps_r": ER}], "pec": ["plates"], "pmc": ["sides"],
            "frequency": {"min_ghz": 0.5, "max_ghz": 10.0, "step_ghz": 0.5}}
    spec_a = dict(base, lumped_ports=[{"index": 1, "tags": ["port1"], "direction": "+Z", "R": 50.0},
                                      {"index": 2, "tags": ["port2"], "direction": "+Z", "R": 50.0}])
    spec_b = dict(base, lumped_ports=[{"index": 1, "tags": ["port1"], "direction": "+Z", "R": 50.0}],
                  lumped_elements=[{"index": 2, "tags": ["port2"], "direction": "+Z", "R": Z0}])
    cost = api.call("palace_estimate_cost", {"unknowns": tags["mesh"]["unknowns"]["2"], "problem": "Driven",
                                             "order": 2, "n_freq": 20, "n_excitations": 2, "ranks": a.ranks})
    print("Z0", Z0, "ohm; w", round(W, 4), "mm; unknowns", tags["mesh"]["unknowns"]["2"],
          "estimate", cost["memory_gb"]["typ"], "GB", cost["wall_min"], "min")
    rec_a, chk_a, sp_a = run_case("a_two_ports", spec_a, a.ranks)
    rec_b, chk_b, sp_b = run_case("b_matched_load", spec_b, a.ranks)
    f = sp_a["f_ghz"]
    s_ref = rf.line_sparams(f, Z0, ER, L * 1e-3, 50.0)
    d_a = np.abs(sp_a["s"] - s_ref)
    s11_b = sp_b["s"][:, 0, 0]
    out = {"example": "lumped_port_line", "z0_ohm": Z0, "eps_r": ER, "h_mm": H, "w_mm": W, "length_mm": L,
           "mesh": {"h_mm": a.h, "tets": tags["mesh"]["tets"], "unknowns_p2": tags["mesh"]["unknowns"]["2"]},
           "a_two_ports": {"max_abs_dS_vs_analytic": float(d_a.max()), "at_ghz": float(f[d_a.max(axis=(1, 2)).argmax()]),
                           "s11_db_at_5ghz": float(20 * np.log10(abs(sp_a["s"][np.argmin(abs(f - 5)), 0, 0]))),
                           "s11_db_analytic_at_5ghz": float(20 * np.log10(abs(s_ref[np.argmin(abs(f - 5)), 0, 0]))),
                           "check": chk_a["overall"], "checks": {c["name"]: c["status"] for c in chk_a["checks"]},
                           "wall_s": rec_a["wall_s"], "peak_mem_gb": rec_a["summary"]["peak_mem_gb"]},
           "b_matched_load": {"s11_db_min": float(20 * np.log10(np.abs(s11_b).min())),
                              "s11_db_max": float(20 * np.log10(np.abs(s11_b).max())),
                              "s11_db_expected": float(20 * np.log10((Z0 - 50) / (Z0 + 50))),
                              "s21_db_at_5ghz": float(20 * np.log10(abs(sp_b["s"][np.argmin(abs(f - 5)), 1, 0]))),
                              # power into the matched 75-ohm load is 1 - |S11|^2; Palace 0.16 reports |S21|^2 too
                              # small by R1/R2 when lumped ports have different R (fixed in 0.17, PR 743)
                              "s21_db_expected_power_wave": float(10 * np.log10(1 - (abs(s11_b) ** 2).mean())),
                              "palace_s21_power_ratio": float((np.abs(sp_b["s"][:, 1, 0]) ** 2 /
                                                               (1 - np.abs(s11_b) ** 2)).mean()),
                              "r1_over_r2": 50.0 / Z0,
                              "check": chk_b["overall"], "wall_s": rec_b["wall_s"]},
           "ranks": a.ranks, "estimate": {"memory_gb": cost["memory_gb"], "wall_min": cost["wall_min"]},
           "palace_version": rec_a["summary"]["palace_version"]}
    (HERE / "results.json").write_text(json.dumps(api._clean(out), indent=1), encoding="utf-8")
    print(json.dumps(api._clean(out), indent=1))


if __name__ == "__main__":
    main()
