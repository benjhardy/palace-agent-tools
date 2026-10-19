"""Microstrip thru line, driven: S11/S21, Z0 and eeff vs Hammerstad-Jensen (Z0, eeff0) + Kirschning-Jansen (eeff(f)).

A 50-ohm (Hammerstad-Jensen) strip on a 0.508 mm, eps_r 3.55 substrate (RO4003C-like but lossless, so the
comparison is clean), zero-thickness PEC strip, PEC ground, inside a PEC shield 9 substrate heights away at the
sides and above (close to an open line).

--ports lumped (default): edge lumped ports, vertical 50-ohm sheets from the ground to the strip, with PMC end
  walls (lesson: absorbing end walls next to edge ports leak power; PEC end walls short the strip). Lumped ports
  are not calibrated: each adds a series inductance (~0.1 nH here), which biases a single-line extraction (eeff
  +5 %, Z0 swinging +-20 % above 5 GHz). Two lengths (10 and 20 mm) fix it: eig(T_20 T_10^-1) gives gamma free of
  the ports, then Z0 and the port inductance are fitted to both lines.
--ports wave: full-face wave ports on the end walls; eeff from the S21 phase (ports de-embedded by construction).

    python examples/microstrip_thru/run_example.py [--ranks 2] [--ports lumped|wave]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from palace_tools import api, rf
from palace_tools.mesh import microstrip_line
from palace_tools.results import read_port_s

HERE = Path(__file__).resolve().parent
H, ER, SIDE, TOP = 0.508, 3.55, 4.5, 4.5
FREQ = {"min_ghz": 0.5, "max_ghz": 10.0, "step_ghz": 0.25, "adaptive_tol": 1e-4}


def case(name: str, length: float, ports: str, w: float, a) -> dict:
    run = HERE / name
    run.mkdir(exist_ok=True)
    tags = microstrip_line(run / "line", w, H, length, ER, side=SIDE, top=TOP, h_mesh=a.h_mesh, h_edge=a.h_edge,
                           ports=ports, f_ghz=10.0, units="mm", verbosity=0)
    spec = {"problem": "driven", "mesh": "line.msh", "units": "mm", "tags": "line.tags.json", "order": 2,
            "materials": [{"tags": ["substrate"], "eps_r": ER}, {"tags": ["air"], "eps_r": 1.0}],
            "pec": ["ground", "strip", "sides", "top"], "frequency": FREQ}
    if ports == "lumped":
        spec["lumped_ports"] = [{"index": 1, "tags": ["port1"], "direction": "+Z", "R": 50.0},
                                {"index": 2, "tags": ["port2"], "direction": "+Z", "R": 50.0}]
        spec["pmc"] = ["ends"]
    else:
        spec["wave_ports"] = [{"index": 1, "tags": ["port1"]}, {"index": 2, "tags": ["port2"]}]
    built = api.call("palace_build_config", {"spec": spec, "out": str(run / "config.json")})
    cost = api.call("palace_estimate_cost", {"config": str(run / "config.json"), "tags": str(run / "line.tags.json"),
                                             "ranks": a.ranks})
    print(name, "lint:", built["lint"]["counts"], [f["id"] for f in built["lint"]["findings"]])
    print(name, "estimate:", cost["unknowns"], "unknowns", cost["memory_gb"]["typ"], "GB",
          cost["wall_min"], "min")
    if a.no_run:
        return {}
    rec = api.call("palace_run", {"case_dir": str(run), "ranks": a.ranks})
    if not rec.get("ok"):
        raise SystemExit(json.dumps(rec, indent=1))
    res = api.call("palace_results", {"case_dir": str(run), "touchstone": str(run / "line")})
    chk = api.call("palace_check", {"case_dir": str(run), "symmetric_pairs": [[1, 2]]})
    sp = read_port_s(run / "postpro" / "port-S.csv")
    return {"sp": sp, "check": chk, "length_mm": length,
            "summary": {"length_mm": length, "tets": tags["mesh"]["tets"], "unknowns": res["log"]["unknowns"],
                        "ranks": a.ranks, "wall_s": rec["wall_s"], "palace_elapsed_s": res["log"]["elapsed_s"],
                        "peak_mem_gb": res["log"]["peak_mem_gb"], "check": chk["overall"],
                        "checks": {c["name"]: c["status"] for c in chk["checks"]},
                        "estimate": {"memory_gb": cost["memory_gb"], "wall_min": cost["wall_min"]}},
            "version": res["log"]["palace_version"]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ranks", type=int, default=2)
    ap.add_argument("--ports", default="lumped", choices=["lumped", "wave"])
    ap.add_argument("--h-mesh", type=float, default=1.6)
    ap.add_argument("--h-edge", type=float, default=0.12)
    ap.add_argument("--no-run", action="store_true")
    a = ap.parse_args()
    w = rf.microstrip_width_for_z0(50.0, H, ER)
    hj = rf.microstrip_hj(w, H, ER)
    print(f"width {w:.4f} mm; Hammerstad-Jensen Z0 {hj['z0']:.3f} ohm, eeff0 {hj['eeff']:.4f}")
    pick = [1.0, 2.0, 5.0, 8.0, 10.0]
    out = {"example": "microstrip_thru", "ports": a.ports, "w_mm": w, "h_mm": H, "eps_r": ER,
           "hammerstad_jensen": hj}
    if a.ports == "wave":
        c = case("run_wave", 20.0, "wave", w, a)
        if a.no_run:
            return
        f, s = c["sp"]["f_ghz"], c["sp"]["s"]
        ph = -np.unwrap(np.angle(s[:, 1, 0]))
        eeff = (ph * rf.C0 / (2 * np.pi * f * 1e9 * 20e-3)) ** 2
        idx = [int(np.argmin(np.abs(f - x))) for x in pick]
        out["cases"] = [c["summary"]]
        out["table"] = [{"f_ghz": f[i], "s11_db": 20 * np.log10(abs(s[i, 0, 0])), "s21_db": 20 * np.log10(abs(s[i, 1, 0])),
                         "eeff_palace": eeff[i], "eeff_kj": rf.microstrip_eeff_kj(f[i], w, H, ER, hj["eeff"])}
                        for i in idx]
    else:
        c1 = case("run_lumped_L10", 10.0, "lumped", w, a)
        c2 = case("run_lumped_L20", 20.0, "lumped", w, a)
        if a.no_run:
            return
        f = c1["sp"]["f_ghz"]
        s1, s2 = c1["sp"]["s"], c2["sp"]["s"]
        g = rf.two_line_gamma(s1, s2, 10e-3, f)
        gamma = g["alpha_np_m"] * 0 + 1j * g["beta_rad_m"]  # lossless model; alpha from the fit is numerical noise
        fit = rf.fit_line_with_series_ports([s1, s2], [10e-3, 20e-3], f, gamma)
        single = rf.extract_line_from_s(s2, f, 20e-3)
        idx = [int(np.argmin(np.abs(f - x))) for x in pick]
        out["cases"] = [c1["summary"], c2["summary"]]
        out["table"] = [{"f_ghz": f[i], "s11_db_L20": 20 * np.log10(abs(s2[i, 0, 0])),
                         "s21_db_L20": 20 * np.log10(abs(s2[i, 1, 0])),
                         "z0_two_line_fit": fit["z0"][i], "z0_hj": hj["z0"],
                         "port_series_l_ph": fit["lp_h"][i] * 1e12,
                         "eeff_two_line": g["eeff"][i], "eeff_kj": rf.microstrip_eeff_kj(f[i], w, H, ER, hj["eeff"]),
                         "z0_single_line_naive": single["z0"][i].real, "eeff_single_line_naive": single["eeff"][i]}
                        for i in idx]
        band = (f >= 0.5) & (f <= 10)
        out["band_summary"] = {
            "z0_two_line_median": float(np.median(fit["z0"][band])),
            "z0_two_line_range": [float(fit["z0"][band].min()), float(fit["z0"][band].max())],
            "eeff_two_line_vs_kj_max_rel": float(np.max(np.abs(
                g["eeff"][band] / np.array([rf.microstrip_eeff_kj(x, w, H, ER, hj["eeff"]) for x in f[band]]) - 1))),
            "port_series_l_ph_median": float(np.median(fit["lp_h"][band]) * 1e12)}
    out["palace_version"] = (c if a.ports == "wave" else c1)["version"]
    out = api._clean(out)
    (HERE / f"results_{a.ports}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    for r in out["table"]:
        print(json.dumps({k: round(v, 4) if isinstance(v, float) else v for k, v in r.items()}))
    if "band_summary" in out:
        print(json.dumps(out["band_summary"]))
    print("checks:", [cs["checks"] for cs in out["cases"]])


if __name__ == "__main__":
    main()
