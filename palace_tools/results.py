"""Read Palace postprocessing CSVs and turn them into numbers an agent can check.

Files handled (Palace 0.16 names, in ``config["Problem"]["Output"]``):
  port-S.csv  (driven: |S| dB and arg deg per column)   eig.csv (eigenmode: f, Q, errors)
  port-V.csv / port-I.csv (complex peak values)          port-EPR.csv / port-Q.csv (eigenmode EPR, coupling Q)
  domain-E.csv (energies)  error-indicators.csv (AMR)   surface-F.csv / surface-Q.csv / probe-*.csv
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

SIJ = re.compile(r"\|S\[(\d+)\]\[(\d+)\]\| \(dB\)")


def read_csv(path: Path) -> tuple[list[str], np.ndarray]:
    """Palace CSV: one header line, comma separated, padded numbers. Returns (header names, 2-D float array)."""
    text = Path(path).read_text(encoding="utf-8", errors="replace").strip().splitlines()
    hdr = [h.strip() for h in text[0].split(",")]
    rows = []
    for ln in text[1:]:
        if not ln.strip():
            continue
        vals = []
        for x in ln.split(","):
            x = x.strip()
            try:
                vals.append(float(x))
            except ValueError:
                vals.append(float("nan"))
        rows.append(vals)
    width = len(hdr)
    data = np.array([r + [np.nan] * (width - len(r)) if len(r) < width else r[:width] for r in rows], float)
    return hdr, data.reshape(-1, width)


def read_port_s(path: Path) -> dict:
    """port-S.csv -> {"f_ghz": (nf,), "s": (nf, n, n) complex with NaN where a column was not driven,
    "driven": [port indices whose column is present], "ports": [all port indices]}"""
    hdr, data = read_csv(path)
    f = data[:, 0]
    idx = {}
    for k, h in enumerate(hdr):
        m = SIJ.match(h)
        if m:
            idx[(int(m.group(1)), int(m.group(2)))] = k
    ports = sorted({i for i, _ in idx} | {j for _, j in idx})
    pos = {p: n for n, p in enumerate(ports)}
    s = np.full((len(f), len(ports), len(ports)), np.nan + 1j * np.nan)
    for (i, j), k in idx.items():
        ph = next((kk for kk, h in enumerate(hdr) if h.startswith(f"arg(S[{i}][{j}])")), None)
        mag = 10 ** (data[:, k] / 20)
        ang = np.deg2rad(data[:, ph]) if ph is not None else 0.0
        s[:, pos[i], pos[j]] = mag * np.exp(1j * ang)
    driven = sorted({j for _, j in idx})
    return {"f_ghz": f, "s": s, "ports": ports, "driven": driven}


def merge_port_s(parts: list[dict]) -> dict:
    """Combine single-excitation runs (one column each) on the same frequency grid into one S-matrix."""
    f = parts[0]["f_ghz"]
    ports = sorted({p for d in parts for p in d["ports"]})
    pos = {p: n for n, p in enumerate(ports)}
    s = np.full((len(f), len(ports), len(ports)), np.nan + 1j * np.nan)
    for d in parts:
        if len(d["f_ghz"]) != len(f) or not np.allclose(d["f_ghz"], f):
            raise ValueError("frequency grids differ; resample before merging")
        for j in d["driven"]:
            for i in d["ports"]:
                s[:, pos[i], pos[j]] = d["s"][:, d["ports"].index(i), d["ports"].index(j)]
    return {"f_ghz": f, "s": s, "ports": ports, "driven": sorted({j for d in parts for j in d["driven"]})}


def port_impedances(config: dict) -> dict:
    """Reference impedance of each port as Palace normalises S: lumped port R (circuit value; per-square Rs is
    converted by Palace with the port's l/w, which needs the mesh) and None for wave ports (normalised to the
    port's own mode, so there is no fixed ohmic reference)."""
    out = {}
    for p in config.get("Boundaries", {}).get("LumpedPort", []):
        out[int(p["Index"])] = float(p["R"]) if p.get("R") else None
    for p in config.get("Boundaries", {}).get("WavePort", []):
        out[int(p["Index"])] = None
    return out


def to_touchstone(sp: dict, out: Path, z0: float | list | None = 50.0) -> dict:
    """Write a Touchstone file with scikit-rf (falls back to a plain writer for 1- and 2-port data).

    Only complete matrices are written; missing columns raise. ``z0`` may be one number or one per port; Palace's
    S is referenced to each lumped port's own R, so pass those values (``port_impedances``).
    """
    s = sp["s"]
    if np.isnan(s).any():
        missing = sorted({sp["ports"][j] for j in range(s.shape[2]) if np.isnan(s[:, :, j]).any()})
        raise ValueError(f"S-matrix incomplete: columns for ports {missing} were not driven; merge more runs")
    n = s.shape[1]
    out = Path(out).with_suffix(f".s{n}p")
    zz = np.broadcast_to(np.asarray(50.0 if z0 is None else z0, float), (n,))
    try:
        import skrf as rf

        net = rf.Network(frequency=rf.Frequency.from_f(sp["f_ghz"], unit="ghz"), s=s, z0=np.tile(zz, (len(s), 1)))
        net.write_touchstone(str(out.with_suffix("")), form="ma")
        writer = "scikit-rf"
    except ImportError:
        if n > 2 or not np.allclose(zz, zz[0]):
            raise
        rows = [f"# GHz S MA R {zz[0]:g}"]
        order = [(0, 0), (1, 0), (0, 1), (1, 1)] if n == 2 else [(0, 0)]
        for k, fk in enumerate(sp["f_ghz"]):
            rows.append(f"{fk:.9g} " + " ".join(f"{abs(s[k, i, j]):.9e} {np.angle(s[k, i, j], deg=True):.6f}"
                                               for i, j in order))
        out.write_text("\n".join(rows) + "\n", encoding="utf-8")
        writer = "builtin"
    return {"path": str(out), "ports": n, "points": int(len(sp["f_ghz"])), "z0": zz.tolist(), "writer": writer}


def read_eig(path: Path) -> list[dict]:
    hdr, data = read_csv(path)
    col = {h: k for k, h in enumerate(hdr)}

    def get(name_start):
        k = next((v for h, v in col.items() if h.startswith(name_start)), None)
        return data[:, k] if k is not None else np.full(len(data), np.nan)

    m, fr, fi, q = get("m"), get("Re{f}"), get("Im{f}"), get("Q")
    eb, ea = get("Error (Bkwd.)"), get("Error (Abs.)")
    return [{"m": int(m[i]), "f_ghz": float(fr[i]), "f_imag_ghz": float(fi[i]), "q": float(q[i]),
             "error_bkwd": float(eb[i]), "error_abs": float(ea[i])} for i in range(len(data))]


def read_table(path: Path) -> dict:
    """Any Palace CSV as {column name: list}."""
    hdr, data = read_csv(path)
    return {h: data[:, k].tolist() for k, h in enumerate(hdr)}


def read_error_indicators(path: Path) -> list[dict]:
    hdr, data = read_csv(path)
    return [dict(zip([h.lower() for h in hdr], map(float, row))) for row in data]


def amr_history(output_dir: Path) -> list[dict]:
    """Per-AMR-iteration S or eigenfrequencies from ``iterationN`` sub-folders plus the final folder."""
    output_dir = Path(output_dir)
    dirs = sorted((d for d in output_dir.glob("iteration*") if d.is_dir()),
                  key=lambda d: int(re.sub(r"\D", "", d.name) or 0)) + [output_dir]
    hist = []
    for d in dirs:
        item: dict = {"folder": d.name}
        if (d / "port-S.csv").exists():
            item["port_s"] = read_port_s(d / "port-S.csv")
        if (d / "eig.csv").exists():
            item["eig"] = read_eig(d / "eig.csv")
        if (d / "palace.json").exists():
            try:
                item["unknowns"] = json.loads((d / "palace.json").read_text())["Problem"]["DegreesOfFreedom"]
            except (KeyError, ValueError):
                pass
        if len(item) > 1:
            hist.append(item)
    return hist


def fields_summary(output_dir: Path) -> dict:
    """What field output exists (ParaView / GLVis), and how big it is."""
    output_dir = Path(output_dir)
    res = {}
    for sub in ("paraview", "gridfunction"):
        d = output_dir / sub
        if d.exists():
            files = [p for p in d.rglob("*") if p.is_file()]
            res[sub] = {"files": len(files), "size_mb": round(sum(p.stat().st_size for p in files) / 2 ** 20, 1),
                        "collections": sorted({p.parent.relative_to(d).as_posix() for p in files if p.suffix == ".pvd"}
                                              | {p.stem for p in d.glob("*.pvd")}),
                        "steps": len(list(d.rglob("Cycle*")))}
    return res


def summarize_output(output_dir: Path, config: dict | None = None) -> dict:
    """Compact JSON summary of everything in a Palace output folder."""
    output_dir = Path(output_dir)
    res: dict = {"output_dir": str(output_dir), "files": sorted(p.name for p in output_dir.glob("*.csv"))}
    ps = output_dir / "port-S.csv"
    if ps.exists():
        sp = read_port_s(ps)
        s = sp["s"]
        res["port_s"] = {"ports": sp["ports"], "driven": sp["driven"], "n_freq": int(len(sp["f_ghz"])),
                         "f_range_ghz": [float(sp["f_ghz"].min()), float(sp["f_ghz"].max())],
                         "max_abs_s": float(np.nanmax(np.abs(s))), "complete": bool(not np.isnan(s).any())}
        if config:
            res["port_s"]["z_ref"] = port_impedances(config)
    elif config and config.get("Problem", {}).get("Type") == "Driven":
        res["port_s"] = {"missing": True,
                         "why": "No port-S.csv. Palace 0.16 writes none when lumped and wave ports are mixed, when one "
                                "excitation drives several ports, or when a SurfaceCurrent source is present."}
    if (output_dir / "eig.csv").exists():
        res["eig"] = read_eig(output_dir / "eig.csv")
    for name in ("port-EPR.csv", "port-Q.csv", "domain-E.csv", "surface-Q.csv", "surface-F.csv"):
        if (output_dir / name).exists():
            t = read_table(output_dir / name)
            res[name.replace(".csv", "")] = {"columns": list(t), "rows": len(next(iter(t.values()), []))}
    if (output_dir / "error-indicators.csv").exists():
        res["error_indicators"] = read_error_indicators(output_dir / "error-indicators.csv")
    hist = amr_history(output_dir)
    if len(hist) > 1:
        res["amr_iterations"] = len(hist)
    res["fields"] = fields_summary(output_dir)
    return res
