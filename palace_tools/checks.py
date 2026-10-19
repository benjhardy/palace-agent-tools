"""Sanity checks on Palace results. Each check returns {"name", "status": pass|warn|fail|skip, "value", "detail"}.

The checks catch the failures seen in practice: gain from a passive model (|S| > 1 from an unresolved via or a
floating port), lossy walls that fake a frequency-flat loss, broken reciprocity from inconsistent port
normalisation, a wave-port mode solve that jumped to a spurious mode, an eigenmode far from where it should be,
and a mesh study that is not converging.
"""
from __future__ import annotations

import math

import numpy as np


def _res(name, status, value=None, detail="", **kw):
    d = {"name": name, "status": status, "value": value, "detail": detail}
    d.update(kw)
    return d


def passivity(s: np.ndarray, f_ghz=None, tol: float = 1e-3) -> dict:
    """Largest eigenvalue of S^H S over frequency must be <= 1 (+tol). Needs the full matrix."""
    s = np.asarray(s)
    if np.isnan(s).any():
        cols = ~np.isnan(s).any(axis=(0, 1))
        if not cols.any():
            return _res("passivity", "skip", detail="no complete columns")
        colpow = np.nansum(np.abs(s[:, :, cols]) ** 2, axis=1)
        worst = float(colpow.max())
        st = "pass" if worst <= 1 + tol else "fail"
        return _res("passivity", st, round(worst, 6),
                    "column power sum_i |S_ij|^2 (partial matrix); > 1 means the model creates energy")
    ev = np.array([np.linalg.eigvalsh(m.conj().T @ m).max() for m in s])
    k = int(ev.argmax())
    st = "pass" if ev[k] <= 1 + tol else "fail"
    det = "max eigenvalue of S^H S; > 1 means gain from a passive model (check vias, floating ports, port normalisation)"
    return _res("passivity", st, round(float(ev[k]), 6), det,
                at_ghz=None if f_ghz is None else float(np.asarray(f_ghz)[k]))


def reciprocity(s: np.ndarray, f_ghz=None, tol: float = 1e-3) -> dict:
    s = np.asarray(s)
    if s.shape[1] < 2:
        return _res("reciprocity", "skip", detail="1-port")
    d = np.abs(s - np.transpose(s, (0, 2, 1)))
    if np.isnan(d).all():
        return _res("reciprocity", "skip", detail="needs at least two driven columns")
    worst = float(np.nanmax(d))
    k = int(np.nanargmax(np.nanmax(d, axis=(1, 2))))
    st = "pass" if worst <= tol else ("warn" if worst <= 10 * tol else "fail")
    return _res("reciprocity", st, worst, "max |S_ij - S_ji| (needs equal reference impedances or renormalised S)",
                at_ghz=None if f_ghz is None else float(np.asarray(f_ghz)[k]))


def energy_balance(s: np.ndarray, f_ghz, flat_loss_warn: float = 0.01) -> dict:
    """Dissipated fraction 1 - sum_i |S_ij|^2 per driven port. Negative is gain. A loss that is already large at the
    lowest frequency and nearly flat over the band is the signature of absorbing walls too close to the structure
    (or of a wave port that does not capture the whole mode)."""
    s = np.asarray(s)
    f = np.asarray(f_ghz, float)
    out = {}
    worst_status = "pass"
    for j in range(s.shape[2]):
        col = s[:, :, j]
        if np.isnan(col).any():
            continue
        diss = 1 - np.sum(np.abs(col) ** 2, axis=1)
        lo = float(diss[0])
        slope_ratio = float(diss[0] / diss[-1]) if abs(diss[-1]) > 1e-12 else math.inf
        notes = []
        st = "pass"
        if diss.min() < -1e-3:
            st, notes = "fail", notes + ["negative dissipation (gain)"]
        elif lo > flat_loss_warn and slope_ratio > 0.5:
            st = "warn"
            notes.append(f"{100 * lo:.1f} % loss already at the lowest frequency and roughly flat or falling: check "
                         "absorbing-wall distance and that wave ports span the whole field")
        out[j + 1] = {"min": float(diss.min()), "max": float(diss.max()), "at_fmin": lo, "at_fmax": float(diss[-1]),
                      "notes": notes, "status": st}
        if st == "fail" or (st == "warn" and worst_status == "pass"):
            worst_status = st
    if not out:
        return _res("energy_balance", "skip", detail="no complete columns")
    return _res("energy_balance", worst_status, out, "per driven column (index = column position)")


def symmetric_ports(s: np.ndarray, pairs: list[tuple[int, int]], tol_db: float = 0.5) -> dict:
    """For a structure symmetric between ports a and b, |S_aa| and |S_bb| should agree (G01: walls too close gave
    S11 = -42.6 dB vs S22 = -26.8 dB on a bare line)."""
    worst = 0.0
    for a, b in pairs:
        da = 20 * np.log10(np.abs(s[:, a - 1, a - 1]) + 1e-30)
        db = 20 * np.log10(np.abs(s[:, b - 1, b - 1]) + 1e-30)
        m = (da > -40) | (db > -40)  # ignore deep nulls
        if m.any():
            worst = max(worst, float(np.nanmax(np.abs(da - db)[m])))
    st = "pass" if worst <= tol_db else "warn"
    return _res("port_symmetry", st, round(worst, 3), "max dB difference of reflection at mirror ports (above -40 dB)")


def compare_s(s: np.ndarray, s_ref: np.ndarray, tol: float = 0.02, name: str = "vs_reference") -> dict:
    d = np.abs(np.asarray(s) - np.asarray(s_ref))
    worst = float(np.nanmax(d))
    return _res(name, "pass" if worst <= tol else "fail", worst, f"max |S - S_ref| (tolerance {tol})")


def eigen_expected(modes: list[dict], expected_ghz: float, rel_tol: float = 0.01, q_min: float | None = None) -> dict:
    """Is there a mode within rel_tol of the expected frequency (and is its Q plausible)?"""
    if not modes:
        return _res("eigen_expected", "fail", None, "no eigenmodes in eig.csv")
    fs = np.array([m["f_ghz"] for m in modes])
    k = int(np.argmin(np.abs(fs - expected_ghz)))
    err = float(fs[k] / expected_ghz - 1)
    st = "pass" if abs(err) <= rel_tol else "fail"
    det = f"closest mode m={modes[k]['m']} at {fs[k]:.6g} GHz ({100 * err:+.3f} % vs {expected_ghz:.6g} GHz)"
    if q_min is not None and modes[k]["q"] < q_min:
        st, det = "fail", det + f"; Q {modes[k]['q']:.3g} < {q_min}"
    below = int(np.sum(fs < fs[k] * (1 - 1e-6)))
    return _res("eigen_expected", st, err, det, mode=modes[k], modes_below=below)


def eigen_errors(modes: list[dict], max_bkwd: float = 1e-6) -> dict:
    eb = [m["error_bkwd"] for m in modes if not math.isnan(m["error_bkwd"])]
    if not eb:
        return _res("eigen_errors", "skip")
    worst = max(eb)
    return _res("eigen_errors", "pass" if worst <= max_bkwd else "warn", worst, "largest backward error in eig.csv")


def convergence(values: list, unknowns: list | None = None, tol: float = 0.01, kind: str = "relative") -> dict:
    """Mesh/AMR study: values per refinement level (scalars, or arrays such as S-matrices).

    Reports the change between successive levels; pass when the last change <= tol and changes shrink. With three or
    more scalar levels, also a Richardson-style extrapolated value from the observed rate.
    """
    vals = [np.asarray(v) for v in values]
    if len(vals) < 2:
        return _res("convergence", "skip", detail="need at least two levels")
    deltas = []
    for a, b in zip(vals[:-1], vals[1:]):
        d = float(np.nanmax(np.abs(b - a)))
        if kind == "relative" and b.ndim == 0 and b != 0:
            d = abs(float(b - a) / float(b))
        deltas.append(d)
    shrinking = all(d2 <= d1 * 1.05 for d1, d2 in zip(deltas[:-1], deltas[1:]))
    st = "pass" if deltas[-1] <= tol and shrinking else ("warn" if deltas[-1] <= tol else "fail")
    res = _res("convergence", st, deltas[-1], f"successive changes ({kind}): " + ", ".join(f"{d:.3g}" for d in deltas),
               deltas=deltas, shrinking=shrinking)
    if len(vals) >= 3 and vals[-1].ndim == 0:
        v0, v1, v2 = (float(x) for x in vals[-3:])
        d1, d2 = v1 - v0, v2 - v1
        if d1 != 0 and 0 < d2 / d1 < 1:
            r = d2 / d1
            res["extrapolated"] = v2 + d2 * r / (1 - r)
            res["observed_ratio"] = r
    if unknowns:
        res["unknowns"] = list(unknowns)
    return res


def wave_port_modes(kn_records: list[dict], f_ghz_key: str = "f_ghz", eps_r_max: float | None = None,
                    jump_tol: float = 0.05) -> dict:
    """Spurious wave-port mode detection from the log's kn values (logparse.parse_log()["wave_port_kn"]).

    The effective index n_eff = kn / k0 of a port should vary smoothly with frequency and between AMR passes. A jump,
    or n_eff equal to 1 (air) or sqrt(eps_r_max) (substrate-filled) on a line that should sit in between, means the
    port eigen-solve switched mode (seen with nonconformal AMR on large meshes in Palace 0.16)."""
    if not kn_records:
        return _res("wave_port_modes", "skip", detail="no wave-port kn lines in the log")
    by_port: dict[int, list] = {}
    for r in kn_records:
        if r.get(f_ghz_key):
            k0 = 2 * math.pi * r[f_ghz_key] * 1e9 / 299_792_458.0
            by_port.setdefault(r["port"], []).append((r[f_ghz_key], r["kn_re"] / k0))
    notes, status = [], "pass"
    summary = {}
    for p, rows in by_port.items():
        rows.sort()
        neff = np.array([n for _, n in rows])
        summary[p] = {"n_eff_min": float(neff.min()), "n_eff_max": float(neff.max())}
        # the same frequency repeated (AMR passes) or neighbours: relative jumps
        if len(neff) > 1:
            rel = np.abs(np.diff(neff)) / neff[:-1]
            if rel.max() > jump_tol:
                status = "fail"
                notes.append(f"port {p}: n_eff jumps by {100 * rel.max():.1f} %")
        if eps_r_max:
            for val in (1.0, math.sqrt(eps_r_max)):
                if np.any(np.abs(neff / val - 1) < 0.002):
                    status = "warn" if status == "pass" else status
                    notes.append(f"port {p}: n_eff = {val:.4f} (air- or substrate-filled mode?)")
    return _res("wave_port_modes", status, summary, "; ".join(notes) or "smooth effective index")


def run_all_sparams(sp: dict, symmetric_pairs=None, ref=None, tol_ref: float = 0.02) -> list[dict]:
    s, f = sp["s"], sp["f_ghz"]
    out = [passivity(s, f), reciprocity(s, f), energy_balance(s, f)]
    if symmetric_pairs:
        out.append(symmetric_ports(s, symmetric_pairs))
    if ref is not None:
        out.append(compare_s(s, ref, tol_ref))
    return out


def overall(results: list[dict]) -> str:
    st = [r["status"] for r in results]
    return "fail" if "fail" in st else ("warn" if "warn" in st else "pass")
