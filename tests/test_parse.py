import json

import numpy as np
import pytest

from palace_tools import logparse, results

SYNTH = """>> /usr/lib64/mpich/bin/mpirun -n 6 /opt/palace/bin/palace-x86_64.bin config.json
Git changeset ID: v0.16.0-34-gea2e7b23
Running with 6 MPI processes, 1 OpenMP thread

--> Warning!
Unknown absorbing boundary attributes!
Solver will just ignore them!

--> Warning!
One or more external boundary attributes has no associated boundary condition!
PMC/ZeroCharge condition is assumed!
 ** ERROR RETURN ** FROM DMUMPS INFO(1)=   -9
 ** INFO(2)=            8723
SuperLUSolver: Found a singular matrix
PCG solver did NOT converge in 30 iterations (avg. reduction factor: 8.767e-01)
It 3/41: ω/2π = 8.100e+00 GHz (total elapsed time = 1.00e+02 s)
"""


def test_parse_synthetic_errors():
    r = logparse.parse_log(SYNTH)
    assert r["ranks"] == 6 and r["palace_version"].startswith("v0.16")
    sev = {m["line"][:20]: m["severity"] for m in r["messages"]}
    assert r["n_errors"] >= 3
    assert any("Unknown absorbing" in m["line"] and m["severity"] == "error" for m in r["messages"])
    assert any("INFO(1)" in m["line"] and "workspace" in m["meaning"] for m in r["messages"])
    assert any(m["severity"] == "info" and "PMC" in m["meaning"] for m in r["messages"])
    assert r["progress"] == {"step": 3, "of": 41, "f_ghz": 8.1}
    assert not r["finished"]
    assert sev


def test_parse_cavity_log(data):
    r = logparse.parse_log((data / "cavity_palace.log").read_text(encoding="utf-8"))
    assert r["problem"] == "Eigenmode" and r["finished"]
    assert r["order"] == 2 and r["unknowns"] > 10_000 and r["ranks"] == 2
    assert r["eigen_shift_ghz"] == pytest.approx(5.0, rel=1e-3)
    assert r["peak_mem_gb"] > 0 and r["elapsed_s"] > 0
    assert r["multigrid_levels"] == 1 and r["n_errors"] == 0


def test_wave_port_kn_frequency_is_converted(data):
    r = logparse.parse_log((data / "ms_wave_palace.log").read_text(encoding="utf-8"))
    kn = r["wave_port_kn"]
    assert kn and kn[0]["f_ghz"] == pytest.approx(0.5, rel=1e-3)  # log prints 2*pi*f
    k0 = 2 * np.pi * 0.5e9 / 299792458.0
    assert (kn[0]["kn_re"] / k0) ** 2 == pytest.approx(2.78, rel=0.01)


def test_time_v(data):
    t = logparse.parse_time_v((data / "cavity_time.log").read_text(encoding="utf-8"))
    assert t["wall_s"] == pytest.approx(7.71) and t["max_rss_one_process_gb"] > 0


def test_read_eig(data):
    modes = results.read_eig(data / "cavity_eig.csv")
    assert modes[0]["f_ghz"] == pytest.approx(8.2439, rel=1e-4)
    assert modes[0]["q"] == pytest.approx(1000, rel=1e-4)


def test_read_port_s_full_and_touchstone(data, tmp_path):
    sp = results.read_port_s(data / "pp_port-S.csv")
    assert sp["ports"] == [1, 2] and sp["driven"] == [1, 2]
    assert not np.isnan(sp["s"]).any()
    cfg = json.loads((data / "pp_config.json").read_text())
    z = results.port_impedances(cfg)
    assert z == {1: 50.0, 2: 50.0}
    pytest.importorskip("skrf")
    info = results.to_touchstone(sp, tmp_path / "pp", [50.0, 50.0])
    import skrf

    net = skrf.Network(info["path"])
    assert np.allclose(net.s, sp["s"], atol=1e-6)


def test_partial_matrix_and_merge(data):
    sp = results.read_port_s(data / "pp_load_port-S.csv")
    assert sp["driven"] == [1] and np.isnan(sp["s"][:, :, 1]).all()
    with pytest.raises(ValueError, match="incomplete"):
        results.to_touchstone(sp, "x", 50)
    full = results.read_port_s(data / "pp_port-S.csv")
    col1 = {"f_ghz": full["f_ghz"], "s": full["s"].copy(), "ports": [1, 2], "driven": [1]}
    col1["s"][:, :, 1] = np.nan
    col2 = {"f_ghz": full["f_ghz"], "s": full["s"].copy(), "ports": [1, 2], "driven": [2]}
    col2["s"][:, :, 0] = np.nan
    m = results.merge_port_s([col1, col2])
    assert np.allclose(m["s"], full["s"])


def test_summarize_output(tmp_path, data):
    (tmp_path / "eig.csv").write_text((data / "cavity_eig.csv").read_text())
    s = results.summarize_output(tmp_path)
    assert len(s["eig"]) >= 4 and s["fields"] == {}
    s = results.summarize_output(tmp_path, {"Problem": {"Type": "Driven"}})
    assert s["port_s"]["missing"]
