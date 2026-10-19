import numpy as np
import pytest

from palace_tools import checks, estimate, rf, results


def test_passivity_reciprocity_on_exact_line(data):
    sp = results.read_port_s(data / "pp_port-S.csv")
    r = checks.run_all_sparams(sp, [(1, 2)])
    assert checks.overall(r) == "pass"
    ref = rf.line_sparams(sp["f_ghz"], 75.0, 2.2, 0.03)
    assert checks.compare_s(sp["s"], ref, 1e-4)["status"] == "pass"


def test_energy_balance_flags_016_unequal_r(data):
    """Palace 0.16 with a 50-ohm port and a 75-ohm load: |S21|^2 is low by 50/75, the check sees missing power."""
    sp = results.read_port_s(data / "pp_load_port-S.csv")
    e = checks.energy_balance(sp["s"], sp["f_ghz"])
    assert e["status"] == "warn"
    assert e["value"][1]["at_fmin"] == pytest.approx(1 - 0.04 - 0.96 * 50 / 75, abs=1e-6)


def test_passivity_fails_on_gain():
    s = np.zeros((3, 2, 2), complex)
    s[:, 1, 0] = s[:, 0, 1] = 1.05
    assert checks.passivity(s)["status"] == "fail"
    assert checks.energy_balance(s, [1, 2, 3])["status"] == "fail"


def test_flat_wall_loss_detected():
    f = np.linspace(1, 10, 10)
    s = np.zeros((10, 2, 2), complex)
    s[:, 1, 0] = s[:, 0, 1] = np.sqrt(0.95)  # 5 % flat loss on a "lossless" line
    assert checks.energy_balance(s, f)["status"] == "warn"


def test_eigen_expected():
    modes = [{"m": 1, "f_ghz": 8.2439, "q": 1000.0, "error_bkwd": 1e-9, "error_abs": 1e-6}]
    assert checks.eigen_expected(modes, 8.2439, 1e-3, q_min=990)["status"] == "pass"
    assert checks.eigen_expected(modes, 9.0, 1e-2)["status"] == "fail"
    assert checks.eigen_errors(modes)["status"] == "pass"


def test_convergence_and_extrapolation():
    c = checks.convergence([112.95, 111.54, 111.08], tol=0.005)
    assert c["status"] == "pass" and c["shrinking"]
    assert 110.0 < c["extrapolated"] < 111.08
    assert checks.convergence([1.0, 1.5, 1.0], tol=0.01)["status"] == "fail"


def test_wave_port_mode_jump():
    k0 = lambda f: 2 * np.pi * f * 1e9 / 299792458.0  # noqa: E731
    recs = [{"f_ghz": 3.0, "port": 1, "kn_re": 1.73 * k0(3.0)}, {"f_ghz": 3.0, "port": 1, "kn_re": 1.733 * k0(3.0)},
            {"f_ghz": 3.0, "port": 1, "kn_re": 2.0976 * k0(3.0)}]
    r = checks.wave_port_modes(recs, eps_r_max=4.4)
    assert r["status"] == "fail" and "jumps" in r["detail"]


# ------------------------------------------------------------------------------------------- estimator

def test_nedelec_count_matches_palace_log():
    # birdcage production mesh: 212,703 edges, 349,706 faces, 169,135 tets -> Palace 'ND (p = 2): 1124818'
    assert estimate.nedelec_unknowns(212703, 349706, 169135, 2) == 1124818
    assert estimate.nedelec_unknowns(212703, 349706, 169135, 1) == 212703


def test_estimate_anchor_eigen_1p12M():
    r = estimate.estimate(1124818, "Eigenmode", 2, "direct", eigen_n=20, ranks=4, ram_gb=25)
    assert 14 <= r["memory_gb"]["typ"] <= 18
    assert r["wall_min"]["low"] <= 14 and r["wall_min"]["high"] >= 19
    assert r["fits_in_ram"]


def test_estimate_driven_multigrid_and_adaptive():
    r = estimate.estimate(1230882, "Driven", 2, "multigrid", n_freq=3, n_excitations=2, ranks=8)
    assert r["wall_min"]["low"] <= 12.9 <= r["wall_min"]["high"]  # recorded Wilkinson fine mesh: 12.9 min
    assert 6 <= r["memory_gb"]["typ"] <= 11
    a = estimate.estimate(258302, "Driven", 2, "direct", n_freq=317, n_excitations=7, adaptive=True, ranks=8)
    assert a["wall_min"]["low"] <= 128 <= a["wall_min"]["high"]  # recorded 7-port sweep: 128 min
    assert a["solves"]["online_evaluations"] == 317 * 7


def test_estimate_monotonic_and_tets_fallback():
    m = [estimate.estimate(n, "Driven", 2, "direct")["memory_gb"]["typ"] for n in (1e4, 1e5, 1e6)]
    assert m[0] < m[1] < m[2]
    r = estimate.estimate(None, tets=10000, order=2)
    assert r["unknowns"] == 67000 and r["notes"]
    small = estimate.estimate(20000, "Driven", 2, "direct", ranks=8)
    assert any("tiny" in a for a in small["advice"])


def test_estimate_from_config():
    from palace_tools.config import build_config

    c = build_config({"problem": "driven", "mesh": "m.msh", "units": "mm", "materials": [{"tags": [1]}],
                      "lumped_ports": [{"index": 1, "tags": [101], "R": 50}, {"index": 2, "tags": [102], "R": 50}],
                      "frequency": {"min_ghz": 1, "max_ghz": 10, "step_ghz": 0.25, "adaptive_tol": 1e-3}})
    r = estimate.estimate_from_config(c, unknowns=100_000, ranks=2)
    assert r["solver_mode"] == "direct" and r["solves"]["online_evaluations"] == 37 * 2


def test_calibration_summary():
    s = estimate.calibration_summary()
    assert s["records"] > 100 and "Eigenmode/direct/p2" in s["memory_gb"]
