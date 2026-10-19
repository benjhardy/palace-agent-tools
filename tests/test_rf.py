import math

import numpy as np
import pytest

from palace_tools import rf


def test_cavity_te101_wr90():
    modes = rf.cavity_modes(22.86e-3, 10.16e-3, 30e-3, n=3)
    assert modes[0]["mode"] == "TE101"
    f = rf.C0 / 2 * math.sqrt((1 / 22.86e-3) ** 2 + (1 / 30e-3) ** 2) / 1e9
    assert modes[0]["f_ghz"] == pytest.approx(f, rel=1e-12)
    assert modes[1]["mode"] == "TE102"


def test_cavity_q_dielectric():
    assert rf.cavity_q_dielectric(1e-3) == pytest.approx(1000.0, rel=1e-6)
    assert rf.cavity_q_dielectric(0.1) == pytest.approx(10.025, rel=1e-3)


def test_hj_against_scikit_rf():
    skrf = pytest.importorskip("skrf")
    from skrf.media import MLine

    for w, h, er in [(1.1366, 0.508, 3.55), (0.3, 0.635, 9.8), (3.0, 1.6, 4.4)]:
        ours = rf.microstrip_hj(w, h, er)
        freq = skrf.Frequency(1, 1, 1, unit="MHz")
        ml = MLine(frequency=freq, w=w * 1e-3, h=h * 1e-3, t=None, ep_r=er, disp="none", diel="frequencyinvariant",
                   z0_port=50, rough=0, tand=0, rho=0)
        assert ours["z0"] == pytest.approx(float(np.real(ml.z0[0])), rel=2e-3)
        assert ours["eeff"] == pytest.approx(float(np.real(ml.ep_reff_f[0])), rel=2e-3)


def test_width_for_z0_inverts_hj():
    w = rf.microstrip_width_for_z0(50.0, 0.508, 3.55)
    assert rf.microstrip_hj(w, 0.508, 3.55)["z0"] == pytest.approx(50.0, abs=1e-6)
    assert 1.0 < w < 1.25


def test_kj_dispersion_monotonic_and_bounded():
    hj = rf.microstrip_hj(1.1366, 0.508, 3.55)
    e = [rf.microstrip_eeff_kj(f, 1.1366, 0.508, 3.55, hj["eeff"]) for f in (0.1, 5, 20, 100)]
    assert e[0] == pytest.approx(hj["eeff"], rel=1e-3)
    assert all(a < b for a, b in zip(e[:-1], e[1:]))
    assert e[-1] < 3.55


def test_parallel_plate_and_line_roundtrip():
    z0 = rf.parallel_plate_z0(1.0, 3.386553673622131, 2.2)
    assert z0 == pytest.approx(75.0, rel=1e-9)
    f = np.linspace(0.5, 10, 20)
    s = rf.line_sparams(f, 75.0, 2.2, 0.03)
    for m in s:  # lossless: unitary
        assert np.allclose(m.conj().T @ m, np.eye(2), atol=1e-12)
    ex = rf.extract_line_from_s(s, f, 0.03)
    assert np.allclose(ex["z0"].real, 75.0, rtol=1e-9)
    assert np.allclose(ex["eeff"], 2.2, rtol=1e-9)


def _series_port_line(f, z0, eeff, length, lp):
    w = 2 * np.pi * f * 1e9
    out = []
    for wk in w:
        gl = 1j * wk * np.sqrt(eeff) / rf.C0 * length
        line = np.array([[np.cosh(gl), z0 * np.sinh(gl)], [np.sinh(gl) / z0, np.cosh(gl)]])
        ser = np.array([[1, 1j * wk * lp], [0, 1]])
        a, b, c, d = (ser @ line @ ser).ravel()
        den = a + b / 50 + c * 50 + d
        out.append([[(a + b / 50 - c * 50 - d) / den, 2 * (a * d - b * c) / den],
                    [2 / den, (-a + b / 50 - c * 50 + d) / den]])
    return np.array(out)


def test_two_line_deembedding_recovers_line():
    pytest.importorskip("scipy")
    f = np.linspace(0.5, 10, 39)
    s1 = _series_port_line(f, 49.0, 2.8, 10e-3, 140e-12)
    s2 = _series_port_line(f, 49.0, 2.8, 20e-3, 140e-12)
    g = rf.two_line_gamma(s1, s2, 10e-3, f)
    assert np.allclose(g["eeff"], 2.8, rtol=1e-8)
    fit = rf.fit_line_with_series_ports([s1, s2], [10e-3, 20e-3], f, 1j * g["beta_rad_m"])
    assert np.allclose(fit["z0"], 49.0, rtol=1e-6)
    assert np.allclose(fit["lp_h"], 140e-12, rtol=1e-4)
    naive = rf.extract_line_from_s(s2, f, 20e-3)
    assert abs(naive["eeff"][0] / 2.8 - 1) > 0.01  # a single-line extraction is biased by the port inductance


def test_skin_depth_copper():
    assert rf.skin_depth(1.0, 5.8e7) == pytest.approx(2.09e-6, rel=0.01)
