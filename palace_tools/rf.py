"""Closed-form RF references used to check Palace results.

Everything here is textbook material (Pozar, Hammerstad & Jensen 1980, Kirschning & Jansen 1982) and is used by
`check` and by the examples to compare a full-wave answer with an analytic one.
SI units unless a name says otherwise (``*_ghz``, ``*_mm``).
"""
from __future__ import annotations

import math
from itertools import product

import numpy as np

C0 = 299_792_458.0
MU0 = 4e-7 * math.pi  # Palace and most solvers use the classical value; the 2019 SI change is 1e-10 relative
EPS0 = 1.0 / (MU0 * C0 * C0)
ETA0 = MU0 * C0

UNITS = {"m": 1.0, "cm": 1e-2, "mm": 1e-3, "um": 1e-6, "mil": 25.4e-6, "in": 25.4e-3, "nm": 1e-9}


def wavelength(f_ghz: float, eps_r: float = 1.0, mu_r: float = 1.0) -> float:
    """Wavelength in m at f (GHz) in a medium."""
    return C0 / (f_ghz * 1e9 * math.sqrt(eps_r * mu_r))


def skin_depth(f_ghz: float, sigma: float, mu_r: float = 1.0) -> float:
    """Skin depth in m."""
    return 1.0 / math.sqrt(math.pi * f_ghz * 1e9 * MU0 * mu_r * sigma)


def surface_resistance(f_ghz: float, sigma: float, mu_r: float = 1.0) -> float:
    return 1.0 / (sigma * skin_depth(f_ghz, sigma, mu_r))


# ----------------------------------------------------------------------------------------------- rectangular cavity

def cavity_modes(a: float, b: float, d: float, eps_r: float = 1.0, mu_r: float = 1.0, n: int = 10,
                 max_index: int = 4) -> list[dict]:
    """Resonances of a PEC rectangular cavity (a along x, b along y, d along z; propagation along z).

    Returns the lowest ``n`` modes (TE_mnp and TM_mnp, degenerate pairs listed separately), sorted by frequency.
    TE_mnp: m, n not both 0, p >= 1. TM_mnp: m, n >= 1, p >= 0.
    """
    out = []
    v = C0 / math.sqrt(eps_r * mu_r)
    for m, nn, p in product(range(max_index + 1), repeat=3):
        k = math.sqrt((m * math.pi / a) ** 2 + (nn * math.pi / b) ** 2 + (p * math.pi / d) ** 2)
        if k == 0:
            continue
        f = v * k / (2 * math.pi) / 1e9
        if (m or nn) and p >= 1:
            out.append({"mode": f"TE{m}{nn}{p}", "f_ghz": f})
        if m >= 1 and nn >= 1:
            out.append({"mode": f"TM{m}{nn}{p}", "f_ghz": f})
    out.sort(key=lambda r: r["f_ghz"])
    return out[:n]


def cavity_q_dielectric(tan_d: float) -> float:
    """Exact Q of a cavity completely filled with a dielectric of loss tangent tan_d (PEC walls).

    eps = eps' (1 - j tan_d) gives omega = omega0 (1 - j tan_d)^(-1/2), so Q = Re(w) / (2 |Im(w)|) = 1 / (2 tan(theta/2))
    with theta = atan(tan_d). For small tan_d this is 1 / tan_d.
    """
    th = math.atan(tan_d)
    return 1.0 / (2.0 * math.tan(th / 2.0))


def cavity_q_conductor_te10l(a: float, b: float, d: float, l: int, sigma: float, eps_r: float = 1.0) -> float:
    """Conductor Q of the TE_10l mode of a rectangular cavity with walls of conductivity sigma (Pozar eq. 6.46)."""
    f = cavity_modes_te10l_f(a, d, l, eps_r)
    k = 2 * math.pi * f * 1e9 * math.sqrt(eps_r) / C0
    eta = ETA0 / math.sqrt(eps_r)
    rs = surface_resistance(f, sigma)
    num = (k * a * d) ** 3 * b * eta
    den = 2 * math.pi ** 2 * rs * (2 * l ** 2 * a ** 3 * b + 2 * b * d ** 3 + l ** 2 * a ** 3 * d + a * d ** 3)
    return num / den


def cavity_modes_te10l_f(a: float, d: float, l: int, eps_r: float = 1.0) -> float:
    return C0 / math.sqrt(eps_r) / 2 * math.sqrt((1 / a) ** 2 + (l / d) ** 2) / 1e9


# ------------------------------------------------------------------------------------------------------- microstrip

def _hj_z01(u: float) -> float:
    fu = 6 + (2 * math.pi - 6) * math.exp(-((30.666 / u) ** 0.7528))
    return ETA0 / (2 * math.pi) * math.log(fu / u + math.sqrt(1 + (2 / u) ** 2))


def _hj_eeff(u: float, er: float) -> float:
    a = 1 + math.log((u ** 4 + (u / 52) ** 2) / (u ** 4 + 0.432)) / 49 + math.log(1 + (u / 18.1) ** 3) / 18.7
    b = 0.564 * ((er - 0.9) / (er + 3)) ** 0.053
    return (er + 1) / 2 + (er - 1) / 2 * (1 + 10 / u) ** (-a * b)


def microstrip_hj(w: float, h: float, eps_r: float, t: float = 0.0) -> dict:
    """Quasi-static microstrip Z0 and effective permittivity, Hammerstad & Jensen (1980), with their thickness
    correction when t > 0. Lengths in any consistent unit. Accuracy ~0.2-1 % for 0.01 <= w/h <= 100, eps_r <= 128."""
    u = w / h
    if t > 0:
        tn = t / h
        du1 = tn / math.pi * math.log(1 + 4 * math.e / (tn * (1 / math.tanh(math.sqrt(6.517 * u))) ** 2))
        dur = 0.5 * (1 + 1 / math.cosh(math.sqrt(eps_r - 1))) * du1
        u1, ur = u + du1, u + dur
        eeff = _hj_eeff(ur, eps_r) * (_hj_z01(u1) / _hj_z01(ur)) ** 2
        z0 = _hj_z01(ur) / math.sqrt(_hj_eeff(ur, eps_r))
    else:
        eeff = _hj_eeff(u, eps_r)
        z0 = _hj_z01(u) / math.sqrt(eeff)
    return {"z0": z0, "eeff": eeff, "w_over_h": u}


def microstrip_eeff_kj(f_ghz: float, w: float, h_mm: float, eps_r: float, eeff0: float) -> float:
    """Kirschning & Jansen (1982) dispersion of eeff. w and h in mm (fn = f[GHz] * h[mm])."""
    u = w / h_mm
    fn = f_ghz * h_mm
    p1 = 0.27488 + (0.6315 + 0.525 / (1 + 0.0157 * fn) ** 20) * u - 0.065683 * math.exp(-8.7513 * u)
    p2 = 0.33622 * (1 - math.exp(-0.03442 * eps_r))
    p3 = 0.0363 * math.exp(-4.6 * u) * (1 - math.exp(-((fn / 38.7) ** 4.97)))
    p4 = 1 + 2.751 * (1 - math.exp(-((eps_r / 15.916) ** 8)))
    p = p1 * p2 * ((0.1844 + p3 * p4) * fn) ** 1.5763
    return eps_r - (eps_r - eeff0) / (1 + p)


def microstrip_width_for_z0(z0: float, h: float, eps_r: float, t: float = 0.0) -> float:
    """Width (same unit as h) that gives Z0 by Hammerstad-Jensen (bisection on w/h)."""
    lo, hi = 0.01 * h, 50 * h
    for _ in range(100):
        mid = math.sqrt(lo * hi)
        if microstrip_hj(mid, h, eps_r, t)["z0"] > z0:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


# -------------------------------------------------------------------------------------- parallel-plate TEM line

def parallel_plate_z0(h: float, w: float, eps_r: float = 1.0, mu_r: float = 1.0) -> float:
    """Exact TEM impedance of a parallel-plate line of height h and width w with PMC side walls."""
    return ETA0 * math.sqrt(mu_r / eps_r) * h / w


# -------------------------------------------------------------------------------------------- line S-parameters

def line_sparams(f_ghz, z0: float, eeff: float, length: float, z_ref: float = 50.0, alpha_np_m: float = 0.0):
    """2-port S-matrix (nf, 2, 2) of a uniform TEM/quasi-TEM line of impedance z0 between z_ref ports."""
    f = np.atleast_1d(np.asarray(f_ghz, float)) * 1e9
    gamma = alpha_np_m + 1j * 2 * np.pi * f * np.sqrt(eeff) / C0
    gl = gamma * length
    a = np.cosh(gl)
    b = z0 * np.sinh(gl)
    c = np.sinh(gl) / z0
    d = a
    den = a + b / z_ref + c * z_ref + d
    s = np.empty((len(f), 2, 2), complex)
    s[:, 0, 0] = (a + b / z_ref - c * z_ref - d) / den
    s[:, 0, 1] = 2 * (a * d - b * c) / den
    s[:, 1, 0] = 2 / den
    s[:, 1, 1] = (-a + b / z_ref - c * z_ref + d) / den
    return s


def extract_line_from_s(s: np.ndarray, f_ghz, length: float, z_ref: float = 50.0) -> dict:
    """Characteristic impedance and effective permittivity of a symmetric 2-port line from its S-matrix.

    S -> ABCD, Z0 = sqrt(B / C), cosh(gamma l) = A. The phase is unwrapped over frequency, so pass a sweep that
    starts below the first half-wave resonance of the line. Port discontinuities are not removed.
    """
    s = np.asarray(s, complex)
    f = np.atleast_1d(np.asarray(f_ghz, float)) * 1e9
    s11, s12, s21, s22 = s[:, 0, 0], s[:, 0, 1], s[:, 1, 0], s[:, 1, 1]
    a = ((1 + s11) * (1 - s22) + s12 * s21) / (2 * s21)
    b = z_ref * ((1 + s11) * (1 + s22) - s12 * s21) / (2 * s21)
    c = ((1 - s11) * (1 - s22) - s12 * s21) / (2 * s21 * z_ref)
    z0 = np.sqrt(b / c)
    z0 = np.where(z0.real < 0, -z0, z0)
    # e^(gamma l) = cosh + sinh = A + C Z0 (with the Re Z0 > 0 branch); unwrap its phase over frequency
    egl = a + c * z0
    beta_l = np.unwrap(np.angle(egl))
    if beta_l[0] < 0:
        beta_l = beta_l + 2 * np.pi * np.ceil(-beta_l[0] / (2 * np.pi))
    eeff = (beta_l * C0 / (2 * np.pi * f * length)) ** 2
    alpha = np.log(np.abs(egl)) / length
    return {"z0": z0, "eeff": eeff, "alpha_np_m": alpha}


# ------------------------------------------------------------------------------- two-line (multiline) de-embedding

def s_to_t(s: np.ndarray) -> np.ndarray:
    """2-port S (nf, 2, 2) -> wave-cascading T with b1,a1 = T (a2, b2) ordering (T = [[-det/S21, S11/S21],
    [-S22/S21, 1/S21]])."""
    s11, s12, s21, s22 = s[:, 0, 0], s[:, 0, 1], s[:, 1, 0], s[:, 1, 1]
    det = s11 * s22 - s12 * s21
    t = np.empty_like(s)
    t[:, 0, 0] = -det / s21
    t[:, 0, 1] = s11 / s21
    t[:, 1, 0] = -s22 / s21
    t[:, 1, 1] = 1 / s21
    return t


def two_line_gamma(s_short: np.ndarray, s_long: np.ndarray, dl: float, f_ghz) -> dict:
    """Propagation constant of a line from two lengths with identical ports (multiline/TRL idea):
    eig(T_long T_short^-1) = exp(+-gamma dl), independent of the port discontinuities. dl in m."""
    f = np.atleast_1d(np.asarray(f_ghz, float))
    t1, t2 = s_to_t(np.asarray(s_short, complex)), s_to_t(np.asarray(s_long, complex))
    rough = np.unwrap(np.angle(s_short[:, 1, 0])) - np.unwrap(np.angle(s_long[:, 1, 0]))
    beta_dl = np.empty(len(f))
    alpha = np.empty(len(f))
    for k in range(len(f)):
        lam = np.linalg.eigvals(t2[k] @ np.linalg.inv(t1[k]))
        th = abs(np.angle(lam[0]))
        n = np.round((rough[k] - th) / (2 * np.pi))
        cands = [th + 2 * np.pi * n, -th + 2 * np.pi * (n + 1), -th + 2 * np.pi * n, th + 2 * np.pi * (n - 1)]
        beta_dl[k] = min(cands, key=lambda c: abs(c - rough[k]))
        alpha[k] = abs(np.log(abs(lam[0]))) / dl
    beta = beta_dl / dl
    eeff = (beta * C0 / (2 * np.pi * f * 1e9)) ** 2
    return {"beta_rad_m": beta, "alpha_np_m": alpha, "eeff": eeff}


def fit_line_with_series_ports(s_list, lengths, f_ghz, gamma, z_ref: float = 50.0) -> dict:
    """Per frequency, fit Z0 (real) and a series port inductance Lp at each end so that
    [jwLp] [line Z0, gamma, L] [jwLp] reproduces the S-matrices of all given lengths (least squares)."""
    from scipy.optimize import least_squares

    f = np.atleast_1d(np.asarray(f_ghz, float)) * 1e9
    z0s, lps, res = [], [], []

    def abcd_s(abcd):
        a, b, c, d = abcd[0, 0], abcd[0, 1], abcd[1, 0], abcd[1, 1]
        den = a + b / z_ref + c * z_ref + d
        return np.array([[(a + b / z_ref - c * z_ref - d) / den, 2 * (a * d - b * c) / den],
                         [2 / den, (-a + b / z_ref - c * z_ref + d) / den]])

    x0 = [50.0, 50.0]
    for k, fk in enumerate(f):
        w = 2 * np.pi * fk

        def model(x, L):
            z0, lp = x[0], x[1] * 1e-12
            gl = gamma[k] * L
            line = np.array([[np.cosh(gl), z0 * np.sinh(gl)], [np.sinh(gl) / z0, np.cosh(gl)]])
            ser = np.array([[1, 1j * w * lp], [0, 1]])
            return abcd_s(ser @ line @ ser)

        def resid(x):
            r = []
            for s, L in zip(s_list, lengths):
                d = (model(x, L) - s[k]).ravel()
                r += list(d.real) + list(d.imag)
            return r

        sol = least_squares(resid, x0)
        x0 = sol.x
        z0s.append(sol.x[0])
        lps.append(sol.x[1] * 1e-12)
        res.append(float(np.sqrt(np.mean(np.square(sol.fun)))))
    return {"z0": np.array(z0s), "lp_h": np.array(lps), "rms_residual": np.array(res)}
