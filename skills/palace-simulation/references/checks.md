# What `palace-tools check` tests, and what to do when it fails

| check | fails when | usual causes | fix |
|---|---|---|---|
| `passivity` | max eigenvalue of S^H S > 1 (or a column power sum > 1) | a single tiny via under a big pad (unresolved current path), a floating port, inconsistent port normalisation | solid via posts, check port grounding, equal port R |
| `reciprocity` | max \|S_ij - S_ji\| > 1e-3 | unequal port references, a non-converged iterative solve | equal R, MUMPS/direct, check `log` |
| `energy_balance` | dissipation < 0 (gain), or > 1 % already at the lowest frequency and flat | absorbing walls in the near field; a wave port smaller than the mode; Palace 0.16 unequal-R lumped ports (\|S21\|^2 low by R1/R2) | move walls out / PEC shield; full-face wave ports; equal R |
| `port_symmetry` | \|S_aa\| and \|S_bb\| of mirror ports differ > 0.5 dB | asymmetric mesh or walls; walls too close (S11 -42.6 vs S22 -26.8 dB seen) | symmetric mesh, walls further out |
| `eigen_expected` | no mode within rel_tol of the expected frequency, or Q too low | Target above the mode, wrong units (L0), wrong lumped C (symmetry cut), mesh too coarse | lower Target, check L0, C/2 rule, refine |
| `eigen_errors` | backward error > 1e-6 | loose Tol, nonlinear problem | tighter Tol, more MaxIts |
| `log` | Palace error lines, or warnings | "Unknown ... boundary attributes ... ignore them" (a typo silently drops a BC), "Linear solver did not converge", MUMPS INFO(1) < 0 | fix tags; stronger solver; fewer ranks / more RAM |
| `wave_port_modes` | n_eff = kn/k0 jumps between frequencies/passes, or equals 1 or sqrt(eps_r max) | the port eigen-solve picked a spurious (air- or substrate-filled) mode, seen with nonconformal AMR in 0.16 | `"Nonconformal": false`, check the port face |
| `amr_convergence_max_dS` | max \|dS\| between the last two AMR passes > tolerance, or not shrinking | AMR marking does not target your figure of merit; resonant structures converge slowly | seed the mesh at resonators, order 3, more passes |
| `convergence` (user values) | last relative change > tol or changes not shrinking | not converged | refine; report the extrapolated value with its rate |

A pass is necessary, not sufficient: always compare with an independent reference (closed form, published or
measured data, another solver) and do at least one refinement step on the quantity you report.
