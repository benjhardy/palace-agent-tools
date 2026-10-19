# Examples

Small, fast cases with closed-form references. Each script builds the mesh (Gmsh), builds and lints the config,
estimates the cost, runs Palace through `palace_tools.run` (with the busy-host guard), checks the result and
writes `results*.json`. Solver output goes to `run*/` (not committed).

| script | what it verifies | time on 2 ranks |
|---|---|---|
| `cavity_te101/run_example.py` | eigenmode frequencies of a PEC box vs analytic TE/TM modes; Q of a lossy fill vs the exact 1/tan d | ~10 s |
| `lumped_port_line/run_example.py` | lumped ports and a lumped load on an exact TEM line (parallel plate with PMC sides): S vs the closed form to 1e-6 | ~1.5 min |
| `microstrip_thru/run_example.py --ports wave` | microstrip eeff(f) from full-face wave ports vs Kirschning-Jansen | ~2 min |
| `microstrip_thru/run_example.py` (lumped) | Z0 and eeff vs Hammerstad-Jensen / Kirschning-Jansen with edge lumped ports and two-line de-embedding | ~2.5 min |

Run from the repository root with the package installed (`pip install -e ".[all]"`), e.g.

```
python examples/cavity_te101/run_example.py --ranks 2
```

Palace must be reachable (see `palace-tools doctor`). All results are simulation only.
