# AGENTS.md: using palace-agent-tools from any coding agent

You are working with AWS Palace (FEM electromagnetics). Use the `palace-tools` CLI (or the `palace_*` MCP tools of
the same names; `python -m palace_tools.mcp_server`). Every tool takes one JSON object and prints one JSON object.
`palace-tools list` shows all tools; `palace-tools spec` shows the config spec.

## Workflow (do not skip steps 4, 5 and 7)

1. **Spec.** Write the figure of merit, its tolerance, the reference you will compare with, the band and the
   materials. Eigenmode for resonances/Q, Driven for S-parameters.
2. **Geometry and mesh.** `palace_tools.mesh.RFModel` or `palace-tools mesh-template`. Named physical groups for
   every volume and every boundary/port; check the mesh report (`verdict: ok`).
3. **Config.** `palace-tools build-config '{"spec": ..., "out": "case/config.json"}'`.
4. **Lint.** `palace-tools lint-config '{"config": "case/config.json", "geometry": "case/<name>.tags.json"}'`.
   Fix all errors; justify or fix each warning.
5. **Estimate and host check.** `palace-tools estimate ...` and `palace-tools status '{"ranks": N}'`. If the
   typical time is > 20 min or memory > 60 % of RAM, show the plan and estimate to the user before running.
   If another Palace job is running, use <= 2 ranks or wait.
6. **Run.** `palace-tools run '{"case_dir": "case", "ranks": N}'` (blocks; refuses on a busy or short-of-RAM
   host unless `force`). Under WSL never detach Palace with nohup/&: run the command itself in the background.
7. **Check.** `palace-tools check '{"case_dir": "case", ...}'`: passivity, reciprocity, energy balance, symmetry,
   expected eigenmode, log errors, spurious port modes, convergence. Then compare with the reference and refine
   the mesh once on the figure of merit.
8. **Report.** `palace-tools results '{"case_dir": "case", "touchstone": "case/out"}'`; state the reference, the
   difference, convergence, cost (unknowns, ranks, GB, minutes) and the simplifications.

## Hard rules (Palace 0.16)

- `Solver.Linear.Type: "MUMPS"`. Never "Default" (it is SuperLU in common builds: rank-dependent wrong answers).
- Set `Model.L0` (Palace assumes micrometres).
- One port family per driven model: all lumped or all wave ports (mixed -> no S-parameters). Resistors next to
  wave ports: Impedance sheet with Rs = R * width / length.
- Each excited port gets `Excitation` equal to its `Index`.
- Lumped R/L/C are circuit values in ohm/H/F on flat rectangular patches with Direction along an edge.
- Lumped ports with different R give wrong |S| in 0.16 (low by R1/R2); use equal R and renormalise.
- Edge lumped ports: PMC end walls (leave them unlisted), not absorbing or PEC.
- Absorbing walls >= 0.3-0.5 free-space wavelengths away, or a PEC shield; check the energy balance.
- Eigenmode: `Target` below the lowest wanted mode; direct solve (`MGMaxLevels: 1`) with lumped elements.
- AMR: set `MaxSize`; with wave ports use `"Nonconformal": false`; adapt at one frequency, then sweep.
- Symmetry planes: PEC (tangential E = 0) or PMC (tangential H = 0); rescale cut lumped elements (C/2 when cut
  along the width, 2C across the length); a driven excitation must respect the symmetry.

All 44 lint rules with their lessons: `palace-tools gotchas '{}'` or
`skills/palace-simulation/references/gotchas.md`.
