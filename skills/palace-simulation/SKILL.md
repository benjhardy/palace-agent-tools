---
name: palace-simulation
description: Set up, lint, cost, run, check and post-process AWS Palace finite-element EM simulations (eigenmode and driven S-parameters) with the palace-agent-tools package. Use whenever a task involves Palace configs, Gmsh meshes for Palace, running Palace (WSL/Apptainer/Docker/native), Palace output CSVs (port-S.csv, eig.csv), or judging whether a Palace result can be trusted.
---

# Palace simulation workflow

Palace is easy to run and easy to run wrong: most failures give a plausible-looking answer and no error. Work
through these steps in order and do not skip the lint, estimate and check steps. Every tool takes one JSON
object and returns one (`palace-tools <tool> '<json>'`, or the `palace_*` MCP tools of the same names).

Target: Palace 0.16 (container builds from early 2026). Version differences that change results are in
[references/version-notes.md](references/version-notes.md).

## 1. Spec

Write down, before any geometry: the quantity you need (S-parameters over a band, a resonance and its Q, a
loss), its acceptance tolerance, the reference you will compare with (closed form, published data, a second
solver, a mesh study), the frequency range, and the materials. Choose the problem type:

- **Eigenmode** for resonances / Q / energy participation. Lumped C/L are LumpedPorts with only C or L set.
- **Driven** for S-parameters. Pick one port family: all lumped or all wave ports (Palace 0.16 writes no S
  when they are mixed; model resistors as an Impedance sheet next to wave ports).

## 2. Geometry and mesh

Use `palace_tools.mesh.RFModel` (Gmsh OCC) or a template (`palace-tools mesh-template`). Rules:

- Every volume is a named 3-D physical group; every boundary, port and lumped patch a named 2-D group.
  Untagged outer faces become PMC silently.
- Remove thick conductors from the mesh and put a Conductivity (or PEC) boundary on their surfaces.
- Lumped patches must be flat rectangles with the Direction along one edge (Palace aborts > 1 degree off).
- Edge lumped ports: PMC end walls (leave them unlisted), not absorbing (20 % leak) or PEC (shorts the strip).
- Wave ports: make each the whole domain end face.
- Absorbing walls at >= 0.3-0.5 wavelengths (free space, highest frequency) from the structure, or a PEC
  shield; walls in the near field fake a flat loss.
- Use symmetry planes when the excitation allows (PEC plane: tangential E = 0; PMC: tangential H = 0).
  A lumped element cut along its width needs C/2 (L x 2, R x 2); cut across its length 2C.

The mesh writer saves `<name>.msh` (MSH 2.2) and `<name>.tags.json` (group names -> tags, geometry used by the
linter, and the exact Nedelec unknown count per order). Check the mesh report: `verdict` must be `ok`.

## 3. Config

Build it from a compact spec (`palace-tools spec` prints every key):

```
palace-tools build-config '{"spec": "spec.json", "out": "case/config.json"}'
```

The builder sets: units (L0, Palace defaults to micrometres!), MUMPS (never 'Default'/SuperLU), a direct
factorisation for eigenmode and small/medium driven problems, ComplexCoarseSolve, the cheap error estimator
when AMR is off, Excitation = Index for every excited port, no ParaView output unless asked.

## 4. Lint

```
palace-tools lint-config '{"config": "case/config.json", "geometry": "case/line.tags.json"}'
```

Fix every `error`. Read every `warning` and either fix it or write down why it does not apply. The findings
carry the lesson and its source. Full list: [references/gotchas.md](references/gotchas.md).

## 5. Cost estimate, then the host

```
palace-tools estimate '{"config": "case/config.json", "tags": "case/line.tags.json", "ranks": 4}'
palace-tools status '{"ranks": 4}'
```

The estimate (calibrated on ~200 recorded Palace 0.16 runs) gives RAM and wall time as quiet / typical / busy.
If the typical time is over ~20 minutes or memory over ~60 % of the machine, tell the user the plan and the
estimate and get a go-ahead before running. If another Palace job is running, use at most 2 ranks or wait.
See [references/running.md](references/running.md).

## 6. Run

```
palace-tools run '{"case_dir": "case", "ranks": 4}'
```

It refuses when another job would be oversubscribed or RAM is short (`force` overrides). It blocks until Palace
exits and writes `palace.log`, `time.log` and `run_record.json`. Under WSL never `nohup`/`&` Palace: WSL kills
detached jobs when the last `wsl.exe` exits. To run in the background, start this command itself as a background
job of your shell/agent tool and poll `palace-tools status '{"case_dir": "case"}'`.

## 7. Check

```
palace-tools check '{"case_dir": "case", "expect": {"eigen_ghz": 8.24, "rel_tol": 0.001}, "symmetric_pairs": [[1, 2]]}'
```

Passivity (|S| <= 1), reciprocity, energy balance (fake wall loss shows up as a flat loss at the lowest
frequency), port symmetry, log errors (Palace only warns about ignored boundary attributes), wave-port spurious
modes, eigenmode at the expected frequency, AMR convergence. Then compare with your reference and do at least
one mesh-convergence step on the figure of merit. Meanings and fixes: [references/checks.md](references/checks.md).

## 8. Results and report

```
palace-tools results '{"case_dir": "case", "touchstone": "case/line"}'
```

Report: the figure of merit with its reference and the difference; mesh/AMR convergence; cost (unknowns, ranks,
peak GB, minutes); the model simplifications; "simulation only" if nothing was measured.

## Recipes

Spec examples for eigenmode with lumped capacitors, driven lumped/wave ports, resistive sheets, symmetry and
AMR: [references/config-recipes.md](references/config-recipes.md). Palace's own docs (0.16):
https://awslabs.github.io/palace/v0.16.0/ (config reference under `config/`, guides under `guide/`).
