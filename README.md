# palace-agent-tools

[![tests](https://github.com/benjhardy/palace-agent-tools/actions/workflows/tests.yml/badge.svg)](https://github.com/benjhardy/palace-agent-tools/actions/workflows/tests.yml)

Part of agent.rf — the RF lab: https://benjhardy.github.io/agent-rf/

Tools that let AI coding agents (Claude Code, Codex, Cursor, or a plain script) set up, lint, run, check and
post-process [AWS Palace](https://github.com/awslabs/palace) finite-element electromagnetic simulations
reliably.

Palace is a capable open-source FEM solver, but most mistakes in a Palace model produce a plausible answer and
no error: a mesh read in micrometres, a boundary tag that Palace silently ignores, absorbing walls that fake a
loss, a solver choice that gives wrong points on some rank counts, S-parameters that are never written. This
package turns the lessons from a series of Palace projects (MRI birdcage coil, Wilkinson divider, lumped-filter
port tuning, an HFSS-vs-Palace comparison, a GDS-to-Palace benchmark, a mmWave PA) into tools an agent calls
with one JSON object and gets one JSON object back:

| step | tool | what it does |
|---|---|---|
| geometry, mesh | `mesh-template`, `palace_tools.mesh.RFModel` | Gmsh OCC parts (boxes, substrates, traces, vias, port rectangles, air box padded in wavelengths), physical-group conventions, wavelength and feature size fields, MSH 2.2 + a `*.tags.json` with group tags, geometry and a mesh report (quality, exact Nedelec unknowns per order) |
| config | `build-config`, `spec` | Palace JSON from a compact spec with names instead of tag numbers and safe defaults |
| lint | `lint-config`, `lint-rules`, `gotchas` | Palace 0.16 schema validation plus 44 rules from real failures, each with the lesson and its source |
| cost | `estimate` | peak RAM and wall time (quiet / typical / busy) calibrated on 190 recorded Palace 0.16 runs |
| run | `status`, `run`, `doctor` | WSL + Apptainer, Apptainer, Docker or native; refuses to oversubscribe a busy machine or run out of RAM; logs wall time and memory |
| results | `results` | port-S / eig / EPR / Q / energy / AMR / field summaries, Touchstone with the right reference impedances |
| check | `check` | passivity, reciprocity, energy balance, port symmetry, expected eigenmode, log errors, spurious wave-port modes, convergence |

It ships as a Python package with a CLI, an MCP server (stdio) and a Claude Code skill.

**Who it is for:** RF/microwave and EM engineers who use Palace with an AI agent, and agent builders who want
an EM solver step their agent cannot quietly get wrong.

**Status:** targets Palace 0.16 (the version the rules and cost data come from), with version switches for the
behaviour that changed in 0.17 and 0.18. Simulation-only examples; nothing here was measured.

## Install

```
python -m venv .venv
. .venv/bin/activate            # Windows: .venv\Scripts\activate
pip install -e ".[all]"          # numpy + gmsh + scikit-rf + jsonschema;  ".[dev]" adds pytest
```

The core (config, lint, estimate, results, checks) needs only numpy. `gmsh` is needed for meshing, `jsonschema`
for schema validation, `scikit-rf` for Touchstone files. Palace itself is not installed by this package: use
the official container or a native build (see [Palace install](https://awslabs.github.io/palace/v0.16.0/install/)),
then tell the runner how to call it (`PALACE_TOOLS_BACKEND`, `PALACE_TOOLS_IMAGE`, `PALACE_TOOLS_CMD`,
`PALACE_TOOLS_DISTRO`; `palace-tools doctor` shows what it found).

## Quick start (human)

```
python examples/cavity_te101/run_example.py --ranks 2       # 10 s: eigenmodes vs analytic TE101 ...
python examples/lumped_port_line/run_example.py --ranks 2   # 1.5 min: lumped ports on an exact TEM line
python examples/microstrip_thru/run_example.py --ranks 2    # 2.5 min: microstrip Z0, eeff vs Hammerstad-Jensen
```

Or step by step with the CLI:

```
palace-tools mesh-template '{"template": "microstrip_line", "out": "case/line", "w": 1.137, "h": 0.508,
    "length": 20, "eps_r": 3.55, "side": 4.5, "top": 4.5, "h_mesh": 1.6, "h_edge": 0.12, "ports": "wave", "f_ghz": 10}'
palace-tools build-config @spec.json          # {"spec": {...}, "out": "case/config.json"}
palace-tools lint-config '{"config": "case/config.json", "geometry": "case/line.tags.json"}'
palace-tools estimate '{"config": "case/config.json", "tags": "case/line.tags.json", "ranks": 4}'
palace-tools status '{"ranks": 4}'
palace-tools run '{"case_dir": "case", "ranks": 4}'
palace-tools check '{"case_dir": "case", "symmetric_pairs": [[1, 2]]}'
palace-tools results '{"case_dir": "case", "touchstone": "case/line"}'
```

`palace-tools list` prints every tool. Arguments can come inline, from `@file.json`, or from stdin (`-`). Exit
code 1 means the tool reported a problem (`ok: false` or `overall: fail`), 2 means bad input.

## Quick start (agents)

**MCP (any MCP client).** The server is dependency-free JSON-RPC over stdio:

```
claude mcp add palace-tools -- python -m palace_tools.mcp_server
```

or the repository's `.mcp.json`. Tools: `palace_spec`, `palace_build_config`, `palace_lint_config`,
`palace_lint_rules`, `palace_gotchas`, `palace_mesh_template`, `palace_mesh_report`, `palace_estimate_cost`,
`palace_status`, `palace_run`, `palace_results`, `palace_check`, `palace_doctor`. Each has a JSON input schema;
results come back as JSON text and `structuredContent`.

**Claude Code skill.** Copy `skills/palace-simulation/` into `~/.claude/skills/` (or a project's
`.claude/skills/`). It teaches the workflow spec -> geometry -> mesh -> config -> lint -> cost estimate -> run
-> check -> report, with the gotchas, recipes, run rules, check meanings and version notes in `references/`.

**Other agents.** [AGENTS.md](AGENTS.md) has the same workflow and rules in plain form.

## The config spec

```json
{"problem": "driven", "mesh": "line.msh", "units": "mm", "tags": "line.tags.json", "order": 2,
 "materials": [{"tags": ["substrate"], "eps_r": 3.55, "tan_d": 0.0027}, {"tags": ["air"]}],
 "pec": ["ground", "strip"], "pmc": ["ends"],
 "lumped_ports": [{"index": 1, "tags": ["port1"], "direction": "+Z", "R": 50},
                  {"index": 2, "tags": ["port2"], "direction": "+Z", "R": 50}],
 "frequency": {"min_ghz": 0.5, "max_ghz": 10, "step_ghz": 0.25, "adaptive_tol": 1e-3}}
```

Also: `absorbing`, `conductivity`, `impedance` (give `R`, `length`, `width` and it sets Rs = R w / l),
`lumped_elements` (with `cut: half-width | half-length` for elements split by a symmetry plane), `wave_ports`,
`wave_port_pec`, `symmetry`, `eigen`, `amr`, `solver`, `probes`, `energy_domains`, `save_fields`, and `extra`
(raw Palace keys merged last). Defaults: MUMPS; a direct factorisation (`MGMaxLevels: 1`) for eigenmode and
for driven problems up to 600k unknowns; `ComplexCoarseSolve`; `EstimatorMG` + `EstimatorTol 1e-4` when AMR is
off; `Excitation` = port index; no ParaView output unless asked. `palace-tools spec` prints every key.

## Examples (Palace 0.16, 8-thread desktop, 2 MPI ranks unless noted)

| example | reference | result | size / cost |
|---|---|---|---|
| [`cavity_te101`](examples/cavity_te101): 22.86 x 10.16 x 30 mm PEC cavity, tan d = 1e-3 | analytic TE/TM resonances; exact dielectric Q | TE101 8.243934 vs 8.243877 GHz (+0.0007 %); first 5 modes within 0.02 %; Q = 1000.00 vs 1000.00 | 17.8k unknowns (order 2), 0.24 GB, 9 s |
| [`lumped_port_line`](examples/lumped_port_line) a: 75-ohm parallel-plate line (PMC sides), 50-ohm lumped ports | closed-form line S-matrix | max \|S - S_ref\| = 8.7e-7 over 0.5-10 GHz | 34.5k unknowns, 0.65 GB, 55 s |
| `lumped_port_line` b: same line terminated in a 75-ohm lumped element | \|S11\| = 0.2 (-13.979 dB), flat | -13.97940 to -13.97942 dB; also exposes the 0.16 unequal-R bug: \|S21\|^2 = 2/3 of the delivered power (= R1/R2), flagged by `check` | 23 s |
| [`microstrip_thru`](examples/microstrip_thru) (wave ports): 50-ohm line, eps_r 3.55, h 0.508 mm, 20 mm, PEC shield at 9h | Kirschning-Jansen eeff(f) | eeff within 0.4 % over 1-10 GHz; \|S11\| < -40 dB | 105k unknowns, 2.5 GB, 78 s on 4 ranks (104 s on 2) |
| `microstrip_thru` (edge lumped ports, 10 + 20 mm lines, two-line de-embedding) | Hammerstad-Jensen Z0 = 50.0 ohm; Kirschning-Jansen eeff | Z0 49.5 ohm median (48.6-51.7) over 0.5-10 GHz; eeff within 0.25 %; fitted port inductance 143 pH. Single-line extraction without de-embedding: eeff +5.6 %, Z0 40-60 ohm above 5 GHz | 70k + 109k unknowns, 1.5 + 2.4 GB, 45 + 84 s |

Every run passed the passivity, reciprocity, energy-balance, port-symmetry and log checks, except case b,
where the energy check flags the 0.16 bug, as it should. The examples write their summaries to
`results*.json` next to each script.

## Gotchas (what the lint rules encode)

All 44 rules with their full lessons: [skills/palace-simulation/references/gotchas.md](skills/palace-simulation/references/gotchas.md)
(generated from the code). The ones that cost the most time:

| gotcha | what happened | where it was learned | rule |
|---|---|---|---|
| `Linear.Type: "Default"` is SuperLU in the 0.16 container | partition-dependent wrong points (S21 ~ -45 dB at 2 of 21 frequencies on 6 ranks) and "singular matrix" aborts; MUMPS identical on any rank count | GDS-to-Palace benchmark | SOL001 |
| p-multigrid on an eigenproblem with lumped C | ~1 min per shift-and-invert solve; direct MUMPS (`MGMaxLevels: 1`) did the whole 484k run in 6.6 min | birdcage | SOL002 |
| default error estimator when AMR is off | half of the run time; `EstimatorMG` + `EstimatorTol 1e-4` gave bit-identical S 5.8x faster | GDS-to-Palace benchmark | SOL003 |
| mixed wave and lumped ports | Palace 0.16 writes no `port-S.csv`; resistor as an Impedance sheet, or all ports lumped | Wilkinson, port tuning | PRT001 |
| lumped ports with different R | 0.16 reports \|S21\|^2 low by R1/R2 (fixed in 0.17: drop the correction after upgrading) | `lumped_port_line` example, benchmark combiner | PRT011 |
| edge lumped port beside an absorbing wall | ~20 % of the port power lost into the wall; PMC end walls are lossless | port tuning | PRT012 |
| curved or tilted lumped patch | Palace aborts (projected-length / 1-degree axis check); make patches flat, Direction along an edge | birdcage | PRT009 |
| element cut by a symmetry plane | needs C/2 along its width (full C: -7.6 % frequency error), 2C across its length | birdcage | BND004 |
| absorbing walls in the near field | 0.03 dB fake loss at 0.2 lambda (Wilkinson); 5 % flat loss and S11 != S22 at 1 mm (benchmark) | Wilkinson, benchmark | BND001 |
| wave port smaller than the mode | 3-6 % loss falling with frequency on a bare line; use the whole end face | benchmark | PRT013 |
| adaptive sweep in a closed box | aborted after the first sample ("singular matrix", SuperLU); direct points worked | benchmark | BND002 |
| adaptive sweep online phase | 2.6 s per frequency per excitation at 258k unknowns: 97 of 128 min | port tuning | DRV003 |
| nonconformal AMR with wave ports | the port mode jumped to a spurious mode (n_eff = sqrt(4.4)); conformal refinement fixed it | HFSS-vs-Palace | PRT014 |
| wave-port de-embedding with lossy conductors | the de-embedded length carries no conductor loss in 0.16 (-0.001 vs -0.30 dB per inch) | HFSS-vs-Palace | PRT015 |
| AMR size cap | checked after refining: a 450k cap ended at 623k unknowns, 14 GB | HFSS-vs-Palace | AMR001 |
| several excitations with field output | only the last excitation's ParaView fields are kept (0.16) | Wilkinson | PRT016 |
| `Model.L0` missing | the mesh is read in micrometres | Palace docs | MOD001 |
| tag typo | Palace prints "Unknown ... boundary attributes! Solver will just ignore them!" and runs on | Palace source | BND005 |
| WSL and `nohup` | WSL kills detached jobs when the last `wsl.exe` exits; keep the run in the foreground of a background job | benchmark, PA | (runner) |
| floating port / unresolved via | \|S21\| ~ -50 dB with no error; S22 = +0.8 dB from one tiny via under a pad | PA, benchmark | (`check`) |

Version-dependent behaviour (0.16 vs 0.17/0.18) is listed in
[references/version-notes.md](skills/palace-simulation/references/version-notes.md).

## Cost estimate

`palace-tools estimate` fits power laws in the number of unknowns to 190 recorded Palace 0.16 runs (eigenmode
and driven, direct and p-multigrid, orders 1-3, 1k to 1.6M unknowns) from the projects above. Anchors: a 1.12 M
unknown order-2 eigenmode solve with MUMPS on 4 ranks used 16.4 GB and 11-19 min (estimate 15.5 GB, 12 min);
a 1.32 M unknown driven sweep with p-multigrid used 10 GB (estimate 10.7 GB). On the recorded runs the median
memory error is 4 % (eigenmode, direct), 10 % (driven, p-multigrid) and 14 % (driven, direct; 80th percentile
41 %). Wall time is only good to 2-3x, because many recorded runs shared the CPU with another job. The estimate is a
go/no-go aid for "ask the user first", not a schedule.

## Limits

- Rules and cost data come from Palace 0.16 on one 8-thread, 32 GB desktop (WSL2 with 25 GB). Times do not
  scale with rank count in the model; clusters will differ.
- The lint rules catch known failure patterns; they do not prove a model right. Always compare with a reference
  and refine the mesh once on the figure of merit.
- Geometry-aware rules need the `*.tags.json` from `palace_tools.mesh`; meshes from other tools get the
  config-only rules.
- The mesh helpers cover axis-aligned planar RF structures (boxes, traces, vias, rectangular ports) and three
  templates; complex 3-D CAD should be built with Gmsh directly, using the same group conventions.
- The MCP server runs tools synchronously; `palace_run` blocks until Palace exits.
- Not affiliated with AWS or the Palace developers.

## Development

```
pip install -e ".[dev]"
pytest                      # 67 tests, ~3 s, no Palace needed
python tools/gen_gotchas.py # regenerate the gotchas reference after editing rules
python tools/harvest_costs.py label=<runs folder> ...   # add your own finished runs to the cost model
```

## Licence

Apache License 2.0, the same licence as Palace. The vendored Palace config schema
(`palace_tools/data/schema_0.16/`) is Copyright Amazon.com, Inc. or its affiliates, Apache-2.0; see
[NOTICE](NOTICE).
