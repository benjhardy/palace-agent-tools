# Config recipes (compact spec for `palace-tools build-config`)

Tags are physical-group names from the mesh (`*.tags.json`). Lengths are in the mesh units given by `units`.
Frequencies are in GHz. Circuit values are SI (ohm, H, F).

## Eigenmode: cavity / resonator

```json
{"problem": "eigenmode", "mesh": "cavity.msh", "units": "mm", "tags": "cavity.tags.json", "order": 2,
 "materials": [{"tags": ["cavity"], "eps_r": 1.0, "tan_d": 1e-3}],
 "pec": ["walls"],
 "eigen": {"target_ghz": 5.0, "n": 4, "tol": 1e-8, "expected_ghz": 8.24}}
```

`target_ghz` must be below the lowest mode you want: Palace finds the `n` modes above it. With a crowded
spectrum (birdcage coils: 20 modes in 128-483 MHz) put it well below (0.5 x f0) and ask for more modes.
Q comes out directly when anything is lossy (here Q = 1/tan_d = 1000).

## Eigenmode with lumped capacitors (birdcage-style)

```json
{"problem": "eigenmode", "mesh": "coil.msh", "units": "m", "tags": "coil.tags.json", "order": 2,
 "materials": [{"tags": ["air"]}], "pec": ["coil", "shield"],
 "lumped_elements": [{"index": 1, "tags": ["gap1"], "direction": "+Z", "C": 1.83e-12},
                     {"index": 2, "tags": ["gap2"], "direction": "+Z", "C": 1.83e-12, "cut": "half-width"}],
 "symmetry": [{"tags": ["sym_x"], "type": "pec"}, {"tags": ["sym_y"], "type": "pmc"}],
 "eigen": {"target_ghz": 0.07, "n": 10}}
```

The builder uses MUMPS with `MGMaxLevels: 1` (direct): p-multigrid took ~1 min per shift-and-invert solve on
such a model. `cut` rescales an element split by a symmetry plane.

## Driven, edge lumped ports (microstrip, PCB)

```json
{"problem": "driven", "mesh": "line.msh", "units": "mm", "tags": "line.tags.json", "order": 2,
 "materials": [{"tags": ["substrate"], "eps_r": 3.55, "tan_d": 0.0027}, {"tags": ["air"]}],
 "pec": ["ground"], "conductivity": [{"tags": ["copper"], "sigma": 5.8e7}],
 "pmc": ["ends"], "absorbing": {"tags": ["sides", "top"], "order": 1},
 "lumped_ports": [{"index": 1, "tags": ["port1"], "direction": "+Z", "R": 50},
                  {"index": 2, "tags": ["port2"], "direction": "+Z", "R": 50}],
 "frequency": {"min_ghz": 0.5, "max_ghz": 10, "step_ghz": 0.05, "adaptive_tol": 1e-3}}
```

Edge lumped ports add their own series inductance (~0.1 nH for a 0.5 mm-tall sheet). For Z0/eeff extraction
use two line lengths (see `examples/microstrip_thru`).

## Driven, wave ports, resistor as a sheet (Wilkinson-style)

```json
{"problem": "driven", "mesh": "wilk.msh", "units": "mm", "tags": "wilk.tags.json", "order": 2,
 "materials": [{"tags": ["substrate"], "eps_r": 3.55, "tan_d": 0.0027}, {"tags": ["air"]}],
 "pec": ["ground"], "conductivity": [{"tags": ["copper"], "sigma": 5.8e7}],
 "absorbing": {"tags": ["walls"]}, "wave_port_pec": ["copper"],
 "impedance": [{"tags": ["r_iso"], "R": 100, "length": 0.5, "width": 0.5}],
 "wave_ports": [{"index": 1, "tags": ["port1"], "offset": 4.7}, {"index": 2, "tags": ["port2"], "offset": 5.0},
                {"index": 3, "tags": ["port3"], "offset": 5.0}],
 "frequency": {"min_ghz": 4, "max_ghz": 12, "step_ghz": 0.05, "adaptive_tol": 1e-3}}
```

`impedance` with `R`, `length` (along the current) and `width` sets Rs = R w / l. Palace 0.16 gives no
S-parameters with a lumped port next to wave ports, hence the sheet.

## Internal ports for port tuning

Every component footprint becomes a 50-ohm lumped port (`lumped_ports`, all excited). The tuner then terminates
the internal ports in circuit simulation. Keep the external ports lumped too (0.16), and expect the online phase
of the adaptive sweep to cost ~2-3 s per frequency per excitation at ~250k unknowns: use a coarse step.

## AMR

```json
"amr": {"max_its": 6, "max_size": 400000, "nonconformal": false, "update_fraction": 0.7}
```

Adapt at one frequency (`"frequency": {"points_ghz": [10.0]}`), then sweep on the adapted mesh. The size cap is
checked after refining (450k cap ended at 623k). Use conformal refinement with wave ports in 0.16.

## Escape hatch

`"extra": {...}` is deep-merged into the generated config last, for any Palace key the spec does not cover.
