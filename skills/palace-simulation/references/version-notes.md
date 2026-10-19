# Palace version notes that change results or configs

The rules target Palace 0.16 (`palace-tools lint-config` takes `"palace_version"` to switch version-dependent
rules). Checked against the Palace CHANGELOG and the 0.16 source.

| topic | 0.16 behaviour | changed in |
|---|---|---|
| Mixed lumped + wave ports | no port-S.csv (S measurement returns early) | 0.17 enables S with both |
| Lumped ports with different R | \|S_ij\|^2 low by R_j/R_i (measured 2/3 for 50 -> 75 ohm) | 0.17 fix (PR 743): drop any manual R-ratio correction after upgrading |
| Wave-port mode with lossy conductors | Conductivity treated as PEC in the 2-D port solve: lossless kn, de-embedded lengths carry no conductor loss | 0.17 changed the wave-port formulation to support impedance boundaries |
| Nonconformal AMR with wave ports | erroneous results / spurious port modes | 0.17 fix |
| `Problem.Output` missing | abort ("Invalid output directory") | 0.17 defaults to "postpro" |
| ParaView output with several excitations | only the last excitation kept | 0.18 fix |
| Excited lumped port with L/C | not allowed (must be purely resistive) | 0.18 allows reactive excited ports |
| `Eigenmode.Target` | required | (0.16 made it always required) |
| `Linear.Type: "Default"` for frequency-domain problems | SuperLU if built with it, else STRUMPACK, else MUMPS | unchanged |
| Driven `Samples` (Linear/Log/Point), integer `Excitation`, unique Index across port kinds, ComplexCoarseSolve | available | introduced in 0.14 |
| `Problem.OutputFormats`, `Eigenmode.NonlinearType`, `Periodic.BoundaryPairs` | available | introduced in 0.15 |
| Schema: PEC/Ground and PMC/ZeroCharge mutually exclusive | enforced in later 0.16.x | 0.16.1 |
| Log line "boundary modes at wave ports for w/2pi = X GHz" | X is 2 pi f (a print bug); `palace_tools.logparse` divides by 2 pi | |

Defaults worth knowing (0.16): `Model.L0` = 1e-6 (micrometres), `Solver.Order` = 1, `MGMaxLevels` = 100
(p-multigrid on), `EstimatorTol` = 1e-6 with Jacobi (`EstimatorMG` false), `Refinement.Nonconformal` = true,
`UpdateFraction` = 0.7, `AdaptiveMaxSamples` = 20 per excitation, ParaView output on.
