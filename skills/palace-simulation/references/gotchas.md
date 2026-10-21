# Palace gotchas (generated from palace_tools.lint; do not edit by hand)

Each lint rule below is checked by `palace-tools lint-config`. Rules marked *geometry* need the `*.tags.json` written by `palace_tools.mesh`. Targets Palace 0.16; version differences are stated in the rule.

## SOL001: Sparse direct solver: use MUMPS, not SuperLU / 'Default'

SuperLU (what 'Default' resolves to when Palace is built with it, as the 0.16 container is) gave partition-dependent wrong points (S21 ~ -45 dB at 2 of 21 frequencies on 6 ranks) and 'singular matrix' aborts in closed boxes. MUMPS gave identical results on any rank count.

Source: gds-bench (G01, validation V4). Docs: https://awslabs.github.io/palace/v0.16.0/config/solver/

## SOL002: Direct factorisation vs p-multigrid

Eigenmode with lumped elements: p-multigrid took ~1 min per shift-and-invert solve at order 2; a direct MUMPS factorisation (MGMaxLevels 1) did the whole 484k-unknown run in 6.6 min. Driven, 102k unknowns: full complex MUMPS gave the same S to 7e-9, 5.6x faster than p-multigrid, at 2x the memory. Above roughly 0.6-1 M unknowns the factorisation memory (~15-25 kB per unknown) decides.

Source: birdcage, port-tuning. Docs: https://awslabs.github.io/palace/v0.16.0/config/solver/

## SOL003: Error-estimator cost when AMR is off

Palace always runs the flux-projection error estimate. With the defaults (Jacobi PCG, EstimatorTol 1e-6) it took 1850 of 3660 s on an inductor; EstimatorMG true + EstimatorTol 1e-4 gave bit-identical S-parameters 5.8x faster overall.

Source: gds-bench (validation V1/V3). Docs: https://awslabs.github.io/palace/v0.16.0/config/solver/

## SOL004: Complex coarse solve for driven problems

With MGMaxLevels 1 the 'coarse' solve is the whole problem. ComplexCoarseSolve true factorises the true complex matrix, so GMRES converges in 1-2 iterations instead of tens.

Source: port-tuning. Docs: https://awslabs.github.io/palace/v0.16.0/config/solver/

## SOL005: Krylov iteration limit

An iterative solve that stops at MaxIts gives a wrong answer at that frequency; Palace only warns in the log.

Source: Palace docs. Docs: https://awslabs.github.io/palace/v0.16.0/config/solver/

## SOL006: Element order

Order 1 shifted a birdcage k=1 mode by 1.25 % against order 2 on the same mesh; order 3 costs ~3x the unknowns per element of order 2 (21 vs 7 per tet).

Source: birdcage, microstrip-amr. Docs: https://awslabs.github.io/palace/v0.16.0/config/solver/

## PRT001: Mixed wave and lumped ports

Palace 0.16 writes no port-S.csv when lumped and wave ports are both present (PostOperator::MeasureSParameter returns early). Fixed in 0.17. Workarounds: model a resistor as an Impedance sheet (Rs = R w / l) and keep wave ports, or make every port lumped.

Source: wilkinson, port-tuning. Docs: https://awslabs.github.io/palace/v0.16.0/guide/boundaries/

## PRT002: Excitation index must equal port index

For one-port-per-excitation runs Palace assumes excitation index == port index; anything else gives no or mislabelled S columns.

Source: Palace docs (guide/boundaries), source. Docs: https://awslabs.github.io/palace/v0.16.0/guide/boundaries/

## PRT003: Several ports in one excitation

Ports sharing an excitation index are driven together: Palace writes port-V/port-I but no port-S. Useful for grouped/combiner drive; compute the metric from port-V.csv and port-I.csv.

Source: gds-bench (G02 combiner). Docs: https://awslabs.github.io/palace/v0.16.0/guide/boundaries/

## PRT004: Surface current source disables S-parameters

Surface currents are added to every excitation, so no excitation is 'simple' and S is not computed.

Source: Palace source (portexcitations). Docs: https://awslabs.github.io/palace/v0.16.0/config/boundaries/

## PRT005: Port indices must be unique across port kinds

Index values must be unique across LumpedPort, WavePort, SurfaceCurrent and Terminal (Palace >= 0.14).

Source: Palace docs. Docs: https://awslabs.github.io/palace/v0.16.0/guide/boundaries/

## PRT006: Driven problem without an excitation

A driven run needs at least one excited port, a SurfaceCurrent or a CurrentDipole.

Source: Palace docs. Docs: https://awslabs.github.io/palace/v0.16.0/guide/problem/

## PRT007: Lumped port values: circuit or per-square, not both, not neither

Palace aborts if a LumpedPort has both R/L/C and Rs/Ls/Cs, or none of them.

Source: Palace source (lumpedportoperator). Docs: https://awslabs.github.io/palace/v0.16.0/config/boundaries/

## PRT008: Excited lumped port must be resistive

Palace <= 0.17 aborts on an excited lumped port with R = 0 or with L/C. Put reactive parts on separate, non-excited lumped elements (0.18 allows reactive excited ports).

Source: Palace source / changelog 0.18. Docs: https://awslabs.github.io/palace/v0.16.0/config/boundaries/

## PRT009: Lumped element direction vs patch shape *(geometry)*

Palace measures a lumped element along its Direction on the patch's bounding box: the direction must be within 1 degree of a box axis (abort) and the patch must be flat along it (projected-length check). Curved ring-gap patches failed until made planar; a direction along the patch normal gives zero length.

Source: birdcage. Docs: https://awslabs.github.io/palace/v0.16.0/config/boundaries/

## PRT010: Lumped element units (H, F, ohm; circuit vs per-square)

R/L/C are circuit values in ohm/H/F; Palace converts them to per-square values with the patch: Rs = R w/l, Ls = L w/l, Cs = C l/w (the log prints the per-square values). A value of 1.8 meant as pF is 1.8 F.

Source: birdcage, port-tuning. Docs: https://awslabs.github.io/palace/v0.16.0/config/boundaries/

## PRT011: Lumped ports with different R in Palace 0.16

Palace 0.16 reports |S_ij|^2 between lumped ports of different R low by exactly R_j/R_i (measured: 50 -> 75 ohm gave |S21|^2 = 2/3 of the delivered power; a 400-ohm node port needed |S21|^2 x 8). Palace 0.17 fixed it ('incorrect S-parameters when lumped port pairs have different R values', PR 743), so the x R_i/R_j correction must NOT be applied to 0.17+ results.

Source: lumped_port_line example, gds-bench (G02), Palace changelog 0.17. Docs: https://awslabs.github.io/palace/v0.16.0/config/boundaries/

## PRT012: Edge lumped port next to an absorbing wall *(geometry)*

With first-order absorbing end walls an edge lumped port lost ~20 % of its power into the wall beside it (eigenvalues of S^H S 0.45-0.54). A PEC end wall shorts the strip; a PMC (natural, unlisted) end wall is lossless and does not touch the strip.

Source: port-tuning. Docs: https://awslabs.github.io/palace/v0.16.0/config/boundaries/

## PRT013: Wave port that does not span the domain face *(geometry)*

A 5w x 6h wave port on an open face next to absorbing walls did not hold the whole microstrip field: 3-6 % loss on a bare line that fell with frequency. Ports spanning the full domain end face (in a shielded box) fixed it; Palace's 2-D port solve treats unlabelled port edges as PMC.

Source: gds-bench (G01). Docs: https://awslabs.github.io/palace/v0.16.0/guide/boundaries/

## PRT014: Wave ports with nonconformal AMR

In Palace 0.16 the port eigen-solve switched to a spurious mode during nonconformal AMR (to the substrate-filled mode n_eff = sqrt(4.4) at 39k unknowns on a CPW; to n_eff = 1 at 295k on a waveguide tee). 'Nonconformal': false ran to 470k unknowns with a steady kn. Palace 0.17 fixed 'nonconformal AMR erroneous results when waveports are present'.

Source: microstrip-amr. Docs: https://awslabs.github.io/palace/v0.16.0/config/model/

## PRT015: Wave-port de-embedding ignores conductor loss in 0.16

Palace 0.16's port-mode solve treats Conductivity boundaries as PEC (lossless kn). A de-embedded (Offset) length therefore carries no conductor loss: a 1-inch CPW read -0.001 dB instead of -0.30 dB, and beta is low by about alpha. 0.17 changed the wave-port formulation to support impedance boundaries.

Source: microstrip-amr. Docs: https://awslabs.github.io/palace/v0.16.0/config/boundaries/

## PRT016: Field output with several excitations

Palace 0.16 kept the ParaView fields of the last excitation only (fixed in 0.18). Run one excitation per field case. ParaView output is ~230 MB per excitation and frequency on a 1.3 M-unknown model.

Source: wilkinson. Docs: https://awslabs.github.io/palace/v0.16.0/config/solver/

## PRT017: Wave ports must be on the outer boundary *(geometry)*

Wave ports are one-sided: mesh elements on one side only. An internal face needs a lumped port.

Source: Palace docs. Docs: https://awslabs.github.io/palace/v0.16.0/guide/boundaries/

## BND001: Absorbing boundary too close (distance in wavelengths) *(geometry)*

First-order absorbing walls in the reactive near field act as a resistive load and fake a nearly frequency-flat loss: Wilkinson walls at 8/6 mm (0.21/0.16 lambda0 at 8 GHz) added 0.03 dB; 12/10 and 16/14 mm agreed within 0.003 dB. Walls 1 mm from a microstrip at 10 GHz gave 5 % flat loss and S11 != S22.

Source: wilkinson, gds-bench (G01). Docs: https://awslabs.github.io/palace/v0.16.0/config/boundaries/

## BND002: Adaptive sweep in a closed lossless box

Palace 0.16's adaptive sweep aborted after the first sample in a closed PEC box ('SuperLUSolver: Found a singular matrix') on 4 and 6 ranks; direct point solves worked. With absorbing walls the adaptive sweep matched direct solves within 0.2 dB.

Source: gds-bench (G01). Docs: https://awslabs.github.io/palace/v0.16.0/config/solver/

## BND003: Finite-conductivity boundary needs thickness >> skin depth

The Conductivity boundary is a surface-impedance model for conductors much thicker than the skin depth; below ~3 skin depths it misestimates the loss.

Source: Palace docs (guide/boundaries). Docs: https://awslabs.github.io/palace/v0.16.0/guide/boundaries/

## BND004: Symmetry planes: PEC/PMC choice and what they keep

A PEC plane forces tangential E = 0 (B1 lies in it; currents odd), a PMC plane tangential H = 0 (B1 crosses it; currents even). PMC is Palace's natural BC, so an unlisted outer face is already PMC. A symmetry model keeps only modes of matching parity, and a lumped element cut by the plane must be rescaled: cut along its width -> C/2 (full C gave -7.6 %), cut across its length -> 2C. A driven port must respect the symmetry.

Source: birdcage. Docs: https://awslabs.github.io/palace/v0.16.0/guide/boundaries/

## BND005: Attributes that are not in the mesh / volumes without material *(geometry)*

Palace only warns 'Unknown ... boundary attributes! Solver will just ignore them!' and runs on: a typo silently removes a boundary condition. Every volume needs a material.

Source: Palace source (iodata). Docs: https://awslabs.github.io/palace/v0.16.0/config/boundaries/

## BND006: Second-order absorbing boundary support

Order-2 absorbing boundaries are only available for frequency-domain problems; transient supports order 1.

Source: Palace docs. Docs: https://awslabs.github.io/palace/v0.16.0/config/boundaries/

## BND007: One attribute, two boundary conditions

An attribute listed under two conditions (e.g. PEC and Absorbing) gives an ill-defined model.

Source: Palace docs. Docs: https://awslabs.github.io/palace/v0.16.0/config/boundaries/

## MOD001: Mesh units (Model.L0)

Palace's default L0 is 1e-6: mesh coordinates are read as micrometres unless L0 is set. A mm mesh without L0 is 1000x too small and every resonance lands 1000x too high.

Source: Palace docs (config/model). Docs: https://awslabs.github.io/palace/v0.16.0/config/model/

## MOD002: Mesh file exists

A relative mesh path is resolved from Palace's working directory, so run Palace from the config folder.

Source: Palace docs (config/model). Docs: https://awslabs.github.io/palace/v0.16.0/config/model/

## MOD003: Output folder

Palace 0.14-0.16 abort with 'Invalid output directory' when Problem.Output is missing (0.17 defaults to 'postpro').

Source: Palace changelog 0.17. Docs: https://awslabs.github.io/palace/v0.16.0/config/problem/

## MOD004: Field output size

ParaView output is written for every saved step and is large (~230 MB per excitation and frequency at 1.3 M unknowns). It is not needed for S-parameters or eigenfrequencies.

Source: wilkinson. Docs: https://awslabs.github.io/palace/v0.16.0/config/problem/

## AMR001: AMR without a size cap

Palace's Dorfler marking (UpdateFraction 0.7) grew meshes 1.3-2.8x per pass; MaxSize is checked after refining, so a 450k cap ended at 623k unknowns (14 GB). Without MaxSize the last pass can exhaust RAM.

Source: microstrip-amr. Docs: https://awslabs.github.io/palace/v0.16.0/config/model/

## AMR002: AMR over a whole frequency sweep

Every AMR pass re-solves every frequency point. The cheaper pattern (HFSS-like): adapt at one or a few frequencies, then sweep on the adapted mesh (load the saved mesh, or reuse the refinement in a second run).

Source: microstrip-amr. Docs: https://awslabs.github.io/palace/v0.16.0/config/model/

## AMR003: AMR convergence is slow on resonant structures

On an 8-rod interdigital filter the band-centre phase of S21 still moved 7 deg between the last two AMR passes at 623k (p2) and 673k (p3) unknowns; seeding the mesh along the resonators and order 3 got closest per unknown. Judge convergence on the quantity you need, between passes.

Source: microstrip-amr. Docs: https://awslabs.github.io/palace/v0.16.0/config/model/

## EIG001: Eigenmode target

Palace 0.16 requires Solver.Eigenmode.Target (GHz) and finds the N modes ABOVE it (shift-and-invert). Put it below the lowest mode you need (the birdcage used 0.55 x f0 for a crowded spectrum).

Source: birdcage, Palace changelog 0.16. Docs: https://awslabs.github.io/palace/v0.16.0/config/solver/

## EIG002: Lossy / nonlinear eigenproblems are slow

A conducting phantom made the birdcage eigenproblem complex: GMRES ~19 iterations per solve and hours for a dozen modes. Wave ports and Conductivity boundaries make it nonlinear (Palace's Hybrid solver).

Source: birdcage. Docs: https://awslabs.github.io/palace/v0.16.0/config/solver/

## EIG003: Excitations are ignored in eigenmode runs

Lumped ports in an eigenmode run are circuit elements (R adds port-Q, L/C enter the eigenproblem); Excitation has no effect. A lumped port with only C set is a lumped capacitor.

Source: birdcage. Docs: https://awslabs.github.io/palace/v0.16.0/guide/problem/

## DRV001: Frequency values and units

Palace frequencies are in GHz. Zero or negative frequencies are invalid for driven runs.

Source: Palace docs. Docs: https://awslabs.github.io/palace/v0.16.0/config/solver/

## DRV002: Adaptive sweep settings

AdaptiveTol 1e-3 kept S within ~1e-3 of direct solves (0.005 dB in a passband) but deep stopbands at -32...-65 dB differed by up to 0.33 dB. Palace needs more than a handful of points, or reverts to a uniform sweep.

Source: port-tuning, wilkinson. Docs: https://awslabs.github.io/palace/v0.16.0/config/solver/

## DRV003: Adaptive sweep online cost

After building the reduced model, Palace 0.16 re-assembles operators at every frequency of every excitation: ~2.6 s each at 258k unknowns (97 of 128 min for 317 points x 7 excitations); with wave ports it also re-solves the port mode per point (a quarter of a 90 min run).

Source: port-tuning, wilkinson. Docs: https://awslabs.github.io/palace/v0.16.0/config/solver/

## DRV004: Many direct frequency points

Each direct point is a full factorisation/solve per excitation; the adaptive sweep needs ~5-12 full solves per excitation for a smooth band.

Source: wilkinson, microstrip-amr. Docs: https://awslabs.github.io/palace/v0.16.0/config/solver/

# Workflow lessons (not lint rules)

## WSL kills detached jobs

When the last wsl.exe exits, WSL stops the VM and kills nohup'd or &-backgrounded Palace runs. Keep a foreground wsl.exe alive for the whole run (palace-tools run blocks; run that command itself in the background).

Source: gds-bench, rfpa-open.

## wsl.exe re-parses bash -lc arguments

Loop variables in 'wsl bash -lc "for ..."' expand early; put loops in a script file.

Source: wilkinson.

## Windows default encoding

Python on Windows writes cp1252 by default and breaks on Greek letters in Palace logs (omega, kn): always pass encoding='utf-8'.

Source: wilkinson.

## Unequal lumped-port R in Palace 0.16

Between lumped ports of different R, 0.16 reports |S21|^2 low by R1/R2 (50 -> 400 ohm: power gain = |S21|^2 x 8). This is a 0.16 bug fixed in 0.17; after the upgrade the x8 correction would be wrong. Prefer equal R and renormalise in scikit-rf.

Source: gds-bench (G02), lumped_port_line example.

## Unresolved vias create gain

A single 0.9 um via under a 70 um pad gave S22 = +0.8 dB; a solid post fixed it. Any |S| > 1 means a modelling error.

Source: gds-bench (G01).

## Floating ports and keep-outs

A via port inside a ground keep-out floats: |S21| ~ -50 dB and no error from the solver. Check port grounding geometrically before a 70-minute solve.

Source: rfpa-open.

## Holes need booleans

Multi-loop addPlaneSurface gave an invalid solid (holes counted as material); cut holes with OCC booleans.

Source: gds-bench.

## Wave-port log frequency

Palace 0.16 prints 'boundary modes at wave ports for w/2pi = X GHz' with X = 2 pi f (rad/ns); divide by 2 pi.

Source: Palace 0.16 source.

## Coarse meshes shift resonances down

A coarse Wilkinson mesh put the minima 1-2 % low; converge on the figure of merit, not only on the AMR indicator.

Source: wilkinson.

## Symmetry pays

Quarter / eighth birdcage models: 4.5x / 7.9x less memory, 5x faster, and the filtered spectrum makes mode identification easier.

Source: birdcage.

## Internal lumped ports are not calibrated

A lumped sheet across a pad gap has its own self-inductance (~0.2 nH for 0.5 mm of 0.6 mm-wide line): port-tuned values carry that error.

Source: port-tuning.
