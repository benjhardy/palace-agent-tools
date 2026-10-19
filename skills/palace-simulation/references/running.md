# Running Palace: backends, guards, cost

## Backends

| backend | when | command run |
|---|---|---|
| `wsl` (default on Windows) | Palace in an Apptainer image inside WSL2 | `run_palace config.json N` if that wrapper exists, else `palace -np N config.json`, or `apptainer run $PALACE_TOOLS_IMAGE -np N config.json` |
| `apptainer` | Linux with a `.sif` | `apptainer run <image> -np N config.json` |
| `docker` | Docker image with `palace` on PATH | `docker run --rm -v <case>:/work -w /work <image> palace -np N config.json` |
| `native` | Palace installed | `palace -np N config.json` |

Environment: `PALACE_TOOLS_BACKEND`, `PALACE_TOOLS_DISTRO` (WSL distribution), `PALACE_TOOLS_IMAGE`,
`PALACE_TOOLS_CMD` (full template with `{config}` and `{np}`). Check with `palace-tools doctor`.

The run is wrapped in `/usr/bin/time -v` (wall time, max RSS of one rank) and Palace's stdout goes to
`palace.log`. Palace also writes `palace.json` (DOFs, ranks, peak memory, timings) into the output folder.

## WSL rules

- WSL stops the VM, and kills every process in it, when the last `wsl.exe` exits. A `nohup ... &` inside
  `wsl -- bash -lc` dies with it. `palace-tools run` blocks, so the `wsl.exe` it starts stays alive; run that
  command as a background job of your own tool if you need to keep working.
- `wsl.exe` re-parses the arguments of `bash -lc`: `$var` and loops expand early. Put loops in a script.
- Paths: `C:\x\y` is `/mnt/c/x/y` inside WSL. Running on `/mnt/c` is fine for Palace-sized I/O.
- Killing a run: stop the process that owns the `wsl.exe` (your background job); `pkill` inside WSL may be
  blocked by sandboxing.

## Guard

`palace-tools status '{"ranks": N}'` lists Palace ranks already running (one `palace-*.bin` process per MPI
rank), CPUs and available RAM, and says whether N ranks may start. With another job running, use <= 2 ranks or
wait. `run` applies the same guard and refuses unless `force` is true.

## Ranks

On an 8-thread desktop, 4 ranks was the sweet spot; 8 ranks was slower than 4 for problems under ~50k unknowns
(MPI overhead plus contention with the host OS), and MUMPS can fail with INFO(1) = -9 when a tiny mesh is split
over many ranks. Memory is summed over ranks; each rank adds ~75 MB of baseline.

## Cost model

`palace-tools estimate` fits `y = a (N/1e5)^b` per class to recorded runs (`palace_tools/data/cost_records.json`):

| class | memory, 1e5 unknowns | exponent b | notes |
|---|---|---|---|
| eigenmode, direct (MUMPS, MGMaxLevels 1), order 2 | ~1.3 GB | ~1.04 | 1.12 M unknowns -> ~16 GB |
| driven, direct, order 2 | ~2.2 GB | ~0.83 | complex factorisation, ~15-25 kB per unknown at 200-600k |
| driven, p-multigrid + MUMPS coarse, order 2 | ~1.2 GB | ~0.85 | ~8 kB per unknown at 1.3 M |

Time: driven = (full solves) x (time per solve) + (online evaluations of the adaptive sweep) x (time per
point); eigenmode = a N^b (modes/10)^c. The spread between a quiet and a busy machine is 2-5x, and the typical
value over-predicted the small examples by about 2x. Use the estimate to decide whether to ask the user first,
not as a promise.

Anchors from the records: birdcage eigenmode 1.12 M unknowns, order 2, 4 ranks: 16.4 GB, 11-19 min (20
modes). Wilkinson driven 1.32 M unknowns, p-multigrid, 8 ranks, 161 points x 3 excitations adaptive: 10 GB,
84 min. Port-tuning 7-port, 258k unknowns, direct, 317 points x 7 excitations: 7 GB, 128 min (97 min of it in
the adaptive sweep's online phase).
