"""Launch Palace (WSL + Apptainer, Apptainer, Docker or a native binary) with guards and cost logging.

Guards (refuse unless ``force``):
  * another Palace job is running and the requested ranks would oversubscribe the CPU
    (with a job running, use at most 2 ranks or wait);
  * the estimated peak memory does not fit in the available RAM.

Logging: Palace stdout -> ``palace.log``, ``/usr/bin/time -v`` -> ``time.log`` (wall time, max RSS of one rank),
and a ``run_record.json`` with the command, ranks, wall time, Palace's peak-memory estimate and parsed warnings.

WSL note: WSL stops the VM, killing ``nohup``/``&`` jobs, when the last ``wsl.exe`` exits. ``run`` blocks until
Palace finishes, so the calling process holds ``wsl.exe`` open. To run in the background, start this command
itself as a background job of your agent/terminal; never detach Palace inside WSL.

Configuration (arguments override environment variables):
  PALACE_TOOLS_BACKEND  wsl | apptainer | docker | native   (default: wsl on Windows, else apptainer if an image
                        is set, else native)
  PALACE_TOOLS_DISTRO   WSL distribution (default: the WSL default distribution)
  PALACE_TOOLS_IMAGE    Apptainer .sif path or Docker image
  PALACE_TOOLS_CMD      full command template, e.g. "run_palace {config} {np}" (overrides the backend command)
"""
from __future__ import annotations

import json
import os
import platform
import shlex
import shutil
import subprocess
import time
from pathlib import Path

from .logparse import parse_log, parse_time_v

PROC_PATTERN = r"palace-[A-Za-z0-9_]+\.bin"


def _wslpath(p: Path) -> str:
    p = Path(p).resolve()
    drive = p.drive.rstrip(":").lower()
    return f"/mnt/{drive}" + p.as_posix()[len(p.drive):]


def _backend(backend: str | None) -> str:
    b = backend or os.environ.get("PALACE_TOOLS_BACKEND")
    if b:
        return b
    if platform.system() == "Windows":
        return "wsl"
    return "apptainer" if os.environ.get("PALACE_TOOLS_IMAGE") else "native"


def _shell(cmd: str, backend: str, distro: str | None, timeout: float = 60) -> subprocess.CompletedProcess:
    """Run a short shell command where Palace runs (inside WSL for the wsl backend)."""
    if backend == "wsl":
        args = ["wsl"] + (["-d", distro] if distro else []) + ["--", "bash", "-lc", cmd]
    else:
        args = ["bash", "-lc", cmd]
    return subprocess.run(args, capture_output=True, text=True, timeout=timeout, encoding="utf-8", errors="replace")


def host_status(backend: str | None = None, distro: str | None = None) -> dict:
    """Running Palace processes (MPI ranks), CPU count and available RAM where Palace runs."""
    backend = _backend(backend)
    distro = distro or os.environ.get("PALACE_TOOLS_DISTRO")
    # no '$' in the command: wsl.exe re-parses its arguments and would expand them early
    probe = f"pgrep -af '{PROC_PATTERN}' | grep -v pgrep; echo '@@'; nproc; echo '@@'; free -b | grep '^Mem:'"
    try:
        r = _shell(probe, "wsl" if backend == "wsl" else "local", distro)
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        return {"backend": backend, "error": f"cannot reach the Palace host: {e}"}
    parts = r.stdout.split("@@")
    lines = [ln.strip() for ln in parts[0].splitlines() if ln.strip()] if parts else []
    # one line per MPI rank: the command itself is the Palace binary (skip mpirun / wrappers that mention it)
    procs = [ln for ln in lines if len(ln.split()) > 1 and ln.split()[1].endswith(".bin")]
    try:
        cpus = int(parts[1].strip())
        mem = parts[2].split()  # Mem: total used free shared buff/cache available
        tot, avail = int(mem[1]), int(mem[6])
    except (IndexError, ValueError):
        cpus, tot, avail = os.cpu_count() or 1, 0, 0
    cases = sorted({ln.split()[-1] for ln in procs})
    return {"backend": backend, "distro": distro, "palace_ranks_running": len(procs), "processes": procs[:16],
            "configs_running": cases, "cpus": cpus, "ram_total_gb": round(tot / 2 ** 30, 2),
            "ram_available_gb": round(avail / 2 ** 30, 2)}


def guard(ranks: int, est_mem_gb: float | None = None, status: dict | None = None, max_ranks_if_busy: int = 2) -> dict:
    """Decide whether a new run may start now. Returns {"ok", "reasons", "suggested_ranks", "status"}."""
    st = status or host_status()
    reasons = []
    busy = st.get("palace_ranks_running", 0)
    cpus = st.get("cpus") or 1
    sugg = ranks
    if st.get("error"):
        return {"ok": False, "reasons": [st["error"]], "suggested_ranks": ranks, "status": st}
    if busy:
        free_cpus = max(0, cpus - busy)
        sugg = max(1, min(ranks, max_ranks_if_busy, free_cpus or 1))
        if ranks > max_ranks_if_busy or ranks > free_cpus:
            reasons.append(f"another Palace job is using {busy} ranks on {cpus} CPUs: use <= {sugg} ranks or wait")
    elif ranks > cpus:
        reasons.append(f"{ranks} ranks > {cpus} CPUs")
        sugg = cpus
    avail = st.get("ram_available_gb") or 0
    if est_mem_gb and avail and est_mem_gb > 0.9 * avail:
        reasons.append(f"estimated peak memory {est_mem_gb:.1f} GB > 90 % of the {avail:.1f} GB available")
    return {"ok": not reasons, "reasons": reasons, "suggested_ranks": sugg, "status": st}


def palace_command(config: str, ranks: int, backend: str, image: str | None = None, cmd: str | None = None,
                   case_dir: Path | None = None) -> str:
    tpl = cmd or os.environ.get("PALACE_TOOLS_CMD")
    image = image or os.environ.get("PALACE_TOOLS_IMAGE")
    if tpl:
        return tpl.format(config=shlex.quote(config), np=ranks)
    if backend in ("wsl", "apptainer"):
        if image:
            return f"apptainer run {shlex.quote(image)} -np {ranks} {shlex.quote(config)}"
        # Palace's own wrapper (`palace -np N config.json`) or a site wrapper such as run_palace
        return (f"if command -v run_palace >/dev/null; then run_palace {shlex.quote(config)} {ranks}; "
                f"else palace -np {ranks} {shlex.quote(config)}; fi")
    if backend == "docker":
        if not image:
            raise ValueError("docker backend needs an image (PALACE_TOOLS_IMAGE)")
        d = str(Path(case_dir).resolve()) if case_dir else "$PWD"
        return (f"docker run --rm -v {shlex.quote(d)}:/work -w /work {shlex.quote(image)} "
                f"palace -np {ranks} {shlex.quote(config)}")
    return f"palace -np {ranks} {shlex.quote(config)}"


def run_case(case_dir: Path, config: str = "config.json", ranks: int = 4, backend: str | None = None,
             distro: str | None = None, image: str | None = None, cmd: str | None = None, force: bool = False,
             est_mem_gb: float | None = None, timeout_s: float | None = None, dry_run: bool = False) -> dict:
    """Run Palace on case_dir/config. Blocks until Palace exits. Returns the run record."""
    case_dir = Path(case_dir)
    if not (case_dir / config).exists():
        return {"ok": False, "error": f"{config} not found in {case_dir}"}
    backend = _backend(backend)
    distro = distro or os.environ.get("PALACE_TOOLS_DISTRO")
    pc = palace_command(config, ranks, backend, image, cmd, case_dir)
    timed = f"( command -v /usr/bin/time >/dev/null && /usr/bin/time -v bash -c {shlex.quote(pc)} " \
            f"|| bash -c {shlex.quote(pc)} ) > palace.log 2> time.log"
    if backend == "wsl":
        inner = f"cd {shlex.quote(_wslpath(case_dir))} && {timed}"
        args = ["wsl"] + (["-d", distro] if distro else []) + ["--", "bash", "-lc", inner]
    else:
        inner = f"cd {shlex.quote(str(case_dir.resolve()))} && {timed}"
        args = ["bash", "-lc", inner]
    rec = {"case": str(case_dir), "config": config, "ranks": ranks, "backend": backend, "command": pc}
    if dry_run:
        rec.update({"ok": True, "dry_run": True, "argv": args})
        return rec
    if not force:
        g = guard(ranks, est_mem_gb, host_status(backend, distro))
        rec["guard"] = g
        if not g["ok"]:
            rec.update({"ok": False, "refused": True,
                        "error": "; ".join(g["reasons"]) + " (pass force=true to override)"})
            return rec
    t0 = time.time()
    rec["started"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    try:
        p = subprocess.run(args, timeout=timeout_s)
        rec["exit_code"] = p.returncode
    except subprocess.TimeoutExpired:
        rec["exit_code"] = None
        rec["timed_out"] = True
    rec["wall_s"] = round(time.time() - t0, 1)
    log = case_dir / "palace.log"
    tl = case_dir / "time.log"
    lg = parse_log(log.read_text(encoding="utf-8", errors="replace")) if log.exists() else {}
    tv = parse_time_v(tl.read_text(encoding="utf-8", errors="replace")) if tl.exists() else {}
    rec["time"] = tv
    rec["summary"] = {k: lg.get(k) for k in ("palace_version", "problem", "order", "unknowns", "ranks",
                                              "peak_mem_gb", "elapsed_s", "finished", "n_errors", "n_warnings")}
    rec["messages"] = lg.get("messages", [])
    rec["ok"] = rec.get("exit_code") == 0 and bool(lg.get("finished")) and not lg.get("n_errors")
    (case_dir / "run_record.json").write_text(json.dumps(rec, indent=1), encoding="utf-8")
    return rec


def progress(case_dir: Path) -> dict:
    """Where a running (or finished) case is: driven step k/N, adaptive samples, AMR pass, elapsed log time."""
    case_dir = Path(case_dir)
    lp = case_dir / "palace.log"
    if not lp.exists():
        return {"case": str(case_dir), "state": "not started"}
    lg = parse_log(lp.read_text(encoding="utf-8", errors="replace"))
    age = time.time() - lp.stat().st_mtime
    st = "finished" if lg.get("finished") else ("failed" if lg.get("n_errors") else "running")
    if st == "running" and age > 3600:
        st = "stalled or killed (log not written for over an hour)"
    return {"case": str(case_dir), "state": st, "log_age_s": round(age), "unknowns": lg.get("unknowns"),
            "progress": lg.get("progress"), "adaptive": lg.get("adaptive"),
            "greedy_errors": lg.get("adaptive_errors", [])[-3:], "amr": lg.get("amr"),
            "messages": lg.get("messages")}


def which_palace() -> dict:
    """Quick environment probe for the README/skill 'doctor' step."""
    out = {"platform": platform.system(), "backend": _backend(None),
           "wsl": bool(shutil.which("wsl")), "docker": bool(shutil.which("docker")),
           "apptainer": bool(shutil.which("apptainer")), "palace": bool(shutil.which("palace"))}
    if out["backend"] == "wsl" and out["wsl"]:
        try:
            r = _shell("command -v run_palace palace apptainer; ls /opt/palace 2>/dev/null | head -5", "wsl",
                       os.environ.get("PALACE_TOOLS_DISTRO"), timeout=60)
            out["inside_wsl"] = r.stdout.split()
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            out["inside_wsl_error"] = str(e)
    return out
