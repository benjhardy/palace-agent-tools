"""Command line: one JSON object in, one JSON object out.

    palace-tools <tool> '<json>'        arguments inline
    palace-tools <tool> @args.json      arguments from a file
    palace-tools <tool> -               arguments from stdin
    palace-tools list                   tool names and descriptions

Tool names are the MCP tool names without the 'palace_' prefix, with '-' or '_' (e.g. lint-config).
Exit code 0 on success, 1 when the tool reports ok = false or overall = fail, 2 on bad input.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from . import api


def _parse_args(arg: str | None) -> dict:
    if arg is None:
        return {}
    if arg == "-":
        return json.loads(sys.stdin.read() or "{}")
    if arg.startswith("@"):
        return json.loads(Path(arg[1:]).read_text(encoding="utf-8"))
    return json.loads(arg)


def main(argv: list[str] | None = None) -> None:
    argv = sys.argv[1:] if argv is None else argv
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except AttributeError:  # pragma: no cover
        pass
    if not argv or argv[0] in ("-h", "--help", "help"):
        print(__doc__)
        sys.exit(0)
    if argv[0] == "list":
        print(json.dumps({n: t["description"] for n, t in api.TOOLS.items()}, indent=1))
        return
    name = "palace_" + argv[0].replace("-", "_")
    try:
        args = _parse_args(argv[1] if len(argv) > 1 else None)
        res = api.call(name, args)
    except (KeyError, ValueError, FileNotFoundError, json.JSONDecodeError) as e:
        print(json.dumps({"ok": False, "error": f"{type(e).__name__}: {e}"}))
        sys.exit(2)
    print(json.dumps(res, indent=1))
    bad = res.get("ok") is False or res.get("overall") == "fail"
    sys.exit(1 if bad else 0)
