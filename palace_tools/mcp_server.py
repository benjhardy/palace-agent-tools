"""Minimal MCP server (JSON-RPC 2.0 over stdio, newline-delimited) exposing palace_tools as tools.

No dependency on an MCP SDK. Register it with Claude Code:

    claude mcp add palace-tools -- python -m palace_tools.mcp_server

or in ``.mcp.json``: {"mcpServers": {"palace-tools": {"command": "python", "args": ["-m", "palace_tools.mcp_server"]}}}

Anything a tool prints (Gmsh writes to the C-level stdout) is redirected to stderr, so the protocol stream on
stdout stays clean. Tools run synchronously; ``palace_run`` blocks until Palace finishes.
"""
from __future__ import annotations

import json
import os
import sys
import traceback

from . import __version__, api

SUPPORTED = ["2025-06-18", "2025-03-26", "2024-11-05"]


def tool_list() -> list[dict]:
    return [{"name": n, "description": t["description"], "inputSchema": t["schema"]} for n, t in api.TOOLS.items()]


def handle(msg: dict) -> dict | None:
    mid = msg.get("id")
    method = msg.get("method")
    params = msg.get("params") or {}
    if mid is None:  # notification (e.g. notifications/initialized): no reply
        return None
    try:
        if method == "initialize":
            want = params.get("protocolVersion")
            res = {"protocolVersion": want if want in SUPPORTED else SUPPORTED[0],
                   "capabilities": {"tools": {"listChanged": False}},
                   "serverInfo": {"name": "palace-agent-tools", "version": __version__},
                   "instructions": "Tools for AWS Palace (FEM EM solver) runs. Workflow: palace_mesh_template or your "
                                   "own Gmsh script -> palace_build_config -> palace_lint_config -> "
                                   "palace_estimate_cost -> palace_status -> palace_run -> palace_check -> "
                                   "palace_results. Read palace_gotchas once."}
        elif method == "ping":
            res = {}
        elif method == "tools/list":
            res = {"tools": tool_list()}
        elif method == "tools/call":
            name = params.get("name")
            try:
                out = api.call(name, params.get("arguments") or {})
                bad = out.get("ok") is False or out.get("overall") == "fail"
                res = {"content": [{"type": "text", "text": json.dumps(out, indent=1)}], "isError": False,
                       "structuredContent": out}
                if bad:
                    res["content"].insert(0, {"type": "text", "text": "Result reports a problem (ok=false or "
                                                                      "overall=fail); see the findings."})
            except Exception as e:  # noqa: BLE001  tool errors are returned to the model, not raised
                res = {"content": [{"type": "text", "text": f"{type(e).__name__}: {e}\n{traceback.format_exc(limit=3)}"}],
                       "isError": True}
        else:
            return {"jsonrpc": "2.0", "id": mid, "error": {"code": -32601, "message": f"method not found: {method}"}}
        return {"jsonrpc": "2.0", "id": mid, "result": res}
    except Exception as e:  # noqa: BLE001
        return {"jsonrpc": "2.0", "id": mid, "error": {"code": -32603, "message": str(e)}}


def main() -> None:
    out = os.fdopen(os.dup(1), "w", encoding="utf-8", newline="\n")
    os.dup2(2, 1)  # stray prints (Python or C) go to stderr
    sys.stdout = sys.stderr
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            msg = json.loads(line)
        except json.JSONDecodeError:
            out.write(json.dumps({"jsonrpc": "2.0", "id": None, "error": {"code": -32700, "message": "parse error"}})
                      + "\n")
            out.flush()
            continue
        reply = handle(msg)
        if reply is not None:
            out.write(json.dumps(reply) + "\n")
            out.flush()


if __name__ == "__main__":
    main()
