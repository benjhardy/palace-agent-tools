import json
import subprocess
import sys

import pytest

from palace_tools import api
from palace_tools.lint import lint


def test_cavity_template_counts(tmp_path):
    pytest.importorskip("gmsh")
    from palace_tools import mesh

    info = mesh.rectangular_cavity(tmp_path / "cav", 22.86, 10.16, 30.0, 4.0, verbosity=0)
    m = info["mesh"]
    assert m["tets"] > 100 and m["verdict"] == "ok"
    assert m["unknowns"]["2"] == 2 * m["edges"] + 2 * m["faces"]
    assert info["surfaces"]["walls"]["outer"] and set(info["surfaces"]["walls"]["sides"]) == set(mesh.SIDES)
    assert (tmp_path / "cav.msh").exists() and (tmp_path / "cav.tags.json").exists()


def test_microstrip_template_geometry_and_lint(tmp_path):
    pytest.importorskip("gmsh")
    from palace_tools import mesh
    from palace_tools.config import build_config

    info = mesh.microstrip_line(tmp_path / "ms", 1.1, 0.508, 10.0, 3.55, side=3, top=3, h_mesh=2.0, h_edge=0.3,
                                verbosity=0)
    s = info["surfaces"]
    assert {"strip", "port1", "port2", "ends", "sides", "ground", "top"} <= set(s)
    assert s["port1"]["sides"] == ["-x"] and s["strip"]["outer"] is False
    assert s["port1"]["lengths"][0] == 0.0  # flat in x
    assert info["structure_bbox"][5] == pytest.approx(0.508, abs=1e-5)
    spec = {"problem": "driven", "mesh": "ms.msh", "units": "mm", "tags": str(tmp_path / "ms.tags.json"),
            "materials": [{"tags": ["substrate"], "eps_r": 3.55}, {"tags": ["air"]}],
            "pec": ["ground", "strip"], "pmc": ["ends"], "absorbing": {"tags": ["sides", "top"]},
            "lumped_ports": [{"index": 1, "tags": ["port1"], "direction": "+Z"},
                             {"index": 2, "tags": ["port2"], "direction": "+Z"}],
            "frequency": {"min_ghz": 1, "max_ghz": 20, "step_ghz": 1}}
    cfg = build_config(spec)
    r = lint(cfg, info)
    assert r["counts"]["error"] == 0
    assert "BND001" in {f["id"] for f in r["findings"]}  # 3 mm walls at 20 GHz (0.2 lambda) -> caution
    cfg["Boundaries"]["LumpedPort"][0]["Direction"] = "+X"
    assert "PRT009" in {f["id"] for f in lint(cfg, info)["findings"] if f["severity"] == "error"}


def test_api_build_lint_estimate_roundtrip(tmp_path):
    spec = dict(api.EXAMPLE_SPEC, tags={"substrate": 1, "air": 2, "ground": 101, "strip": 102, "port1": 103,
                                        "port2": 104})
    r = api.call("palace_build_config", {"spec": spec, "out": str(tmp_path / "config.json")})
    assert r["schema"]["valid"]
    assert {f["id"] for f in r["lint"]["findings"] if f["severity"] == "error"} == {"MOD002"}  # no mesh written yet
    r2 = api.call("palace_lint_config", {"config": str(tmp_path / "config.json")})
    assert any(f["id"] == "MOD002" for f in r2["findings"])  # the mesh file does not exist here
    e = api.call("palace_estimate_cost", {"config": str(tmp_path / "config.json"), "unknowns": 80000})
    assert e["memory_gb"]["typ"] > 0
    assert api.call("palace_gotchas", {})["gotchas"]


def _mcp(lines):
    p = subprocess.run([sys.executable, "-m", "palace_tools.mcp_server"], input="\n".join(map(json.dumps, lines)) + "\n",
                       capture_output=True, text=True, timeout=120)
    return [json.loads(x) for x in p.stdout.splitlines() if x.strip()]


def test_mcp_server_protocol():
    out = _mcp([{"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {"protocolVersion": "2025-06-18"}},
                {"jsonrpc": "2.0", "method": "notifications/initialized"},
                {"jsonrpc": "2.0", "id": 2, "method": "tools/list"},
                {"jsonrpc": "2.0", "id": 3, "method": "tools/call",
                 "params": {"name": "palace_estimate_cost", "arguments": {"unknowns": 500000, "problem": "Eigenmode"}}},
                {"jsonrpc": "2.0", "id": 4, "method": "tools/call", "params": {"name": "nope", "arguments": {}}},
                {"jsonrpc": "2.0", "id": 5, "method": "bogus"}])
    by = {m["id"]: m for m in out}
    assert by[1]["result"]["protocolVersion"] == "2025-06-18"
    names = {t["name"] for t in by[2]["result"]["tools"]}
    assert {"palace_build_config", "palace_lint_config", "palace_run", "palace_check", "palace_results"} <= names
    assert all(t["inputSchema"]["type"] == "object" for t in by[2]["result"]["tools"])
    assert by[3]["result"]["structuredContent"]["memory_gb"]["typ"] > 0
    assert by[4]["result"]["isError"] is True
    assert by[5]["error"]["code"] == -32601
    assert len(out) == 5  # no reply to the notification


def test_cli_json_io(tmp_path, data):
    p = subprocess.run([sys.executable, "-m", "palace_tools", "lint-config",
                        json.dumps({"config": str(data / "pp_config.json")})], capture_output=True, text=True)
    out = json.loads(p.stdout)
    assert out["schema"]["valid"] is True
    # the fixture's mesh is not shipped, so the only error is MOD002 and the exit code is 1
    assert {f["id"] for f in out["findings"] if f["severity"] == "error"} == {"MOD002"} and p.returncode == 1
    p = subprocess.run([sys.executable, "-m", "palace_tools", "lint-config", "{bad json"], capture_output=True,
                       text=True)
    assert p.returncode == 2 and json.loads(p.stdout)["ok"] is False
    p = subprocess.run([sys.executable, "-m", "palace_tools", "run",
                        json.dumps({"case_dir": str(tmp_path), "dry_run": True})], capture_output=True, text=True)
    assert json.loads(p.stdout)["ok"] is False  # no config.json there
