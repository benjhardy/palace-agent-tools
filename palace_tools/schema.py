"""Validate a config against the Palace 0.16 JSON schema (vendored from awslabs/palace, Apache-2.0)."""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

SCHEMA_DIR = Path(__file__).resolve().parent / "data" / "schema_0.16"


@lru_cache(maxsize=1)
def _load():
    root = json.loads((SCHEMA_DIR / "config-schema.json").read_text(encoding="utf-8"))
    subs = {}
    for p in (SCHEMA_DIR / "config").glob("*.json"):
        subs[f"./config/{p.name}"] = json.loads(p.read_text(encoding="utf-8"))
    return root, subs


def _validator():
    import jsonschema
    from referencing import Registry, Resource
    from referencing.jsonschema import DRAFT7

    root, subs = _load()
    base = "https://palace.invalid/schema/"
    root = dict(root, **{"$id": base + "config-schema.json"})
    reg = Registry()
    for rel, sch in subs.items():
        uri = base + rel[2:]
        reg = reg.with_resource(uri, Resource.from_contents(dict(sch, **{"$id": uri}), default_specification=DRAFT7))
    return jsonschema.Draft7Validator(root, registry=reg)


def validate(cfg: dict) -> dict:
    """{"valid": bool, "errors": [{"path", "message"}], "validator": "jsonschema"|"unavailable"}"""
    try:
        import jsonschema
    except ImportError:
        return {"valid": None, "errors": [], "validator": "unavailable",
                "note": "pip install jsonschema to validate against the Palace schema"}
    v = _validator()
    errs = []
    for e in sorted(v.iter_errors(cfg), key=lambda e: list(e.absolute_path)):
        path = "/".join(str(x) for x in e.absolute_path) or "(root)"
        msg = e.message
        if len(msg) > 300:
            msg = msg[:300] + "..."
        errs.append({"path": path, "message": msg})
    return {"valid": not errs, "errors": errs[:50], "validator": "jsonschema", "schema": "Palace 0.16"}
