"""Regenerate skills/palace-simulation/references/gotchas.md from the lint rules and the workflow lessons.

    python tools/gen_gotchas.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from palace_tools.api import OTHER_GOTCHAS  # noqa: E402
from palace_tools.lint import rules_table  # noqa: E402


def main() -> None:
    lines = ["# Palace gotchas (generated from palace_tools.lint; do not edit by hand)", "",
             "Each lint rule below is checked by `palace-tools lint-config`. Rules marked *geometry* need the "
             "`*.tags.json` written by `palace_tools.mesh`. Targets Palace 0.16; version differences are stated in "
             "the rule.", ""]
    for r in rules_table():
        lines.append(f"## {r['id']}: {r['title']}" + (" *(geometry)*" if r["needs"] == "geometry" else ""))
        lines.append("")
        lines.append(r["lesson"])
        lines.append("")
        lines.append(f"Source: {r['source']}." + (f" Docs: {r['docs']}" if r["docs"] else ""))
        lines.append("")
    lines += ["# Workflow lessons (not lint rules)", ""]
    for g in OTHER_GOTCHAS:
        lines += [f"## {g['title']}", "", g["lesson"], "", f"Source: {g['source']}.", ""]
    out = ROOT / "skills" / "palace-simulation" / "references" / "gotchas.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines), encoding="utf-8")
    print(out, len(rules_table()), "rules")


if __name__ == "__main__":
    main()
