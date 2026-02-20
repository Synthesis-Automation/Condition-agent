"""
Build knowledge_base/notes/_index.json from all note subfolders.

Run after adding new notes or when the index is stale:
    python knowledge_base/notes/build_index.py

Or call build_index() from Python:
    from knowledge_base.notes.build_index import build_index
    build_index()
"""
from __future__ import annotations

import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Set


# Note type subfolders to scan
NOTE_SUBDIRS = ["reactions", "mechanisms", "substrates", "protocols", "routes"]


# ---------------------------------------------------------------------------
# Minimal YAML-like front matter parser (no pyyaml dependency)
# ---------------------------------------------------------------------------

def _parse_front_matter(text: str) -> Dict[str, Any]:
    """
    Parse YAML front matter between opening and closing --- lines.

    Supports:
      key: scalar value
      key: [item1, item2]   (inline list)
    """
    lines = text.splitlines()
    if not lines or lines[0].strip() != "---":
        return {}

    fm_lines = []
    for line in lines[1:]:
        if line.strip() == "---":
            break
        fm_lines.append(line)

    result: Dict[str, Any] = {}
    for line in fm_lines:
        m = re.match(r"^(\w[\w_]*)\s*:\s*(.*)$", line.strip())
        if not m:
            continue
        key, val = m.group(1), m.group(2).strip()

        # Inline list: [a, b, c]
        if val.startswith("[") and val.endswith("]"):
            inner = val[1:-1]
            items = [s.strip().strip("'\"") for s in inner.split(",")]
            result[key] = [i for i in items if i]
        else:
            # Scalar — strip optional quotes
            result[key] = val.strip("'\"")

    return result


# ---------------------------------------------------------------------------
# Index builder
# ---------------------------------------------------------------------------

def build_index(notes_root: Optional[Path] = None) -> Path:
    """
    Scan all note subfolders, parse front matter, and write _index.json.

    Returns the path to the written index file.
    """
    if notes_root is None:
        # Resolve relative to THIS file's location (notes/)
        notes_root = Path(__file__).parent

    entries: List[Dict[str, Any]] = []
    tag_index: Dict[str, List[str]] = {}
    bond_index: Dict[str, List[str]] = {}
    metal_index: Dict[str, List[str]] = {}
    type_index: Dict[str, List[str]] = {}

    def _register(key: str, value_or_list: Any, index: Dict[str, List[str]]) -> None:
        vals = value_or_list if isinstance(value_or_list, list) else [value_or_list]
        for v in vals:
            v = str(v).strip()
            if not v:
                continue
            index.setdefault(v, [])
            if key not in index[v]:
                index[v].append(key)

    for subdir_name in NOTE_SUBDIRS:
        subdir = notes_root / subdir_name
        if not subdir.exists():
            continue

        for md_file in sorted(subdir.glob("*.md")):
            if md_file.name.startswith(("README", "_", ".")):
                continue

            try:
                text = md_file.read_text(encoding="utf-8")
            except Exception:
                continue

            fm = _parse_front_matter(text)
            note_id = fm.get("id") or md_file.stem
            note_type = fm.get("type") or subdir_name.rstrip("s")  # reactions→reaction
            title = fm.get("title") or note_id.replace("_", " ").title()
            aliases = fm.get("aliases", [])
            tags = fm.get("tags", [])
            bond_formed = fm.get("bond_formed", [])
            metal = fm.get("metal", [])
            related = (
                fm.get("related_reactions", [])
                + fm.get("related_mechanisms", [])
                + fm.get("related_substrates", [])
                + fm.get("related", [])
            )

            entry: Dict[str, Any] = {
                "id": note_id,
                "type": note_type,
                "title": title,
                "file": f"{subdir_name}/{md_file.name}",
            }
            if aliases:
                entry["aliases"] = aliases
            if tags:
                entry["tags"] = tags
            if bond_formed:
                entry["bond_formed"] = bond_formed
            if metal:
                entry["metal"] = metal
            if related:
                entry["related"] = related

            # Add optional fields if present (includes route-specific fields)
            for extra in ("mechanism", "substrates", "applies_to", "used_in",
                          "variants", "applies_to_reactions", "bond_broken",
                          "target_smiles", "steps", "overall_yield"):
                if extra in fm:
                    entry[extra] = fm[extra]

            entries.append(entry)

            # Build inverted indices
            _register(note_id, tags, tag_index)
            _register(note_id, bond_formed, bond_index)
            _register(note_id, metal, metal_index)
            type_index.setdefault(note_type, [])
            if note_id not in type_index[note_type]:
                type_index[note_type].append(note_id)

    index = {
        "version": 1,
        "last_built": datetime.now(timezone.utc).isoformat(),
        "entries": entries,
        "tag_index": tag_index,
        "bond_index": bond_index,
        "metal_index": metal_index,
        "type_index": type_index,
    }

    out_path = notes_root / "_index.json"
    out_path.write_text(json.dumps(index, indent=2, ensure_ascii=False), encoding="utf-8")
    return out_path


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    path = build_index(root)
    entries = json.loads(path.read_text())["entries"]
    print(f"Built {path} — {len(entries)} entries")
    for e in entries:
        print(f"  [{e['type']:12s}] {e['id']}")
