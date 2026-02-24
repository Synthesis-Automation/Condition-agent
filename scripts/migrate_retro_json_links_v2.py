#!/usr/bin/env python3
"""
Prepare v2 retro JSON payloads with explicit canonical taxonomy-family links.

v1 retro JSONs (`retron_patterns.json`, `hte_templates.json`) use `taxonomy_id`
as a link field, but in practice it has mixed semantics (canonical family IDs,
aliases, and subtype labels). This script generates v2-style previews that keep
backward-compatible fields while adding explicit canonical link fields:

- retrons: `taxonomy_family_id`, `retro_transform_id`
- templates: `taxonomy_family_id`

By default, writes preview files next to the originals with `.v2.preview.json`.
Use `--in-place` only when you are ready to migrate the source schema.
"""

from __future__ import annotations

import argparse
import json
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.taxonomy import loader as taxonomy_loader
from chemtools.taxonomy.reaction_catalog import get_reaction_type, resolve_reaction_type


AMBIGUOUS_RETRO_OVERRIDES: Dict[str, str] = {
    # Explicit chemistry-aware overrides (avoid ambiguous alias resolution like "SNAr")
    "aryl_fluorine_snarf": "C_N_Coupling",
    "snar_amination": "C_N_Coupling",
    "snar_co": "C_O_Coupling",
    "cuaac_triazole": "Click_azide_alkyne_cycloaddition",
}


def _canonical_taxonomy_id(raw_label: Any, *, entry_name: str = "") -> Tuple[Optional[str], str]:
    text = str(raw_label or "").strip()
    if not text:
        return None, "missing"

    direct = get_reaction_type(text)
    if direct:
        # `get_reaction_type` resolves aliases too; classify the path for reporting.
        if direct.id == text:
            return direct.id, "canonical"

    # Entry-specific overrides next (prevents ambiguous alias collapse, e.g. "SNAr").
    if entry_name and entry_name in AMBIGUOUS_RETRO_OVERRIDES:
        return AMBIGUOUS_RETRO_OVERRIDES[entry_name], "override"

    if direct:
        return direct.id, "alias"

    resolved = resolve_reaction_type(text)
    if resolved:
        return resolved, "alias"
    return None, "orphan"


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _migrate_retrons(payload: Mapping[str, Any]) -> Tuple[Dict[str, Any], Dict[str, int]]:
    out = deepcopy(dict(payload))
    entries = out.get("retrons")
    counts: Dict[str, int] = {"total": 0, "canonical": 0, "alias": 0, "override": 0, "orphan": 0, "missing": 0}
    if not isinstance(entries, list):
        return out, counts

    for entry in entries:
        if not isinstance(entry, dict):
            continue
        counts["total"] += 1
        name = str(entry.get("name") or "").strip()
        canonical, status = _canonical_taxonomy_id(entry.get("taxonomy_id"), entry_name=name)
        counts[status] = counts.get(status, 0) + 1
        if canonical:
            entry["taxonomy_family_id"] = canonical
        if entry.get("reaction_name") and not entry.get("retro_transform_id"):
            entry["retro_transform_id"] = str(entry.get("reaction_name"))

    out.setdefault("schema", {})
    if isinstance(out["schema"], dict):
        out["schema"]["retro_link_version"] = "v2-preview"
    return out, counts


def _migrate_templates(payload: Mapping[str, Any]) -> Tuple[Dict[str, Any], Dict[str, int]]:
    out = deepcopy(dict(payload))
    entries = out.get("templates")
    counts: Dict[str, int] = {"total": 0, "canonical": 0, "alias": 0, "override": 0, "orphan": 0, "missing": 0}
    if not isinstance(entries, list):
        return out, counts

    for entry in entries:
        if not isinstance(entry, dict):
            continue
        counts["total"] += 1
        name = str(entry.get("name") or "").strip()
        canonical, status = _canonical_taxonomy_id(entry.get("taxonomy_id"), entry_name=name)
        counts[status] = counts.get(status, 0) + 1
        if canonical:
            entry["taxonomy_family_id"] = canonical

    out.setdefault("schema", {})
    if isinstance(out["schema"], dict):
        out["schema"]["retro_link_version"] = "v2-preview"
    return out, counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Overwrite source JSONs (otherwise writes *.v2.preview.json files).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional directory for preview files (ignored with --in-place).",
    )
    args = parser.parse_args()

    retron_file = taxonomy_loader.RETRON_PATTERNS_FILE
    hte_file = taxonomy_loader.HTE_TEMPLATES_FILE
    retron_payload = taxonomy_loader.load_retron_patterns()
    hte_payload = taxonomy_loader.load_hte_templates()

    retron_v2, retron_counts = _migrate_retrons(retron_payload)
    hte_v2, hte_counts = _migrate_templates(hte_payload)

    if args.in_place:
        _write_json(retron_file, retron_v2)
        _write_json(hte_file, hte_v2)
        print(f"Wrote in-place: {retron_file}")
        print(f"Wrote in-place: {hte_file}")
    else:
        out_dir = args.output_dir
        if out_dir:
            out_dir.mkdir(parents=True, exist_ok=True)
            retron_out = out_dir / retron_file.name.replace(".json", ".v2.preview.json")
            hte_out = out_dir / hte_file.name.replace(".json", ".v2.preview.json")
        else:
            retron_out = retron_file.with_name(retron_file.stem + ".v2.preview.json")
            hte_out = hte_file.with_name(hte_file.stem + ".v2.preview.json")
        _write_json(retron_out, retron_v2)
        _write_json(hte_out, hte_v2)
        print(f"Wrote preview: {retron_out}")
        print(f"Wrote preview: {hte_out}")

    print("\nMigration summary")
    print("-" * 18)
    print(f"Retrons   : {retron_counts}")
    print(f"Templates : {hte_counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
