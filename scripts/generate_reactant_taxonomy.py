"""
Legacy: reactant_types.json is retired.

Reactant type SMARTS are now derived from organic_compounds.v1.3.json via
chemtools.taxonomy.calculable_spec. This script is kept for reference only.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional


ROOT = Path(__file__).resolve().parents[1]
SPEC_PATH = ROOT / "chemtools" / "taxonomy" / "data" / "calculable_features.json"
OUTPUT_PATH = ROOT / "chemtools" / "taxonomy" / "data" / "reactant_types.json"

if not OUTPUT_PATH.exists():
    raise SystemExit(
        "reactant_types.json has been removed; use organic_compounds.v1.3.json instead."
    )


def _load_existing() -> Dict[str, dict]:
    if not OUTPUT_PATH.exists():
        return {}
    existing = json.loads(OUTPUT_PATH.read_text(encoding="utf-8"))
    return {entry["id"]: entry for entry in existing if isinstance(entry, dict) and "id" in entry}


def _member_metadata(meta: dict, existing: Optional[dict]) -> dict:
    """Merge existing metadata with calculable reactant metadata."""
    merged = dict(existing.get("metadata", {})) if existing else {}
    for key in ("compound_type", "coupling_role", "group"):
        if key in meta and key not in merged:
            merged[key] = meta[key]
    return merged


def _member_name(member_id: str, meta: dict, existing: Optional[dict]) -> str:
    if existing and existing.get("name"):
        return existing["name"]
    return meta.get("member_name") or meta.get("compound_type") or member_id


def build_reactant_taxonomy() -> List[dict]:
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8"))
    existing = _load_existing()

    categories: Dict[str, dict] = {}

    for feature in spec.get("features", []):
        meta = feature.get("reactant_metadata")
        if not meta:
            continue

        cat_id = meta.get("reactant_category")
        member_id = meta.get("reactant_member")
        if not cat_id or not member_id:
            continue

        cat_entry = categories.get(cat_id)
        if cat_entry is None:
            existing_cat = existing.get(cat_id, {})
            cat_entry = {
                "id": cat_id,
                "name": existing_cat.get("name") or cat_id,
                "description": existing_cat.get("description"),
                "category": existing_cat.get("category"),
                "aliases": list(existing_cat.get("aliases", [])),
                "metadata": dict(existing_cat.get("metadata", {})),
                "members": [],
            }
            categories[cat_id] = cat_entry

        members_map = {member["id"]: member for member in cat_entry["members"] if "id" in member}
        if member_id in members_map:
            continue

        existing_member: Optional[dict] = None
        if cat_id in existing:
            for member in existing[cat_id].get("members", []):
                if member.get("id") == member_id:
                    existing_member = member
                    break

        cat_entry["members"].append(
            {
                "id": member_id,
                "name": _member_name(member_id, meta, existing_member),
                "aliases": list(existing_member.get("aliases", [])) if existing_member else [],
                "metadata": _member_metadata(meta, existing_member),
            }
        )

    # Deterministic ordering: categories by id, members by id
    ordered_categories: List[dict] = []
    for cat_id in sorted(categories):
        entry = categories[cat_id]
        entry["members"] = sorted(entry["members"], key=lambda m: m.get("id", ""))
        ordered_categories.append(entry)

    return ordered_categories


def main() -> int:
    taxonomy = build_reactant_taxonomy()
    OUTPUT_PATH.write_text(json.dumps(taxonomy, indent=2), encoding="utf-8")
    print(f"Wrote {len(taxonomy)} reactant type categories to {OUTPUT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
