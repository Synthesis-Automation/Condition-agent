from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

from chemtools.taxonomy.compound_catalog import build_documented_compound_catalog
from chemtools.taxonomy import loader as taxonomy_loader


def _norm_id(entry: Dict[str, Any]) -> str:
    cid = str(entry.get("id") or "").strip()
    if cid:
        return cid
    a = str(entry.get("A") or "").strip()
    b = str(entry.get("B") or "").strip()
    if not a or not b:
        return ""
    return f"{a}{'' if b.startswith('-') else '-'}{b}"


def _norm_template(entry: Dict[str, Any]) -> str:
    template = str(entry.get("template") or "").strip()
    if template:
        return template
    if entry.get("smarts") or entry.get("smarts_any"):
        return ""
    return "single_bond"


def _canonical_entry(entry: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "id": _norm_id(entry),
        "A": entry.get("A"),
        "B": entry.get("B"),
        "template": _norm_template(entry),
        "smarts": entry.get("smarts"),
        "smarts_any": entry.get("smarts_any"),
        "description": entry.get("description"),
        "priority": entry.get("priority"),
        "reactivity_weight": entry.get("reactivity_weight"),
    }


def _canonical_map(payload: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    for raw in payload.get("compounds", []) or []:
        if not isinstance(raw, dict):
            continue
        item = _canonical_entry(raw)
        cid = item["id"]
        assert cid, f"compound entry missing canonical id: {raw}"
        out[cid] = item
    return out


def test_generated_compound_catalog_matches_legacy_snapshot() -> None:
    generated = build_documented_compound_catalog()
    assert isinstance(generated, dict)
    assert generated.get("compounds"), "generated compound catalog is empty"

    # Loader should expose the same canonical generated payload.
    loaded = taxonomy_loader.load_organic_compounds()
    assert _canonical_map(loaded) == _canonical_map(generated)

    legacy_path = (
        Path(__file__).resolve().parents[1]
        / "chemtools"
        / "taxonomy"
        / "data"
        / "organic_compounds.v1.3.json"
    )
    if not legacy_path.exists():
        return
    with legacy_path.open("r", encoding="utf-8") as handle:
        legacy = json.load(handle)

    gen_map = _canonical_map(generated)
    legacy_map = _canonical_map(legacy)

    assert set(gen_map) == set(legacy_map)
    assert gen_map == legacy_map
