from __future__ import annotations

import json
from pathlib import Path

import pytest

from chemtools.featurizers.motifs.registry import _default_registry_paths, build_compound_registry
from chemtools.featurizers.nearby_groups import _resolve_group_id
from chemtools.taxonomy import loader as taxonomy_loader
from chemtools.taxonomy.substituent_composer import (
    compose_groups_from_fragments,
    load_organic_groups_with_compositions,
    validate_substituent_fragments_payload,
)
from chemtools.util.rdkit_helpers import rdkit_available
from chemtools.util.smarts_cache import compile_smarts


def _taxonomy_data_dir() -> Path:
    return Path(__file__).resolve().parents[1] / "chemtools" / "taxonomy" / "data"


def _load_compounds() -> list[dict]:
    payload = taxonomy_loader.load_organic_compounds()
    return [entry for entry in (payload.get("compounds") or []) if isinstance(entry, dict)]


def _load_compound_logic_sets() -> dict[str, set[str]]:
    data_dir = _taxonomy_data_dir()
    logic_path = data_dir / "compound_logic.json"
    payload = json.loads(logic_path.read_text(encoding="utf-8"))
    motif_sets = payload.get("motif_sets", {}) or {}
    out: dict[str, set[str]] = {}
    for set_name, entry in motif_sets.items():
        if not isinstance(entry, dict):
            continue
        members = entry.get("members") or []
        out[str(set_name)] = {str(m) for m in members if str(m).strip()}
    return out


def _load_raw_groups() -> list[dict]:
    data_dir = _taxonomy_data_dir()
    groups_path = data_dir / "organic_groups.v1.3.json"
    payload = json.loads(groups_path.read_text(encoding="utf-8"))
    return [entry for entry in (payload.get("groups") or []) if isinstance(entry, dict)]


def _load_substituent_fragments() -> dict:
    data_dir = _taxonomy_data_dir()
    path = data_dir / "substituent_fragments.v1.json"
    return json.loads(path.read_text(encoding="utf-8"))


def test_composed_groups_are_merged_without_duplicates() -> None:
    data_dir = _taxonomy_data_dir()
    groups_path = data_dir / "organic_groups.v1.3.json"
    payload = load_organic_groups_with_compositions(groups_path)
    groups = payload.get("groups", []) or []
    ids = [str(entry.get("id") or "") for entry in groups if isinstance(entry, dict)]

    assert "-CON3" in ids
    assert "-CONHNH2" in ids
    assert "-COOH" in ids
    assert "-COOR" in ids
    assert "-OCONH2" in ids
    assert "-SO2NHNH2" in ids
    assert "-PO2OH" in ids
    assert "-PO2OR" in ids
    assert "-PO2NH2" in ids

    assert ids.count("-CONH2") == 1
    assert ids.count("-COOH") == 1
    assert ids.count("-COOR") == 1
    assert ids.count("-SO2NH2") == 1


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_composed_groups_smarts_compile() -> None:
    data_dir = _taxonomy_data_dir()
    payload = load_organic_groups_with_compositions(data_dir / "organic_groups.v1.3.json")
    group_map = {
        str(entry.get("id") or ""): entry
        for entry in (payload.get("groups", []) or [])
        if isinstance(entry, dict) and entry.get("id")
    }
    for group_id in (
        "-CON3",
        "-CONHNH2",
        "-COOH",
        "-COOR",
        "-OCONH2",
        "-SO2NHNH2",
        "-PO2OH",
        "-PO2OR",
        "-PO2NH2",
    ):
        smarts = str(group_map[group_id].get("smarts") or "")
        assert smarts
        assert compile_smarts(smarts, validate=True) is not None


def test_registry_includes_composed_groups() -> None:
    registry = build_compound_registry(_default_registry_paths())
    groups = registry.get("groups", {}) or {}
    assert "-CON3" in groups
    assert "-CONHNH2" in groups
    assert "-COOH" in groups
    assert "-COOR" in groups
    assert "-OCONH2" in groups
    assert "-SO2NHNH2" in groups
    assert "-PO2OH" in groups
    assert "-PO2OR" in groups
    assert "-PO2NH2" in groups


def test_registry_has_curated_composed_compounds() -> None:
    registry = build_compound_registry(_default_registry_paths())
    compound_map = registry.get("compound_map", {}) or {}
    assert "Ar-CON3" in compound_map
    assert "Alkyl-CONHNH2" in compound_map
    assert "Alkyl-OCONH2" in compound_map
    assert "Ar-COOH" in compound_map
    assert "Ar-COOR" in compound_map
    assert "Bn-SO2NHNH2" in compound_map
    assert "Ar-PO2OH" in compound_map
    assert "Alkyl-PO2OR" in compound_map
    assert "RCH2-PO2NH2" in compound_map


def test_compound_expansion_added_for_composed_groups() -> None:
    compounds = _load_compounds()
    pairs = {(str(entry.get("A") or ""), str(entry.get("B") or "")) for entry in compounds}
    expected = {
        ("Ar", "-CON3"),
        ("Ar", "-CONHNH2"),
        ("Ar", "-SO2NHNH2"),
        ("Alkyl", "-CON3"),
        ("Alkyl", "-CONHNH2"),
        ("Alkyl", "-OCONH2"),
        ("Alkyl", "-SO2NHNH2"),
        ("RCH2", "-CON3"),
        ("RCH2", "-CONHNH2"),
        ("RCH2", "-SO2NHNH2"),
        ("R2CH", "-CON3"),
        ("R2CH", "-CONHNH2"),
        ("R2CH", "-SO2NHNH2"),
        ("R3C", "-CON3"),
        ("R3C", "-CONHNH2"),
        ("R3C", "-SO2NHNH2"),
        ("Bn", "-CON3"),
        ("Bn", "-CONHNH2"),
        ("Bn", "-SO2NHNH2"),
        ("Allyl", "-CON3"),
        ("Allyl", "-CONHNH2"),
        ("Allyl", "-SO2NHNH2"),
        ("Alkenyl", "-CON3"),
        ("Alkenyl", "-CONHNH2"),
        ("Alkenyl", "-SO2NHNH2"),
        ("Alkynyl", "-CON3"),
        ("Alkynyl", "-CONHNH2"),
        ("Alkynyl", "-SO2NHNH2"),
        ("Ar", "-PO2OH"),
        ("Ar", "-PO2OR"),
        ("Ar", "-PO2NH2"),
        ("Ar", "-PO2NHR"),
        ("Ar", "-PO2NR2"),
        ("Alkyl", "-PO2OH"),
        ("Alkyl", "-PO2OR"),
        ("Alkyl", "-PO2NH2"),
        ("Alkyl", "-PO2NHR"),
        ("Alkyl", "-PO2NR2"),
        ("RCH2", "-PO2OH"),
        ("RCH2", "-PO2OR"),
        ("RCH2", "-PO2NH2"),
        ("RCH2", "-PO2NHR"),
        ("RCH2", "-PO2NR2"),
        ("R2CH", "-PO2OH"),
        ("R2CH", "-PO2OR"),
        ("R2CH", "-PO2NH2"),
        ("R2CH", "-PO2NHR"),
        ("R2CH", "-PO2NR2"),
        ("R3C", "-PO2OH"),
        ("R3C", "-PO2OR"),
        ("R3C", "-PO2NH2"),
        ("R3C", "-PO2NHR"),
        ("R3C", "-PO2NR2"),
    }
    missing = sorted(expected - pairs)
    assert not missing, f"Missing composed compound pairs: {missing}"


def test_compound_logic_sets_include_composed_compounds() -> None:
    sets = _load_compound_logic_sets()

    assert "Ar-CONHNH2" in sets["amines_nh"]
    assert "Alkyl-CONHNH2" in sets["amines_nh"]
    assert "Alkyl-OCONH2" in sets["amines_nh"]
    assert "H-NH2" in sets["amines_nh"]
    assert "H-CONH2" in sets["amines_nh"]
    assert "Ar-SO2NHNH2" in sets["amines_nh"]
    assert "Alkyl-SO2NHNH2" in sets["amines_nh"]

    assert "Ar-CONHNH2" in sets["amides"]
    assert "RCH2-CONHNH2" in sets["amides"]
    assert "H-CONH2" in sets["amides"]

    assert "Ar-SO2NHNH2" in sets["aryl_sulfonamides"]
    assert "Alkyl-SO2NHNH2" in sets["aryl_sulfonamides"]

    assert "Ar-CON3" in sets["azides"]
    assert "Alkyl-CON3" in sets["azides"]


def test_generated_only_groups_not_defined_in_base_organic_groups() -> None:
    ids = {str(entry.get("id") or "") for entry in _load_raw_groups()}
    forbidden = {
        "-COOH",
        "-COOR",
        "-CONH2",
        "-CONHR",
        "-CONR2",
        "-CON3",
        "-CONHNH2",
        "-OCONH2",
        "-SO2NH2",
        "-SO2NHR",
        "-SO2NR2",
        "-SO2NHNH2",
        "-PO2OH",
        "-PO2OR",
        "-PO2NH2",
        "-PO2NHR",
        "-PO2NR2",
    }
    overlap = sorted(ids & forbidden)
    assert not overlap, f"Generated-only groups present in base organic_groups: {overlap}"


def test_substituent_fragments_uses_terminal_group_and_no_terminals_section() -> None:
    payload = _load_substituent_fragments()
    assert "terminals" not in payload
    groups = payload.get("groups") or []
    assert groups
    for entry in groups:
        assert isinstance(entry, dict)
        assert entry.get("terminal_group")


def test_substituent_fragments_schema_validation_passes_current_payload() -> None:
    payload = _load_substituent_fragments()
    errors = validate_substituent_fragments_payload(payload)
    assert errors == []


def test_substituent_fragments_schema_rejects_deprecated_terminal_and_lowercase_label() -> None:
    payload = {
        "linkers": [{"label": "co", "smarts_template": "[CX3:2](=O){TAIL}"}],
        "groups": [{"id": "-CONH2", "label": "co", "terminal": "-NH2"}],
    }
    errors = validate_substituent_fragments_payload(payload)
    assert any("must be uppercase" in msg for msg in errors)
    assert any("deprecated key 'terminal'" in msg for msg in errors)


def test_substituent_fragments_schema_rejects_policy_overlap() -> None:
    payload = {
        "linkers": [
            {
                "label": "CO",
                "smarts_template": "[CX3:2](=O){TAIL}",
                "allowed_terminal_groups": ["-NH2", "-OH"],
                "blocked_terminal_groups": ["-OH"],
            }
        ],
        "groups": [{"id": "-CONH2", "label": "CO", "terminal_group": "-NH2"}],
    }
    errors = validate_substituent_fragments_payload(payload)
    assert any("policy overlap" in msg for msg in errors)


def test_compose_groups_respects_linker_terminal_policy() -> None:
    base_groups = [
        {"id": "-NH2", "smarts": "[NX3H2:2]", "kind": "substituent", "priority": 1},
        {"id": "-OH", "smarts": "[OX2H:2]", "kind": "substituent", "priority": 1},
    ]
    fragments = {
        "linkers": [
            {
                "label": "CO",
                "smarts_template": "[CX3:2](=O){TAIL}",
                "allowed_terminal_groups": ["-NH2"],
                "blocked_terminal_groups": ["-OH"],
            }
        ],
        "groups": [
            {"id": "-CONH2", "label": "CO", "terminal_group": "-NH2"},
            {"id": "-COOH", "label": "CO", "terminal_group": "-OH"},
        ],
    }
    generated, errors = compose_groups_from_fragments(
        base_groups=base_groups,
        fragments_payload=fragments,
    )
    ids = {str(entry.get("id") or "") for entry in generated}
    assert "-CONH2" in ids
    assert "-COOH" not in ids
    assert any("is not allowed for linker 'CO'" in msg for msg in errors)


def test_resolve_group_id_preserves_full_suffix_when_unmapped() -> None:
    assert _resolve_group_id("Ar-CONH2", None) == "-CONH2"
    assert _resolve_group_id("Ar-CO-NH2", None) == "-CO-NH2"

