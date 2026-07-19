import json
from pathlib import Path

import pytest
from reactive_taxonomy import (
    REACTION_SIGNATURE_SCHEMA_VERSION,
    reaction_signature_definition_versions,
)

from condition_recommender.evaluation import (
    evaluate_generic_index,
    grouped_holdout_split,
)
from condition_recommender.generic_indexing import (
    build_generic_index,
    load_generic_index,
    load_persisted_generic_index,
    save_generic_index,
)


def _signature(index: int) -> dict:
    return {
        "schema_version": REACTION_SIGNATURE_SCHEMA_VERSION,
        "definition_versions": reaction_signature_definition_versions(),
        "signature_id": f"RS1:{index}",
        "exact_signature_key": f"L0:{index}",
        "handle_signature_key": "L1:shared",
        "transformation_signature_key": "L2:shared",
        "bond_edit_signature_key": "L3:shared",
        "environment_signature_key": f"L4:{index % 2}",
        "formed_bond_types": ["C-N:SINGLE"],
        "broken_bond_types": ["B-C:SINGLE"],
        "order_changes": [],
        "partners": [
            {
                "handle_tokens": ["B(OH)2", "N-H"],
                "anchor_contexts": ["Ar"],
                "steric": {"class": "open"},
                "electronic": {"class": "neutral"},
                "nearby_groups": [],
            }
        ],
        "spectator_groups": [],
        "transformation_class": "generic_c_n_coupling",
        "named_family": None,
        "family_confidence": 0.0,
        "topology": {
            "reaction_scope": "intermolecular",
            "formed_bond_scopes": ["intermolecular"],
            "formed_ring_sizes": [],
            "ring_count_delta": 0,
        },
    }


def _record(index: int, *, canonical_group: str | None = None) -> dict:
    recipe_id = f"RCR1:{index % 2}"
    return {
        "schema_version": "1.7",
        "converter_definition_version": "generic_conversion.v1.2",
        "admission_tier": "verified",
        "reaction_id": f"reaction-{index}",
        "observation_id": f"observation-{index}",
        "canonical_reaction_id": canonical_group or f"CRX1:{index}",
        "reaction_smiles": "C.N>>CN",
        "yield_pct": 50.0 + index,
        "source_dataset": f"dataset-{index % 2}",
        "reaction_signature": _signature(index),
        "resolved_recipe_id": recipe_id,
        "resolved_recipe": {"recipe_id": recipe_id},
        "condition_resolution": {"has_uncertainty": False},
    }


def test_persisted_index_round_trip_is_deterministic(tmp_path: Path) -> None:
    index = build_generic_index([_record(1), _record(2)])
    first_path = tmp_path / "first.json"
    second_path = tmp_path / "second.json"
    first = save_generic_index(index, first_path)
    second = save_generic_index(index, second_path)
    loaded = load_persisted_generic_index(first_path)
    assert first["index_id"] == second["index_id"]
    assert first_path.read_bytes() == second_path.read_bytes()
    assert loaded == index
    assert load_generic_index(first_path) == index
    payload = json.loads(first_path.read_text(encoding="utf-8"))
    assert payload["schema_version"] == "1.1"
    assert payload["reaction_signature_schema_version"] == "1.2"
    assert payload["record_schema_versions"] == ["1.7"]


def test_index_rejects_stale_converted_records() -> None:
    stale = _record(1)
    stale["reaction_signature"] = {
        **stale["reaction_signature"],
        "schema_version": "1.0",
    }
    with pytest.raises(ValueError, match="regenerate converted records"):
        build_generic_index([stale])


def test_loader_rejects_stale_index_schema(tmp_path: Path) -> None:
    path = tmp_path / "index.json"
    save_generic_index(build_generic_index([_record(1)]), path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["schema_version"] = "1.0"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="rebuild the index"):
        load_persisted_generic_index(path)


def test_persisted_index_detects_tampering(tmp_path: Path) -> None:
    path = tmp_path / "index.json"
    save_generic_index(build_generic_index([_record(1)]), path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["rows"][0]["yield_pct"] = 99.0
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="integrity"):
        load_persisted_generic_index(path)


def test_grouped_split_keeps_duplicate_reactions_together() -> None:
    index = build_generic_index(
        [
            _record(1, canonical_group="CRX1:duplicate"),
            _record(2, canonical_group="CRX1:duplicate"),
            _record(3),
            _record(4),
        ]
    )
    first = grouped_holdout_split(index.rows, test_fraction=0.5, seed=7)
    second = grouped_holdout_split(index.rows, test_fraction=0.5, seed=7)
    assert first == second
    assert not set(first.train_group_ids).intersection(first.test_group_ids)
    duplicate_locations = {
        "train" if row in first.train_rows else "test"
        for row in index.rows
        if row.canonical_reaction_id == "CRX1:duplicate"
    }
    assert len(duplicate_locations) == 1


def test_grouped_evaluation_writes_leakage_safe_metrics(tmp_path: Path) -> None:
    index_path = tmp_path / "index.json"
    output = tmp_path / "evaluation"
    save_generic_index(
        build_generic_index([_record(index) for index in range(10)]),
        index_path,
    )
    report = evaluate_generic_index(
        index_path,
        output,
        test_fraction=0.3,
        seed=11,
        top_k=2,
        minimum_pool_size=1,
    )
    assert report["split"]["leakage_group_count"] == 0
    assert report["split"]["test_group_count"] == 3
    assert report["metrics"]["query_count"] == 3
    assert report["metrics"]["coverage_rate"] == 1.0
    assert report["metrics"]["hard_incompatible_recommendation_count"] == 0
    assert (output / "evaluation_report.json").is_file()
    assert len((output / "evaluation_cases.jsonl").read_text().splitlines()) == 3


def test_grouped_split_validates_inputs() -> None:
    one_group = build_generic_index([_record(1)]).rows
    with pytest.raises(ValueError, match="two canonical"):
        grouped_holdout_split(one_group)
    with pytest.raises(ValueError, match="between zero and one"):
        grouped_holdout_split(one_group, test_fraction=1.0)
