"""Tests for the isolated weak-label recommendation fallback."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

from condition_recommender import (
    generate_weak_label_screening_array,
    load_weak_label_index,
    recommend_weak_label_conditions,
    validate_weak_label_retrieval_rules,
)


QUERY = "Brc1ccccc1.CN>>CNc1ccccc1"
FIELDS = (
    "source_reaction_type",
    "reactive_site_1_display_label",
    "reactive_site_1_signature",
    "reactive_site_1_center_class",
    "reactive_site_1_attachment_class",
    "reactive_site_1_alpha_branched",
    "reactive_site_2_display_label",
    "reactive_site_2_signature",
    "reactive_site_2_center_class",
    "reactive_site_2_attachment_class",
    "reactive_site_2_alpha_branched",
    "yield_pct",
    "z_score",
    "condition_recipe_id",
)


def _recipe(recipe_id: str, solvent_id: str) -> dict:
    recipe = {
        bucket: []
        for bucket in (
            "catalysts",
            "ligands",
            "bases",
            "acids",
            "oxidants",
            "reductants",
            "condensation_agents",
            "solvents",
            "additives",
            "other_components",
        )
    }
    recipe.update(
        {
            "recipe_id": recipe_id,
            "schema_version": "1.1",
            "temperature_c": 80.0,
            "time_h": 12.0,
            "concentration_m": None,
            "atmosphere": None,
            "declared_absences": [],
            "stages": [],
            "warnings": [],
        }
    )
    recipe["solvents"] = [
        {
            "substance_id": solvent_id,
            "canonical_name": solvent_id,
            "raw_identifier": solvent_id,
            "identity_status": "resolved",
            "role_status": "assigned",
            "primary_role": "solvent",
        }
    ]
    return recipe


def _dataset(tmp_path: Path) -> Path:
    path = tmp_path / "weak.csv"
    recipes = {
        "RCR1:a": _recipe("RCR1:a", "solvent:a"),
        "RCR1:b": _recipe("RCR1:b", "solvent:b"),
        "RCR1:wrong": _recipe("RCR1:wrong", "solvent:c"),
    }
    rows = [
        {
            "source_reaction_type": "Buchwald-Hartwig",
            "reactive_site_1_display_label": "Ar–Br",
            "reactive_site_1_signature": "LG|Ar|Br",
            "reactive_site_2_display_label": "Alk–NH₂",
            "reactive_site_2_signature": "XH|N|H2|Alkyl",
            "yield_pct": "82",
            "z_score": "1.2",
            "condition_recipe_id": "RCR1:a",
        },
        {
            "source_reaction_type": "CN-Coupling",
            "reactive_site_1_display_label": "Alk–NH₂",
            "reactive_site_1_signature": "XH|N|H2|Alkyl",
            "reactive_site_2_display_label": "Ar–Br",
            "reactive_site_2_signature": "LG|Ar|Br",
            "yield_pct": "90",
            "z_score": "1.5",
            "condition_recipe_id": "RCR1:a",
        },
        {
            "source_reaction_type": "SNAr",
            "reactive_site_1_display_label": "Ar–Br",
            "reactive_site_1_signature": "LG|Ar|Br",
            "reactive_site_2_display_label": "Alk–NH₂",
            "reactive_site_2_signature": "XH|N|H2|Alkyl",
            "yield_pct": "75",
            "z_score": "0.4",
            "condition_recipe_id": "RCR1:b",
        },
        {
            "source_reaction_type": "Buchwald-Hartwig",
            "reactive_site_1_display_label": "Ar–Br",
            "reactive_site_1_signature": "LG|Ar|Br",
            "reactive_site_2_display_label": "Alk–OH",
            "reactive_site_2_signature": "XH|O|H1|Alkyl",
            "yield_pct": "99",
            "z_score": "4.0",
            "condition_recipe_id": "RCR1:wrong",
        },
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    catalog = path.with_name("weak.condition_recipes.jsonl.gz")
    with gzip.open(catalog, "wt", encoding="utf-8") as handle:
        for recipe in recipes.values():
            handle.write(json.dumps(recipe, sort_keys=True) + "\n")
    return path


def test_weak_label_rules_and_index_are_valid(tmp_path: Path) -> None:
    path = _dataset(tmp_path)
    index = load_weak_label_index(path)

    assert validate_weak_label_retrieval_rules() == []
    assert len(index.rows) == 4
    assert len(index.select_types(("Buchwald-Hartwig",))) == 2
    assert index.rows[0].participants[0].display_label == "Ar–Br"


def test_fallback_requires_graph_hint_and_matches_unordered_sites(
    tmp_path: Path,
) -> None:
    result = recommend_weak_label_conditions(
        QUERY,
        records_path=_dataset(tmp_path),
        top_k=5,
    )

    assert result.valid is True
    assert result.reaction_type_id == "sp2_c_n_substitution"
    assert result.source_reaction_type_candidates == (
        "Buchwald-Hartwig",
        "CN-Coupling",
        "SNAr",
    )
    assert result.candidate_count == 4
    assert result.compatible_candidate_count == 3
    assert result.excluded_candidate_count == 1
    assert [item.recipe_id for item in result.recommendations] == [
        "RCR1:a",
        "RCR1:b",
    ]
    assert result.recommendations[0].support == 2
    assert result.schema_version == "1.1"
    assert result.recommendations[0].source_matches[0].source_reaction_type == (
        "CN-Coupling"
    )
    assert result.recommendations[0].source_matches[0].participant_roles == (
        "electrophile",
        "nitrogen_partner",
    )
    assert (
        result.recommendations[0].source_matches[0].participant_display_labels
        == ("Ar–Br", "Alk–NH₂")
    )
    assert result.recommendations[0].source_matches[0].participant_signatures == (
        "LG|Ar|Br",
        "XH|N|H2|Alkyl",
    )
    assert "WEAK_LABEL_PRECEDENTS_NOT_STRUCTURE_VERIFIED" in result.warnings


def test_source_hint_can_narrow_but_cannot_contradict_graph(tmp_path: Path) -> None:
    path = _dataset(tmp_path)
    narrowed = recommend_weak_label_conditions(
        QUERY,
        records_path=path,
        source_reaction_type_hint="SNAr",
    )
    contradicted = recommend_weak_label_conditions(
        QUERY,
        records_path=path,
        source_reaction_type_hint="Suzuki-Miyaura",
    )

    assert narrowed.valid is True
    assert narrowed.source_reaction_type_candidates == ("SNAr",)
    assert [item.recipe_id for item in narrowed.recommendations] == ["RCR1:b"]
    assert contradicted.valid is False
    assert contradicted.error == "INCOMPATIBLE_SOURCE_REACTION_TYPE_HINT"


def test_screening_array_is_deterministic_and_diversified(tmp_path: Path) -> None:
    path = _dataset(tmp_path)
    first = generate_weak_label_screening_array(
        QUERY,
        records_path=path,
        array_size=2,
    )
    second = generate_weak_label_screening_array(
        QUERY,
        records_path=path,
        array_size=2,
    )

    assert first.valid is True
    assert first.recommendation_mode == "weak_label_screening"
    assert first.to_dict() == second.to_dict()
    assert {item.recipe_id for item in first.recommendations} == {
        "RCR1:a",
        "RCR1:b",
    }


def test_fallback_rejects_unverified_or_unsupported_graph_query(
    tmp_path: Path,
) -> None:
    result = recommend_weak_label_conditions(
        "CCO>>CCO",
        records_path=_dataset(tmp_path),
    )

    assert result.valid is False
    assert result.error in {
        "QUERY_HAS_NO_VERIFIED_REACTION_SIGNATURE",
        "QUERY_NOT_SUPPORTED_BY_WEAK_LABEL_DATASET",
    }
