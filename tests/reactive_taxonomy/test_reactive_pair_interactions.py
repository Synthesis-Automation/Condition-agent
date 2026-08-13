"""Maintainable intramolecular reactive-pair rule regressions."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

from reactive_taxonomy import (
    assess_reactive_pair_interactions,
    load_reactive_pair_interaction_definition,
    validate_reactive_pair_interaction_definition,
    validate_taxonomy,
)


REPORTED_ACYL_CHLORIDE = (
    "COc1cc(O)c(C(=O)Nc2nc(C(=O)Cl)cs2)cc1OC"
)


def test_reactive_pair_definition_is_versioned_and_taxonomy_valid() -> None:
    definition = load_reactive_pair_interaction_definition()

    assert definition["definition_id"] == "reactive_pair_interactions.v1"
    assert definition["schema_version"] == "1.0"
    assert definition["allowed_scopes"] == ["same_component"]
    assert not validate_reactive_pair_interaction_definition(definition)
    assert not validate_taxonomy()


def test_warning_rules_do_not_invalidate_signature_or_condition_indexes() -> None:
    manifest_path = (
        Path(__file__).resolve().parents[2]
        / "reactive_taxonomy"
        / "definitions"
        / "taxonomy_manifest.v4.json"
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert "reactive_pair_interactions.v1.json" not in {
        *manifest["identity_definitions"],
        *manifest["annotation_definitions"],
    }


def test_reactive_pair_definition_rejects_unknown_selector_operator() -> None:
    definition = deepcopy(load_reactive_pair_interaction_definition())
    definition["rules"][0]["left_site"]["constraints"][0]["operator"] = "python"

    errors = validate_reactive_pair_interaction_definition(definition)

    assert any("unknown_operator" in error for error in errors)


def test_grignard_and_alcohol_warn_only_in_one_component() -> None:
    same_component = assess_reactive_pair_interactions("C([Mg]Br)CCO")
    separate_components = assess_reactive_pair_interactions("C[Mg]Br.CCO")

    assert [item.rule_id for item in same_component] == [
        "RPI001_ORGANOMETALLIC_PROTIC_SELF_QUENCH"
    ]
    assert same_component[0].scope == "same_component"
    assert same_component[0].warning_strength == "strong"
    assert same_component[0].intrinsic_severity == "critical"
    assert separate_components == ()


def test_acyl_halide_and_amine_warn_only_in_one_component() -> None:
    same_component = assess_reactive_pair_interactions("NCCC(=O)Cl")
    separate_components = assess_reactive_pair_interactions("CN.CC(=O)Cl")

    assert [item.rule_id for item in same_component] == [
        "RPI003_ACYL_HALIDE_FREE_AMINE"
    ]
    assert same_component[0].potential_closure_ring_size == 4
    assert separate_components == ()


def test_reported_acyl_chloride_has_strong_nine_membered_phenol_alert() -> None:
    assessments = assess_reactive_pair_interactions(REPORTED_ACYL_CHLORIDE)

    assert [item.rule_id for item in assessments] == [
        "RPI005_ACYL_HALIDE_FREE_HYDROXYL"
    ]
    assessment = assessments[0]
    assert assessment.left_site.chemist_label == "HetAr–C(O)Cl"
    assert assessment.right_site.chemist_label == "Ar–OH"
    assert assessment.graph_distance == 8
    assert assessment.potential_closure_ring_size == 9
    assert assessment.warning_strength == "strong"


def test_grignard_and_ketone_intramolecular_pair_is_detected() -> None:
    assessments = assess_reactive_pair_interactions("C([Mg]Br)CCC(=O)C")

    assert "RPI002_ORGANOMETALLIC_CARBONYL_SELF_REACTION" in {
        item.rule_id for item in assessments
    }
