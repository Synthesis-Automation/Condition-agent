"""Deterministic strategic-complexity regressions."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

import pytest

from reactive_taxonomy import (
    assess_molecule_complexity,
    assess_retrosynthetic_complexity_reduction,
    complex_target_requires_strategic_candidate,
    load_strategic_complexity_definition,
    validate_strategic_complexity_definition,
    validate_taxonomy,
)


FUSED_ESTERIFICATION = (
    "[CH3:1][I:900]."
    "[OH:2][C:3](=[O:4])[c:5]1[cH:6][n:7][c:8]2[cH:9][n:10]"
    "[c:11]3[nH:12][cH:13][cH:14][c:15]3[n:16]12>>"
    "[CH3:1][O:2][C:3](=[O:4])[c:5]1[cH:6][n:7][c:8]2[cH:9]"
    "[n:10][c:11]3[nH:12][cH:13][cH:14][c:15]3[n:16]12"
)
BALANCED_SCAFFOLD_COUPLING = (
    "[CH3:1][CH2:2][CH2:3][Br:10]."
    "[CH3:4][CH2:5][CH2:6][Cl:11]>>"
    "[CH3:1][CH2:2][CH2:3][CH2:6][CH2:5][CH3:4]"
)
RING_CONSTRUCTION = (
    "[NH2:1][CH2:2][CH2:3][CH2:4][CH2:5][CH2:6][Br:9]>>"
    "[NH:1]1[CH2:2][CH2:3][CH2:4][CH2:5][CH2:6]1"
)


def test_strategic_complexity_definition_is_versioned_and_valid() -> None:
    definition = load_strategic_complexity_definition()

    assert definition["definition_id"] == "strategic_complexity.v1"
    assert definition["schema_version"] == "1.0"
    assert not validate_strategic_complexity_definition(definition)
    assert not validate_taxonomy()


def test_definition_rejects_non_normalized_component_weights() -> None:
    definition = deepcopy(load_strategic_complexity_definition())
    definition["retrosynthetic_reduction"]["component_weights"][
        "core_fragmentation"
    ] = 0.9

    errors = validate_strategic_complexity_definition(definition)

    assert "strategic_complexity:weights_do_not_sum_to_one" in errors


def test_strategic_ranking_policy_does_not_change_signature_identity() -> None:
    manifest_path = (
        Path(__file__).resolve().parents[2]
        / "reactive_taxonomy"
        / "definitions"
        / "taxonomy_manifest.v4.json"
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert "strategic_complexity.v1.json" not in {
        *manifest["identity_definitions"],
        *manifest["annotation_definitions"],
    }


def test_fused_target_esterification_is_tactical_not_strategic() -> None:
    assessment = assess_retrosynthetic_complexity_reduction(
        FUSED_ESTERIFICATION
    )

    assert assessment.target.heavy_atom_count == 16
    assert assessment.target.cycle_rank == 3
    assert assessment.product_derived_component_heavy_atom_counts == (15, 1)
    assert assessment.largest_retained_core_fraction == pytest.approx(15 / 16)
    assert assessment.ring_topology_reduction_score == 0.0
    assert assessment.strategic_class == "functional_group_interconversion"
    assert assessment.score == 0.0
    assert assessment.is_strategic is False
    assert complex_target_requires_strategic_candidate(assessment) is True


def test_balanced_scaffold_split_is_strategic_and_order_invariant() -> None:
    assessment = assess_retrosynthetic_complexity_reduction(
        BALANCED_SCAFFOLD_COUPLING
    )
    reversed_components = assess_retrosynthetic_complexity_reduction(
        BALANCED_SCAFFOLD_COUPLING.replace(
            "[CH3:1][CH2:2][CH2:3][Br:10]."
            "[CH3:4][CH2:5][CH2:6][Cl:11]",
            "[CH3:4][CH2:5][CH2:6][Cl:11]."
            "[CH3:1][CH2:2][CH2:3][Br:10]",
        )
    )

    assert assessment.product_derived_component_heavy_atom_counts == (3, 3)
    assert assessment.core_fragmentation_score == 1.0
    assert assessment.convergency_score == 1.0
    assert assessment.strategic_class == "scaffold_split"
    assert assessment.is_strategic is True
    assert assessment.score == reversed_components.score
    assert assessment.to_dict() == reversed_components.to_dict()


def test_ring_construction_is_strategic_without_fragmentation() -> None:
    assessment = assess_retrosynthetic_complexity_reduction(RING_CONSTRUCTION)

    assert assessment.product_derived_component_heavy_atom_counts == (6,)
    assert assessment.core_fragmentation_score == 0.0
    assert assessment.formed_product_ring_bond_count == 1
    assert assessment.ring_topology_reduction_score == 1.0
    assert assessment.strategic_class == "ring_construction"
    assert assessment.is_strategic is True


def test_molecule_trace_rejects_invalid_input_and_clears_maps() -> None:
    trace = assess_molecule_complexity("[CH3:1][CH2:2]O")

    assert trace.canonical_smiles == "CCO"
    assert trace.heavy_atom_count == 3
    with pytest.raises(ValueError, match="invalid molecule"):
        assess_molecule_complexity("not smiles")
