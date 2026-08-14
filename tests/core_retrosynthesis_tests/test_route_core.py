"""Route-core projection and cross-step lineage regressions."""

from __future__ import annotations

from dataclasses import replace

from core_retrosynthesis.route_conversion import build_observed_route_tree
from core_retrosynthesis.route_core import (
    ROUTE_CORE_SCHEMA_VERSION,
    RouteCoreProjection,
    build_route_core_projection,
    validate_route_core_projection,
)


def _two_step_record(*, symmetric_intermediate: bool = False) -> dict:
    if symmetric_intermediate:
        target = "CCBr"
        steps = [
            {
                "source_reaction_id": "consumer",
                "product_smiles": target,
                "precursor_smiles": ["BrBr", "CC"],
                "reactants_smiles": "[CH3:10][CH3:11].[Br:12][Br:13]",
                "reagents_smiles": "",
                "product_smiles_mapped": "[CH3:10][CH2:11][Br:12]",
                "reaction_smiles": (
                    "[CH3:10][CH3:11].[Br:12][Br:13]"
                    ">>[CH3:10][CH2:11][Br:12]"
                ),
                "abstracted_reaction_smiles": "",
            },
            {
                "source_reaction_id": "producer",
                "product_smiles": "CC",
                "precursor_smiles": ["C=C"],
                "reactants_smiles": "[CH2:1]=[CH2:2]",
                "reagents_smiles": "",
                "product_smiles_mapped": "[CH3:1][CH3:2]",
                "reaction_smiles": "[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]",
                "abstracted_reaction_smiles": "",
            },
        ]
    else:
        target = "CNC"
        steps = [
            {
                "source_reaction_id": "consumer",
                "product_smiles": target,
                "precursor_smiles": ["CBr", "CN"],
                "reactants_smiles": "[CH3:10][NH2:11].[CH3:12][Br:13]",
                "reagents_smiles": "",
                "product_smiles_mapped": "[CH3:10][NH:11][CH3:12]",
                "reaction_smiles": (
                    "[CH3:10][NH2:11].[CH3:12][Br:13]"
                    ">>[CH3:10][NH:11][CH3:12]"
                ),
                "abstracted_reaction_smiles": "",
            },
            {
                "source_reaction_id": "producer",
                "product_smiles": "CN",
                "precursor_smiles": ["CBr", "N"],
                "reactants_smiles": "[CH3:1][Br:2].[NH3:3]",
                "reagents_smiles": "",
                "product_smiles_mapped": "[CH3:1][NH2:3]",
                "reaction_smiles": (
                    "[CH3:1][Br:2].[NH3:3]>>[CH3:1][NH2:3]"
                ),
                "abstracted_reaction_smiles": "",
            },
        ]
    return {
        "schema_version": "1.0",
        "route_id": "route-a",
        "patent_id": "patent-a",
        "split": "train",
        "target_smiles": target,
        "original_reaction_count": 2,
        "higher_level_reaction_count": 2,
        "higher_level_depth": 2,
        "abstraction_reduction": 0,
        "steps": steps,
    }


def test_route_core_projection_contains_step_chemistry_and_unique_lineage() -> None:
    tree = build_observed_route_tree(_two_step_record())

    projection = build_route_core_projection(tree)

    assert projection.schema_version == ROUTE_CORE_SCHEMA_VERSION
    assert projection.fully_chemistry_resolved
    assert projection.fully_lineage_connected
    assert len(projection.steps) == 2
    assert all(step.reaction_signature_id for step in projection.steps)
    assert all(step.reaction_core_id for step in projection.steps)
    assert all(step.minimum_reaction_smiles for step in projection.steps)
    assert len(projection.motifs) == 2
    assert {motif.abstraction_level for motif in projection.motifs} == {
        "typed",
        "shape",
    }
    assert all(len(motif.core_keys) == 2 for motif in projection.motifs)
    assert projection.lineage_links[0].status == "unique"
    assert projection.lineage_links[0].candidates[0].atom_map_pairs == (
        (1, 10),
        (3, 11),
    )
    restored = RouteCoreProjection.from_dict(projection.to_dict())
    assert restored == projection


def test_route_core_projection_preserves_symmetric_lineage_candidates() -> None:
    tree = build_observed_route_tree(
        _two_step_record(symmetric_intermediate=True)
    )

    projection = build_route_core_projection(tree)

    link = projection.lineage_links[0]
    assert link.status == "symmetry_ambiguous"
    assert link.candidate_count == 2
    assert len(link.candidates) == 2
    assert "ROUTE_CORE_LINEAGE_HAS_SYMMETRY_AMBIGUITY" in projection.warnings


def test_route_core_validator_rejects_wrong_declared_reaction_count() -> None:
    projection = build_route_core_projection(
        build_observed_route_tree(_two_step_record())
    )

    report = validate_route_core_projection(
        replace(projection, reaction_count=3)
    )

    assert not report.valid
    assert "reaction_count_mismatch" in report.issues


def test_route_core_motifs_exclude_incomplete_chemistry_sequences() -> None:
    record = _two_step_record()
    producer = record["steps"][1]
    producer["reactants_smiles"] = "[CH3:1][Cl:3]"
    producer["reaction_smiles"] = (
        "[CH3:1][Cl:3]>>[CH3:1][NH2:3]"
    )

    projection = build_route_core_projection(
        build_observed_route_tree(record)
    )

    assert any(step.quality_status == "unavailable" for step in projection.steps)
    assert projection.motifs == ()
    assert projection.typed_motif_keys == ()
    assert projection.shape_motif_keys == ()
