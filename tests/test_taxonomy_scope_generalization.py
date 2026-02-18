from __future__ import annotations

from chemtools.featurizers.formatters.reaction import _project_formed_motifs_by_taxonomy
from chemtools.reaction_inference import _taxonomy_consistency_check


def test_taxonomy_consistency_check_accepts_scope_child_reactant() -> None:
    ok, reason = _taxonomy_consistency_check(
        reaction_type="Epoxidation",
        reacted_motifs=["RCH2-Alkenyl"],
        formed_motifs=["Epoxide"],
        formed_bonds=["C-O"],
        broken_bonds=[],
    )
    assert ok is True, reason


def test_product_projection_accepts_scope_child_product_motif() -> None:
    projected = _project_formed_motifs_by_taxonomy(
        reaction_type="Oxidation_primary_alcohol_to_aldehyde",
        formed_in_product={"RCH2-CHO"},
        inferred_in_product=[],
    )
    assert "RCH2-CHO" in projected
