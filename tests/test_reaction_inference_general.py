from __future__ import annotations

import pytest

from chemtools.reaction_inference import analyze_reaction_general
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_general_inference_detects_amide_formation() -> None:
    rxn = "O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1"
    result = analyze_reaction_general(rxn)

    assert result.decision.reaction_type == "Amide_formation"
    assert result.decision.source == "deterministic"
    assert result.decision.confidence >= 0.9
    assert result.validation.passed is True
    assert result.taxonomy_candidates


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_general_inference_detects_oxidation() -> None:
    rxn = "CCO>>CC=O"
    result = analyze_reaction_general(rxn)

    assert result.decision.reaction_type == "Oxidation_primary_alcohol_to_aldehyde"
    assert result.decision.mechanism_family == "oxidation"
    assert result.validation.passed is True


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_general_inference_flags_taxonomy_gap_without_forcing_type() -> None:
    rxn = "Clc1ncc(-c2ccccc2)cn1.NN>>NN=c1ncc(-c2ccccc2)c[nH]1"
    result = analyze_reaction_general(rxn)

    assert result.decision.reaction_type == "unknown"
    assert result.taxonomy_candidates == []
    assert result.formula_delta["Cl"] == -1
    assert result.formula_delta["N"] == 2
    assert result.features["core_reactant_aromatic_ring_n_count"] >= 2
    assert result.features["side_reactants_have_nn"] is True
    assert "tautomer" in result.summary.lower()


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_general_inference_is_not_specific_to_single_reaction_family() -> None:
    rxn = "CN.CBr>>CNC"
    result = analyze_reaction_general(rxn)

    assert result.decision.reaction_type == "Alkyl_Nucleophilic_Substitution"
    assert result.validation.passed is True
    assert result.taxonomy_candidates[0]["reaction_type"] == "Alkyl_Nucleophilic_Substitution"


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_general_inference_can_emit_mechanism_class_for_taxonomy_gap_annulation() -> None:
    rxn = (
        "Cc1ccc(Nc2ccc(C)cc2)cc1.O=C(O)C=Cc1ccoc1>>"
        "Cc1ccc(-n2c(=O)cc(-c3ccoc3)c3cc(C)ccc32)cc1"
    )
    result = analyze_reaction_general(rxn)

    assert result.decision.reaction_type == "unknown"
    assert result.decision.mechanism_family == "annulation_cyclization"
    assert "Ann" in (result.evidence.get("certain") or {}).get("event_tokens", [])


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_general_inference_detects_stille_with_terminal_alkenyl_product_motif() -> None:
    rxn = "Ic1ccncc1.C=C[Sn](C)(C)C>>C=Cc1ccncc1"
    result = analyze_reaction_general(rxn)

    assert result.decision.reaction_type.lower() == "stille"
    assert result.decision.mechanism_family == "cross_coupling"
    assert any(
        str(row.get("reaction_type", "")).lower() == "stille"
        for row in result.taxonomy_candidates
    )
