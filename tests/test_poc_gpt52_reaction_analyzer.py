"""Tests for general reaction analysis evidence output."""

from chemtools.reaction_inference import analyze_reaction_general


def test_new_md_example_flags_taxonomy_gap_with_snar_like_evidence() -> None:
    reaction = "Clc1ncc(-c2ccccc2)cn1>>NN=c1ncc(-c2ccccc2)c[nH]1"
    result = analyze_reaction_general(reaction, use_llm=False)

    assert result.decision.reaction_type == "unknown"
    assert result.decision.mechanism_family == "nucleophilic_substitution"
    assert result.formula_delta["Cl"] == -1
    assert result.formula_delta["N"] == 2
    assert result.features["core_reactant_aromatic_ring_n_count"] >= 2
    assert "tautomer" in result.summary.lower()


def test_multi_reactant_selects_principal_pair_and_keeps_core_features() -> None:
    reaction = "Clc1ncc(-c2ccccc2)cn1.NN>>NN=c1ncc(-c2ccccc2)c[nH]1"
    result = analyze_reaction_general(reaction, use_llm=False)

    assert len(result.reactant_smiles_list) == 2
    assert result.features["reactant_count"] == 2
    assert result.principal_pair["reactant_smiles"] == "Clc1ncc(-c2ccccc2)cn1"
    assert result.formula_delta["Cl"] == -1
    assert result.formula_delta["N"] == 2
    assert result.features["side_reactants_have_nn"] is True
