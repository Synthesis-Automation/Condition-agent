"""Tests for standalone GPT-5.2 reaction analysis PoC."""

from poc_gpt52_reaction import analyze_reaction


def test_new_md_example_prefers_snar_hydrazinolysis() -> None:
    reaction = "Clc1ncc(-c2ccccc2)cn1>>NN=c1ncc(-c2ccccc2)c[nH]1"
    result = analyze_reaction(reaction, use_llm=False)

    assert result.best_hypothesis["reaction_class"] == "SNAr_hydrazinolysis"
    assert result.confidence >= 0.8
    assert result.formula_delta["Cl"] == -1
    assert result.formula_delta["N"] == 2
    assert result.features["core_reactant_aromatic_ring_n_count"] >= 2
    assert "tautomer" in result.summary.lower()


def test_multi_reactant_selects_principal_pair_and_keeps_classification() -> None:
    reaction = "Clc1ncc(-c2ccccc2)cn1.NN>>NN=c1ncc(-c2ccccc2)c[nH]1"
    result = analyze_reaction(reaction, use_llm=False)

    assert result.best_hypothesis["reaction_class"] == "SNAr_hydrazinolysis"
    assert len(result.reactant_smiles_list) == 2
    assert result.features["reactant_count"] == 2
    assert result.selected_pair["reactant_smiles"] == "Clc1ncc(-c2ccccc2)cn1"
    assert result.formula_delta["Cl"] == -1
    assert result.formula_delta["N"] == 2
    assert result.features["side_reactants_have_nn"] is True
