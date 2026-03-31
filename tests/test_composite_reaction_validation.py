from chem_coworker.tools.composite import (
    _forward_synthesis_step,
    _recommend_reaction_conditions,
    _retrosynthesis_step,
)


def test_recommend_reaction_conditions_rejects_missing_product_side() -> None:
    result = _recommend_reaction_conditions(reaction_smiles="CCO>>")
    assert result["success"] is False
    assert "missing product" in result["error"].lower()


def test_forward_synthesis_step_skips_conditions_when_no_top_product(monkeypatch) -> None:
    monkeypatch.setattr(
        "chem_coworker.tools.forward_synthesis._inspect_reactants",
        lambda **kwargs: {"success": True},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.forward_synthesis._identify_reactions",
        lambda **kwargs: {"success": True, "total_matched": 0},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.forward_synthesis._generate_products",
        lambda **kwargs: {"success": True, "total_products": 0, "products": []},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.forward_synthesis._rank_products",
        lambda **kwargs: {"success": True, "ranked_products": []},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.forward_synthesis._find_forward_precedent",
        lambda **kwargs: {"success": True},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.forward_synthesis._search_reactant_precedent",
        lambda **kwargs: {"success": True},  # noqa: ARG005
    )

    result = _forward_synthesis_step(smiles_a="CCO", include_conditions=True)
    assert result["success"] is True
    cond = result.get("conditions_for_top_product", {})
    assert cond.get("success") is True
    assert "skipped condition recommendation" in str(cond.get("note", "")).lower()


def test_forward_synthesis_step_flags_invalid_top_product_reaction_smiles(monkeypatch) -> None:
    monkeypatch.setattr(
        "chem_coworker.tools.forward_synthesis._inspect_reactants",
        lambda **kwargs: {"success": True},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.forward_synthesis._identify_reactions",
        lambda **kwargs: {"success": True, "total_matched": 1},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.forward_synthesis._generate_products",
        lambda **kwargs: {
            "success": True,
            "total_products": 1,
            "products": [{"product_smiles": "not_a_smiles", "template_name": "x", "taxonomy_id": "x"}],
        },  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.forward_synthesis._rank_products",
        lambda **kwargs: {"success": True, "ranked_products": []},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.forward_synthesis._find_forward_precedent",
        lambda **kwargs: {"success": True},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.forward_synthesis._search_reactant_precedent",
        lambda **kwargs: {"success": True},  # noqa: ARG005
    )

    result = _forward_synthesis_step(smiles_a="CCO", include_conditions=True)
    assert result["success"] is True
    cond = result.get("conditions_for_top_product", {})
    assert cond.get("success") is False
    assert "invalid top-product reaction smiles" in str(cond.get("error", "")).lower()


def test_retrosynthesis_step_promotes_concrete_fallback_and_surfaces_reaction_smiles(monkeypatch) -> None:
    monkeypatch.setattr(
        "chem_coworker.tools.retrosynthesis._inspect_target",
        lambda **kwargs: {"success": True},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.retrosynthesis._identify_retrons",
        lambda **kwargs: {"success": True, "total_matched": 1},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.retrosynthesis._generate_disconnections",
        lambda **kwargs: {
            "success": True,
            "disconnections": [
                {
                    "rank": 1,
                    "reaction_name": "conceptual_only",
                    "precursor_1": "",
                    "precursor_2": "",
                    "overall_score": 0.99,
                }
            ],
        },  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.retrosynthesis._apply_hte_templates",
        lambda **kwargs: {
            "success": True,
            "disconnections": [
                {
                    "template_name": "buchwald_hartwig",
                    "precursor_1": "ArBr",
                    "precursor_2": "HNR2",
                    "reaction_smiles": "ArBr.HNR2>>Target",
                    "difficulty": 0.3,
                }
            ],
        },  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.retrosynthesis._search_by_product_similarity",
        lambda **kwargs: {"success": True, "precedents": []},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.reaction_eval._check_retro_consistency",
        lambda **kwargs: {"success": True, "verdict": "PASS"},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.composite._evaluate_synthesis_proposal",
        lambda **kwargs: {
            "success": True,
            "verdict": "PASS",
            "evaluation_score": 91.0,
            "evaluation_grade": "A",
        },  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.composite._recommend_conditions_with_strategy",
        lambda reaction_smiles, **kwargs: {"success": True, "reaction_smiles": reaction_smiles},  # noqa: ARG005
    )

    result = _retrosynthesis_step(target_smiles="Target", include_conditions=True)

    assert result["success"] is True
    assert result["summary"]["fallback_used"] is True
    assert result["top_disconnection"]["reaction_name"] == "buchwald_hartwig"
    assert result["top_disconnection"]["reaction_smiles"] == "ArBr.HNR2>>Target"
    assert result["disconnections"]["disconnections"][0]["reaction_smiles"] == "ArBr.HNR2>>Target"
    assert result["top_reaction_smiles_evaluation"]["success"] is True
    assert result["conditions_for_top_disconnection"]["reaction_smiles"] == "ArBr.HNR2>>Target"


def test_retrosynthesis_step_adds_reaction_smiles_to_standard_disconnections(monkeypatch) -> None:
    monkeypatch.setattr(
        "chem_coworker.tools.retrosynthesis._inspect_target",
        lambda **kwargs: {"success": True},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.retrosynthesis._identify_retrons",
        lambda **kwargs: {"success": True, "total_matched": 1},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.retrosynthesis._generate_disconnections",
        lambda **kwargs: {
            "success": True,
            "disconnections": [
                {
                    "rank": 1,
                    "reaction_name": "buchwald_hartwig",
                    "precursor_1": "ArBr",
                    "precursor_2": "HNR2",
                    "overall_score": 0.8,
                }
            ],
        },  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.retrosynthesis._apply_hte_templates",
        lambda **kwargs: {"success": True, "disconnections": []},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.reaction_eval._check_retro_consistency",
        lambda **kwargs: {"success": True, "verdict": "PASS"},  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.composite._evaluate_synthesis_proposal",
        lambda **kwargs: {
            "success": True,
            "verdict": "PASS",
            "evaluation_score": 90.0,
            "evaluation_grade": "A",
        },  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.composite._recommend_conditions_with_strategy",
        lambda reaction_smiles, **kwargs: {"success": True, "reaction_smiles": reaction_smiles},  # noqa: ARG005
    )

    result = _retrosynthesis_step(target_smiles="Target", include_conditions=True)

    assert result["success"] is True
    assert result["summary"]["fallback_used"] is False
    assert result["disconnections"]["disconnections"][0]["reaction_smiles"] == "ArBr.HNR2>>Target"
    assert result["top_disconnection"]["reaction_smiles"] == "ArBr.HNR2>>Target"
    assert result["top_reaction_smiles_evaluation"]["success"] is True
