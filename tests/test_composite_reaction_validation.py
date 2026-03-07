from chem_coworker.tools.composite import (
    _forward_synthesis_step,
    _recommend_reaction_conditions,
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
