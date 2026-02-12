from __future__ import annotations

import pytest

from chemtools.featurizers.formatters import detection_validation as dv
from chemtools.taxonomy.reaction_catalog import ReactionTypeDefinition, SlotRequirement


def _mk_defn(
    rid: str,
    *,
    reactants: dict[str, list[str]],
    products: dict[str, list[str]],
) -> ReactionTypeDefinition:
    return ReactionTypeDefinition(
        id=rid,
        name=rid,
        category="test",
        aliases=[],
        description=None,
        reactants={k: SlotRequirement(allowed=v) for k, v in reactants.items()},
        products={k: SlotRequirement(allowed=v) for k, v in products.items()},
        catalysts=[],
        conditions=None,
        metadata={},
        reference_reactions=[],
        notes=None,
        constraints={
            "include_reacted": [],
            "exclude_reacted": [],
            "include_formed": [],
            "exclude_formed": [],
            "include_bond_formed": [],
            "exclude_bond_formed": [],
            "include_bond_broken": [],
            "exclude_bond_broken": [],
            "min_reactant_slot_matches": 0,
            "min_product_slot_matches": 0,
        },
    )


def test_validate_detection_with_crk_key_includes_evidence_payload() -> None:
    result = dv.validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key="|Ar-NH2 -> Ar-Cl",
    )
    evidence = result.get("evidence")
    assert isinstance(evidence, dict)
    assert evidence.get("matcher") == "taxonomy_specificity_v2"
    assert isinstance(evidence.get("top_candidates"), list)


def test_specificity_matcher_can_outrank_legacy_first_hit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    broad = _mk_defn(
        "Broad_class",
        reactants={"substrate": ["Ar-NH2"]},
        products={"product": ["HeteroAr-Cl", "Pyrimidine"]},
    )
    specific = _mk_defn(
        "Specific_class",
        reactants={
            "amine": ["Ar-NH2"],
            "amide": ["Ar-NHCOR"],
            "c1_source": ["CH3-OH"],
        },
        products={"heterocycle": ["Pyrimidine"]},
    )
    definitions = {"Broad_class": broad, "Specific_class": specific}

    monkeypatch.setattr(dv, "_get_catalog", lambda: (definitions, {}))
    crk = "|Ar-NH2|Ar-NHCOR|CH3-OH -> HeteroAr-Cl|Pyrimidine"

    legacy = dv.validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key=crk,
        use_legacy=True,
    )
    improved = dv.validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key=crk,
        use_legacy=False,
    )

    assert legacy["reaction_type"] == "Broad_class"
    assert improved["reaction_type"] == "Specific_class"
    assert improved["evidence"]["matcher"] == "taxonomy_specificity_v2"


def test_crk_validation_can_use_bond_change_constraints(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    broad = _mk_defn(
        "EventLike",
        reactants={"substrate": ["Alkenyl-CHO"]},
        products={"product": ["Alkenyl-Alkenyl"]},
    )
    wittig_like = _mk_defn(
        "Wittig_like",
        reactants={"substrate": ["Alkenyl-CHO"]},
        products={"product": ["Alkenyl-Alkenyl"]},
    )
    broad.constraints["exclude_bond_broken"] = ["C-P"]
    wittig_like.constraints["include_bond_broken"] = ["C-P"]
    definitions = {"EventLike": broad, "Wittig_like": wittig_like}

    monkeypatch.setattr(dv, "_get_catalog", lambda: (definitions, {}))
    crk = "|Alkenyl-CHO -> Alkenyl-Alkenyl | bond_formed: C-C | bond_broken: C-P; C-O"
    result = dv.validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key=crk,
        use_legacy=False,
    )
    assert result["reaction_type"] == "Wittig_like"
