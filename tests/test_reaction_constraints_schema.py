from __future__ import annotations

from chemtools.taxonomy.reaction_catalog import (
    REACTION_CONSTRAINT_KEYS,
    load_reaction_catalog,
    normalize_reaction_constraints,
)


def test_normalize_reaction_constraints_provides_full_schema() -> None:
    normalized = normalize_reaction_constraints(None)
    for key in REACTION_CONSTRAINT_KEYS:
        assert key in normalized
    assert normalized["include_reacted"] == []
    assert normalized["exclude_reacted"] == []
    assert normalized["include_formed"] == []
    assert normalized["exclude_formed"] == []
    assert normalized["include_bond_formed"] == []
    assert normalized["exclude_bond_formed"] == []
    assert normalized["include_bond_broken"] == []
    assert normalized["exclude_bond_broken"] == []
    assert normalized["min_reactant_slot_matches"] == 0
    assert normalized["min_product_slot_matches"] == 0


def test_loaded_reaction_definitions_have_normalized_constraints() -> None:
    definitions, _ = load_reaction_catalog()
    assert definitions
    required = set(REACTION_CONSTRAINT_KEYS)
    for defn in definitions.values():
        assert required.issubset(set(defn.constraints.keys()))


def test_loaded_reaction_definitions_expose_redox_neutral_flags() -> None:
    definitions, _ = load_reaction_catalog()
    assert definitions["Suzuki_miyaura"].redox_neutral is True
    assert definitions["Reductive_amination"].redox_neutral is False
