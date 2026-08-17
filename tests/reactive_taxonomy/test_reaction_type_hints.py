"""Graph-derived reaction-type hint regressions."""

from reactive_taxonomy import (
    featurize_reaction,
    validate_reaction_type_hint_definitions,
)


SUZUKI = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"


def _primary(reaction_smiles: str):
    analysis = featurize_reaction(reaction_smiles)
    assert analysis.valid is True
    assert analysis.interpretation is not None
    hint_id = analysis.interpretation.primary_reaction_type_hint_id
    return next(
        hint
        for hint in analysis.interpretation.reaction_type_hints
        if hint.hint_id == hint_id
    )


def test_reaction_type_hint_definitions_are_valid() -> None:
    assert validate_reaction_type_hint_definitions() == []


def test_suzuki_hint_binds_roles_to_graph_sites() -> None:
    hint = _primary(SUZUKI)

    assert hint.type_id == "c_c_boron_transfer_coupling"
    assert hint.pattern_id == "organoboron_c_c_coupling_like"
    assert [(item.role, item.canonical_signature) for item in hint.participants] == [
        ("electrophile", "LG|Ar|Br"),
        ("transfer_partner", "TM|Ar|B(OH)2"),
    ]
    assert hint.warnings == ()


def test_hint_identity_is_invariant_to_reactant_order() -> None:
    forward = _primary(SUZUKI)
    reversed_order = _primary(
        "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )

    assert forward.hint_id == reversed_order.hint_id
    assert [item.role for item in forward.participants] == [
        item.role for item in reversed_order.participants
    ]


def test_ambiguous_sp2_c_n_hint_preserves_condition_dependence() -> None:
    hint = _primary("Brc1ccccc1.CN>>CNc1ccccc1")

    assert hint.type_id == "sp2_c_n_substitution"
    assert hint.requires_condition_evidence is True
    assert {item.role for item in hint.participants} == {
        "electrophile",
        "nitrogen_partner",
    }


def test_specific_sulfur_and_aromatic_ch_hints_precede_generic_coupling() -> None:
    sulfur = _primary("Brc1ccccc1.CS>>CSc1ccccc1")
    aromatic_ch = _primary(
        "Brc1ccccc1.c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )

    assert sulfur.type_id == "sp2_c_s_substitution"
    assert aromatic_ch.type_id == "sp2_c_aromatic_c_h_coupling"
    assert aromatic_ch.participants[1].canonical_signature == "CH|ArH"


def test_terminal_alkyne_hint_selects_the_hydrogen_bearing_site() -> None:
    hint = _primary("Brc1ccccc1.C#C>>C#Cc1ccccc1")

    assert hint.type_id == "c_c_terminal_alkyne_coupling"
    assert hint.participants[1].canonical_signature.startswith("XH|Csp|")
    assert hint.participants[1].role == "alkyne_partner"
