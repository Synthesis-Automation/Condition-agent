"""Selected-bond evidence never guesses a target site or atom correspondence."""

from dataclasses import replace

import pytest

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.disconnection_focus import (
    check_disconnection_bond, disconnection_focus_policy, prepare_disconnection_focus,
)


REACTION = "[CH3:1][CH2:2]Br.[NH3:3]>>[CH3:1][CH2:2][NH2:3]"


def test_selected_bond_has_typed_forward_formation_witness():
    focus = prepare_disconnection_focus("NCC", [2, 1], "CCN")
    assert focus.atom_ids == (1, 2)
    assert focus == prepare_disconnection_focus("CCN", [1, 2], "NCC")
    assert disconnection_focus_policy()["required_forward_edit"] == "formed"
    assert check_disconnection_bond(focus, REACTION, (1, 2, 3)).status == "matched"
    check = check_disconnection_bond(focus, REACTION, (1, 2, 3),
                                   observation=featurize_reaction(REACTION).observation)
    assert check.status == "verified"
    assert check.target_atom_maps == (2, 3)
    assert check.formed_edit.edit_type == "formed"


def test_unselected_bond_and_order_change_do_not_qualify():
    focus = prepare_disconnection_focus("CCN", [0, 1], "CCN")
    assert check_disconnection_bond(focus, REACTION, (1, 2, 3)).status == "not_matched"
    focus = prepare_disconnection_focus("CO", [0, 1], "CO")
    reduction = "[CH2:1]=[O:2]>>[CH3:1][OH:2]"
    assert check_disconnection_bond(focus, reduction, (1, 2)).status == "not_matched"


@pytest.mark.parametrize("bond,reference", [
    ([0, 0], "CCN"), ([0, 2], "CCN"), ([0, 8], "CCN"), ([-1, 1], "CCN"),
    ([True, 2], "CCN"), ([1], "CCN"), ([1, 2], "CCO"), ([1, 2], None),
])
def test_invalid_or_stale_target_selection_is_rejected(bond, reference):
    with pytest.raises(ValueError):
        prepare_disconnection_focus("CCN", bond, reference)


@pytest.mark.parametrize("reaction,maps", [
    ("[CH3:1][CH2:2]Br.[NH3:2]>>[CH3:1][CH2:2][NH2:3]", (1, 2, 3)),
    ("[CH3:1][CH2:2]Br>>[CH3:1][CH2:2][NH2:3]", (1, 2, 3)),
    (REACTION, (1, 3, 2)), (REACTION, (1, 2, 2)),
])
def test_invalid_correspondence_never_passes(reaction, maps):
    focus = prepare_disconnection_focus("CCN", [1, 2], "CCN")
    assert check_disconnection_bond(focus, reaction, maps).status == "unresolved"


@pytest.mark.parametrize("changes", [
    {"evidence_quality": "ambiguous"}, {"warnings": ("MAPPED_OPERATOR_CONFLICT",)},
    {"input_reaction_smiles": "CCBr.N>>CCN"},
])
def test_final_ambiguity_or_conflict_overrides_preliminary_match(changes):
    focus = prepare_disconnection_focus("CCN", [1, 2], "CCN")
    observation = replace(featurize_reaction(REACTION).observation, **changes)
    result = check_disconnection_bond(focus, REACTION, (1, 2, 3), observation=observation)
    assert result.status == "unresolved"


def test_intramolecular_ring_closure_does_not_require_two_precursors():
    reaction = "[NH2:1][CH2:2][CH2:3][CH2:4][C:5](=[O:6])O>>[NH:1]1[CH2:2][CH2:3][CH2:4][C:5]1=[O:6]"
    focus = prepare_disconnection_focus("O=C1CCCN1", [1, 5], "O=C1CCCN1")
    assert check_disconnection_bond(focus, reaction, (6, 5, 4, 3, 2, 1)).status == "matched"


def test_target_stereo_and_isotopes_cannot_be_changed_by_the_binding():
    focus = prepare_disconnection_focus("C[C@H](F)N", [1, 3], "C[C@H](F)N")
    reaction = "[CH3:1][C@@H:2]([F:3])Br.[NH3:4]>>[CH3:1][C@@H:2]([F:3])[NH2:4]"
    assert check_disconnection_bond(focus, reaction, (1, 2, 3, 4)).status == "unresolved"
    isotopic = prepare_disconnection_focus("[13CH3]N", [0, 1], "[13CH3]N")
    unlabelled = "[CH3:1]Br.[NH3:2]>>[CH3:1][NH2:2]"
    assert check_disconnection_bond(isotopic, unlabelled, (1, 2)).status == "unresolved"
