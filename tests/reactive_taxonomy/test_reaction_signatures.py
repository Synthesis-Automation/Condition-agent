"""Generic reaction-signature regressions for the graph-only pipeline."""

import pytest

from reactive_taxonomy import featurize_reaction


@pytest.mark.parametrize(
    "reaction",
    [
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        "Brc1ccccc1.CN>>CNc1ccccc1",
        "Brc1ccccc1.CO>>COc1ccccc1",
        "Brc1ccccc1.CS>>CSc1ccccc1",
    ],
)
def test_common_edit_chemistries_receive_generic_signatures(reaction: str) -> None:
    result = featurize_reaction(reaction)
    assert result.valid
    assert result.observation is not None
    assert result.observation.edits
    assert result.reaction_signature is not None
    assert result.reaction_signature.named_family is None


def test_partner_order_does_not_change_signature() -> None:
    forward = featurize_reaction(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )
    reverse = featurize_reaction(
        "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )
    assert forward.reaction_signature is not None
    assert reverse.reaction_signature is not None
    assert forward.reaction_signature.signature_id == reverse.reaction_signature.signature_id


def test_mapped_unknown_family_receives_deterministic_signature() -> None:
    reaction = "[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]"
    first = featurize_reaction(reaction)
    second = featurize_reaction(reaction)
    assert first.reaction_signature is not None
    assert first.named_family is None
    assert first.reaction_signature.signature_id == second.reaction_signature.signature_id
    assert first.reaction_signature.signature_id.startswith("RS3:")


def test_invalid_mapping_is_not_promoted_to_a_signature() -> None:
    result = featurize_reaction("[CH3:1][OH:1]>>[CH3:1]Cl")
    assert result.valid
    assert result.evidence_quality == "unresolved"
    assert result.reaction_signature is None


def test_signature_uses_versioned_hierarchical_keys() -> None:
    signature = featurize_reaction("CC=CC>>CCCC").reaction_signature
    assert signature is not None
    assert signature.exact_signature_key.startswith("L0:")
    assert signature.handle_signature_key.startswith("L1:")
    assert signature.transformation_signature_key.startswith("L2:")
    assert signature.bond_edit_signature_key.startswith("L3:")
    assert signature.environment_signature_key.startswith("L4:")
