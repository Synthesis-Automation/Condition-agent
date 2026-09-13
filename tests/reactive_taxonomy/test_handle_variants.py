"""Positive, ambiguous and conflicting evidence for virtual handle lookups."""

from copy import deepcopy
from dataclasses import asdict

import pytest

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.handle_variants import aromatic_leaving_handle_variant


def variant(smiles, signature=None):
    if signature is None:
        analysis = featurize_reaction(smiles)
        signature = asdict(analysis.reaction_signature) if analysis.reaction_signature else {}
    return aromatic_leaving_handle_variant(
        smiles, signature, query_element="I", precedent_element="Br",
        formed_partner_elements=("C", "N", "O", "S"),
    )


@pytest.mark.parametrize("smiles", [
    "Ic1ccccc1.CN>>CNc1ccccc1", "CN.Ic1ccccc1>>CNc1ccccc1",
    "Ic1ccccc1.CO>>COc1ccccc1", "Ic1ccccc1.CS>>CSc1ccccc1",
    "Ic1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
    "Ic1ccc(Br)cc1.CN>>CNc1ccc(Br)cc1",
    "[I:1][c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1.[CH3:8][NH2:9]>>[CH3:8][NH:9][c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1",
])
def test_one_observed_aromatic_leaving_handle_generates_search_hypothesis(smiles):
    result = variant(smiles)
    assert result is not None
    assert result.query_element == "I" and result.precedent_element == "Br"
    assert result.reaction_smiles.split(">")[2] == smiles.split(">")[2]
    assert featurize_reaction(result.reaction_smiles).reaction_signature is not None


@pytest.mark.parametrize("smiles", [
    "CCI.CN>>CCNC", "Fc1ccccc1.CN>>CNc1ccccc1",
    "Ic1ccc(Br)cc1.CN>>CNc1ccc(I)cc1",
    "Ic1ccc(I)cc1.CN.CN>>CNc1ccc(NC)cc1",
    "invalid>>invalid",
])
def test_retained_unsupported_ambiguous_or_multiple_handles_do_not_generate_variants(smiles):
    assert variant(smiles) is None


def test_conflicting_atom_reference_is_rejected_without_mutation():
    smiles = "Ic1ccccc1.CN>>CNc1ccccc1"
    signature = asdict(featurize_reaction(smiles).reaction_signature)
    signature["edits"][0]["atom_1"]["atom_index"] = 2
    before = deepcopy(signature)
    assert variant(smiles, signature) is None
    assert signature == before


def test_ambiguous_duplicate_broken_edit_is_rejected():
    smiles = "Ic1ccccc1.CN>>CNc1ccccc1"
    signature = asdict(featurize_reaction(smiles).reaction_signature)
    signature["edits"] = (*signature["edits"], signature["edits"][0])
    assert variant(smiles, signature) is None


def test_conflicting_map_reference_is_rejected():
    smiles = "Ic1ccccc1.CN>>CNc1ccccc1"
    signature = asdict(featurize_reaction(smiles).reaction_signature)
    signature["edits"][0]["atom_1"]["atom_map_number"] = 999
    assert variant(smiles, signature) is None
