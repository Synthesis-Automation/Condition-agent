"""Cross-step and symmetry controls for projection of retained constructed bonds."""

from copy import deepcopy

import pytest

from core_retrosynthesis.chemistry import canonical_smiles
from core_retrosynthesis.external_route_admission import ExternalRouteProposal, assess_external_route_proposal
from core_retrosynthesis.generic_library import build_generic_library
from core_retrosynthesis.route_constructions import construction_inventory


METHYLATION = "[CH3:1][CH2:2][OH:3].[CH3:4]Br>>[CH3:1][CH2:2][O:3][CH3:4]"
ETHYLATION = "[CH3:1][CH2:2]Br.[OH:3][CH3:4]>>[CH3:1][CH2:2][O:3][CH3:4]"
ETHANOL = "[CH3:8][CH2:5]Cl.[OH2:7]>>[CH3:8][CH2:5][OH:7]"
PROPANE = "[CH3:8][CH2:3]Br.[CH3:7]I>>[CH3:8][CH2:3][CH3:7]"
PROPANOL = "[CH3:1][CH2:2][CH3:3].[OH2:4]>>[CH3:1][CH2:2][CH2:3][OH:4]"


@pytest.fixture(scope="module")
def library():
    return build_generic_library(())


def _record(library, *reactions):
    steps = []
    for step_id, rxn in reactions:
        left, right = rxn.split(">>")
        steps.append({
            "external_step_id": step_id, "precursor_smiles": canonical_smiles(left),
            "target_smiles": canonical_smiles(right), "mapped_reaction_smiles": rxn,
        })
    route = ExternalRouteProposal.from_dict({"target_smiles": steps[-1]["target_smiles"], "steps": steps})
    return assess_external_route_proposal(route, library).to_dict()


def test_same_target_site_with_different_handles_and_different_depths(library):
    route = _record(library, ("ethyl", ETHANOL), ("methyl", METHYLATION))
    inventory = construction_inventory(route)
    assert inventory["unresolved_construction_step_ids"] == []
    assert set(inventory["leaf_smiles"]) == {"CCCl", "O", "CBr"}
    cuts = {row["step_id"]: row["bond_id"] for row in inventory["key_constructions"]}
    other = construction_inventory(_record(library, ("direct", ETHYLATION)))
    assert cuts["ethyl"] == other["key_constructions"][0]["bond_id"]
    assert cuts["methyl"] != cuts["ethyl"]
    iodomethane = construction_inventory(_record(library, ("iodide", METHYLATION.replace("Br", "I"))))
    assert cuts["methyl"] == iodomethane["key_constructions"][0]["bond_id"]


def test_symmetry_is_reported_and_ambiguous_downstream_correspondence_abstains(library):
    one_step = construction_inventory(_record(library, ("assemble", PROPANE)))
    assert one_step["key_constructions"][0]["equivalent_target_atom_pairs"] == [[0, 1], [1, 2]]
    later = construction_inventory(_record(library, ("assemble", PROPANE), ("differentiate", PROPANOL)))
    assert later["unresolved_construction_step_ids"] == ["assemble"]
    assert [c["step_id"] for c in later["key_constructions"]] == ["differentiate"]


def test_lost_bonds_and_order_changes_do_not_count_as_constructions(library):
    reduction = "[CH3:1][CH:2]=[O:3]>>[CH3:1][CH2:2][OH:3]"
    assert construction_inventory(_record(library, ("reduce", reduction)))["key_constructions"] == []
    # The intermediate's C-O bond is actually cleaved as it is transformed into a chloride.
    replace_oxygen = "[CH3:1][CH2:2][OH:3].[ClH:4]>>[CH3:1][CH2:2][Cl:4]"
    record = _record(library, ("make_ether_bond", ETHANOL.replace("Cl", "Br")),
                     ("substitute", replace_oxygen))
    result = construction_inventory(record)
    assert result["unresolved_construction_step_ids"] == []
    assert [c["step_id"] for c in result["key_constructions"]] == ["substitute"]


def test_missing_legacy_maps_and_contradictory_saved_sides_remain_unresolved(library):
    record = _record(library, ("ethyl", ETHANOL), ("methyl", METHYLATION))
    legacy = deepcopy(record)
    del legacy["step_assessments"][1]["assessment"]["normalized_mapped_reaction_smiles"]
    contradictory = deepcopy(record)
    contradictory["step_assessments"][1]["assessment"]["normalized_mapped_reaction_smiles"] = ETHANOL
    for value in (legacy, contradictory):
        result = construction_inventory(value)
        assert result["key_constructions"] == []
        assert result["unresolved_construction_step_ids"] == ["ethyl", "methyl"]
