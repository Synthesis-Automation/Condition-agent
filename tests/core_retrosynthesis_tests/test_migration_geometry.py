"""Graph and stereochemistry regressions for the cyclic oxime migration advisory."""

import pytest
from rdkit import Chem

from core_retrosynthesis.chemistry import canonical_smiles
from core_retrosynthesis.external_proposal_assessment import (
    ExternalRetrosynthesisProposal,
    assess_external_retrosynthesis_proposal,
)
from core_retrosynthesis.generic_models import GenericTemplateLibrary
from core_retrosynthesis.migration_geometry import assess_migration_geometry
from core_retrosynthesis.retro_validity import assess_retro_validity


OXIME = "[CH3:1][CH:2]1[CH2:3][CH2:4][CH2:5][CH2:6]{left}[C:7]1=[N:8]{right}[OH:9]"
LACTAM = "[CH3:1][CH:2]1[CH2:3][CH2:4][CH2:5][CH2:6][C:7](=[O:9])[NH:8]1"
OTHER = "[CH3:1][CH:2]1[CH2:3][CH2:4][CH2:5][CH2:6][NH:8][C:7]1=[O:9]"


@pytest.mark.parametrize("left,right,status", [
    ("", "", "unspecified"), ("/", "/", "syn_aligned"), ("/", "\\", "anti_aligned"),
])
def test_geometry_and_the_other_constitutional_product(left, right, status):
    precursor = OXIME.format(left=left, right=right)
    result = assess_migration_geometry(f"{precursor}>>{LACTAM}")
    assert result.status == status
    assert result.alternative_product_smiles == "CC1CCCCNC1=O"
    alternate = assess_migration_geometry(f"{precursor}>>{OTHER}")
    assert alternate.alternative_product_smiles == "CC1CCCCC(=O)N1"
    assert alternate.status == {
        "anti_aligned": "syn_aligned", "syn_aligned": "anti_aligned", "unspecified": "unspecified",
    }[status]


def test_stereo_atom_order_maps_and_oxygen_provenance():
    oxime = OXIME.format(left="/", right="/")
    reaction = f"{oxime}.[OH2:10]>>{LACTAM.replace('[O:9]', '[O:10]')}"
    expected = assess_migration_geometry(reaction)
    sides = []
    for side in reaction.split(">>"):
        mol = Chem.MolFromSmiles(side)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(30 - atom.GetAtomMapNum())
        reversed_mol = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
        sides.append(Chem.MolToSmiles(reversed_mol, canonical=False))
    assert expected.status == "syn_aligned"
    assert assess_migration_geometry(">>".join(sides)) == expected


@pytest.mark.parametrize("reaction", [
    "ON=C1CCCCC1>>O=C1CCCCCN1",
    "[OH:1][N:2]=[C:3]1[CH2:4][CH2:5][CH2:6][CH2:7][CH2:8]1>>"
    "[O:1]=[C:3]1[NH:2][CH2:4][CH2:5][CH2:6][CH2:7][CH2:8]1",
    OXIME.format(left="", right="") + ">>" + LACTAM.replace("NH:8", "NH:99"),
    OXIME.format(left="", right="").replace("CH2:3", "CH2:2") + ">>" + LACTAM,
    "C1CCC(=NO)CC1.CCO>>CC=O",
    "CC(=NO)c1ccccc1>>CC(=O)Nc1ccccc1",
    "invalid",
])
def test_equivalent_migrations_or_unsupported_correspondences(reaction):
    assert assess_migration_geometry(reaction) is None


def test_external_proposal_preserves_structural_status_and_geometry_evidence():
    oxime = OXIME.format(left="", right="")
    proposal = ExternalRetrosynthesisProposal(
        precursor_smiles=canonical_smiles(oxime), target_smiles=canonical_smiles(LACTAM),
        mapped_reaction_smiles=f"{oxime}>>{LACTAM}",
    )
    result = assess_external_retrosynthesis_proposal(proposal, GenericTemplateLibrary((), 0, 0, {}, {}))
    assert result.status == "structurally_verified"
    assert result.to_dict()["compatibility_evidence"]["migration_geometry"]["status"] == "unspecified"


@pytest.mark.parametrize("left,right,expect_caution", [("", "", True), ("/", "/", True), ("/", "\\", False)])
def test_validity_cautions_follow_the_drawn_geometry(left, right, expect_caution):
    oxime = OXIME.format(left=left, right=right)
    proposal = ExternalRetrosynthesisProposal(
        precursor_smiles=canonical_smiles(oxime), target_smiles=canonical_smiles(LACTAM),
        mapped_reaction_smiles=f"{oxime}>>{LACTAM}",
    )
    result = assess_retro_validity(proposal, GenericTemplateLibrary((), 0, 0, {}, {}))
    geometry = result.structural_assessment.compatibility_evidence.migration_geometry
    assert (geometry.summary in result.cautions) is expect_caution
