"""Scoped parser reuse preserves atoms and cannot leak molecule mutations."""

from dataclasses import asdict

import pytest

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.chemistry.rdkit_utils import MoleculeParseCache, molecule_parse_scope, parse_smiles


def test_cached_molecules_are_fresh_copies_and_scope_is_reset():
    cache = MoleculeParseCache()
    with molecule_parse_scope(cache):
        mol = parse_smiles("[13CH3:7][C@H:8](O)Cl")
        mol.GetAtomWithIdx(0).SetAtomMapNum(99)
        another = parse_smiles("[13CH3:7][C@H:8](O)Cl")
        assert another.GetAtomWithIdx(0).GetAtomMapNum() == 7
        assert another.GetAtomWithIdx(0).GetIsotope() == 13
        assert parse_smiles("not a molecule") is None
    before = len(cache.entries)
    parse_smiles("CCC")
    assert len(cache.entries) == before


@pytest.mark.parametrize("reaction", [
    "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
    "[CH3:1][Br:2].[OH:3][CH3:4]>>[CH3:1][O:3][CH3:4]",
    "[CH3:1][Br:2].[NH2:3][CH3:4]>>[CH3:1][NH:3][CH3:4]",
    "[CH3:1][Br:2].[SH:3][CH3:4]>>[CH3:1][S:3][CH3:4]",
    "[CH3:1][OH:2]>>[CH2:1]=[O:2]",
    "[CH3:1][OH:1]>>[CH2:1]=[O:2]",
    "CC.O>>CCO", "invalid>>CO",
])
def test_cached_and_uncached_chemistry_are_identical(reaction):
    baseline = asdict(featurize_reaction(reaction))
    with molecule_parse_scope(MoleculeParseCache()):
        assert asdict(featurize_reaction(reaction)) == baseline
