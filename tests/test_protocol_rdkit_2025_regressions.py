import pytest

from chemtools.util.rdkit_helpers import rdkit_available


pytestmark = pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")


def test_reaction_smiles_parser_and_smarts_matching_recommend():
    from chemtools.protocol.recommend import match_reaction_smarts

    reaction_smiles = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    patterns = ["[c,C,n,o,s]Br.OB(O)[c,C,n,o,s]>>[c,C,n,o,s][c,C,n,o,s]"]

    assert match_reaction_smarts(reaction_smiles, patterns) is True


def test_reaction_smiles_parser_and_smarts_matching_validate_protocols():
    from chemtools.protocol.validate_protocols import match_reaction_smarts

    reaction_smiles = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    patterns = ["[c,C,n,o,s]Br.OB(O)[c,C,n,o,s]>>[c,C,n,o,s][c,C,n,o,s]"]

    matched, errors = match_reaction_smarts(reaction_smiles, patterns)
    assert matched is True
    assert errors == []


def test_aromatic_smarts_requires_sanitize_for_reaction_components():
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions

    from chemtools.protocol.recommend import match_reaction_smarts

    reaction_smiles = "BrC1=CC=CC=C1.CC([Si](C)(C)C)=O>>CC(C2=CC=CC=C2)=O"
    aromatic_pattern = Chem.MolFromSmarts("[Br][c]")
    assert aromatic_pattern is not None

    rxn = rdChemReactions.ReactionFromSmiles(reaction_smiles)
    assert rxn is not None
    first_reactant = rxn.GetReactants()[0]

    # RDKit reaction components need sanitization to recover aromatic matching.
    assert first_reactant.HasSubstructMatch(aromatic_pattern) is False

    patterns = ["[Br][c].C[Si]>>C[c]"]
    assert match_reaction_smarts(reaction_smiles, patterns) is True
