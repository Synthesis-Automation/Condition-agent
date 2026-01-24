import pytest

from chemtools.featurizers import reaction_detection as rd
from chemtools.util.rdkit_helpers import rdkit_available


def test_detects_amide_formation_from_benzoic_acid() -> None:
    if not rdkit_available():
        pytest.skip("rdkit not available")

    rxn = "O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1"
    result = rd.detect_reaction_types(rxn)
    assert any(match.reaction_type == "amide_formation" for match in result.matches)
