import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pytest

from chemtools.reagent import (
    build_reaction_lookup,
    build_reactant_lookup,
    classify_reactant_smiles,
    normalize_reaction_type,
    required_reactant_categories,
)


def test_build_reactant_lookup_contains_core_alias():
    alias_map, id_to_category = build_reactant_lookup()
    assert alias_map["arbr"] == "ArBr"
    assert id_to_category["ArBr"] == "ArX*"


def test_build_reaction_lookup_aliases():
    id_to_meta, alias_map = build_reaction_lookup()
    assert alias_map["suzuki"] == "suzuki_miyaura"
    assert alias_map["buchwald-hartwig"] == "c_n_cross_coupling"
    assert id_to_meta["suzuki_miyaura"]["category"] == "c_c_coupling_reactions"


def test_normalize_reaction_type_alias():
    assert normalize_reaction_type("Suzuki-Miyaura, in situ") is None


def test_required_reactant_categories_structure():
    required = required_reactant_categories("suzuki_miyaura")
    assert required is not None
    assert "electrophiles" in required
    assert "nucleophiles" in required
    assert required["electrophiles"][0].startswith("Ar-")


@pytest.mark.parametrize("smiles", ["c1ccccc1Br", "CCBr"])
def test_classify_reactant_smiles_rdkit_disabled(monkeypatch, smiles):
    monkeypatch.setenv("CHEMTOOLS_DISABLE_RDKIT", "1")
    try:
        assert classify_reactant_smiles(smiles) is None
    finally:
        monkeypatch.delenv("CHEMTOOLS_DISABLE_RDKIT", raising=False)
