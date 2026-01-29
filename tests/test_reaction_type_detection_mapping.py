import pytest

from chemtools.featurizers import detection as rd
from chemtools.taxonomy.reaction_catalog import ReactionTypeDefinition


def _make_definition(reaction_id: str, reference: str) -> ReactionTypeDefinition:
    return ReactionTypeDefinition(
        id=reaction_id,
        name=reaction_id,
        category="test",
        aliases=[],
        description=None,
        reactants={},
        products={},
        catalysts=[],
        conditions=None,
        metadata={},
        reference_reactions=[reference],
        notes=None,
        constraints={},
    )


def test_bond_change_detection_maps_cc_coupling(monkeypatch: pytest.MonkeyPatch) -> None:
    reference_cc = "[CH3:1][Br].[CH3:2]B(O)O>>[CH3:1][CH3:2]"
    reference_cn = "[CH3:1][Cl].[NH3:2]>>[CH3:1][NH2:2]"

    definitions = {
        "mock_cc": _make_definition("mock_cc", reference_cc),
        "mock_cn": _make_definition("mock_cn", reference_cn),
    }

    def fake_load_reaction_catalog(path=None):
        return definitions, {key: key for key in definitions}

    monkeypatch.setattr(rd, "load_reaction_catalog", fake_load_reaction_catalog)
    rd.clear_bond_change_cache()

    result = rd.detect_reaction_types(reference_cc, use_bond_changes=True)
    assert result.matches
    assert result.matches[0].reaction_type == "mock_cc"
    assert result.matches[0].slot_evidence.get("bond_changes")


def test_bond_change_detection_maps_cn_coupling(monkeypatch: pytest.MonkeyPatch) -> None:
    reference_cc = "[CH3:1][Br].[CH3:2]B(O)O>>[CH3:1][CH3:2]"
    reference_cn = "[CH3:1][Cl].[NH3:2]>>[CH3:1][NH2:2]"

    definitions = {
        "mock_cc": _make_definition("mock_cc", reference_cc),
        "mock_cn": _make_definition("mock_cn", reference_cn),
    }

    def fake_load_reaction_catalog(path=None):
        return definitions, {key: key for key in definitions}

    monkeypatch.setattr(rd, "load_reaction_catalog", fake_load_reaction_catalog)
    rd.clear_bond_change_cache()

    result = rd.detect_reaction_types(reference_cn, use_bond_changes=True)
    assert result.matches
    assert result.matches[0].reaction_type == "mock_cn"
    assert result.matches[0].slot_evidence.get("bond_changes")
