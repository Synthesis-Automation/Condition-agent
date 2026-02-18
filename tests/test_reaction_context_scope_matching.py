from __future__ import annotations

import pytest

from chemtools.featurizers.analysis import reaction_context as rc
from chemtools.taxonomy import reaction_catalog
from chemtools.taxonomy.reaction_catalog import ReactionTypeDefinition, SlotRequirement


def test_infer_reactant_role_uses_scope_compatibility(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    definition = ReactionTypeDefinition(
        id="Epoxidation_like",
        name="Epoxidation_like",
        category="test",
        aliases=[],
        description=None,
        reactants={"substrate": SlotRequirement(allowed=["Alkyl-Alkenyl"])},
        products={"product": SlotRequirement(allowed=["Epoxide"])},
        catalysts=[],
        conditions=None,
        metadata={},
        reference_reactions=[],
        notes=None,
        constraints={},
    )

    monkeypatch.setattr(
        reaction_catalog,
        "load_reaction_catalog",
        lambda: ({"Epoxidation_like": definition}, {}),
    )

    role = rc._infer_reactant_role("RCH2-Alkenyl", "Epoxidation_like")
    assert role == "substrate"
