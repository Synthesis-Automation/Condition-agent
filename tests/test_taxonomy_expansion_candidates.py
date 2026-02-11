from __future__ import annotations

from chemtools.taxonomy.reaction_catalog import ReactionTypeDefinition, SlotRequirement
from scripts import suggest_taxonomy_expansion_candidates as sugg


def _make_def(
    rid: str,
    reactants: list[str],
    products: list[str],
) -> ReactionTypeDefinition:
    return ReactionTypeDefinition(
        id=rid,
        name=rid,
        category="test",
        aliases=[],
        description=None,
        reactants={"electrophile": SlotRequirement(allowed=reactants)},
        products={"product": SlotRequirement(allowed=products)},
        catalysts=[],
        conditions=None,
        metadata={},
        reference_reactions=[],
        notes=None,
        constraints={},
    )


def test_suggest_candidates_returns_ranked_matches(monkeypatch) -> None:
    defs = {
        "C_N_Coupling": _make_def("C_N_Coupling", ["Ar-I", "Ar-Br"], ["Ar-NR2"]),
        "C_O_Coupling": _make_def("C_O_Coupling", ["Ar-I", "Ar-Br"], ["Ar-OR"]),
    }
    monkeypatch.setattr(sugg, "load_reaction_catalog", lambda: (defs, {}))

    payload = {
        "clusters": [
            {
                "cluster_id": "k1",
                "count": 7,
                "motif_signature": "Ar-I -> Ar-NR2",
                "event_signature": "c_n_bond_formation",
                "source_reaction_labels": [["UnknownX", 7]],
            }
        ]
    }
    out = sugg.suggest_candidates(payload, top_candidates=2)
    assert out["input_cluster_count"] == 1
    cluster = out["clusters"][0]
    top = cluster["top_taxonomy_candidates"][0]
    assert top["reaction_id"] == "C_N_Coupling"
    assert top["score"] >= 0.8
