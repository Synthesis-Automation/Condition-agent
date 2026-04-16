from __future__ import annotations

from types import SimpleNamespace

from chem_coworker.tools.composite import _recommend_reaction_conditions
from chem_coworker.tools.conditions import (
    _build_condition_reaction_context,
    _compose_condition_candidates,
    _get_literature_condition_evidence,
)


_VALID_RXN = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"


def test_build_condition_reaction_context_returns_atomic_context(monkeypatch) -> None:
    monkeypatch.setattr(
        "chemtools.recommend.query_analysis.analyze_recommendation_query",
        lambda request: SimpleNamespace(
            reaction_smiles_normalized=_VALID_RXN,
            reactants=("Brc1ccccc1", "OB(O)c1ccccc1"),
            agents=(),
            products=("c1ccc(-c2ccccc2)cc1",),
            reactant_a_smiles="Brc1ccccc1",
            reactant_b_smiles="OB(O)c1ccccc1",
            product_smiles="c1ccc(-c2ccccc2)cc1",
            detected_reaction_type="suzuki_miyaura",
            detected_reaction_type_id="suzuki_miyaura",
            detected_reaction_type_name="Suzuki-Miyaura",
            detected_reaction_type_category="cross_coupling",
            reaction_type_confidence=0.91,
            reaction_key="rk-demo",
            reacted_motifs=("aryl_halide", "aryl_boron"),
            formed_motifs=("biaryl",),
            spectator_motifs=("ether",),
            spectator_groups=("aryl_ether",),
            requested_reaction_type_filter=None,
            requested_reaction_type_filter_canonical=None,
            raw_feature_summary={"reacted_motif_count": 2},
            warnings=(),
        ),
    )

    result = _build_condition_reaction_context(_VALID_RXN)

    assert result["success"] is True
    assert result["reaction_type"] == "suzuki_miyaura"
    assert result["reaction_key"] == "rk-demo"
    assert result["reactant_a_smiles"] == "Brc1ccccc1"
    assert result["reactant_b_smiles"] == "OB(O)c1ccccc1"
    assert result["formed_motifs"] == ["biaryl"]


def test_source_specific_condition_evidence_returns_common_schema(monkeypatch) -> None:
    monkeypatch.setattr(
        "chem_coworker.tools.conditions._recommend_conditions",
        lambda **kwargs: {
            "success": True,
            "reaction_smiles": kwargs["reaction_smiles"],
            "recommendations": [
                {
                    "rank": 1,
                    "catalyst": "Pd(OAc)2",
                    "ligand": "SPhos",
                    "base": "K3PO4",
                    "solvent": "dioxane",
                    "temperature": "90 C",
                    "confidence": 0.88,
                    "success_rate": 0.81,
                    "avg_yield": 82.0,
                    "median_yield": 84.0,
                    "match_score": 0.93,
                    "num_experiments": 6,
                    "reaction_type": "Suzuki_miyaura",
                    "reaction_category": "cross_coupling",
                    "reactant_types": ["aryl_halide", "aryl_boron"],
                    "precedent_ids": "lit-1",
                }
            ],
            "detected_reaction_type": "Suzuki_miyaura",
            "reaction_type_confidence": 0.88,
            "evidence": {
                "matched_transformation": "RXNEVT|form=C-C",
                "database_coverage": 0.23,
            },
            "_warnings": ["Low support warning"],
            "selection_mode": kwargs.get("selection_mode", "best"),
            "total_available": 1,
            "hte_timing_ms": {"total_ms": 10.0},
            "hte_processing_time_ms": 10.0,
            "hte_recommender_stage_timing_ms": {"match_retrieval_ms": 7.0},
        },
    )

    result = _get_literature_condition_evidence(_VALID_RXN, top_k=3)

    assert result["success"] is True
    assert result["source"] == "literature"
    assert len(result["condition_evidence"]) == 1
    hit = result["condition_evidence"][0]
    assert hit["source"] == "literature"
    assert hit["proposed_fragments"]["catalyst"] == "Pd(OAc)2"
    assert hit["matched_features"]["matched_transformation"] == "RXNEVT|form=C-C"
    assert hit["provenance"]["precedent_ids"] == "lit-1"
    assert "Low support warning" in result["_warnings"]


def test_compose_condition_candidates_merges_atomic_source_outputs(monkeypatch) -> None:
    def _fake_get_condition_evidence(reaction_smiles: str, *, source: str, **kwargs):
        assert reaction_smiles == _VALID_RXN
        base = {
            "success": True,
            "reaction_smiles": reaction_smiles,
            "source": source,
            "condition_evidence": [{"source": source, "score": 0.9}],
            "recommendations": [
                {
                    "rank": 1,
                    "catalyst": "Pd(OAc)2",
                    "ligand": "SPhos",
                    "base": "K3PO4",
                    "solvent": "dioxane",
                    "confidence": 0.9 if source == "literature" else 0.8,
                    "success_rate": 0.85,
                    "avg_yield": 86.0,
                    "median_yield": 87.0,
                    "match_score": 0.95,
                    "num_experiments": 6 if source == "literature" else 4,
                    "source": source,
                }
            ],
        }
        return base

    monkeypatch.setattr(
        "chem_coworker.tools.conditions._get_condition_evidence",
        _fake_get_condition_evidence,
    )

    result = _compose_condition_candidates(
        _VALID_RXN,
        top_k=3,
        sources=["literature", "motif"],
    )

    assert result["success"] is True
    assert len(result["recommendations"]) == 1
    top = result["recommendations"][0]
    assert top["consensus_count"] == 2
    assert set(top["source_consensus"]) == {"literature", "motif"}
    assert result["source_counts"]["literature"] == 1
    assert result["source_counts"]["motif"] == 1


def test_legacy_full_strategy_delegates_to_atomic_composer(monkeypatch) -> None:
    seen: dict[str, object] = {}

    def _fake_compose_condition_candidates(reaction_smiles: str, top_k: int, sources, selection_mode: str):
        seen["reaction_smiles"] = reaction_smiles
        seen["top_k"] = top_k
        seen["sources"] = list(sources)
        seen["selection_mode"] = selection_mode
        return {
            "success": True,
            "reaction_smiles": reaction_smiles,
            "recommendations": [],
            "source_counts": {},
            "consensus_summary": {"unique_condition_sets": 0, "failed_sources": [], "sources_run": list(sources)},
        }

    monkeypatch.setattr(
        "chem_coworker.tools.conditions._compose_condition_candidates",
        _fake_compose_condition_candidates,
    )

    result = _recommend_reaction_conditions(
        reaction_smiles=_VALID_RXN,
        top_k=4,
        condition_strategy="full",
        condition_selection_mode="diverse",
    )

    assert result["success"] is True
    assert seen["reaction_smiles"] == _VALID_RXN
    assert seen["top_k"] == 4
    assert seen["sources"] == ["literature", "motif", "similarity", "rules"]
    assert seen["selection_mode"] == "diverse"
    assert result["condition_strategy"] == "full"
