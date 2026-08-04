"""Regression tests for exploratory structural-analogue discovery."""

from dataclasses import asdict

from reactive_taxonomy import featurize_reaction

from condition_recommender import ReactionDiscoveryExplorer, load_discovery_rules
from condition_recommender.generic_indexing import build_generic_index


REACTION = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"


def _record(index: int, *, yield_pct: float | None) -> dict:
    analysis = featurize_reaction(REACTION)
    assert analysis.valid and analysis.reaction_signature is not None
    recipe_id = f"RCR1:discovery-{index}"
    recipe_core_id = f"RCORE1:discovery-{index}"
    return {
        "schema_version": "10.0",
        "converter_definition_version": "generic_conversion.v10.0",
        "admission_tier": "verified",
        "index_eligibility": "eligible",
        "precedent_tier": "trusted",
        "core_eligibility": "trusted_core",
        "core_eligibility_definition_version": "core_eligibility.v1@1.0",
        "chemistry_status": "verified",
        "condition_status": "resolved_complete",
        "condition_stage_status": "single_stage",
        "outcome_status": "usable" if yield_pct is not None else "missing",
        "reaction_id": f"discovery-reaction-{index}",
        "observation_id": f"discovery-observation-{index}",
        "canonical_reaction_id": f"discovery-canonical-{index}",
        "reaction_smiles": REACTION,
        "reaction_label": asdict(analysis.reaction_label),
        "yield_pct": yield_pct,
        "source_dataset": "discovery-test",
        "reference_id": f"REF1:discovery-{index}",
        "reaction_signature": asdict(analysis.reaction_signature),
        "reaction_core": (
            asdict(analysis.reaction_core)
            if analysis.reaction_core is not None
            else None
        ),
        "reference_condition_series_id": f"RCS1:discovery-{index}",
        "resolved_recipe_id": recipe_id,
        "resolved_recipe_core_id": recipe_core_id,
        "resolved_recipe": {
            "recipe_id": recipe_id,
            "recipe_core_id": recipe_core_id,
            "catalysts": [{"canonical_name": f"Catalyst {index}"}],
            "temperature_c": 40.0 + index,
        },
        "condition_resolution": {"has_uncertainty": False},
    }


def _explorer() -> ReactionDiscoveryExplorer:
    return ReactionDiscoveryExplorer(
        build_generic_index(
            [
                _record(1, yield_pct=88.0),
                _record(2, yield_pct=5.0),
                _record(3, yield_pct=None),
            ]
        )
    )


def test_discovery_returns_individual_precedents_and_observed_conditions() -> None:
    result = _explorer().discover(REACTION, top_k=3)

    assert result.valid
    assert len(result.hits) == 3
    assert all(hit.relation_class == "direct_signature_analogue" for hit in result.hits)
    assert all(
        hit.score_trace.components["edit_similarity"] == 1.0 for hit in result.hits
    )
    assert result.hits[0].resolved_recipe["catalysts"]
    assert "DISCOVERY_CONDITIONS_ARE_OBSERVED_NOT_RECOMMENDED" in result.warnings
    assert result.to_dict()["schema_version"] == "1.0"


def test_discovery_keeps_failure_evidence_by_default_and_can_filter_it() -> None:
    explorer = _explorer()
    inclusive = explorer.discover(REACTION, top_k=10)
    filtered = explorer.discover(
        REACTION,
        top_k=10,
        include_low_yield=False,
        include_unreported_outcomes=False,
    )

    assert {hit.yield_pct for hit in inclusive.hits} == {88.0, 5.0, None}
    assert [hit.yield_pct for hit in filtered.hits] == [88.0]
    low = next(hit for hit in inclusive.hits if hit.yield_pct == 5.0)
    missing = next(hit for hit in inclusive.hits if hit.yield_pct is None)
    assert any("failure boundary" in value for value in low.cautions)
    assert any("not reported" in value for value in missing.cautions)


def test_discovery_order_and_score_trace_are_deterministic() -> None:
    explorer = _explorer()
    first = explorer.discover(REACTION, top_k=3).to_dict()
    second = explorer.discover(REACTION, top_k=3).to_dict()

    assert first == second
    weights = load_discovery_rules()["weights"]
    assert weights["edit_similarity"] > weights["reactive_scaffold"]
    assert sum(weights.values()) == 1.0


def test_discovery_rejects_unknown_view_without_touching_recommender() -> None:
    result = _explorer().discover(REACTION, view="opaque_magic")

    assert not result.valid
    assert result.error == "UNSUPPORTED_DISCOVERY_VIEW"
