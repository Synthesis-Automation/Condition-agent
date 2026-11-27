"""Contracts for transferring source-paired, model-mapped reaction evidence."""

import pytest

from condition_recommender.conversion.generic import GenericConversionCache, convert_record
from condition_recommender.conversion.intermediate import raw_record_from_intermediate
from condition_recommender.generic_indexing import build_generic_index
from condition_recommender.ingestion.models import INTERMEDIATE_OBSERVATION_SCHEMA_VERSION
from reactive_taxonomy.chemistry.source_mapping import assess_source_mapping


ETHANOL = "[CH3:1][CH2:2][OH:3]"
ACETALDEHYDE = "[CH3:1][CH:2]=[O:3]"


def claim(key, name, smiles):
    return {
        "component_key": key, "source_slot": "materials", "source_role_hint": "reagent",
        "identifiers": [
            {"identifier_type": "name", "value": name},
            {"identifier_type": "smiles", "value": smiles},
        ],
        "provenance": {"requires_participation_check": True},
    }


def observation(components=(), agents="", *, warnings=(), raw=None, outcomes=()):
    return {
        "schema_version": INTERMEDIATE_OBSERVATION_SCHEMA_VERSION,
        "observation_id": "source:route:step:1",
        "source": {"corpus_id": "test", "source_record_id": "step-1"},
        "reaction": {
            "reaction_smiles": f"{ETHANOL}>>{ACETALDEHYDE}",
            "source_labels": {
                "source_reaction_smiles": f"{ETHANOL}>{agents}>{ACETALDEHYDE}",
                "source_unmapped_reaction_smiles": raw or (
                    f"CCO.{agents}>>CC=O" if agents else "CCO>>CC=O"
                ),
                "mapping_origin": "RXNMapper",
            },
        },
        "conditions": {"components": components, "warnings": warnings},
        "outcomes": outcomes,
    }


def test_paired_registry_conflict_blocks_transfer_and_preserves_both_channels():
    isotope_labelled = "[2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]"
    item = claim("s1", "DMSO", isotope_labelled)
    payload = observation([item], isotope_labelled)
    record = convert_record(raw_record_from_intermediate(payload))
    evidence = record.source["source_transfer_assessment"]
    assert evidence["claims"][0]["identifier_status"] == "conflict"
    assert "source_material_identity_conflict" in record.admission_reasons
    assert record.condition_status.value == "resolved_partial"
    assert record.source["source_transfer_blocked"]
    component, = record.resolved_recipe["other_components"]
    assert component["substance_id"] is None
    pair = component["provenance"]["paired_identifier_evidence"]
    assert pair["name_candidates"] == ("cas:67-68-5",)
    assert pair["structure_candidates"] == ("cas:2206-27-1",)
    assert component["provenance"]["source_identifiers"] == tuple(item["identifiers"])
    assert not build_generic_index([record.to_dict()]).rows
    assert record.precedent_tier is None


def test_ordinary_dmso_name_and_structure_corroborate_without_isotopes():
    payload = observation([claim("s1", "DMSO", "CS(C)=O")], "CS(C)=O")
    record = convert_record(raw_record_from_intermediate(payload))
    assert record.source["source_transfer_assessment"]["claims"][0]["identifier_status"] == "corroborated"
    assert "source_material_identity_conflict" not in record.admission_reasons


def test_donor_removed_with_auditable_claim_and_unlisted_middle_retained():
    payload = observation([claim("s1", "ethanol", "CCO"), claim("s2", "water", "O")], "O.CC#N")
    record = convert_record(raw_record_from_intermediate(payload))
    evidence = record.source["source_transfer_assessment"]
    assert [c["participation"] for c in evidence["claims"]] == [
        "mapped_atom_contributor", "no_mapped_product_atoms",
    ]
    assert evidence["claims"][1]["identifier_status"] == "corroborated"
    assert evidence["unlisted_agent_graphs"] == ("CC#N",)
    assert "source_middle_agents_unaccounted" in record.admission_reasons
    assert "source_model_mapping_requires_review" in record.admission_reasons
    assert [c["raw_identifier"] for c in record.condition_resolution["components"]] == ["water"]
    assert record.source["condition_component_inputs"][0]["raw_identifier"] == "ethanol"
    assert record.reaction_signature is not None


@pytest.mark.parametrize(("name", "smiles", "agents", "reason"), [
    ("salt", "CCO.[Na+]", "[Na+]", "source_material_contributor_overlap"),
    ("ethanol", "CCO", "CCO", "source_material_contributor_overlap"),
    ("water", "O", "", "source_material_not_accounted_in_inputs"),
    ("sodium hydroxide", "[OH-].[Na+]", "[Na+]", "source_material_not_accounted_in_inputs"),
])
def test_ambiguous_participation_keeps_entire_claim_for_review(name, smiles, agents, reason):
    record = convert_record(raw_record_from_intermediate(observation([claim("s", name, smiles)], agents)))
    assert reason in record.admission_reasons
    assert record.source["source_transfer_blocked"]
    assert len(record.condition_resolution["components"]) == 1
    assert record.source["source_transfer_assessment"]["unlisted_agent_graphs"] == (
        (agents,) if agents else ()
    )


def test_source_graph_conflict_prevents_donor_exclusion_even_when_analysis_is_cached():
    cache = GenericConversionCache()
    inputs = [claim("s", "ethanol", "CCO")]
    first = convert_record(raw_record_from_intermediate(observation(inputs)), cache=cache)
    outcome = {"outcome_type": "reported_yield_pct", "unit": "%", "value": 60}
    payload = observation(inputs, raw="CO>>CO", warnings=["unpaired_condition_evidence:media"],
                          outcomes=[outcome])
    second = convert_record(raw_record_from_intermediate(payload), cache=cache)
    assert first.condition_resolution["component_count"] == 0
    assert second.condition_resolution["component_count"] == 1
    assert second.yield_pct is None
    assert {"source_mapped_inputs_disagree", "source_mapped_products_disagree",
            "unpaired_condition_evidence:media"}.issubset(second.admission_reasons)
    assert second.source["raw_fields"]["source_observation"] == payload


def test_source_mapping_checks_element_isotope_multiplicity_stereo_and_dative_bonds():
    assert "source_product_atom_origin_unverified" in assess_source_mapping(
        "CO>>[13CH2]=O", "[CH3:1][OH:2]>>[13CH2:1]=[O:2]",
    ).issues
    assert "source_duplicate_atom_maps" in assess_source_mapping(
        "CO>>C=O", "[CH3:1][OH:1]>>[CH2:1]=[O:2]",
    ).issues
    assert "source_mapped_products_disagree" in assess_source_mapping(
        "F[C@H](Cl)Br>>F[C@@H](Cl)Br", "F[C@H](Cl)Br>>F[C@H](Cl)Br",
    ).issues
    assert "source_mapped_inputs_disagree" in assess_source_mapping(
        "CCO.O.O>>CC=O", f"{ETHANOL}>O>{ACETALDEHYDE}",
    ).issues
    assert not assess_source_mapping(
        "CCO.N->[Cu+2]>>CC=O", f"{ETHANOL}>N->[Cu+2]>{ACETALDEHYDE}",
    ).issues


@pytest.mark.parametrize(("kind", "metadata", "expected"), [
    ("reported_yield_pct", {}, 41), ("reported_product_yield_pct", {}, 41),
    ("uv_area_yield_pct", {}, 41), ("sequence_yield_pct", {}, None),
    ("reported_yield_pct", {"scope": "sequence"}, None),
    ("reported_yield_pct", {"primary_eligible": False}, None),
])
def test_only_recognized_step_assays_project_to_primary_yield(kind, metadata, expected):
    outcome = {"outcome_type": kind, "unit": "%", "value": 41, "metadata": metadata}
    raw = raw_record_from_intermediate(observation(outcomes=[outcome]))
    assert raw.yield_pct == expected
    assert raw.raw_fields["source_observation"]["outcomes"] == [outcome]
