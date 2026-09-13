"""Counterfactual diagnosis never authorizes rejected reaction matches."""

from dataclasses import asdict, replace

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.shared_core_diagnostics import diagnose_shared_core_difference
from reactive_taxonomy.shared_reaction_core import build_shared_reaction_core, compare_reaction_cores


def project(reaction):
    analysis = featurize_reaction(reaction)
    return build_shared_reaction_core(
        reaction, asdict(analysis.reaction_signature) if analysis.reaction_signature else {},
        asdict(analysis.reaction_core) if analysis.reaction_core else {},
    )


def test_same_core_is_reported_as_qualified():
    diagnostic = diagnose_shared_core_difference(project("CC=O>>CCO"), project("CCC=O>>CCCO"))
    assert diagnostic.eligible
    assert not diagnostic.matching_counterfactuals


def test_hydrogen_substitution_is_a_review_candidate_only():
    query, precedent = project("CC=O>>CCO"), project("CC(C)=O>>CC(C)O")
    diagnostic = diagnose_shared_core_difference(query, precedent)
    assert "hydrogen_substitution_with_same_delta" in diagnostic.matching_counterfactuals
    assert diagnostic.requires_review and not diagnostic.eligible
    assert not compare_reaction_cores(query, precedent).eligible


def test_opposite_edits_are_not_erased_by_diagnostics():
    diagnostic = diagnose_shared_core_difference(project("CC=O>>CCO"), project("CCO>>CC=O"))
    assert not diagnostic.eligible
    assert not diagnostic.matching_counterfactuals


def test_unavailable_and_conflicting_evidence_stay_unavailable():
    core = project("CC=O>>CCO")
    missing = replace(core, levels=(), unavailable_reasons=("AMBIGUOUS_CORRESPONDENCE",))
    assert diagnose_shared_core_difference(core, missing).classification == "observation_unavailable"
    conflicting = replace(core, definition_hash="other")
    assert diagnose_shared_core_difference(core, conflicting).classification == "definition_mismatch"
