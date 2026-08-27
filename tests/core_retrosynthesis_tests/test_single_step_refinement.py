"""Deterministic single-step repair and verification regressions."""

from __future__ import annotations

from dataclasses import replace

from core_retrosynthesis import (
    GenericDisconnectionCandidate,
    collect_single_step_refinement_issues,
    detect_functional_group_competition,
    enumerate_single_step_repair_proposals,
    group_strategy_candidates,
    verify_single_step_strategy,
)


REACTION = "Cc1[nH]cnc1CCl.NCCS>>SCCNCc1c(C)[nH]cn1"


def _candidate(name: str, *, warning: bool) -> GenericDisconnectionCandidate:
    value = GenericDisconnectionCandidate(
        target_smiles="CCN",
        precursor_smiles="CCBr.N",
        proposed_reaction_smiles="CCBr.N>>CCN",
        transformation_kind="substitution",
        abstraction_level="L2",
        compiler_engine="test",
        template_id=f"template:{name}",
        score=0.9,
        context_similarity=0.8,
        product_similarity=0.9,
        precursor_similarity=0.8,
        template_specificity=0.8,
        independent_reference_support=2,
        forward_validation_status="verified_signature",
        center_transition_key=f"center:{name}",
        disconnection_site_key=f"SITE1:{name}",
        precedent_reaction_ids=(f"precedent:{name}",),
        operator_id=f"OP1:{name}",
        realization_id=f"REAL1:{name}",
        operator_signature=f"operator:{name}",
        synthon_signature=f"SYN1:{name}",
    )
    if not warning:
        return value
    competition = detect_functional_group_competition(REACTION)
    assert competition is not None
    return replace(value, selectivity_warnings=(competition,))


def _strategy(candidate: GenericDisconnectionCandidate):
    return group_strategy_candidates((candidate,), top_k_strategies=1)[0]


def test_existing_lower_issue_strategy_is_an_actionable_repair() -> None:
    source = _strategy(_candidate("source", warning=True))
    target = _strategy(_candidate("target", warning=False))
    issue = collect_single_step_refinement_issues(source)[0]

    proposals = enumerate_single_step_repair_proposals(
        source,
        issue,
        strategies=(source, target),
        condition_evidence_by_strategy={},
    )
    alternate = next(
        item for item in proposals if item.repair_kind == "alternate_strategy"
    )

    assert issue.kind == "selectivity"
    assert alternate.status == "actionable"
    assert alternate.target_strategy_id == target.strategy_id
    assert "SMILES" not in alternate.to_dict()


def test_single_step_verification_rechecks_reaction_and_target_identity() -> None:
    strategy = _strategy(_candidate("valid", warning=False))
    verified = verify_single_step_strategy(strategy)
    invalid_candidate = replace(
        strategy.representative,
        proposed_reaction_smiles="CCBr.N>>CCC",
    )
    invalid_strategy = replace(strategy, representative=invalid_candidate)
    failed = verify_single_step_strategy(invalid_strategy)

    assert verified.status == "verified"
    assert failed.status == "failed"
    assert any(
        gate.gate == "target_identity" and gate.status == "fail"
        for gate in failed.gates
    )


def test_no_issue_produces_no_repair_surface() -> None:
    strategy = _strategy(_candidate("clean", warning=False))

    assert collect_single_step_refinement_issues(strategy) == ()
