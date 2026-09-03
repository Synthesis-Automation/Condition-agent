"""Strict-versus-guarded root operator coverage regressions."""

from __future__ import annotations

import json

from core_retrosynthesis import (
    GenericDisconnectionCandidate,
    GenericTemplateLibrary,
    MoleculeOccurrenceNode,
    ReactionRouteTree,
    RouteActionEvaluationConfig,
    RouteReactionNode,
    RouteStepEvidence,
    compare_root_operator_coverage,
    render_operator_coverage_comparison_html,
    write_operator_coverage_comparison,
)
from core_retrosynthesis.observed_route_action import (
    build_observed_route_action_label,
)
from core_retrosynthesis.cli import _parser


REACTION = (
    "[CH3:1][C:2](=[O:3])[OH:6].[NH2:4][CH2:5][CH3:7]>O>"
    "[CH3:1][C:2](=[O:3])[NH:4][CH2:5][CH3:7]"
)
TARGET = "CCNC(C)=O"
PATENT = "US01234567B2"


def _tree() -> ReactionRouteTree:
    children = (
        MoleculeOccurrenceNode(
            occurrence_id="acid",
            smiles="CC(=O)O",
            depth=1,
            terminal=True,
            terminal_evidence="observed_leaf",
            unresolved_reason=None,
        ),
        MoleculeOccurrenceNode(
            occurrence_id="amine",
            smiles="CCN",
            depth=1,
            terminal=True,
            terminal_evidence="observed_leaf",
            unresolved_reason=None,
        ),
    )
    reaction = RouteReactionNode(
        reaction_node_id="reaction",
        step_id="step",
        depth=1,
        reaction_smiles=REACTION,
        evidence=RouteStepEvidence(
            evidence_kind="observed",
            source_reaction_id="source-reaction",
        ),
        children=children,
    )
    return ReactionRouteTree(
        tree_id="tree",
        route_kind="observed",
        target_smiles=TARGET,
        root=MoleculeOccurrenceNode(
            occurrence_id="target",
            smiles=TARGET,
            depth=0,
            terminal=False,
            terminal_evidence="expanded",
            unresolved_reason=None,
            reaction=reaction,
        ),
        reaction_count=1,
        maximum_depth=1,
        fingerprint_tokens=(),
        source_route_id="route",
        patent_id=PATENT,
        split="test",
    )


def _library(policy: str) -> GenericTemplateLibrary:
    return GenericTemplateLibrary(
        templates=(),
        source_row_count=100,
        accepted_observation_count=(95 if policy == "validated_departures" else 6),
        rejection_counts={},
        definition={"core_admission_policy": policy},
    )


def _candidate(exact: bool) -> GenericDisconnectionCandidate:
    observed = build_observed_route_action_label(
        REACTION,
        route_product_smiles=TARGET,
        reaction_id="source-reaction",
        reference_id=PATENT,
    )
    assert observed.expected_precursor_smiles is not None
    assert observed.retained_operator_id is not None
    assert observed.disconnection_site_key is not None
    assert observed.synthon_signature is not None
    precursors = observed.expected_precursor_smiles if exact else "CCN.C=O"
    return GenericDisconnectionCandidate(
        target_smiles=TARGET,
        precursor_smiles=precursors,
        proposed_reaction_smiles=f"{precursors}>>{TARGET}",
        transformation_kind=None,
        abstraction_level="L2",
        compiler_engine="test",
        template_id="template:exact" if exact else "template:alternative",
        score=0.9,
        context_similarity=0.8,
        product_similarity=0.9,
        precursor_similarity=0.9,
        template_specificity=1.0,
        independent_reference_support=2,
        forward_validation_status="verified_signature",
        center_transition_key="center",
        disconnection_site_key=(
            observed.disconnection_site_key if exact else "site:alternative"
        ),
        precedent_reaction_ids=(
            (f"{PATENT}:source",) if exact else ("US07654321B2:source",)
        ),
        operator_id=(observed.retained_operator_id if exact else "operator:other"),
        realization_id="realization",
        operator_signature="signature",
        synthon_signature=(
            observed.synthon_signature if exact else "synthon:alternative"
        ),
    )


def test_compares_root_recovery_and_visible_patent_overlap(tmp_path) -> None:
    review = {
        "cases": [
            {
                "case_id": "case",
                "selection_rank": 1,
                "route_tree": _tree().to_dict(),
            }
        ]
    }

    def searcher(target, library, **kwargs):
        del target, kwargs
        exact = library.definition["core_admission_policy"] == "validated_departures"
        return (_candidate(exact),), {}

    comparison = compare_root_operator_coverage(
        review,
        _library("pass_only"),
        _library("validated_departures"),
        config=RouteActionEvaluationConfig(
            top_k=1,
            max_templates_to_apply=1,
            max_candidates_to_validate=1,
            use_hierarchical_ranking=False,
        ),
        searcher=searcher,
    )

    assert comparison.strict_summary.exact_recovery_count == 0
    assert comparison.guarded_summary.exact_recovery_count == 1
    assert comparison.guarded_summary.exact_top1_count == 1
    assert comparison.guarded_summary.exact_with_source_patent_overlap_count == 1
    assert comparison.guarded_summary.exact_without_source_patent_overlap_count == 0
    assert comparison.cases[0].guarded.outcome == "exact_top1"
    assert "SOURCE_PATENT_PRECEDENT_OVERLAP_REPORTED_NOT_EXCLUDED" in (
        comparison.warnings
    )

    html = render_operator_coverage_comparison_html(comparison)
    assert "Strict versus guarded operator coverage" in html
    assert "exact observed precursors" in html
    assert html.count("<svg") >= 4

    output_json = tmp_path / "comparison.json"
    output_html = tmp_path / "comparison.html"
    write_operator_coverage_comparison(comparison, output_json, output_html)
    assert json.loads(output_json.read_text("utf-8")) == comparison.to_dict()
    assert "operator-coverage-comparison-data" in output_html.read_text("utf-8")


def test_operator_coverage_comparison_cli_contract() -> None:
    arguments = _parser().parse_args(
        [
            "compare-root-operator-coverage",
            "review.json",
            "strict.json.gz",
            "guarded.json.gz",
            "comparison.json",
            "comparison.html",
            "--top-k",
            "7",
            "--lazy-validation",
        ]
    )

    assert arguments.command == "compare-root-operator-coverage"
    assert arguments.top_k == 7
    assert arguments.lazy_validation
