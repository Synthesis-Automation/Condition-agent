"""Tests for deterministic chemistry-tool contracts used by assistance."""

from __future__ import annotations

from dataclasses import replace

from core_retrosynthesis import (
    GenericCoreTemplate,
    GenericDisconnectionCandidate,
    GenericTemplateLibrary,
    TemplateContext,
    lookup_step_precedents,
    plan_multistep_routes,
    verify_planned_route,
)
from core_retrosynthesis.generic_models import GenericTemplatePrecedent


class _NoMatches:
    def lookup(self, identity, *, provenance_limit: int = 5):
        return None


def _candidate() -> GenericDisconnectionCandidate:
    return GenericDisconnectionCandidate(
        target_smiles="C1=CCCCC1",
        precursor_smiles="C=C.C=CC=C",
        proposed_reaction_smiles="C=C.C=CC=C>>C1=CCCCC1",
        transformation_kind=None,
        abstraction_level="L2",
        compiler_engine="reaction_core",
        template_id="template:tool-test",
        score=0.9,
        context_similarity=0.8,
        product_similarity=1.0,
        precursor_similarity=1.0,
        template_specificity=1.0,
        independent_reference_support=2,
        forward_validation_status="verified_signature",
        center_transition_key="center",
        disconnection_site_key="site",
        precedent_reaction_ids=("reaction:exact", "reaction:analog"),
        operator_id="operator:tool-test",
        realization_id="realization:tool-test",
        operator_signature="signature",
        synthon_signature="synthon",
    )


def _route():
    candidate = _candidate()
    result = plan_multistep_routes(
        candidate.target_smiles,
        object(),
        _NoMatches(),
        max_depth=2,
        molecular_weight_threshold=60.0,
        expander=lambda product, top_k: (candidate,) if product == candidate.target_smiles else (),
    )
    assert result.routes
    return result.routes[0]


def _library() -> GenericTemplateLibrary:
    """Two atom-balanced cycloadditions with an exactly replayable shared disconnection."""
    context = TemplateContext((), (), ())
    precedents = (
        GenericTemplatePrecedent(
            reaction_id="reaction:analog",
            reference_id="reference:B",
            product_smiles="CC1CC=CCC1",
            precursor_smiles="C=CC.C=CC=C",
            mapped_reaction_smiles=(
                "[CH2:1]=[CH:2][CH:3]=[CH2:4].[CH:5]([CH3:7])=[CH2:6]>>"
                "[CH2:1]1[CH:2]=[CH:3][CH2:4][CH:5]([CH3:7])[CH2:6]1"
            ),
            context=context,
        ),
        GenericTemplatePrecedent(
            reaction_id="reaction:exact",
            reference_id="reference:A",
            product_smiles="C1=CCCCC1",
            precursor_smiles="C=C.C=CC=C",
            mapped_reaction_smiles=(
                "[CH2:1]=[CH:2][CH:3]=[CH2:4].[CH2:5]=[CH2:6]>>"
                "[CH2:1]1[CH:2]=[CH:3][CH2:4][CH2:5][CH2:6]1"
            ),
            context=context,
        ),
    )
    product_smarts = "[#6:1]1-[#6:2]=[#6:3]-[#6:4]-[#6:5]-[#6:6]-1"
    precursor_smarts = "[#6:1]=[#6:2]-[#6:3]=[#6:4].[#6:5]=[#6:6]"
    template = GenericCoreTemplate(
        template_id="template:tool-test",
        operator_id="operator:tool-test",
        transformation_kind=None,
        abstraction_level="L2",
        compiler_engine="reaction_core",
        reaction_smarts=f"{product_smarts}>>{precursor_smarts}",
        product_smarts=product_smarts,
        precursor_smarts=precursor_smarts,
        edit_tokens=("formed:C-C",),
        handle_signature="test",
        stereo_policy="exact",
        observation_support=2,
        independent_reference_support=2,
        precedents=precedents,
    )
    return GenericTemplateLibrary(
        templates=(template,),
        source_row_count=2,
        accepted_observation_count=2,
        rejection_counts={},
        definition={"definition_id": "fixture"},
    )


def test_step_precedent_lookup_returns_only_admitted_ranked_evidence() -> None:
    route = _route()

    result = lookup_step_precedents(route.steps[0], _library())

    assert result.available_precedent_count == 2
    assert [item.reaction_id for item in result.matches] == [
        "reaction:exact",
        "reaction:analog",
    ]
    assert all(
        item.admission_status == "operator_library_admitted"
        for item in result.matches
    )
    assert result.matches[0].product_similarity == 1.0
    assert result.matches[0].precursor_similarity == 1.0


def test_route_verification_accepts_graph_consistency_with_condition_caution() -> None:
    report = verify_planned_route(_route())

    assert report.status == "verified_with_cautions"
    statuses = {item.gate: item.status for item in report.gates}
    assert statuses["route_tree_integrity"] == "pass"
    assert statuses["step_graph_consistency"] == "pass"
    assert statuses["forward_validation"] == "pass"
    assert statuses["terminal_completion"] == "pass"
    assert statuses["condition_evidence"] == "warn"


def test_route_verification_fails_partial_or_internally_changed_route() -> None:
    route = _route()
    changed_candidate = replace(route.steps[0].candidate, target_smiles="CCCCC")
    changed_step = replace(route.steps[0], candidate=changed_candidate)
    changed_route = replace(route, solved=False, steps=(changed_step,))

    report = verify_planned_route(changed_route)

    assert report.status == "failed"
    statuses = {item.gate: item.status for item in report.gates}
    assert statuses["step_graph_consistency"] == "fail"
    assert statuses["terminal_completion"] == "fail"
