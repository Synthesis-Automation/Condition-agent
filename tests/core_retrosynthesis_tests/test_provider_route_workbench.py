"""Provider-backed multistep workbench regressions."""

from __future__ import annotations

from core_retrosynthesis import (
    CallableTransitionProvider,
    GenericCoreTemplate,
    GenericDisconnectionCandidate,
    GenericSearchDiagnostics,
    GenericTemplateLibrary,
    ProviderBackedOneStepExpander,
    RouteWorkbenchSettings,
    TemplateContext,
    TransitionProviderBatch,
    TransitionProviderMetadata,
    TransitionProviderOrchestrator,
    expansion_state_from_route,
    run_provider_route_workbench,
)
from core_retrosynthesis.generic_models import GenericTemplatePrecedent


class _NoMatches:
    def lookup(self, identity, *, provenance_limit: int = 5):
        return None


def _candidate() -> GenericDisconnectionCandidate:
    return GenericDisconnectionCandidate(
        target_smiles="CCN",
        precursor_smiles="CCBr.N",
        proposed_reaction_smiles="CCBr.N>>CCN",
        transformation_kind=None,
        abstraction_level="L2",
        compiler_engine="reaction_core",
        template_id="template:test",
        score=0.9,
        context_similarity=0.8,
        product_similarity=1.0,
        precursor_similarity=1.0,
        template_specificity=1.0,
        independent_reference_support=2,
        forward_validation_status="verified_signature",
        center_transition_key="center",
        disconnection_site_key="site",
        precedent_reaction_ids=("reaction:test",),
        operator_id="operator:test",
        realization_id="realization:test",
        operator_signature="signature:test",
        synthon_signature="synthon:test",
        strategic_class="scaffold_split",
        strategic_complexity_score=0.5,
    )


def _library() -> GenericTemplateLibrary:
    context = TemplateContext((), (), ())
    precedent = GenericTemplatePrecedent(
        reaction_id="reaction:test",
        reference_id="reference:test",
        product_smiles="CCN",
        precursor_smiles="CCBr.N",
        mapped_reaction_smiles="[CH3:1][CH2:2][Br:3].[NH3:4]>>[CH3:1][CH2:2][NH2:4]",
        context=context,
    )
    template = GenericCoreTemplate(
        template_id="template:test",
        operator_id="operator:test",
        transformation_kind=None,
        abstraction_level="L2",
        compiler_engine="reaction_core",
        reaction_smarts="[*:1]>>[*:1]",
        product_smarts="[*:1]",
        precursor_smarts="[*:1]",
        edit_tokens=("formed:C-N",),
        handle_signature="test",
        stereo_policy="exact",
        observation_support=2,
        independent_reference_support=2,
        precedents=(precedent,),
    )
    return GenericTemplateLibrary(
        templates=(template,),
        source_row_count=2,
        accepted_observation_count=2,
        rejection_counts={},
        definition={"definition_id": "test"},
    )


def _orchestrator() -> TransitionProviderOrchestrator:
    candidate = _candidate()
    diagnostics = GenericSearchDiagnostics(
        library_template_count=1,
        indexed_template_count=1,
        metadata_filtered_template_count=1,
        product_query_match_count=1,
        applied_template_count=1,
        generated_precursor_count=1,
        validation_attempt_count=1,
        valid_candidate_count=1,
    )

    def expand(target: str, budget: int) -> TransitionProviderBatch:
        candidates = (candidate,) if target == "CCN" else ()
        return TransitionProviderBatch(
            candidates=candidates[:budget],
            diagnostics={
                "levels_attempted": ["L2"],
                "level_diagnostics": {"L2": diagnostics.to_dict()},
                "proposed_action_count": len(candidates),
                "validation_attempt_count": len(candidates),
                "valid_action_count": len(candidates),
            },
        )

    provider = CallableTransitionProvider(
        metadata=TransitionProviderMetadata(
            provider_id="provider:test",
            display_name="Test provider",
            maximum_proposal_budget=5,
        ),
        expansion_function=expand,
    )
    return TransitionProviderOrchestrator((provider,))


def test_provider_expander_preserves_admission_and_ladder_diagnostics() -> None:
    expander = ProviderBackedOneStepExpander(
        _orchestrator(),
        "provider:test",
    )

    batch = expander("CCN", 3)

    assert batch.candidates == (_candidate(),)
    assert batch.diagnostics is not None
    assert batch.diagnostics.levels_attempted == ("L2",)
    assert batch.diagnostics.valid_action_count == 1
    attribution = expander.attribution_for_candidate(batch.candidates[0])
    assert attribution is not None
    assert attribution.provider_id == "provider:test"
    assert attribution.provider_rank == 1


def test_workbench_composes_provider_precedent_and_verification_evidence() -> None:
    result = run_provider_route_workbench(
        "CCN",
        _library(),
        _NoMatches(),
        _orchestrator(),
        "provider:test",
        settings=RouteWorkbenchSettings(
            max_depth=2,
            molecular_weight_threshold=150.0,
            top_k_routes=2,
            per_step_top_k=3,
            beam_width=4,
            max_expansions=4,
        ),
    )

    assert result.route_kind == "solved"
    assert len(result.routes) == 1
    report = result.routes[0]
    assert report.route.solved is True
    assert report.verification.status == "verified_with_cautions"
    assert report.weakest_issue_id is None
    assert report.repair_proposals == ()
    assert len(report.step_evidence) == 1
    evidence = report.step_evidence[0]
    assert evidence.provider_id == "provider:test"
    assert evidence.provider_rank == 1
    assert evidence.transition_id
    assert evidence.precedent_evidence is not None
    assert evidence.precedent_evidence.available_precedent_count == 1
    state = expansion_state_from_route(report.route)
    assert state.leaves
    assert all(not leaf.expandable for leaf in state.leaves)
    serialized = result.to_dict()
    assert serialized["routes"][0]["expansion_state"]["state_id"]


def test_partial_route_exposes_unresolved_leaf_as_weakest_issue() -> None:
    result = run_provider_route_workbench(
        "CCN",
        _library(),
        _NoMatches(),
        _orchestrator(),
        "provider:test",
        settings=RouteWorkbenchSettings(
            max_depth=2,
            molecular_weight_threshold=10.0,
            top_k_routes=1,
            per_step_top_k=3,
            beam_width=4,
            max_expansions=4,
        ),
    )

    assert result.route_kind == "partial"
    report = result.routes[0]
    assert report.route.solved is False
    weakest = next(
        item for item in report.issues if item.issue_id == report.weakest_issue_id
    )
    assert weakest.kind == "unresolved_leaf"
    assert weakest.severity == "strong"
    state = expansion_state_from_route(report.route)
    assert any(leaf.expandable for leaf in state.leaves)
