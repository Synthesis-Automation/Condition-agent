"""Role-neutral synthetic partition and landscape regressions."""

from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

import core_retrosynthesis.cli as cli_module
from cas_tools import MoleculeIndexMatch, molecule_identity
from core_retrosynthesis import (
    PARTITION_ASSESSMENT_POLICY_PATH,
    PRECURSOR_STATE_FEASIBILITY_POLICY_PATH,
    GenericDisconnectionCandidate,
    GenericCoreTemplate,
    GenericTemplateLibrary,
    ModuleAnnotation,
    MoleculeOccurrenceNode,
    ReactionRouteTree,
    RouteReactionNode,
    RouteStepEvidence,
    PARTITION_REALIZATION_POLICY_PATH,
    SYNTHETIC_PARTITION_POLICY_PATH,
    PartitionRealizationConfig,
    PartitionRealizationResult,
    PartitionAssessmentResult,
    RetrosynthesisConditionEvidence,
    SyntheticPartition,
    SyntheticPartitionLandscape,
    TemplateContext,
    analyze_partition_target,
    assess_partition_realizations,
    build_partition_blind_review_packet,
    build_module_id,
    build_operator_partition_landscape,
    create_synthetic_partition,
    load_synthetic_partition_policy,
    plan_multistep_routes,
    project_reaction_to_target,
    project_route_partitions,
    realize_synthetic_partition,
    render_partition_assessment_html,
    render_partition_landscape_html,
    validate_partition_assessment_policy,
    validate_precursor_state_feasibility_policy,
    validate_partition_realization_policy,
    validate_synthetic_partition_policy,
    write_partition_landscape_review,
)
from core_retrosynthesis.generic_models import GenericTemplatePrecedent


TARGET = "CNC(C)=O"
PROJECT_ROOT = Path(__file__).resolve().parents[2]
AMIDE_FORMATION = (
    "[CH3:1][C:2](=[O:3])[OH:6].[NH2:4][CH3:5]"
    ">>[CH3:1][C:2](=[O:3])[NH:4][CH3:5]"
)
N_METHYLATION = (
    "[CH3:1][C:2](=[O:3])[NH2:4].[CH3:5][Br:6]"
    ">>[CH3:1][C:2](=[O:3])[NH:4][CH3:5]"
)


def _candidate(
    name: str,
    reaction: str,
    *,
    score: float = 0.9,
    support: int = 2,
    context: float = 0.7,
    status: str = "verified_signature",
) -> GenericDisconnectionCandidate:
    precursors, _ = reaction.split(">>")
    return GenericDisconnectionCandidate(
        target_smiles=TARGET,
        precursor_smiles=precursors,
        proposed_reaction_smiles=reaction,
        transformation_kind=None,
        abstraction_level="L2",
        compiler_engine="test",
        template_id=f"template:{name}",
        score=score,
        context_similarity=context,
        product_similarity=score,
        precursor_similarity=score,
        template_specificity=1.0,
        independent_reference_support=support,
        forward_validation_status=status,
        center_transition_key=f"center:{name}",
        disconnection_site_key=f"site:{name}",
        precedent_reaction_ids=(f"precedent:{name}",),
        operator_id=f"operator:{name}",
        realization_id=f"realization:{name}",
        operator_signature=f"signature:{name}",
        synthon_signature=f"synthon:{name}",
        condition_query_reaction_smiles=reaction,
    )


def _route_candidate(
    product: str,
    precursors: str,
    reaction: str,
    name: str,
) -> GenericDisconnectionCandidate:
    return GenericDisconnectionCandidate(
        target_smiles=product,
        precursor_smiles=precursors,
        proposed_reaction_smiles=reaction,
        transformation_kind=None,
        abstraction_level="L2",
        compiler_engine="test",
        template_id=f"template:{name}",
        score=0.9,
        context_similarity=0.8,
        product_similarity=0.9,
        precursor_similarity=0.9,
        template_specificity=1.0,
        independent_reference_support=2,
        forward_validation_status="verified_signature",
        center_transition_key=f"center:{name}",
        disconnection_site_key=f"site:{name}",
        precedent_reaction_ids=(f"precedent:{name}",),
        operator_id=f"operator:{name}",
        realization_id=f"realization:{name}",
        operator_signature=f"signature:{name}",
        synthon_signature=f"synthon:{name}",
        condition_query_reaction_smiles=reaction,
    )


class _PartitionLiteratureIndex:
    def __init__(self, known: tuple[str, ...]) -> None:
        self.known = {
            molecule_identity(value).canonical_smiles for value in known
        }

    def lookup(self, identity, *, provenance_limit: int = 5):
        if identity.canonical_smiles not in self.known:
            return None
        return MoleculeIndexMatch(
            canonical_smiles=identity.canonical_smiles,
            inchi_key=identity.inchi_key,
            occurrence_count=1,
            source_records=(
                {"reaction_id": "known", "source_role": "reactant"},
            ),
        )


def _partition_expander(include_second_step: bool = True):
    root = _route_candidate(
        TARGET,
        "CC(=O)O.CN",
        AMIDE_FORMATION,
        "amide-route",
    )
    methylamine_reaction = (
        "[NH3:11].[CH3:12][Br:13]>>[NH2:11][CH3:12]"
    )
    methylamine = _route_candidate(
        "CN",
        "CBr.N",
        methylamine_reaction,
        "methylamine-route",
    )
    expansions = {TARGET: (root,)}
    if include_second_step:
        expansions["CN"] = (methylamine,)

    def expand(product: str, top_k: int):
        return expansions.get(product, ())[:top_k]

    return expand


def _assessment_library() -> GenericTemplateLibrary:
    context = TemplateContext((), (), ())

    def template(
        name: str,
        product: str,
        precursors: str,
        reaction: str,
    ) -> GenericCoreTemplate:
        return GenericCoreTemplate(
            template_id=f"template:{name}",
            operator_id=f"operator:{name}",
            transformation_kind=None,
            abstraction_level="L2",
            compiler_engine="reaction_core",
            reaction_smarts="[*:1]>>[*:1]",
            product_smarts="[*:1]",
            precursor_smarts="[*:1]",
            edit_tokens=("test",),
            handle_signature="test",
            stereo_policy="exact",
            observation_support=1,
            independent_reference_support=1,
            precedents=(
                GenericTemplatePrecedent(
                    reaction_id=f"precedent:{name}",
                    reference_id=f"reference:{name}",
                    product_smiles=product,
                    precursor_smiles=precursors,
                    mapped_reaction_smiles=reaction,
                    context=context,
                ),
            ),
        )

    return GenericTemplateLibrary(
        templates=(
            template(
                "amide-route",
                TARGET,
                "CC(=O)O.CN",
                AMIDE_FORMATION,
            ),
            template(
                "methylamine-route",
                "CN",
                "CBr.N",
                "[NH3:11].[CH3:12][Br:13]>>[NH2:11][CH3:12]",
            ),
        ),
        source_row_count=2,
        accepted_observation_count=2,
        rejection_counts={},
        definition={"definition_id": "partition-assessment-test"},
    )


def _direct_condition(reaction_smiles: str) -> RetrosynthesisConditionEvidence:
    return RetrosynthesisConditionEvidence(
        status="recommended_direct",
        query_reaction_smiles=reaction_smiles,
        recommender_valid=True,
        recommendation_mode="verified_signature",
        retrieval_level="exact_signature",
        uses_type_agnostic_fallback=False,
        candidate_count=1,
        independent_candidate_count=1,
        compatible_candidate_count=1,
        independent_compatible_candidate_count=1,
        excluded_candidate_count=0,
        best_recipe_score=0.9,
        best_recipe_compatibility_score=1.0,
        best_recipe_reference_support=1,
        recommendations=(
            {
                "recipe_id": "recipe:test",
                "recipe_core_id": "core:test",
                "resolved_recipe": {"components": []},
                "precedent_reference_ids": ["reference:test"],
            },
        ),
        warnings=(),
        error=None,
    )


def _observed_two_step_tree() -> ReactionRouteTree:
    ammonia = MoleculeOccurrenceNode(
        occurrence_id="ammonia",
        smiles="N",
        depth=2,
        terminal=True,
        terminal_evidence="observed_leaf",
        unresolved_reason=None,
    )
    methyl_bromide = MoleculeOccurrenceNode(
        occurrence_id="methyl_bromide",
        smiles="CBr",
        depth=2,
        terminal=True,
        terminal_evidence="observed_leaf",
        unresolved_reason=None,
    )
    methylamine_reaction = RouteReactionNode(
        reaction_node_id="reaction:methylamine",
        step_id="step:methylamine",
        depth=2,
        reaction_smiles="[NH3:11].[CH3:12][Br:13]>>[NH2:11][CH3:12]",
        evidence=RouteStepEvidence(evidence_kind="observed"),
        children=(ammonia, methyl_bromide),
    )
    methylamine = MoleculeOccurrenceNode(
        occurrence_id="methylamine",
        smiles="CN",
        depth=1,
        terminal=False,
        terminal_evidence="expanded",
        unresolved_reason=None,
        reaction=methylamine_reaction,
    )
    acetic_acid = MoleculeOccurrenceNode(
        occurrence_id="acetic_acid",
        smiles="CC(=O)O",
        depth=1,
        terminal=True,
        terminal_evidence="observed_leaf",
        unresolved_reason=None,
    )
    root_reaction = RouteReactionNode(
        reaction_node_id="reaction:amide",
        step_id="step:amide",
        depth=1,
        reaction_smiles=AMIDE_FORMATION,
        evidence=RouteStepEvidence(evidence_kind="observed"),
        children=(acetic_acid, methylamine),
    )
    root = MoleculeOccurrenceNode(
        occurrence_id="target",
        smiles=TARGET,
        depth=0,
        terminal=False,
        terminal_evidence="expanded",
        unresolved_reason=None,
        reaction=root_reaction,
    )
    return ReactionRouteTree(
        tree_id="tree:partition-test",
        route_kind="observed",
        target_smiles=TARGET,
        root=root,
        reaction_count=2,
        maximum_depth=2,
        fingerprint_tokens=("amide", "methylamine"),
    )


def _full_realization_result() -> PartitionRealizationResult:
    desired = project_route_partitions(
        _observed_two_step_tree()
    ).frontiers[-1].partition
    return realize_synthetic_partition(
        desired,
        object(),
        _PartitionLiteratureIndex(("CC(=O)O", "N", "CBr")),
        config=PartitionRealizationConfig(
            max_depth=2,
            molecular_weight_threshold=1.0,
            maximum_realizations=1,
            route_candidate_limit=2,
            per_step_top_k=2,
            beam_width=8,
            maximum_expansions=10,
            maximum_templates_to_apply=10,
            maximum_candidates_to_validate=4,
        ),
        expander=_partition_expander(),
    )


def test_partition_identity_is_symmetric_and_role_neutral() -> None:
    canonical, target_id, atoms, _ = analyze_partition_target(TARGET)
    target_maps = tuple(atom.atom_map for atom in atoms)
    left = target_maps[:2]
    right = target_maps[2:]
    plain = create_synthetic_partition(
        canonical,
        (left, right),
        source_kind="structural_baseline",
        evidence_level="E0",
    )
    annotated = create_synthetic_partition(
        canonical,
        (right, left),
        source_kind="structural_baseline",
        evidence_level="E0",
        annotations=(
            ModuleAnnotation(
                module_id=build_module_id(target_id, left),
                label="optional label",
                proposed_role="complex_fragment",
                confidence=0.6,
            ),
        ),
    )

    assert plain.target_atom_maps == target_maps
    assert plain.partition_id == annotated.partition_id
    assert tuple(module.target_atom_maps for module in plain.modules) == (
        left,
        right,
    )
    assert SyntheticPartition.from_dict(annotated.to_dict()) == annotated


def test_partition_rejects_duplicate_or_incomplete_atom_ownership() -> None:
    with pytest.raises(ValueError, match="duplicates target atom ownership"):
        create_synthetic_partition(
            TARGET,
            ((1, 2, 3), (3, 4, 5)),
            source_kind="structural_baseline",
            evidence_level="E0",
        )

    with pytest.raises(ValueError, match="cover every target heavy atom"):
        create_synthetic_partition(
            TARGET,
            ((1, 2), (3, 4)),
            source_kind="structural_baseline",
            evidence_level="E0",
        )


def test_mapped_reaction_projection_separates_target_and_tactical_atoms() -> None:
    projection = project_reaction_to_target(TARGET, AMIDE_FORMATION)

    assert projection.k == 2
    assert sorted(map(len, projection.module_atom_sets)) == [2, 3]
    assert len(projection.target_boundary_bonds) == 1
    assert 6 not in {
        atom_map
        for module in projection.module_atom_sets
        for atom_map in module
    }
    assert projection.mapping_evidence == "supplied_atom_mapping"


def test_reaction_projection_reports_inferred_mapping_and_rejects_missing_source() -> None:
    inferred = project_reaction_to_target(
        TARGET,
        "CC(=O)O.CN>>CNC(C)=O",
    )

    assert inferred.k == 2
    assert inferred.mapping_evidence == "global_atom_correspondence"
    assert inferred.mapping_confidence < 1.0
    with pytest.raises(ValueError, match="do not cover current target ownership"):
        project_reaction_to_target(
            "CO",
            "[CH3:1]>>[CH3:1][OH:2]",
        )
    with pytest.raises(ValueError, match="changes target element"):
        project_reaction_to_target("C", "[Br:1]>>[CH4:1]")


def test_observed_route_projection_returns_k1_k2_and_k3_frontiers() -> None:
    projection = project_route_partitions(_observed_two_step_tree())

    assert projection.unresolved_occurrence_ids == ()
    assert [frontier.depth for frontier in projection.frontiers] == [0, 1, 2]
    assert [frontier.partition.k for frontier in projection.frontiers] == [1, 2, 3]
    assert all(
        frontier.partition.evidence_level == "E4"
        for frontier in projection.frontiers
    )
    assert len(projection.frontiers[-1].latent_states) == 3
    assert projection.to_landscape().searched_k_values == (1, 2, 3)


def test_partition_realization_finds_a_fully_validated_route() -> None:
    desired = project_route_partitions(_observed_two_step_tree()).frontiers[-1].partition
    result = realize_synthetic_partition(
        desired,
        object(),
        _PartitionLiteratureIndex(("CC(=O)O", "N", "CBr")),
        config=PartitionRealizationConfig(
            max_depth=2,
            molecular_weight_threshold=1.0,
            maximum_realizations=2,
            route_candidate_limit=4,
            per_step_top_k=2,
            beam_width=8,
            maximum_expansions=10,
            maximum_templates_to_apply=10,
            maximum_candidates_to_validate=4,
        ),
        expander=_partition_expander(),
    )

    assert result.status == "fully_realized"
    assert result.partition.realization_status == "fully_realized"
    assert result.diagnostics.fully_realized_count == 1
    realization = result.realizations[0]
    assert realization.best_frontier.exact_partition_match is True
    assert realization.best_frontier.frontier_k == 3
    assert realization.evidence_summary.validated_interface_coverage == 1.0
    assert realization.evidence_summary.validated_step_count == 2
    assert realization.evidence_summary.unsupported_step_count == 0
    assert realization.evidence_summary.dependency_graph_valid is True
    assert realization.evidence_summary.forward_order_exists is True
    assert len(realization.dependency_edges) == 1
    assert len(realization.frontier_states) == 3
    assert all(
        state.source_occurrence_id is not None
        for state in realization.frontier_states
    )
    graph = realization.latent_realization_graph
    assert graph is not None
    assert graph.validation_status == "validated"
    assert len(graph.transitions) == 2
    assert [len(stage.transition_ids) for stage in graph.forward_stages] == [1, 1]
    transitions = {transition.transition_id: transition for transition in graph.transitions}
    assert [
        transitions[stage.transition_ids[0]].retrosynthetic_depth
        for stage in graph.forward_stages
    ] == [2, 1]
    assert realization.evidence_summary.latent_transition_count == 2
    assert realization.evidence_summary.forward_stage_count == 2
    assert realization.evidence_summary.parallel_forward_stage_count == 0
    assert PartitionRealizationResult.from_dict(result.to_dict()) == result

    with pytest.raises(ValueError, match="forward schedule stages must be contiguous"):
        replace(graph, forward_stages=tuple(reversed(graph.forward_stages)))
    invalid_transition = replace(
        graph.transitions[0],
        output_state_id="LATENT_STATE:missing",
    )
    with pytest.raises(ValueError, match="unknown state"):
        replace(
            graph,
            transitions=(invalid_transition, *graph.transitions[1:]),
        )
    assembly = next(
        transition
        for transition in graph.transitions
        if len(transition.input_state_ids) > 1
    )
    invalid_ownership = replace(
        assembly,
        input_state_ids=assembly.input_state_ids[:1],
    )
    with pytest.raises(ValueError, match="changes target-atom ownership"):
        replace(
            graph,
            transitions=tuple(
                invalid_ownership
                if transition.transition_id == assembly.transition_id
                else transition
                for transition in graph.transitions
            ),
        )
    with pytest.raises(ValueError, match="dependencies do not match"):
        replace(graph, dependency_edges=())


def test_partition_realization_retains_partial_progress_and_failure_reason() -> None:
    desired = project_route_partitions(_observed_two_step_tree()).frontiers[-1].partition
    result = realize_synthetic_partition(
        desired,
        object(),
        _PartitionLiteratureIndex(("CC(=O)O", "N", "CBr")),
        config=PartitionRealizationConfig(
            max_depth=2,
            molecular_weight_threshold=1.0,
            maximum_realizations=1,
            route_candidate_limit=2,
            per_step_top_k=1,
            beam_width=4,
            maximum_expansions=4,
            maximum_templates_to_apply=4,
            maximum_candidates_to_validate=2,
        ),
        expander=_partition_expander(include_second_step=False),
    )

    assert result.status == "partially_realized"
    assert result.diagnostics.dead_end_states == 1
    assert "OPERATOR_COVERAGE_INCOMPLETE" in result.warnings
    realization = result.realizations[0]
    assert realization.best_frontier.exact_partition_match is False
    assert realization.best_frontier.boundary_recall == 0.5
    assert "ROUTE_TERMINALS_INCOMPLETE" in realization.warnings


def test_partition_realization_abstains_when_no_route_action_exists() -> None:
    desired = project_route_partitions(_observed_two_step_tree()).frontiers[-1].partition
    result = realize_synthetic_partition(
        desired,
        object(),
        _PartitionLiteratureIndex(()),
        config=PartitionRealizationConfig(
            max_depth=2,
            molecular_weight_threshold=1.0,
            maximum_realizations=1,
            route_candidate_limit=1,
            per_step_top_k=1,
            beam_width=2,
            maximum_expansions=2,
            maximum_templates_to_apply=2,
            maximum_candidates_to_validate=1,
        ),
        expander=lambda product, top_k: (),
    )

    assert result.status == "unrealized_but_plausible"
    assert result.realizations[0].reaction_count == 0
    assert result.diagnostics.dead_end_states == 1
    assert "OPERATOR_COVERAGE_INCOMPLETE" in result.warnings


def test_partition_realization_policy_is_validated() -> None:
    raw = json.loads(PARTITION_REALIZATION_POLICY_PATH.read_text("utf-8"))
    assert validate_partition_realization_policy(raw) == []
    raw["search"]["max_depth"] = 0
    assert "invalid_search" in validate_partition_realization_policy(raw)


def test_partition_phase4_attaches_precedents_conditions_and_round_trips() -> None:
    realization = _full_realization_result()

    result = assess_partition_realizations(
        realization,
        _assessment_library(),
        condition_evaluator=_direct_condition,
    )

    assert result.status == "supported"
    assert result.ranking_influence == "none_review_only"
    assert result.source.realizations == realization.realizations
    route = result.route_assessments[0]
    assert route.precedent_supported_step_count == 2
    assert route.condition_supported_step_count == 2
    assert route.hard_incompatible_step_count == 0
    assert all(
        step.precedent_evidence_level == "E4"
        for step in route.step_assessments
    )
    assert all(
        interface.precedent_match_ids
        for interface in route.interface_assessments
    )
    assert route.precursor_state_feasibility is not None
    assert route.precursor_state_feasibility.status == "supported"
    assert (
        route.precursor_state_feasibility.promotion_recommendation
        == "eligible_for_route_review"
    )
    assert route.precursor_state_feasibility.supported_step_count == 2
    assert all(
        step.precursor_state_feasibility is not None
        and step.precursor_state_feasibility.evidence_level == "E4"
        and step.precursor_state_feasibility.reactant_state_support
        == "exact_reactant_state"
        for step in route.step_assessments
    )
    assert PartitionAssessmentResult.from_dict(result.to_dict()) == result


def test_partition_precursor_state_template_only_support_is_held() -> None:
    library = _assessment_library()
    distant_templates = tuple(
        replace(
            template,
            precedents=tuple(
                replace(
                    precedent,
                    product_smiles="c1ccccc1",
                    precursor_smiles="CCCC.CCCC",
                )
                for precedent in template.precedents
            ),
        )
        for template in library.templates
    )
    result = assess_partition_realizations(
        _full_realization_result(),
        replace(library, templates=distant_templates),
        condition_evaluator=_direct_condition,
    )

    assert result.status == "insufficient_evidence"
    route = result.route_assessments[0]
    assert route.precursor_state_feasibility is not None
    assert route.precursor_state_feasibility.status == "insufficient_evidence"
    assert (
        route.precursor_state_feasibility.promotion_recommendation
        == "hold_for_evidence"
    )
    assert route.precursor_state_feasibility.held_step_count == 2
    for step in route.step_assessments:
        evidence = step.precursor_state_feasibility
        assert evidence is not None
        assert evidence.evidence_level == "E2"
        assert evidence.reactant_state_support == "template_only"
        assert "TEMPLATE_ONLY_PRECURSOR_STATE_SUPPORT" in evidence.warnings


def test_partition_phase4_keeps_missing_conditions_as_review_caution() -> None:
    result = assess_partition_realizations(
        _full_realization_result(),
        _assessment_library(),
    )

    assert result.status == "supported_with_cautions"
    assert "CONDITION_EVIDENCE_NOT_REQUESTED" in result.warnings
    assert all(
        step.condition_evidence.status == "insufficient_evidence"
        for step in result.route_assessments[0].step_assessments
    )
    html = render_partition_assessment_html(result)
    assert "Structure-first route review" in html
    assert "Chemist review" in html
    assert "RETROSYNTHETIC FRONTIER" in html
    assert "FORWARD REALIZATION ORDER" in html
    assert "Forward stage 1" in html
    assert "Precursor-state feasibility" in html
    assert "eligible for route review" in html
    assert "Depicted precedents" in html
    assert html.count("<svg") >= 8
    assert "data-route-tab='1'" in html


def test_partition_phase4_blocks_conditions_after_structural_failure() -> None:
    calls = []
    empty_library = replace(_assessment_library(), templates=())

    def conditions(reaction_smiles: str) -> RetrosynthesisConditionEvidence:
        calls.append(reaction_smiles)
        return _direct_condition(reaction_smiles)

    result = assess_partition_realizations(
        _full_realization_result(),
        empty_library,
        condition_evaluator=conditions,
    )

    assert result.status == "hard_incompatible"
    assert calls == []
    assert all(
        step.structural_status == "invalid"
        for step in result.route_assessments[0].step_assessments
    )
    assert all(
        "CONDITION_QUERY_BLOCKED_BY_STRUCTURAL_VALIDATION"
        in step.condition_evidence.warnings
        for step in result.route_assessments[0].step_assessments
    )
    assert all(
        step.precursor_state_feasibility is not None
        and step.precursor_state_feasibility.evidence_level == "E0"
        and step.precursor_state_feasibility.promotion_recommendation
        == "not_promotable"
        for step in result.route_assessments[0].step_assessments
    )


def test_partition_phase4_policy_and_blind_packet_separate_answer_key() -> None:
    raw = json.loads(PARTITION_ASSESSMENT_POLICY_PATH.read_text("utf-8"))
    assert validate_partition_assessment_policy(raw) == []
    raw["ranking_influence"] = "rerank"
    assert "invalid_ranking_influence" in validate_partition_assessment_policy(raw)

    supported = assess_partition_realizations(
        _full_realization_result(),
        _assessment_library(),
        condition_evaluator=_direct_condition,
    )
    supported_route = supported.route_assessments[0]
    negative_route = replace(
        supported_route,
        status="hard_incompatible",
        hard_incompatible_step_count=1,
    )
    negative = replace(
        supported,
        status="hard_incompatible",
        route_assessments=(negative_route,),
    )
    abstention_route = replace(
        supported_route,
        source_realization_status="partially_realized",
    )
    abstention = replace(
        supported,
        status="supported_with_cautions",
        route_assessments=(abstention_route,),
    )

    packet, answer_key = build_partition_blind_review_packet(
        (supported, negative, abstention),
        seed=9,
    )

    assert packet["warnings"] == []
    assert "source_kind" not in json.dumps(packet)
    assert "system_assessment_status" not in json.dumps(packet)
    assert set(answer_key["source_kinds"]) == {
        "realization",
        "negative_control",
        "abstention",
    }
    assert all(
        "precursor_state_feasibility" in answer
        for answer in answer_key["answers"]
    )


def test_precursor_state_feasibility_policy_is_validated() -> None:
    raw = json.loads(
        PRECURSOR_STATE_FEASIBILITY_POLICY_PATH.read_text("utf-8")
    )
    assert validate_precursor_state_feasibility_policy(raw) == []
    raw["ranking_influence"] = "production"
    assert "invalid_ranking_influence" in (
        validate_precursor_state_feasibility_policy(raw)
    )


def test_partition_realization_distinguishes_search_limit_from_dead_end() -> None:
    desired = project_route_partitions(_observed_two_step_tree()).frontiers[-1].partition
    result = realize_synthetic_partition(
        desired,
        object(),
        _PartitionLiteratureIndex(("CC(=O)O", "N", "CBr")),
        config=PartitionRealizationConfig(
            max_depth=2,
            molecular_weight_threshold=1.0,
            maximum_realizations=1,
            route_candidate_limit=1,
            per_step_top_k=1,
            beam_width=2,
            maximum_expansions=1,
            maximum_templates_to_apply=2,
            maximum_candidates_to_validate=1,
        ),
        expander=_partition_expander(),
    )

    assert result.status == "partially_realized"
    assert result.diagnostics.stopped_by_expansion_limit is True
    assert "REALIZATION_SEARCH_LIMIT_REACHED" in result.warnings
    assert "OPERATOR_COVERAGE_INCOMPLETE" not in result.warnings


def test_partition_guidance_expands_multiple_frontier_components() -> None:
    target = "CCNCCOCC"
    root = _route_candidate(
        target,
        "CCN.CCOCCBr",
        (
            "[CH3:1][CH2:2][NH2:3]."
            "[Br:9][CH2:4][CH2:5][O:6][CH2:7][CH3:8]>>"
            "[CH3:1][CH2:2][NH:3][CH2:4][CH2:5][O:6][CH2:7][CH3:8]"
        ),
        "convergent-root",
    )
    left = _route_candidate(
        "CCN",
        "CBr.CN",
        "[CH3:1][NH2:3].[CH3:10][Br:11]>>[CH3:1][CH2:10][NH2:3]",
        "convergent-left",
    )
    right = _route_candidate(
        "CCOCCBr",
        "CCBr.OCCBr",
        (
            "[Br:9][CH2:4][CH2:5][OH:6]."
            "[Br:12][CH2:7][CH3:8]>>"
            "[Br:9][CH2:4][CH2:5][O:6][CH2:7][CH3:8]"
        ),
        "convergent-right",
    )
    expansions = {
        target: (root,),
        "CCN": (left,),
        "CCOCCBr": (right,),
    }

    def expand(product: str, top_k: int):
        return expansions.get(product, ())[:top_k]

    stock = _PartitionLiteratureIndex(("CBr", "CN", "CCBr", "OCCBr"))
    ordinary = plan_multistep_routes(
        target,
        object(),
        stock,
        max_depth=2,
        molecular_weight_threshold=1.0,
        top_k_routes=1,
        per_step_top_k=1,
        beam_width=8,
        max_expansions=10,
        expander=expand,
    )
    assert ordinary.routes
    desired = project_route_partitions(
        ordinary.routes[0].route_tree
    ).frontiers[-1].partition
    result = realize_synthetic_partition(
        desired,
        object(),
        stock,
        config=PartitionRealizationConfig(
            max_depth=2,
            molecular_weight_threshold=1.0,
            maximum_realizations=1,
            route_candidate_limit=2,
            per_step_top_k=1,
            beam_width=8,
            maximum_expansions=10,
            maximum_templates_to_apply=10,
            maximum_candidates_to_validate=2,
        ),
        expander=expand,
    )

    assert result.status == "fully_realized"
    tree = result.realizations[0].route_tree
    assert tree.root.reaction is not None
    assert sum(child.reaction is not None for child in tree.root.reaction.children) == 2
    assert result.realizations[0].reaction_count == 3
    assert len(result.realizations[0].dependency_edges) == 2
    graph = result.realizations[0].latent_realization_graph
    assert graph is not None
    assert [len(stage.transition_ids) for stage in graph.forward_stages] == [2, 1]
    transitions = {transition.transition_id: transition for transition in graph.transitions}
    assert {
        transitions[transition_id].retrosynthetic_depth
        for transition_id in graph.forward_stages[0].transition_ids
    } == {2}
    assert transitions[graph.forward_stages[1].transition_ids[0]].retrosynthetic_depth == 1
    assert result.realizations[0].evidence_summary.parallel_forward_stage_count == 1


def test_partition_realization_cli_writes_a_result(
    monkeypatch,
    tmp_path,
    capsys,
) -> None:
    projection = project_route_partitions(_observed_two_step_tree())
    landscape_path = tmp_path / "landscape.json"
    output_path = tmp_path / "realization.json"
    landscape_path.write_text(
        json.dumps(projection.to_landscape().to_dict()),
        encoding="utf-8",
    )
    selected = projection.frontiers[-1].partition
    captured = {}

    class _StockContext:
        def __enter__(self):
            return object()

        def __exit__(self, exc_type, exc_value, traceback):
            return False

    class _Result:
        status = "unrealized_but_plausible"
        realizations = ()
        diagnostics = SimpleNamespace(fully_realized_count=0)

        def to_dict(self):
            return {"status": self.status, "realizations": []}

    def realize(partition, library, stock, *, config, route_action_policy):
        captured["partition_id"] = partition.partition_id
        captured["max_depth"] = config.max_depth
        return _Result()

    monkeypatch.setattr(cli_module, "load_generic_library", lambda _: object())
    monkeypatch.setattr(cli_module, "open_stock_lookup", lambda _: _StockContext())
    monkeypatch.setattr(cli_module, "realize_synthetic_partition", realize)

    status = cli_module.main(
        (
            "realize-partition",
            "library.json.gz",
            "stock.sqlite",
            str(landscape_path),
            selected.partition_id,
            str(output_path),
            "--max-depth",
            "2",
        )
    )

    assert status == 0
    assert captured == {"partition_id": selected.partition_id, "max_depth": 2}
    assert json.loads(output_path.read_text("utf-8"))["status"] == (
        "unrealized_but_plausible"
    )
    assert json.loads(capsys.readouterr().out)["output_json"] == str(
        output_path.resolve()
    )


def test_partition_assessment_cli_writes_review_and_blind_artifacts(
    monkeypatch,
    tmp_path,
    capsys,
) -> None:
    realization_path = tmp_path / "realization.json"
    realization_path.write_text(
        json.dumps(_full_realization_result().to_dict()),
        encoding="utf-8",
    )
    output_json = tmp_path / "assessment.json"
    output_html = tmp_path / "assessment.html"
    packet_path = tmp_path / "blind_packet.json"
    key_path = tmp_path / "answer_key.json"
    monkeypatch.setattr(
        cli_module,
        "load_generic_library",
        lambda _: _assessment_library(),
    )

    status = cli_module.main(
        (
            "assess-partition",
            "library.json.gz",
            str(realization_path),
            str(output_json),
            str(output_html),
            "--blind-review-packet",
            str(packet_path),
            "--blind-answer-key",
            str(key_path),
        )
    )

    assert status == 0
    assert PartitionAssessmentResult.from_dict(
        json.loads(output_json.read_text("utf-8"))
    )
    assert "Synthetic partition" in output_html.read_text("utf-8")
    assert json.loads(packet_path.read_text("utf-8"))["case_count"] == 1
    assert json.loads(key_path.read_text("utf-8"))["case_count"] == 1
    assert json.loads(capsys.readouterr().out)["ranking_influence"] == (
        "none_review_only"
    )


def test_operator_landscape_combines_independent_views_without_realizing_them() -> None:
    amide = _candidate("amide", AMIDE_FORMATION)
    methylation = _candidate("methylation", N_METHYLATION, score=0.85)
    tactical_variant = _candidate("amide-variant", AMIDE_FORMATION, score=0.8)

    landscape = build_operator_partition_landscape(
        TARGET,
        (amide, methylation, tactical_variant),
    )

    assert landscape.abstained is False
    assert {partition.k for partition in landscape.partitions} == {1, 2, 3}
    combined = next(partition for partition in landscape.partitions if partition.k == 3)
    assert combined.source_kind == "operator_combination_unrealized"
    assert combined.evidence_level == "E1"
    assert "INTERFACE_COMBINATION_NOT_JOINTLY_REALIZED" in combined.warnings
    amide_partition = next(
        partition
        for partition in landscape.partitions
        if any(
            "operator:amide" in interface.candidate_operator_ids
            for interface in partition.interfaces
        )
        and partition.k == 2
    )
    assert any(
        set(interface.candidate_operator_ids)
        == {"operator:amide", "operator:amide-variant"}
        for interface in amide_partition.interfaces
    )
    assert landscape.diagnostics.accepted_seed_count == 2


def test_landscape_abstains_without_a_partition_changing_interface() -> None:
    invalid = _candidate(
        "unverified",
        AMIDE_FORMATION,
        status="unresolved",
    )

    landscape = build_operator_partition_landscape(TARGET, (invalid,))

    assert landscape.abstained is True
    assert landscape.abstention_reasons == (
        "NO_VALIDATED_PARTITION_CHANGING_INTERFACE",
    )
    assert [partition.k for partition in landscape.partitions] == [1]
    assert landscape.diagnostics.rejection_counts == (
        ("candidate_not_forward_verified", 1),
    )


def test_policy_and_static_review_outputs_are_valid(tmp_path) -> None:
    policy = load_synthetic_partition_policy()
    assert policy.minimum_k == 1
    raw_policy = json.loads(
        SYNTHETIC_PARTITION_POLICY_PATH.read_text(encoding="utf-8")
    )
    assert validate_synthetic_partition_policy(raw_policy) == []
    raw_policy["default_k_range"]["minimum"] = 0
    assert "invalid_k_range" in validate_synthetic_partition_policy(raw_policy)

    landscape = build_operator_partition_landscape(
        TARGET,
        (_candidate("amide", AMIDE_FORMATION),),
    )
    assert SyntheticPartitionLandscape.from_dict(landscape.to_dict()) == landscape
    html = render_partition_landscape_html(landscape)
    assert "Synthetic partition landscape" in html
    assert "Role-neutral modules" in html
    json_output = tmp_path / "landscape.json"
    html_output = tmp_path / "landscape.html"
    write_partition_landscape_review(landscape, json_output, html_output)
    assert json.loads(json_output.read_text(encoding="utf-8"))["target_id"]
    assert "<svg" in html_output.read_text(encoding="utf-8")


def test_phase0_baseline_keeps_panels_hash_bound() -> None:
    baseline = json.loads(
        (
            PROJECT_ROOT
            / "docs/new/synthetic_partition_phase0_baseline.v1.json"
        ).read_text(encoding="utf-8")
    )

    assert baseline["test_baseline"]["passed"] == 1193
    assert baseline["test_baseline"]["failed"] == 0
    for panel in baseline["panels"]:
        payload = (PROJECT_ROOT / panel["path"]).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == panel["sha256"]


def test_partition_landscape_cli_writes_review_artifacts(
    monkeypatch,
    tmp_path,
    capsys,
) -> None:
    monkeypatch.setattr(cli_module, "load_generic_library", lambda _: object())
    monkeypatch.setattr(
        cli_module,
        "disconnect_operator_ladder",
        lambda *args, **kwargs: (_candidate("amide", AMIDE_FORMATION),),
    )
    json_output = tmp_path / "cli-landscape.json"
    html_output = tmp_path / "cli-landscape.html"

    status = cli_module.main(
        (
            "partition-landscape",
            "library.json.gz",
            TARGET,
            str(json_output),
            str(html_output),
            "--maximum-k",
            "3",
        )
    )

    assert status == 0
    assert json_output.is_file()
    assert html_output.is_file()
    summary = json.loads(capsys.readouterr().out)
    assert summary["partition_count"] == 2
    assert summary["abstained"] is False
