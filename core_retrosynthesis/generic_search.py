"""Apply and rank structurally diverse generic retrosynthesis templates."""

from __future__ import annotations

import io
import math
from contextlib import redirect_stdout
from dataclasses import replace
from functools import lru_cache
from typing import Callable, Iterable

from rdkit import Chem
from rdchiral.main import rdchiralReactants, rdchiralReaction, rdchiralRun

from cas_tools import (
    PrecursorRealismAssessment,
    aggregate_precursor_realism_trace,
)
from reactive_taxonomy.chemistry.smarts_cache import compile_smarts
from reactive_taxonomy.strategic_complexity import (
    RetrosyntheticComplexityReduction,
    assess_retrosynthetic_complexity_reduction,
)
from .chemistry import canonical_smiles, maximum_similarity

from .context import context_similarity
from .generic_compiler import analyze_generic_reaction
from .generic_models import (
    GenericDisconnectionCandidate,
    GenericSearchDiagnostics,
    GenericTemplateLibrary,
    OperatorLadderDiagnostics,
)
from .hierarchical_ranking import (
    CompletionPriorIndex,
    build_completion_prior_index,
    rank_hierarchical_candidates,
)
from .precursor_compatibility import assess_precursor_compatibility
from .reaction_compatibility import assess_candidate_reaction_compatibility
from .retrieval_index import indexed_template_ids
from .ranking_policy import (
    RetrosynthesisRankingPolicy,
    diversity_group_key,
    load_retrosynthesis_ranking_policy,
    structural_score_bands,
)
from .search import _forward_analysis
from .selectivity_poc import (
    FunctionalGroupCompetitionWarning,
    detect_functional_group_competition,
)


@lru_cache(maxsize=20_000)
def _reaction(smarts: str) -> object:
    return rdchiralReaction(smarts)


@lru_cache(maxsize=50_000)
def _strategic_complexity(
    mapped_reaction_smiles: str,
) -> RetrosyntheticComplexityReduction | None:
    """Return a cached audit without weakening mandatory graph validation."""

    try:
        return assess_retrosynthetic_complexity_reduction(
            mapped_reaction_smiles
        )
    except ValueError:
        return None


def _apply(smarts: str, target_smiles: str) -> tuple[tuple[str, str], ...]:
    """Return canonical precursors and RDChiral-preserved mapped reactions."""

    try:
        reactants = rdchiralReactants(target_smiles)
        with redirect_stdout(io.StringIO()):
            outcomes, mapped_outcomes = rdchiralRun(
                _reaction(smarts),
                reactants,
                combine_enantiomers=False,
                return_mapped=True,
            )
    except Exception:
        return ()
    mapped_target = Chem.MolToSmiles(
        reactants.reactants,
        canonical=True,
        isomericSmiles=True,
    )
    values = {}
    for outcome in outcomes:
        canonical = canonical_smiles(str(outcome))
        mapped = mapped_outcomes.get(outcome)
        if canonical is None or mapped is None:
            continue
        mapped_precursors = str(mapped[0])
        mapped_reaction = f"{mapped_precursors}>>{mapped_target}"
        current = values.get(canonical)
        if current is None or mapped_reaction < current:
            values[canonical] = mapped_reaction
    return tuple(sorted(values.items()))


def _selectivity_warnings(
    reaction_smiles: str,
) -> tuple[FunctionalGroupCompetitionWarning, ...]:
    """Run the optional audit without changing candidate admission or ranking."""

    try:
        warning = detect_functional_group_competition(reaction_smiles)
    except Exception:
        # This review-only pass must fail open if an unanticipated structure is
        # outside the POC's supported graph-edit topology.
        return ()
    return (warning,) if warning is not None else ()


def disconnect_generic_target_detailed(
    target_smiles: str,
    library: GenericTemplateLibrary,
    *,
    transformations: Iterable[str] = (),
    operator_ids: Iterable[str] = (),
    levels: Iterable[str] = (),
    top_k: int = 20,
    max_templates_to_apply: int = 300,
    max_candidates_to_validate: int = 50,
    use_context: bool = True,
    diversify_sites: bool = False,
    stop_after_valid_candidates: int | None = None,
) -> tuple[
    tuple[GenericDisconnectionCandidate, ...],
    GenericSearchDiagnostics,
]:
    """Generate candidates and expose each retrieval/validation stage."""

    if min(top_k, max_templates_to_apply, max_candidates_to_validate) < 1:
        raise ValueError("search limits must be positive")
    if stop_after_valid_candidates is not None and stop_after_valid_candidates < 1:
        raise ValueError("valid-candidate acceptance target must be positive")
    canonical_target = canonical_smiles(target_smiles)
    if canonical_target is None or "." in canonical_target:
        raise ValueError("target must be one valid molecule")
    target = Chem.MolFromSmiles(canonical_target)
    if target is None:
        raise ValueError("target could not be parsed")
    allowed_transformations = set(transformations)
    allowed_operators = set(operator_ids)
    allowed_levels = set(levels)
    indexed_ids = (
        indexed_template_ids(canonical_target, library.retrieval_index)
        if library.retrieval_index is not None
        else frozenset(template.template_id for template in library.templates)
    )
    applicable = []
    metadata_filtered_count = 0
    for template in library.templates:
        if template.template_id not in indexed_ids:
            continue
        if allowed_operators and template.operator_id not in allowed_operators:
            continue
        if allowed_transformations and (
            template.transformation_kind not in allowed_transformations
        ):
            continue
        if allowed_levels and template.abstraction_level not in allowed_levels:
            continue
        metadata_filtered_count += 1
        pattern = compile_smarts(template.product_smarts, validate=False)
        if pattern is None or not target.HasSubstructMatch(pattern):
            continue
        product_similarity = maximum_similarity(
            canonical_target,
            (precedent.product_smiles for precedent in template.precedents),
        )
        specificity = min(1.0, pattern.GetNumAtoms() / target.GetNumHeavyAtoms())
        applicable.append((product_similarity, specificity, template))
    applicable.sort(
        key=lambda item: (
            -item[0],
            -item[1],
            -item[2].independent_reference_support,
            item[2].template_id,
        )
    )
    seeds = []
    templates_to_apply = applicable[:max_templates_to_apply]
    for product_similarity, specificity, template in templates_to_apply:
        for precursors, mapped_proposed in _apply(
            template.reaction_smarts,
            canonical_target,
        ):
            precursor_similarity = maximum_similarity(
                precursors,
                (precedent.precursor_smiles for precedent in template.precedents),
            )
            support = min(
                1.0,
                math.log1p(template.independent_reference_support) / math.log(11),
            )
            preliminary = (
                0.50 * product_similarity
                + 0.28 * precursor_similarity
                + 0.12 * support
                + 0.10 * specificity
            )
            seeds.append(
                (
                    preliminary,
                    precursors,
                    precursor_similarity,
                    product_similarity,
                    specificity,
                    template,
                    mapped_proposed,
                )
            )
    seeds.sort(key=lambda item: (-item[0], item[1], item[5].template_id))
    candidates: dict[str, GenericDisconnectionCandidate] = {}
    invalid_forward_count = 0
    unresolved_identity_count = 0
    operator_mismatch_count = 0
    validation_seeds = seeds[:max_candidates_to_validate]
    validation_attempt_count = 0
    for (
        preliminary,
        precursors,
        precursor_similarity,
        product_similarity,
        specificity,
        template,
        mapped_proposed,
    ) in validation_seeds:
        validation_attempt_count += 1
        proposed = f"{precursors}>>{canonical_target}"
        status, _, query_context, center_key = _forward_analysis(
            mapped_proposed,
            enabled=True,
        )
        if status in {"invalid", "unresolved"}:
            invalid_forward_count += 1
            continue
        identity = analyze_generic_reaction(mapped_proposed)
        if identity is None:
            unresolved_identity_count += 1
            continue
        if template.operator_signature:
            if identity.operator_signature != template.operator_signature:
                operator_mismatch_count += 1
                continue
        elif identity.named_annotation != template.transformation_kind:
            operator_mismatch_count += 1
            continue
        context_score = max(
            (
                context_similarity(query_context, precedent.context)
                for precedent in template.precedents
                if query_context is not None
            ),
            default=0.0,
        )
        score = (
            0.85 * preliminary + 0.15 * context_score if use_context else preliminary
        )
        compatibility = assess_precursor_compatibility(precursors)
        reaction_compatibility = assess_candidate_reaction_compatibility(
            mapped_proposed
        )
        strategic_complexity = _strategic_complexity(mapped_proposed)
        candidate = GenericDisconnectionCandidate(
            target_smiles=canonical_target,
            precursor_smiles=precursors,
            proposed_reaction_smiles=proposed,
            transformation_kind=template.transformation_kind,
            abstraction_level=template.abstraction_level,
            compiler_engine=template.compiler_engine,
            template_id=template.template_id,
            score=round(min(1.0, score), 8),
            context_similarity=round(context_score, 8),
            product_similarity=round(product_similarity, 8),
            precursor_similarity=round(precursor_similarity, 8),
            template_specificity=round(specificity, 8),
            independent_reference_support=template.independent_reference_support,
            forward_validation_status=status,
            center_transition_key=center_key,
            disconnection_site_key=identity.disconnection_site_key,
            precedent_reaction_ids=tuple(
                sorted(
                    {
                        precedent.reaction_id
                        for precedent in template.precedents
                        if precedent.reaction_id
                    }
                )
            ),
            operator_id=template.operator_id,
            realization_id=template.realization_id,
            operator_signature=identity.operator_signature,
            synthon_signature=identity.synthon_signature,
            condition_query_reaction_smiles=mapped_proposed,
            selectivity_warnings=_selectivity_warnings(mapped_proposed),
            precursor_compatibility_assessments=compatibility.assessments,
            precursor_compatibility_disposition=compatibility.disposition,
            precursor_compatibility_warning_strength=compatibility.warning_strength,
            precursor_compatibility_band_penalty=(
                compatibility.structural_band_penalty
            ),
            precursor_compatibility_policy_definition_id=(
                compatibility.policy_definition_id
            ),
            reaction_compatibility_assessments=(
                reaction_compatibility.assessments
            ),
            reaction_compatibility_disposition=(
                reaction_compatibility.disposition
            ),
            reaction_compatibility_warning_strength=(
                reaction_compatibility.warning_strength
            ),
            reaction_compatibility_band_penalty=(
                reaction_compatibility.structural_band_penalty
            ),
            reaction_compatibility_policy_definition_id=(
                reaction_compatibility.policy_definition_id
            ),
            strategic_complexity=strategic_complexity,
            strategic_complexity_score=(
                strategic_complexity.score
                if strategic_complexity is not None
                else 0.0
            ),
            strategic_class=(
                strategic_complexity.strategic_class
                if strategic_complexity is not None
                else "unresolved"
            ),
            strategic_candidate=(
                strategic_complexity.is_strategic
                if strategic_complexity is not None
                else False
            ),
        )
        current = candidates.get(precursors)
        if current is None or candidate.score > current.score:
            candidates[precursors] = candidate
        if (
            stop_after_valid_candidates is not None
            and len(candidates) >= stop_after_valid_candidates
        ):
            break
    ranked = tuple(
        sorted(
            candidates.values(),
            key=lambda candidate: (
                -candidate.score,
                -candidate.independent_reference_support,
                candidate.precursor_smiles,
            ),
        )
    )
    if diversify_sites:
        ranked = rank_site_diverse(ranked)
    selected = ranked[:top_k]
    return selected, GenericSearchDiagnostics(
        library_template_count=len(library.templates),
        indexed_template_count=len(indexed_ids),
        metadata_filtered_template_count=metadata_filtered_count,
        product_query_match_count=len(applicable),
        applied_template_count=len(templates_to_apply),
        generated_precursor_count=len(seeds),
        validation_attempt_count=validation_attempt_count,
        valid_candidate_count=len(candidates),
        invalid_forward_count=invalid_forward_count,
        unresolved_identity_count=unresolved_identity_count,
        operator_mismatch_count=operator_mismatch_count,
    )


def disconnect_generic_target(
    target_smiles: str,
    library: GenericTemplateLibrary,
    *,
    transformations: Iterable[str] = (),
    operator_ids: Iterable[str] = (),
    levels: Iterable[str] = (),
    top_k: int = 20,
    max_templates_to_apply: int = 300,
    max_candidates_to_validate: int = 50,
    use_context: bool = True,
    diversify_sites: bool = False,
) -> tuple[GenericDisconnectionCandidate, ...]:
    """Generate structurally validated candidates using the product index."""

    candidates, _ = disconnect_generic_target_detailed(
        target_smiles,
        library,
        transformations=transformations,
        operator_ids=operator_ids,
        levels=levels,
        top_k=top_k,
        max_templates_to_apply=max_templates_to_apply,
        max_candidates_to_validate=max_candidates_to_validate,
        use_context=use_context,
        diversify_sites=diversify_sites,
    )
    return candidates


def rank_site_diverse(
    candidates: Iterable[GenericDisconnectionCandidate],
) -> tuple[GenericDisconnectionCandidate, ...]:
    """Interleave ranked precursor forms across distinct product edit sites."""

    groups: dict[str, list[GenericDisconnectionCandidate]] = {}
    order = []
    for candidate in candidates:
        key = candidate.disconnection_site_key or candidate.precursor_smiles
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(candidate)
    values = []
    depth = 0
    while len(values) < sum(len(group) for group in groups.values()):
        for key in order:
            group = groups[key]
            if depth < len(group):
                values.append(group[depth])
        depth += 1
    return tuple(values)


def rank_operator_site_diverse(
    candidates: Iterable[GenericDisconnectionCandidate],
    *,
    policy: RetrosynthesisRankingPolicy | None = None,
) -> tuple[GenericDisconnectionCandidate, ...]:
    """Diversify candidates within structural quality and specificity bands."""

    resolved = policy or load_retrosynthesis_ranking_policy()
    structurally_ranked = tuple(
        sorted(
            candidates,
            key=lambda candidate: (
                (
                    resolved.level_rank(candidate.abstraction_level)
                    if resolved.preserve_abstraction_level_order
                    else 0
                ),
                -candidate.score,
                -candidate.independent_reference_support,
                candidate.precursor_smiles,
                candidate.template_id,
            ),
        )
    )
    bands = structural_score_bands(
        structurally_ranked,
        width=resolved.diversity_score_band_width,
        separate_by_level=resolved.preserve_abstraction_level_order,
    )
    pre_ranks = {
        id(candidate): rank
        for rank, candidate in enumerate(structurally_ranked, start=1)
    }
    realism_enabled = any(
        candidate.precursor_realism_score is not None
        for candidate in structurally_ranked
    )
    realism_penalties = {
        id(candidate): (
            resolved.precursor_realism_band_penalty(
                float(candidate.precursor_realism_score)
            )
            if candidate.precursor_realism_score is not None
            else 0
        )
        for candidate in structurally_ranked
    }
    compatibility_penalties = {
        id(candidate): (
            candidate.precursor_compatibility_band_penalty
            + candidate.reaction_compatibility_band_penalty
        )
        for candidate in structurally_ranked
    }
    effective_bands = {
        id(candidate): (
            bands[id(candidate)]
            + realism_penalties[id(candidate)]
            + compatibility_penalties[id(candidate)]
        )
        for candidate in structurally_ranked
    }
    partitions: dict[
        tuple[int, int],
        dict[tuple[str, ...], list[GenericDisconnectionCandidate]],
    ] = {}
    for candidate in structurally_ranked:
        partition_key = (
            (
                resolved.level_rank(candidate.abstraction_level)
                if resolved.preserve_abstraction_level_order
                else 0
            ),
            effective_bands[id(candidate)],
        )
        groups = partitions.setdefault(partition_key, {})
        groups.setdefault(
            diversity_group_key(candidate, resolved),
            [],
        ).append(candidate)

    realism_ranks: dict[int, int] = {}
    if realism_enabled:
        realism_order = sorted(
            structurally_ranked,
            key=lambda candidate: (
                (
                    resolved.level_rank(candidate.abstraction_level)
                    if resolved.preserve_abstraction_level_order
                    else 0
                ),
                effective_bands[id(candidate)],
                -float(candidate.precursor_realism_score or 0.0),
                pre_ranks[id(candidate)],
            ),
        )
        realism_ranks = {
            id(candidate): rank
            for rank, candidate in enumerate(realism_order, start=1)
        }
    diversified = []
    for partition_key in sorted(partitions):
        group_values = list(partitions[partition_key].values())
        if realism_enabled:
            for group in group_values:
                group.sort(
                    key=lambda candidate: (
                        -float(candidate.precursor_realism_score or 0.0),
                        pre_ranks[id(candidate)],
                    )
                )
            group_values.sort(
                key=lambda group: (
                    -float(group[0].precursor_realism_score or 0.0),
                    min(pre_ranks[id(candidate)] for candidate in group),
                )
            )
        depth = 0
        while True:
            added = False
            for group in group_values:
                if depth < len(group):
                    diversified.append(group[depth])
                    added = True
            if not added:
                break
            depth += 1
    return tuple(
        replace(
            candidate,
            pre_diversity_rank=pre_ranks[id(candidate)],
            diversity_rank=rank,
            diversity_group_key=diversity_group_key(candidate, resolved),
            structural_score_band=bands[id(candidate)],
            ranking_policy_definition_id=resolved.definition_id,
            pre_realism_rank=(
                pre_ranks[id(candidate)] if realism_enabled else 0
            ),
            precursor_realism_rank=(
                realism_ranks[id(candidate)] if realism_enabled else 0
            ),
            precursor_realism_band_penalty=realism_penalties[id(candidate)],
            effective_structural_score_band=effective_bands[id(candidate)],
        )
        for rank, candidate in enumerate(diversified, start=1)
    )


def _attach_precursor_realism(
    candidates: Iterable[GenericDisconnectionCandidate],
    scorer: Callable[[str], tuple[PrecursorRealismAssessment, ...]],
) -> tuple[GenericDisconnectionCandidate, ...]:
    """Attach auditable component assessments before optional reranking."""

    values = []
    for candidate in candidates:
        assessments = scorer(candidate.precursor_smiles)
        aggregation = aggregate_precursor_realism_trace(assessments)
        values.append(
            replace(
                candidate,
                precursor_realism_score=aggregation.score,
                precursor_realism_assessments=assessments,
                precursor_realism_aggregation=aggregation,
            )
        )
    return tuple(values)


def rank_precursor_realism(
    candidates: Iterable[GenericDisconnectionCandidate],
    *,
    policy: RetrosynthesisRankingPolicy | None = None,
) -> tuple[GenericDisconnectionCandidate, ...]:
    """Demote assessed candidates through versioned effective score bands."""

    resolved = policy or load_retrosynthesis_ranking_policy()
    values = tuple(candidates)
    chemistry_order = tuple(
        sorted(
            values,
            key=lambda candidate: (
                (
                    resolved.level_rank(candidate.abstraction_level)
                    if resolved.preserve_abstraction_level_order
                    else 0
                ),
                -candidate.score,
                -candidate.independent_reference_support,
                candidate.precursor_smiles,
                candidate.template_id,
            ),
        )
    )
    bands = structural_score_bands(
        chemistry_order,
        width=resolved.diversity_score_band_width,
        separate_by_level=resolved.preserve_abstraction_level_order,
    )
    pre_ranks = {
        id(candidate): rank
        for rank, candidate in enumerate(chemistry_order, start=1)
    }
    realism_penalties = {
        id(candidate): (
            resolved.precursor_realism_band_penalty(
                float(candidate.precursor_realism_score)
            )
            if candidate.precursor_realism_score is not None
            else 0
        )
        for candidate in chemistry_order
    }
    compatibility_penalties = {
        id(candidate): (
            candidate.precursor_compatibility_band_penalty
            + candidate.reaction_compatibility_band_penalty
        )
        for candidate in chemistry_order
    }
    effective_bands = {
        id(candidate): (
            bands[id(candidate)]
            + realism_penalties[id(candidate)]
            + compatibility_penalties[id(candidate)]
        )
        for candidate in chemistry_order
    }
    ranked = sorted(
        chemistry_order,
        key=lambda candidate: (
            (
                resolved.level_rank(candidate.abstraction_level)
                if resolved.preserve_abstraction_level_order
                else 0
            ),
            effective_bands[id(candidate)],
            -float(candidate.precursor_realism_score or 0.0),
            pre_ranks[id(candidate)],
        ),
    )
    return tuple(
        replace(
            candidate,
            structural_score_band=bands[id(candidate)],
            ranking_policy_definition_id=resolved.definition_id,
            pre_realism_rank=pre_ranks[id(candidate)],
            precursor_realism_rank=rank,
            precursor_realism_band_penalty=realism_penalties[id(candidate)],
            effective_structural_score_band=effective_bands[id(candidate)],
        )
        for rank, candidate in enumerate(ranked, start=1)
    )


def _apply_strategic_candidate_reserve(
    selected: Iterable[GenericDisconnectionCandidate],
    candidates_by_level: dict[
        str,
        tuple[GenericDisconnectionCandidate, ...],
    ],
    *,
    top_k: int,
    policy: RetrosynthesisRankingPolicy,
) -> tuple[GenericDisconnectionCandidate, ...]:
    """Retain bounded scaffold-level alternatives when already generated."""

    values = list(selected)
    if top_k < policy.strategic_reserve_minimum_output_size or not values:
        return tuple(values)
    strategic_count = sum(candidate.strategic_candidate for candidate in values)
    needed = max(0, policy.strategic_reserved_candidates - strategic_count)
    if needed == 0:
        return tuple(values)
    selected_precursors = {candidate.precursor_smiles for candidate in values}
    reserve_pool = sorted(
        (
            candidate
            for level_candidates in candidates_by_level.values()
            for candidate in level_candidates
            if candidate.strategic_candidate
            and candidate.precursor_smiles not in selected_precursors
        ),
        key=lambda candidate: (
            policy.level_rank(candidate.abstraction_level),
            candidate.effective_structural_score_band,
            -candidate.strategic_complexity_score,
            -candidate.score,
            candidate.precursor_smiles,
        ),
    )
    for strategic_candidate in reserve_pool:
        replacement_index = next(
            (
                index
                for index in range(len(values) - 1, -1, -1)
                if not values[index].strategic_candidate
                and values[index].abstraction_level
                == strategic_candidate.abstraction_level
                and strategic_candidate.effective_structural_score_band
                <= (
                    values[index].effective_structural_score_band
                    + policy.strategic_maximum_band_displacement
                )
            ),
            None,
        )
        if replacement_index is None:
            continue
        removed = values[replacement_index]
        selected_precursors.discard(removed.precursor_smiles)
        selected_precursors.add(strategic_candidate.precursor_smiles)
        values[replacement_index] = replace(
            strategic_candidate,
            strategic_reserve_selected=True,
        )
        needed -= 1
        if needed == 0:
            break
    return tuple(values)


def disconnect_operator_ladder(
    target_smiles: str,
    library: GenericTemplateLibrary,
    *,
    top_k: int = 20,
    max_templates_to_apply: int = 500,
    max_candidates_to_validate: int = 100,
    use_context: bool = True,
    include_l0: bool = True,
    diversify: bool = True,
    use_hierarchical_ranking: bool = True,
    minimum_candidates_per_level: int = 0,
    precursor_realism_scorer: (
        Callable[[str], tuple[PrecursorRealismAssessment, ...]] | None
    ) = None,
    completion_prior_index: CompletionPriorIndex | None = None,
) -> tuple[GenericDisconnectionCandidate, ...]:
    """Fill specificity tiers with general operator/site-diverse candidates.

    ``minimum_candidates_per_level`` reserves bounded fallback coverage for
    broader tiers. The default preserves the strict specificity-first one-step
    behavior; multistep search uses a positive reserve so a full L2 tranche
    cannot prevent downstream exploration of L1/L0 alternatives.
    """

    if top_k < 1:
        raise ValueError("top-k must be positive")
    if minimum_candidates_per_level < 0:
        raise ValueError("minimum candidates per level cannot be negative")
    selected = []
    seen = set()
    policy = load_retrosynthesis_ranking_policy()
    hierarchical_prior_index = completion_prior_index
    if diversify and use_hierarchical_ranking and hierarchical_prior_index is None:
        hierarchical_prior_index = build_completion_prior_index(library)
    candidate_pool_size = min(
        max_candidates_to_validate,
        max(top_k, top_k * policy.candidate_pool_multiplier),
    )
    levels = ("L2", "L1", "L0") if include_l0 else ("L2", "L1")
    candidates_by_level: dict[
        str,
        tuple[GenericDisconnectionCandidate, ...],
    ] = {}
    for level in levels:
        if len(selected) >= top_k and minimum_candidates_per_level == 0:
            break
        candidates = disconnect_generic_target(
            target_smiles,
            library,
            levels=(level,),
            top_k=candidate_pool_size,
            max_templates_to_apply=max_templates_to_apply,
            max_candidates_to_validate=max_candidates_to_validate,
            use_context=use_context,
        )
        if precursor_realism_scorer is not None:
            candidates = _attach_precursor_realism(
                candidates,
                precursor_realism_scorer,
            )
        if diversify:
            candidates = rank_operator_site_diverse(
                candidates,
                policy=policy,
            )
            if use_hierarchical_ranking:
                candidates = rank_hierarchical_candidates(
                    candidates,
                    library,
                    structural_policy=policy,
                    prior_index=hierarchical_prior_index,
                )
        elif precursor_realism_scorer is not None:
            candidates = rank_precursor_realism(candidates, policy=policy)
        candidates_by_level[level] = candidates
        if minimum_candidates_per_level > 0:
            continue
        for candidate in candidates:
            if candidate.precursor_smiles in seen:
                continue
            selected.append(candidate)
            seen.add(candidate.precursor_smiles)
            if len(selected) >= top_k:
                break

    if minimum_candidates_per_level > 0:
        for level in levels:
            added = 0
            for candidate in candidates_by_level[level]:
                if candidate.precursor_smiles in seen:
                    continue
                selected.append(candidate)
                seen.add(candidate.precursor_smiles)
                added += 1
                if added >= minimum_candidates_per_level or len(selected) >= top_k:
                    break
            if len(selected) >= top_k:
                break
        for level in levels:
            if len(selected) >= top_k:
                break
            for candidate in candidates_by_level[level]:
                if candidate.precursor_smiles in seen:
                    continue
                selected.append(candidate)
                seen.add(candidate.precursor_smiles)
                if len(selected) >= top_k:
                    break
    if not diversify:
        return tuple(selected)
    return _apply_strategic_candidate_reserve(
        selected,
        candidates_by_level,
        top_k=top_k,
        policy=policy,
    )


def disconnect_operator_ladder_detailed(
    target_smiles: str,
    library: GenericTemplateLibrary,
    *,
    top_k: int = 20,
    max_templates_to_apply: int = 500,
    max_candidates_to_validate: int = 100,
    use_context: bool = True,
    include_l0: bool = True,
    diversify: bool = True,
    use_hierarchical_ranking: bool = True,
    minimum_candidates_per_level: int = 0,
    lazy_validation: bool = False,
    precursor_realism_scorer: (
        Callable[[str], tuple[PrecursorRealismAssessment, ...]] | None
    ) = None,
    completion_prior_index: CompletionPriorIndex | None = None,
) -> tuple[
    tuple[GenericDisconnectionCandidate, ...],
    OperatorLadderDiagnostics,
]:
    """Run the specificity ladder with inspectable staged-validation counts.

    Lazy validation changes only how many ranked proposals receive expensive
    graph validation. Every returned action still passes the same mandatory
    forward and operator-signature checks as exhaustive one-step search.
    """

    if top_k < 1:
        raise ValueError("top-k must be positive")
    if minimum_candidates_per_level < 0:
        raise ValueError("minimum candidates per level cannot be negative")
    selected: list[GenericDisconnectionCandidate] = []
    seen: set[str] = set()
    policy = load_retrosynthesis_ranking_policy()
    hierarchical_prior_index = completion_prior_index
    if diversify and use_hierarchical_ranking and hierarchical_prior_index is None:
        hierarchical_prior_index = build_completion_prior_index(library)
    candidate_pool_size = min(
        max_candidates_to_validate,
        max(top_k, top_k * policy.candidate_pool_multiplier),
    )
    levels = ("L2", "L1", "L0") if include_l0 else ("L2", "L1")
    candidates_by_level: dict[
        str,
        tuple[GenericDisconnectionCandidate, ...],
    ] = {}
    diagnostics_by_level: list[tuple[str, GenericSearchDiagnostics]] = []
    for level in levels:
        if len(selected) >= top_k and minimum_candidates_per_level == 0:
            break
        candidates, diagnostics = disconnect_generic_target_detailed(
            target_smiles,
            library,
            levels=(level,),
            top_k=candidate_pool_size,
            max_templates_to_apply=max_templates_to_apply,
            max_candidates_to_validate=max_candidates_to_validate,
            use_context=use_context,
            # Diversification needs a small verified reservoir; stopping at
            # the final top-k can erase distinct operator/site alternatives.
            # The bounded pool still avoids exhausting large validation
            # budgets (100 by default) once enough actions are available.
            stop_after_valid_candidates=(
                candidate_pool_size if lazy_validation else None
            ),
        )
        diagnostics_by_level.append((level, diagnostics))
        if precursor_realism_scorer is not None:
            candidates = _attach_precursor_realism(
                candidates,
                precursor_realism_scorer,
            )
        if diversify:
            candidates = rank_operator_site_diverse(candidates, policy=policy)
            if use_hierarchical_ranking:
                candidates = rank_hierarchical_candidates(
                    candidates,
                    library,
                    structural_policy=policy,
                    prior_index=hierarchical_prior_index,
                )
        elif precursor_realism_scorer is not None:
            candidates = rank_precursor_realism(candidates, policy=policy)
        candidates_by_level[level] = candidates
        if minimum_candidates_per_level > 0:
            continue
        for candidate in candidates:
            if candidate.precursor_smiles in seen:
                continue
            selected.append(candidate)
            seen.add(candidate.precursor_smiles)
            if len(selected) >= top_k:
                break

    if minimum_candidates_per_level > 0:
        for level in levels:
            added = 0
            for candidate in candidates_by_level.get(level, ()):
                if candidate.precursor_smiles in seen:
                    continue
                selected.append(candidate)
                seen.add(candidate.precursor_smiles)
                added += 1
                if added >= minimum_candidates_per_level or len(selected) >= top_k:
                    break
            if len(selected) >= top_k:
                break
        for level in levels:
            if len(selected) >= top_k:
                break
            for candidate in candidates_by_level.get(level, ()):
                if candidate.precursor_smiles in seen:
                    continue
                selected.append(candidate)
                seen.add(candidate.precursor_smiles)
                if len(selected) >= top_k:
                    break
    diagnostics = OperatorLadderDiagnostics(
        levels_attempted=tuple(level for level, _ in diagnostics_by_level),
        level_diagnostics=tuple(diagnostics_by_level),
    )
    if not diversify:
        return tuple(selected), diagnostics
    return (
        _apply_strategic_candidate_reserve(
            selected,
            candidates_by_level,
            top_k=top_k,
            policy=policy,
        ),
        diagnostics,
    )


__all__ = [
    "disconnect_generic_target",
    "disconnect_generic_target_detailed",
    "disconnect_operator_ladder",
    "disconnect_operator_ladder_detailed",
    "rank_operator_site_diverse",
    "rank_hierarchical_candidates",
    "rank_precursor_realism",
    "rank_site_diverse",
]
