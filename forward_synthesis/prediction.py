"""Retrieve operators, generate products, and audit planned forward steps."""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from functools import lru_cache
from typing import Any, Mapping, Optional, Protocol

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, rdChemReactions

from reactive_taxonomy import (
    BidirectionalReactionOperator,
    REACTION_SIGNATURE_SCHEMA_VERSION,
    apply_forward_operator,
    canonical_molecule_collection,
    featurize_reaction,
    reverse_recovers_precursors,
)

from .library import indexed_forward_operators
from .models import (
    ForwardCompetitionGroup,
    ForwardOperatorLibrary,
    ForwardPredictionResult,
    ForwardProductCandidate,
    ForwardRecipeEvidence,
    ForwardSearchDiagnostics,
    RouteStepForwardAssessment,
)
from .ranking import load_forward_ranking_policy, weighted_forward_score


class RecipeAssessorProtocol(Protocol):
    """Narrow optional condition-compatibility interface."""

    def assess(self, reaction_smiles: str, recipe: Mapping[str, Any]) -> Any: ...


def _digest(namespace: str, *values: str) -> str:
    payload = "\0".join(values).encode("utf-8")
    return f"{namespace}:{hashlib.sha256(payload).hexdigest()}"


@lru_cache(maxsize=1)
def _fingerprint_generator() -> Any:
    return rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


@lru_cache(maxsize=50_000)
def _fingerprint(smiles: str) -> Any | None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    return _fingerprint_generator().GetFingerprint(molecule)


def _maximum_similarity(query: str, references: tuple[str, ...]) -> float:
    query_fp = _fingerprint(query)
    if query_fp is None:
        return 0.0
    scores = []
    for reference in references:
        reference_fp = _fingerprint(reference)
        if reference_fp is not None:
            scores.append(float(DataStructs.TanimotoSimilarity(query_fp, reference_fp)))
    return max(scores, default=0.0)


def _template_specificity(
    operator: BidirectionalReactionOperator,
    canonical_starting_materials: str,
) -> float:
    reaction = rdChemReactions.ReactionFromSmarts(operator.forward_smarts)
    molecule = Chem.MolFromSmiles(canonical_starting_materials)
    if reaction is None or molecule is None:
        return 0.0
    template_atoms = sum(
        reaction.GetReactantTemplate(index).GetNumAtoms()
        for index in range(reaction.GetNumReactantTemplates())
    )
    heavy_atoms = molecule.GetNumHeavyAtoms()
    return min(1.0, template_atoms / heavy_atoms) if heavy_atoms else 0.0


def _support_score(count: int) -> float:
    return min(1.0, math.log1p(max(0, count)) / math.log(11))


def _recipe_evidence(
    reaction_smiles: str,
    recipe: Mapping[str, Any] | None,
    assessor: RecipeAssessorProtocol | None,
) -> ForwardRecipeEvidence:
    if recipe is None:
        return ForwardRecipeEvidence(False, None, None)
    try:
        if assessor is not None:
            raw = assessor.assess(reaction_smiles, recipe)
        else:
            from condition_recommender import assess_reaction_recipe

            raw = assess_reaction_recipe(reaction_smiles, recipe)
    except Exception as exc:
        return ForwardRecipeEvidence(
            evaluated=False,
            compatible=None,
            score=None,
            cautions=(f"RECIPE_ASSESSMENT_FAILED:{type(exc).__name__}",),
        )
    compatible = bool(getattr(raw, "compatible", False))
    score = getattr(raw, "score", None)
    evidence = tuple(str(item) for item in getattr(raw, "evidence", ()) or ())
    conflicts = tuple(str(item) for item in getattr(raw, "hard_conflicts", ()) or ())
    return ForwardRecipeEvidence(
        evaluated=True,
        compatible=compatible,
        score=float(score) if score is not None else None,
        hard_conflicts=conflicts,
        cautions=evidence,
        definition_id=str(getattr(raw, "definition_id", "") or ""),
        definition_version=str(getattr(raw, "definition_version", "") or ""),
    )


@dataclass(frozen=True)
class _CandidateSeed:
    candidate: ForwardProductCandidate
    structural_score: float


def _validated_seed(
    operator: BidirectionalReactionOperator,
    outcome: Any,
    recipe: Mapping[str, Any] | None,
    assessor: RecipeAssessorProtocol | None,
) -> tuple[_CandidateSeed | None, str | None]:
    if not reverse_recovers_precursors(operator, outcome):
        return None, "reverse_round_trip"
    reaction_smiles = (
        f"{outcome.participating_precursor_smiles}>>{outcome.product_smiles}"
    )
    mapped_reaction_smiles = (
        f"{outcome.mapped_participating_precursor_smiles}"
        f">>{outcome.mapped_product_smiles}"
    )
    analysis = featurize_reaction(mapped_reaction_smiles)
    if not analysis.valid or analysis.reaction_core is None:
        return None, "invalid_reaction"
    if analysis.reaction_signature is None:
        return None, "missing_signature"
    observed_tokens = tuple(analysis.reaction_core.edit_tokens)
    agreement = bool(operator.edit_tokens) and set(operator.edit_tokens).issubset(
        observed_tokens
    )
    if not agreement:
        return None, "operator_edit_mismatch"
    recipe_evidence = _recipe_evidence(
        mapped_reaction_smiles,
        recipe,
        assessor,
    )
    if recipe_evidence.evaluated and recipe_evidence.compatible is False:
        return None, "recipe_conflict"
    precursor_similarity = _maximum_similarity(
        outcome.participating_precursor_smiles,
        tuple(precedent.precursor_smiles for precedent in operator.precedents),
    )
    specificity = _template_specificity(
        operator,
        outcome.participating_precursor_smiles,
    )
    raw_components: dict[str, float | None] = {
        "precursor_similarity": precursor_similarity,
        "template_specificity": specificity,
        "independent_support": _support_score(operator.independent_reference_support),
        "operator_edit_agreement": 1.0,
        "recipe_compatibility": recipe_evidence.score,
    }
    policy = load_forward_ranking_policy()
    score, contributions = weighted_forward_score(raw_components, policy)
    structural_components = dict(raw_components)
    structural_components["recipe_compatibility"] = None
    structural_score, _ = weighted_forward_score(structural_components, policy)
    signature = analysis.reaction_signature
    pathway_id = _digest(
        "FWP1",
        operator.forward_operator_id,
        outcome.application_id,
        signature.signature_id,
    )
    warnings = tuple(sorted(set(analysis.warnings)))
    candidate = ForwardProductCandidate(
        rank=0,
        product_smiles=outcome.product_smiles,
        reaction_smiles=reaction_smiles,
        mapped_reaction_smiles=mapped_reaction_smiles,
        pathway_id=pathway_id,
        operator_id=operator.operator_id,
        forward_operator_id=operator.forward_operator_id,
        realization_id=operator.realization_id,
        template_id=operator.template_id,
        abstraction_level=operator.abstraction_level,
        participating_component_indices=outcome.participating_component_indices,
        participating_precursor_smiles=outcome.participating_precursor_smiles,
        assignment=outcome.assignment,
        atom_correspondence=outcome.atom_correspondence,
        score=score,
        score_components=contributions,
        structural_score_band=int(
            max(0.0, 1.0 - structural_score) / policy.score_band_width
        ),
        reverse_round_trip=True,
        reaction_signature_id=signature.signature_id,
        reaction_signature_schema_version=REACTION_SIGNATURE_SCHEMA_VERSION,
        operator_edit_agreement=True,
        observed_edit_tokens=observed_tokens,
        independent_reference_support=operator.independent_reference_support,
        observation_support=operator.observation_support,
        precedent_reaction_ids=tuple(
            sorted(
                {
                    precedent.reaction_id
                    for precedent in operator.precedents
                    if precedent.reaction_id
                }
            )
        ),
        precedent_reference_ids=tuple(
            sorted(
                {
                    precedent.reference_id
                    for precedent in operator.precedents
                    if precedent.reference_id
                }
            )
        ),
        recipe_evidence=recipe_evidence,
        named_annotations=operator.named_annotations,
        warnings=warnings,
    )
    return _CandidateSeed(candidate, structural_score), None


def _competition_groups(
    candidates: tuple[ForwardProductCandidate, ...],
) -> tuple[ForwardCompetitionGroup, ...]:
    operator_groups: dict[str, list[ForwardProductCandidate]] = defaultdict(list)
    site_groups: dict[str, list[ForwardProductCandidate]] = defaultdict(list)
    for candidate in candidates:
        all_operators = (candidate.operator_id, *candidate.alternative_operator_ids)
        for operator_id in all_operators:
            operator_groups[operator_id].append(candidate)
        site_key = _digest(
            "FWS1",
            candidate.operator_id,
            json.dumps(candidate.participating_component_indices),
            json.dumps(candidate.assignment),
        )
        site_groups[site_key].append(candidate)
    groups = []
    for level, values in (("operator", operator_groups), ("site", site_groups)):
        for key, members in sorted(values.items()):
            groups.append(
                ForwardCompetitionGroup(
                    competition_level=level,
                    group_key=key,
                    candidate_ranks=tuple(sorted({item.rank for item in members})),
                    product_smiles=tuple(
                        sorted({item.product_smiles for item in members})
                    ),
                    operator_ids=tuple(
                        sorted(
                            {
                                operator_id
                                for item in members
                                for operator_id in (
                                    item.operator_id,
                                    *item.alternative_operator_ids,
                                )
                            }
                        )
                    ),
                )
            )
    groups.append(
        ForwardCompetitionGroup(
            competition_level="product",
            group_key=_digest("FWC1", *(item.product_smiles for item in candidates)),
            candidate_ranks=tuple(item.rank for item in candidates),
            product_smiles=tuple(item.product_smiles for item in candidates),
            operator_ids=tuple(
                sorted(
                    {
                        operator_id
                        for item in candidates
                        for operator_id in (
                            item.operator_id,
                            *item.alternative_operator_ids,
                        )
                    }
                )
            ),
        )
    )
    return tuple(groups)


def predict_products(
    starting_materials: str,
    library: ForwardOperatorLibrary,
    *,
    recipe: Mapping[str, Any] | None = None,
    recipe_assessor: RecipeAssessorProtocol | None = None,
    operator_ids: tuple[str, ...] = (),
    levels: tuple[str, ...] = (),
    top_k: int = 20,
    max_operators_to_apply: int = 300,
    max_assignments_per_operator: int = 128,
    max_outcomes_per_operator: int = 256,
) -> ForwardPredictionResult:
    """Blindly enumerate and rank products from precursor-observable evidence."""

    if (
        min(
            top_k,
            max_operators_to_apply,
            max_assignments_per_operator,
            max_outcomes_per_operator,
        )
        < 1
    ):
        raise ValueError("forward search limits must be positive")
    canonical = canonical_molecule_collection(starting_materials)
    if canonical is None:
        return ForwardPredictionResult(
            query_starting_materials=starting_materials,
            canonical_starting_materials="",
            conditions_supplied=recipe is not None,
            valid=False,
            status="invalid_query",
            candidates=(),
            competition_groups=(),
            diagnostics=ForwardSearchDiagnostics(
                library_operator_count=len(library.operators)
            ),
            error="INVALID_STARTING_MATERIALS",
        )
    allowed_operators = set(operator_ids)
    allowed_levels = set(levels)
    retrieved = tuple(
        operator
        for operator in indexed_forward_operators(canonical, library)
        if (
            not allowed_operators
            or operator.operator_id in allowed_operators
            or operator.forward_operator_id in allowed_operators
        )
        and (not allowed_levels or operator.abstraction_level in allowed_levels)
    )
    policy = load_forward_ranking_policy()
    ordered_operators = sorted(
        retrieved,
        key=lambda item: (
            policy.level_rank(item.abstraction_level),
            -item.independent_reference_support,
            item.forward_operator_id,
        ),
    )[:max_operators_to_apply]
    counters = Counter()
    seeds = []
    for operator in ordered_operators:
        outcomes = apply_forward_operator(
            operator,
            canonical,
            max_assignments=max_assignments_per_operator,
            max_outcomes=max_outcomes_per_operator,
        )
        if outcomes:
            counters["applied"] += 1
        counters["generated"] += len(outcomes)
        for outcome in outcomes:
            seed, rejection = _validated_seed(
                operator,
                outcome,
                recipe,
                recipe_assessor,
            )
            if rejection:
                counters[rejection] += 1
            elif seed is not None:
                seeds.append(seed)
    seeds.sort(
        key=lambda item: (
            policy.level_rank(item.candidate.abstraction_level),
            item.candidate.structural_score_band,
            -(
                item.candidate.recipe_evidence.score
                if item.candidate.recipe_evidence.score is not None
                else -1.0
            ),
            -item.structural_score,
            -item.candidate.score,
            item.candidate.product_smiles,
            item.candidate.pathway_id,
        )
    )
    by_product: dict[str, list[_CandidateSeed]] = defaultdict(list)
    for seed in seeds:
        by_product[seed.candidate.product_smiles].append(seed)
    representatives = []
    for product_smiles, product_seeds in by_product.items():
        representative = product_seeds[0].candidate
        alternatives = tuple(seed.candidate for seed in product_seeds[1:])
        representatives.append(
            replace(
                representative,
                alternative_pathway_ids=tuple(item.pathway_id for item in alternatives),
                alternative_operator_ids=tuple(
                    sorted(
                        {
                            item.operator_id
                            for item in alternatives
                            if item.operator_id != representative.operator_id
                        }
                    )
                ),
                alternative_template_ids=tuple(
                    sorted(
                        {
                            item.template_id
                            for item in alternatives
                            if item.template_id != representative.template_id
                        }
                    )
                ),
            )
        )
    representatives.sort(
        key=lambda item: (
            policy.level_rank(item.abstraction_level),
            item.structural_score_band,
            -(
                item.recipe_evidence.score
                if item.recipe_evidence.score is not None
                else -1.0
            ),
            -item.score,
            item.product_smiles,
        )
    )
    candidates = tuple(
        replace(candidate, rank=rank)
        for rank, candidate in enumerate(representatives[:top_k], start=1)
    )
    diagnostics = ForwardSearchDiagnostics(
        library_operator_count=len(library.operators),
        indexed_operator_count=len(retrieved),
        applied_operator_count=counters["applied"],
        generated_outcome_count=counters["generated"],
        reverse_round_trip_failure_count=counters["reverse_round_trip"],
        invalid_reaction_count=counters["invalid_reaction"],
        missing_signature_count=counters["missing_signature"],
        operator_edit_mismatch_count=counters["operator_edit_mismatch"],
        recipe_conflict_count=counters["recipe_conflict"],
        valid_pathway_count=len(seeds),
        unique_product_count=len(by_product),
    )
    warnings = []
    if recipe is None:
        warnings.append("CONDITIONS_NOT_SUPPLIED_PRODUCTS_ARE_POSSIBILITIES")
    if not candidates:
        warnings.append("NO_FORWARD_PRODUCT_GENERATED")
    return ForwardPredictionResult(
        query_starting_materials=starting_materials,
        canonical_starting_materials=canonical,
        conditions_supplied=recipe is not None,
        valid=True,
        status="predicted" if candidates else "no_supported_product",
        candidates=candidates,
        competition_groups=_competition_groups(candidates) if candidates else (),
        diagnostics=diagnostics,
        warnings=tuple(warnings),
        ranking_definition_id=policy.definition_id,
    )


def _without_stereo(smiles: str) -> Optional[str]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    Chem.RemoveStereochemistry(molecule)
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)


def _match_kind(expected: str, candidate: str) -> str:
    if expected == candidate:
        return "exact"
    if _without_stereo(expected) == _without_stereo(candidate):
        return "stereo_relaxed"
    return "absent"


def assess_proposed_step(
    starting_materials: str,
    intended_product: str,
    library: ForwardOperatorLibrary,
    *,
    recipe: Mapping[str, Any] | None = None,
    recipe_assessor: RecipeAssessorProtocol | None = None,
    operator_hint: str | None = None,
    levels: tuple[str, ...] = (),
    top_k: int = 20,
) -> RouteStepForwardAssessment:
    """Audit a route step using separate targeted and target-blind passes."""

    canonical_target = canonical_molecule_collection(intended_product)
    blind = predict_products(
        starting_materials,
        library,
        recipe=recipe,
        recipe_assessor=recipe_assessor,
        levels=levels,
        top_k=top_k,
    )
    if canonical_target is None or "." in canonical_target:
        return RouteStepForwardAssessment(
            starting_materials=starting_materials,
            intended_product=intended_product,
            intended_match="invalid",
            targeted_replay_status="invalid_intended_product",
            intended_product_rank=None,
            best_competitor_product=(
                blind.candidates[0].product_smiles if blind.candidates else None
            ),
            score_margin=None,
            disposition="out_of_scope",
            blind_prediction=blind,
            operator_hint=operator_hint,
            warnings=("INVALID_INTENDED_PRODUCT",),
        )
    hinted = tuple(
        operator
        for operator in library.operators
        if operator_hint
        and operator_hint in {operator.operator_id, operator.forward_operator_id}
    )
    if operator_hint and not hinted:
        targeted_status = "operator_hint_not_found"
    elif hinted:
        targeted_matches = []
        for operator in hinted:
            for outcome in apply_forward_operator(operator, starting_materials):
                if _match_kind(canonical_target, outcome.product_smiles) != "absent":
                    targeted_matches.append(outcome)
        targeted_status = (
            "structurally_reproduced" if targeted_matches else "not_reproduced"
        )
    else:
        targeted_status = "not_requested"
    matched = tuple(
        candidate
        for candidate in blind.candidates
        if _match_kind(canonical_target, candidate.product_smiles) != "absent"
    )
    match_kind = (
        _match_kind(canonical_target, matched[0].product_smiles)
        if matched
        else "absent"
    )
    intended_rank = matched[0].rank if matched else None
    competitors = tuple(
        candidate
        for candidate in blind.candidates
        if _match_kind(canonical_target, candidate.product_smiles) == "absent"
    )
    best_competitor = competitors[0] if competitors else None
    margin = (
        round(matched[0].score - best_competitor.score, 8)
        if matched and best_competitor is not None
        else None
    )
    policy = load_forward_ranking_policy()
    if targeted_status in {"not_reproduced", "operator_hint_not_found"}:
        disposition = "structurally_inconsistent"
    elif not blind.candidates:
        disposition = "out_of_scope"
    elif not matched:
        disposition = "unsupported"
    elif intended_rank == 1 and (
        best_competitor is None
        or margin is not None
        and margin >= policy.score_band_width
    ):
        disposition = "clear"
    else:
        disposition = "competitive"
    warnings = []
    if disposition == "competitive":
        warnings.append("COMPETING_FORWARD_PRODUCTS")
    elif disposition == "unsupported":
        warnings.append("INTENDED_PRODUCT_NOT_FOUND_BY_BLIND_FORWARD_SEARCH")
    elif disposition == "structurally_inconsistent":
        warnings.append("RETROSYNTHESIS_OPERATOR_FAILED_TARGETED_FORWARD_REPLAY")
    return RouteStepForwardAssessment(
        starting_materials=blind.canonical_starting_materials or starting_materials,
        intended_product=canonical_target,
        intended_match=match_kind,
        targeted_replay_status=targeted_status,
        intended_product_rank=intended_rank,
        best_competitor_product=(
            best_competitor.product_smiles if best_competitor else None
        ),
        score_margin=margin,
        disposition=disposition,
        blind_prediction=blind,
        operator_hint=operator_hint,
        warnings=tuple(warnings),
    )


__all__ = [
    "RecipeAssessorProtocol",
    "assess_proposed_step",
    "predict_products",
]
