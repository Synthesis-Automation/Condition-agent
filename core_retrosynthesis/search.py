"""Apply, validate, and rank core-derived retrosynthetic templates."""

from __future__ import annotations

import io
import math
from contextlib import redirect_stdout
from dataclasses import dataclass
from functools import lru_cache
from typing import Iterable

from rdkit import Chem
from rdchiral.main import rdchiralReactants, rdchiralReaction, rdchiralRun

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.chemistry.smarts_cache import compile_smarts
from .chemistry import canonical_smiles, digest, maximum_similarity

from .context import context_from_analysis, context_similarity
from .models import (
    CoreDisconnectionCandidate,
    CoreTemplate,
    CoreTemplateLibrary,
    TemplateContext,
)


@dataclass(frozen=True)
class _CandidateSeed:
    precursor_smiles: str
    template: CoreTemplate
    product_similarity: float
    template_specificity: float
    preliminary_score: float


@lru_cache(maxsize=20_000)
def _compiled_reaction(reaction_smarts: str) -> object:
    return rdchiralReaction(reaction_smarts)


def _matches_product(target: object, template: CoreTemplate) -> bool:
    pattern = compile_smarts(template.product_smarts, validate=False)
    return bool(pattern is not None and target.HasSubstructMatch(pattern))


def _specificity(target: object, template: CoreTemplate) -> float:
    pattern = compile_smarts(template.product_smarts, validate=False)
    if pattern is None or target.GetNumHeavyAtoms() == 0:
        return 0.0
    return min(1.0, pattern.GetNumAtoms() / target.GetNumHeavyAtoms())


def _apply(template: CoreTemplate, target_smiles: str) -> tuple[str, ...]:
    try:
        with redirect_stdout(io.StringIO()):
            outcomes = rdchiralRun(
                _compiled_reaction(template.reaction_smarts),
                rdchiralReactants(target_smiles),
            )
    except Exception:
        return ()
    return tuple(
        sorted(
            {
                canonical
                for outcome in outcomes
                if (canonical := canonical_smiles(str(outcome))) is not None
            }
        )
    )


@lru_cache(maxsize=50_000)
def _forward_analysis(
    reaction_smiles: str,
    *,
    enabled: bool,
) -> tuple[str, tuple[str, ...], TemplateContext | None, str]:
    if not enabled:
        return "not_run", (), None, ""
    analysis = featurize_reaction(reaction_smiles)
    if not analysis.valid:
        return "invalid", tuple(analysis.warnings), None, ""
    context = context_from_analysis(analysis, reaction_smiles)
    center_key = (
        analysis.reaction_core.center_transition_key
        if analysis.reaction_core is not None
        else ""
    )
    if analysis.reaction_signature is not None:
        return "verified_signature", tuple(analysis.warnings), context, center_key
    if analysis.reaction_core is not None:
        return "core_only", tuple(analysis.warnings), context, center_key
    return "unresolved", tuple(analysis.warnings), context, center_key


def _support_score(template: CoreTemplate) -> float:
    return min(1.0, math.log1p(template.operator_reference_support) / math.log(11))


def _maximum_context_similarity(
    query: TemplateContext | None,
    template: CoreTemplate,
) -> float:
    if query is None:
        return 0.0
    return max(
        (
            context_similarity(query, precedent.context)
            for precedent in template.precedents
        ),
        default=0.0,
    )


def _candidate(
    *,
    target_smiles: str,
    precursor_smiles: str,
    template: CoreTemplate,
    product_similarity: float,
    template_specificity: float,
    validate_forward: bool,
    use_context: bool,
) -> CoreDisconnectionCandidate:
    precursor_similarity = maximum_similarity(
        precursor_smiles,
        (precedent.precursor_smiles for precedent in template.precedents),
    )
    reaction_smiles = f"{precursor_smiles}>>{target_smiles}"
    status, warnings, query_context, center_key = _forward_analysis(
        reaction_smiles,
        enabled=validate_forward,
    )
    validation_bonus = {
        "verified_signature": 1.0,
        "core_only": 0.5,
        "unresolved": 0.0,
        "invalid": 0.0,
        "not_run": 0.0,
    }[status]
    base_score = (
        0.45 * product_similarity
        + 0.25 * precursor_similarity
        + 0.15 * _support_score(template)
        + 0.10 * template_specificity
        + 0.05 * validation_bonus
    )
    context_score = _maximum_context_similarity(query_context, template)
    score = 0.85 * base_score + 0.15 * context_score if use_context else base_score
    candidate_id = digest(
        "CRD1",
        target_smiles,
        precursor_smiles,
        template.template_id,
        "context" if use_context else "neutral",
    )
    return CoreDisconnectionCandidate(
        candidate_id=candidate_id,
        target_smiles=target_smiles,
        precursor_smiles=precursor_smiles,
        proposed_reaction_smiles=reaction_smiles,
        bond_kind=template.bond_kind,
        abstraction_level=template.abstraction_level,
        template_id=template.template_id,
        operator_id=template.operator_id,
        score=round(min(1.0, score), 8),
        base_score=round(min(1.0, base_score), 8),
        context_similarity=round(context_score, 8),
        product_similarity=round(product_similarity, 8),
        precursor_similarity=round(precursor_similarity, 8),
        template_specificity=round(template_specificity, 8),
        observation_support=template.operator_observation_support,
        independent_reference_support=template.operator_reference_support,
        forward_validation_status=status,
        center_transition_key=center_key,
        precedent_reaction_ids=tuple(
            sorted(
                {
                    precedent.reaction_id
                    for precedent in template.precedents
                    if precedent.reaction_id
                }
            )
        ),
        precedent_reference_ids=tuple(
            sorted(
                {
                    precedent.reference_id
                    for precedent in template.precedents
                    if precedent.reference_id
                }
            )
        ),
        warnings=tuple(sorted(set(warnings))),
    )


def _preliminary_score(
    precursor_smiles: str,
    template: CoreTemplate,
    product_similarity: float,
    template_specificity: float,
) -> float:
    precursor_similarity = maximum_similarity(
        precursor_smiles,
        (precedent.precursor_smiles for precedent in template.precedents),
    )
    return (
        0.45 * product_similarity
        + 0.25 * precursor_similarity
        + 0.15 * _support_score(template)
        + 0.10 * template_specificity
    )


def disconnect_target(
    target_smiles: str,
    library: CoreTemplateLibrary,
    *,
    allowed_bonds: Iterable[str] = ("C-N", "C-O", "C-S"),
    levels: Iterable[str] = ("L1", "L2"),
    max_templates_to_apply: int = 250,
    top_k: int = 20,
    max_candidates_to_validate: int = 50,
    validate_forward: bool = True,
    use_context: bool = True,
) -> tuple[CoreDisconnectionCandidate, ...]:
    """Generate ranked precursor candidates from core-derived SMARTS."""

    if (
        max_templates_to_apply < 1
        or top_k < 1
        or max_candidates_to_validate < top_k
    ):
        raise ValueError("template and result limits must be positive")
    canonical_target = canonical_smiles(target_smiles)
    if canonical_target is None or "." in canonical_target:
        raise ValueError("target must be one valid molecule")
    target = Chem.MolFromSmiles(canonical_target)
    if target is None:
        raise ValueError("target could not be parsed")
    allowed = set(allowed_bonds)
    selected_levels = set(levels)
    applicable = []
    for template in library.templates:
        if (
            template.bond_kind not in allowed
            or template.abstraction_level not in selected_levels
            or not _matches_product(target, template)
        ):
            continue
        product_similarity = maximum_similarity(
            canonical_target,
            (precedent.product_smiles for precedent in template.precedents),
        )
        specificity = _specificity(target, template)
        applicable.append((product_similarity, specificity, template))
    applicable.sort(
        key=lambda value: (
            -value[0],
            -value[1],
            -value[2].operator_reference_support,
            value[2].template_id,
        )
    )

    seeds: dict[tuple[str, str], _CandidateSeed] = {}
    for product_similarity, specificity, template in applicable[
        :max_templates_to_apply
    ]:
        for precursors in _apply(template, canonical_target):
            seed = _CandidateSeed(
                precursor_smiles=precursors,
                template=template,
                product_similarity=product_similarity,
                template_specificity=specificity,
                preliminary_score=_preliminary_score(
                    precursors,
                    template,
                    product_similarity,
                    specificity,
                ),
            )
            key = (precursors, template.template_id)
            current = seeds.get(key)
            if current is None or (
                seed.preliminary_score,
                seed.template.operator_reference_support,
            ) > (
                current.preliminary_score,
                current.template.operator_reference_support,
            ):
                seeds[key] = seed
    shortlist = sorted(
        seeds.values(),
        key=lambda value: (
            -value.preliminary_score,
            -value.template.operator_reference_support,
            value.precursor_smiles,
            value.template.template_id,
        ),
    )[:max_candidates_to_validate]
    candidates: dict[str, CoreDisconnectionCandidate] = {}
    for seed in shortlist:
        candidate = _candidate(
            target_smiles=canonical_target,
            precursor_smiles=seed.precursor_smiles,
            template=seed.template,
            product_similarity=seed.product_similarity,
            template_specificity=seed.template_specificity,
            validate_forward=validate_forward,
            use_context=use_context,
        )
        if validate_forward and candidate.forward_validation_status in {
            "invalid",
            "unresolved",
        }:
            continue
        current = candidates.get(seed.precursor_smiles)
        if current is None or (
            candidate.score,
            candidate.independent_reference_support,
            candidate.template_id,
        ) > (
            current.score,
            current.independent_reference_support,
            current.template_id,
        ):
            candidates[seed.precursor_smiles] = candidate
    return tuple(
        sorted(
            candidates.values(),
            key=lambda value: (
                -value.score,
                -value.independent_reference_support,
                value.precursor_smiles,
                value.template_id,
            ),
        )[:top_k]
    )


__all__ = ["disconnect_target"]
