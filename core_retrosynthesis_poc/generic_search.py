"""Apply and rank structurally diverse generic retrosynthesis templates."""

from __future__ import annotations

import io
import math
from contextlib import redirect_stdout
from functools import lru_cache
from typing import Iterable

from rdkit import Chem
from rdchiral.main import rdchiralReactants, rdchiralReaction, rdchiralRun

from reactive_taxonomy.chemistry.smarts_cache import compile_smarts
from retrosynthesis_poc.chemistry import canonical_smiles, maximum_similarity

from .context import context_similarity
from .generic_compiler import classify_reaction_with_site
from .generic_models import (
    GenericDisconnectionCandidate,
    GenericTemplateLibrary,
)
from .search import _forward_analysis


@lru_cache(maxsize=20_000)
def _reaction(smarts: str) -> object:
    return rdchiralReaction(smarts)


def _apply(smarts: str, target_smiles: str) -> tuple[str, ...]:
    try:
        with redirect_stdout(io.StringIO()):
            outcomes = rdchiralRun(
                _reaction(smarts),
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


def disconnect_generic_target(
    target_smiles: str,
    library: GenericTemplateLibrary,
    *,
    transformations: Iterable[str] = (),
    levels: Iterable[str] = (),
    top_k: int = 20,
    max_templates_to_apply: int = 300,
    max_candidates_to_validate: int = 50,
    use_context: bool = True,
    diversify_sites: bool = False,
) -> tuple[GenericDisconnectionCandidate, ...]:
    """Generate candidates with archetype and forward-chemistry hard filters."""

    if min(top_k, max_templates_to_apply, max_candidates_to_validate) < 1:
        raise ValueError("search limits must be positive")
    canonical_target = canonical_smiles(target_smiles)
    if canonical_target is None or "." in canonical_target:
        raise ValueError("target must be one valid molecule")
    target = Chem.MolFromSmiles(canonical_target)
    if target is None:
        raise ValueError("target could not be parsed")
    allowed_transformations = set(transformations)
    allowed_levels = set(levels)
    applicable = []
    for template in library.templates:
        if allowed_transformations and (
            template.transformation_kind not in allowed_transformations
        ):
            continue
        if allowed_levels and template.abstraction_level not in allowed_levels:
            continue
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
    for product_similarity, specificity, template in applicable[
        :max_templates_to_apply
    ]:
        for precursors in _apply(template.reaction_smarts, canonical_target):
            precursor_similarity = maximum_similarity(
                precursors,
                (
                    precedent.precursor_smiles
                    for precedent in template.precedents
                ),
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
                )
            )
    seeds.sort(key=lambda item: (-item[0], item[1], item[5].template_id))
    candidates: dict[str, GenericDisconnectionCandidate] = {}
    for (
        preliminary,
        precursors,
        precursor_similarity,
        product_similarity,
        specificity,
        template,
    ) in seeds[:max_candidates_to_validate]:
        proposed = f"{precursors}>>{canonical_target}"
        status, _, query_context, center_key = _forward_analysis(
            proposed,
            enabled=True,
        )
        if status in {"invalid", "unresolved"}:
            continue
        classified, site_key = classify_reaction_with_site(proposed)
        if classified != template.transformation_kind:
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
            0.85 * preliminary + 0.15 * context_score
            if use_context
            else preliminary
        )
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
            disconnection_site_key=site_key,
            precedent_reaction_ids=tuple(
                sorted(
                    {
                        precedent.reaction_id
                        for precedent in template.precedents
                        if precedent.reaction_id
                    }
                )
            ),
        )
        current = candidates.get(precursors)
        if current is None or candidate.score > current.score:
            candidates[precursors] = candidate
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
    return ranked[:top_k]


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


__all__ = ["disconnect_generic_target", "rank_site_diverse"]
