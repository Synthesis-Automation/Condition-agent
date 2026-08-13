"""Target-side matching, inverse application, validation, and ranking."""

from __future__ import annotations

import io
import math
from contextlib import redirect_stdout
from functools import lru_cache
from typing import Iterable

from rdkit import Chem
from rdchiral.main import rdchiralReactants, rdchiralReaction, rdchiralRun

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.chemistry.smarts_cache import compile_smarts

from ...chemistry import canonical_smiles, digest, maximum_similarity
from .models import CxTemplate, DisconnectionCandidate, RetrosynthesisLibrary


@lru_cache(maxsize=20_000)
def _compiled_reaction(reaction_smarts: str) -> object:
    return rdchiralReaction(reaction_smarts)


def _matches_product(target: object, template: CxTemplate) -> bool:
    pattern = compile_smarts(template.product_smarts, validate=False)
    return bool(pattern is not None and target.HasSubstructMatch(pattern))


def _specificity(target: object, template: CxTemplate) -> float:
    pattern = compile_smarts(template.product_smarts, validate=False)
    if pattern is None or target.GetNumHeavyAtoms() == 0:
        return 0.0
    return min(1.0, pattern.GetNumAtoms() / target.GetNumHeavyAtoms())


def _apply(template: CxTemplate, target_smiles: str) -> tuple[str, ...]:
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


def _forward_status(reaction_smiles: str, *, enabled: bool) -> tuple[str, tuple[str, ...]]:
    if not enabled:
        return "not_run", ()
    analysis = featurize_reaction(reaction_smiles)
    if not analysis.valid:
        return "invalid", tuple(analysis.warnings)
    if analysis.reaction_signature is not None:
        return "verified_signature", tuple(analysis.warnings)
    if analysis.reaction_core is not None:
        return "core_only", tuple(analysis.warnings)
    return "unresolved", tuple(analysis.warnings)


def _support_score(template: CxTemplate) -> float:
    return min(1.0, math.log1p(template.independent_reference_support) / math.log(11))


def _candidate(
    *,
    target_smiles: str,
    precursor_smiles: str,
    template: CxTemplate,
    product_similarity: float,
    template_specificity: float,
    validate_forward: bool,
) -> DisconnectionCandidate:
    precursor_similarity = maximum_similarity(
        precursor_smiles,
        (precedent.precursor_smiles for precedent in template.precedents),
    )
    status, warnings = _forward_status(
        f"{precursor_smiles}>>{target_smiles}",
        enabled=validate_forward,
    )
    validation_bonus = {
        "verified_signature": 1.0,
        "core_only": 0.5,
        "unresolved": 0.0,
        "invalid": 0.0,
        "not_run": 0.0,
    }[status]
    score = (
        0.45 * product_similarity
        + 0.25 * precursor_similarity
        + 0.15 * _support_score(template)
        + 0.10 * template_specificity
        + 0.05 * validation_bonus
    )
    candidate_id = digest(
        "CXD1", target_smiles, precursor_smiles, template.template_id
    )
    return DisconnectionCandidate(
        candidate_id=candidate_id,
        target_smiles=target_smiles,
        precursor_smiles=precursor_smiles,
        proposed_reaction_smiles=f"{precursor_smiles}>>{target_smiles}",
        bond_kind=template.bond_kind,
        template_id=template.template_id,
        score=round(min(1.0, score), 8),
        product_similarity=round(product_similarity, 8),
        precursor_similarity=round(precursor_similarity, 8),
        template_specificity=round(template_specificity, 8),
        observation_support=template.observation_support,
        independent_reference_support=template.independent_reference_support,
        forward_validation_status=status,  # type: ignore[arg-type]
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
        precedent_support_unit_ids=tuple(
            sorted({precedent.support_unit_id for precedent in template.precedents})
        ),
        warnings=tuple(sorted(set(warnings))),
    )


def disconnect_target(
    target_smiles: str,
    library: RetrosynthesisLibrary,
    *,
    allowed_bonds: Iterable[str] = ("C-N", "C-O", "C-S"),
    max_templates_to_apply: int = 250,
    top_k: int = 20,
    validate_forward: bool = True,
) -> tuple[DisconnectionCandidate, ...]:
    """Generate ranked, precedent-backed one-step precursor candidates."""

    if max_templates_to_apply < 1 or top_k < 1:
        raise ValueError("template and result limits must be positive")
    canonical_target = canonical_smiles(target_smiles)
    if canonical_target is None or "." in canonical_target:
        raise ValueError("target must be one valid molecule")
    target = Chem.MolFromSmiles(canonical_target)
    if target is None:
        raise ValueError("target could not be parsed")
    allowed = set(allowed_bonds)
    applicable = []
    for template in library.templates:
        if template.bond_kind not in allowed or not _matches_product(target, template):
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
            -value[2].independent_reference_support,
            value[2].template_id,
        )
    )

    candidates: dict[str, DisconnectionCandidate] = {}
    for product_similarity, specificity, template in applicable[
        :max_templates_to_apply
    ]:
        for precursors in _apply(template, canonical_target):
            candidate = _candidate(
                target_smiles=canonical_target,
                precursor_smiles=precursors,
                template=template,
                product_similarity=product_similarity,
                template_specificity=specificity,
                validate_forward=validate_forward,
            )
            current = candidates.get(precursors)
            if current is None or (
                candidate.score,
                candidate.independent_reference_support,
                candidate.template_id,
            ) > (
                current.score,
                current.independent_reference_support,
                current.template_id,
            ):
                candidates[precursors] = candidate
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
