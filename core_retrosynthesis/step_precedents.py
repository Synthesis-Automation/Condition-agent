"""Deterministic lookup of admitted precedents for a planned route step."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, TYPE_CHECKING

from .chemistry import digest, maximum_similarity
from .generic_models import GenericTemplateLibrary

if TYPE_CHECKING:
    from .multistep import RetrosynthesisRouteStep


STEP_PRECEDENT_LOOKUP_SCHEMA_VERSION = "step_precedent_lookup.v1"


@dataclass(frozen=True)
class StepPrecedentMatch:
    """One source-round-tripped precedent already admitted to the library."""

    match_id: str
    reaction_id: str
    reference_id: str
    template_id: str
    operator_id: str
    product_smiles: str
    precursor_smiles: str
    mapped_reaction_smiles: str
    product_similarity: float
    precursor_similarity: float
    admission_status: str = "operator_library_admitted"
    schema_version: str = STEP_PRECEDENT_LOOKUP_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible source evidence."""

        return asdict(self)


@dataclass(frozen=True)
class StepPrecedentLookupResult:
    """Ranked admitted precedents supporting one selected route step."""

    step_id: str
    template_id: str
    operator_id: str
    matches: tuple[StepPrecedentMatch, ...]
    available_precedent_count: int
    schema_version: str = STEP_PRECEDENT_LOOKUP_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible lookup result."""

        return {
            **asdict(self),
            "matches": [item.to_dict() for item in self.matches],
        }


def lookup_step_precedents(
    step: "RetrosynthesisRouteStep",
    library: GenericTemplateLibrary,
    *,
    limit: int = 5,
) -> StepPrecedentLookupResult:
    """Return only source precedents attached to the step's admitted template."""

    if limit < 1 or limit > 20:
        raise ValueError("step precedent limit must be between one and twenty")
    template = next(
        (
            item
            for item in library.templates
            if item.template_id == step.candidate.template_id
        ),
        None,
    )
    if template is None:
        raise ValueError("route step template is unavailable in the loaded library")
    admitted_reaction_ids = set(step.candidate.precedent_reaction_ids)
    precedents = tuple(
        item
        for item in template.precedents
        if not admitted_reaction_ids or item.reaction_id in admitted_reaction_ids
    )
    ranked = sorted(
        (
            (
                maximum_similarity(step.product_smiles, (item.product_smiles,)),
                maximum_similarity(
                    step.candidate.precursor_smiles,
                    (item.precursor_smiles,),
                ),
                item,
            )
            for item in precedents
        ),
        key=lambda value: (
            -value[0],
            -value[1],
            value[2].reference_id,
            value[2].reaction_id,
        ),
    )
    matches = tuple(
        StepPrecedentMatch(
            match_id=digest(
                "SPREC1",
                step.step_id,
                precedent.reaction_id,
                precedent.reference_id,
                precedent.mapped_reaction_smiles,
            ),
            reaction_id=precedent.reaction_id,
            reference_id=precedent.reference_id,
            template_id=template.template_id,
            operator_id=template.operator_id,
            product_smiles=precedent.product_smiles,
            precursor_smiles=precedent.precursor_smiles,
            mapped_reaction_smiles=precedent.mapped_reaction_smiles,
            product_similarity=round(product_similarity, 8),
            precursor_similarity=round(precursor_similarity, 8),
        )
        for product_similarity, precursor_similarity, precedent in ranked[:limit]
    )
    return StepPrecedentLookupResult(
        step_id=step.step_id,
        template_id=template.template_id,
        operator_id=template.operator_id,
        matches=matches,
        available_precedent_count=len(precedents),
    )


__all__ = [
    "STEP_PRECEDENT_LOOKUP_SCHEMA_VERSION",
    "StepPrecedentLookupResult",
    "StepPrecedentMatch",
    "lookup_step_precedents",
]
