"""Contracts for structurally diverse core-derived retrosynthesis templates."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Literal, Optional, Tuple

from .models import CenterReactivityContext, TemplateContext


GenericLevel = Literal["RDCHIRAL", "L0", "L1", "L2"]


@dataclass(frozen=True)
class GenericTemplatePrecedent:
    """One source-round-tripped precedent for a generic graph operator."""

    reaction_id: str
    reference_id: str
    product_smiles: str
    precursor_smiles: str
    mapped_reaction_smiles: str
    context: TemplateContext


@dataclass(frozen=True)
class GenericCoreTemplate:
    """One executable template classified only by structural transformation."""

    template_id: str
    operator_id: str
    transformation_kind: Optional[str]
    abstraction_level: GenericLevel
    compiler_engine: Literal["rdchiral", "reaction_core"]
    reaction_smarts: str
    product_smarts: str
    precursor_smarts: str
    edit_tokens: Tuple[str, ...]
    handle_signature: str
    stereo_policy: Literal["exact", "relaxed"]
    observation_support: int
    independent_reference_support: int
    precedents: Tuple[GenericTemplatePrecedent, ...]
    realization_id: str = ""
    operator_signature: str = ""
    synthon_signature: str = ""
    named_annotations: Tuple[str, ...] = ()


@dataclass(frozen=True)
class GenericGraphOperator:
    """One handle-independent graph transformation mined from observations."""

    operator_id: str
    operator_signature: str
    edit_tokens: Tuple[str, ...]
    realization_ids: Tuple[str, ...]
    abstraction_levels: Tuple[str, ...]
    observation_support: int
    independent_reference_support: int
    named_annotations: Tuple[str, ...] = ()


@dataclass(frozen=True)
class GenericTemplateLibrary:
    """Serializable generic template collection."""

    templates: Tuple[GenericCoreTemplate, ...]
    source_row_count: int
    accepted_observation_count: int
    rejection_counts: Dict[str, int]
    definition: Dict[str, Any]
    operators: Tuple[GenericGraphOperator, ...] = ()
    schema_version: str = "2.0"

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-compatible representation."""

        return asdict(self)

    @classmethod
    def from_dict(cls, value: Dict[str, Any]) -> "GenericTemplateLibrary":
        """Load a generic library and its nested contexts."""

        if value.get("schema_version") not in {"1.0", "2.0"}:
            raise ValueError("unsupported generic template library schema")
        templates = []
        for raw_template in value.get("templates") or ():
            precedents = []
            for raw_precedent in raw_template.get("precedents") or ():
                raw_context = raw_precedent.get("context") or {}
                context = TemplateContext(
                    spectator_group_ids=tuple(
                        raw_context.get("spectator_group_ids") or ()
                    ),
                    substituent_feature_tokens=tuple(
                        raw_context.get("substituent_feature_tokens") or ()
                    ),
                    center_profiles=tuple(
                        CenterReactivityContext(**profile)
                        for profile in raw_context.get("center_profiles") or ()
                    ),
                )
                precedents.append(
                    GenericTemplatePrecedent(
                        **{
                            key: item
                            for key, item in raw_precedent.items()
                            if key != "context"
                        },
                        context=context,
                    )
                )
            templates.append(
                GenericCoreTemplate(
                    **{
                        key: item
                        for key, item in raw_template.items()
                        if key
                        not in {
                            "edit_tokens",
                            "precedents",
                            "realization_id",
                            "operator_signature",
                            "synthon_signature",
                            "named_annotations",
                        }
                    },
                    edit_tokens=tuple(raw_template.get("edit_tokens") or ()),
                    precedents=tuple(precedents),
                    realization_id=str(raw_template.get("realization_id") or ""),
                    operator_signature=str(
                        raw_template.get("operator_signature") or ""
                    ),
                    synthon_signature=str(
                        raw_template.get("synthon_signature") or ""
                    ),
                    named_annotations=tuple(
                        raw_template.get("named_annotations") or ()
                    ),
                )
            )
        operators = tuple(
            GenericGraphOperator(
                **{
                    key: item
                    for key, item in raw.items()
                    if key
                    not in {
                        "edit_tokens",
                        "realization_ids",
                        "abstraction_levels",
                        "named_annotations",
                    }
                },
                edit_tokens=tuple(raw.get("edit_tokens") or ()),
                realization_ids=tuple(raw.get("realization_ids") or ()),
                abstraction_levels=tuple(raw.get("abstraction_levels") or ()),
                named_annotations=tuple(raw.get("named_annotations") or ()),
            )
            for raw in value.get("operators") or ()
        )
        return cls(
            templates=tuple(templates),
            source_row_count=int(value.get("source_row_count") or 0),
            accepted_observation_count=int(
                value.get("accepted_observation_count") or 0
            ),
            rejection_counts={
                str(key): int(count)
                for key, count in (value.get("rejection_counts") or {}).items()
            },
            definition=dict(value.get("definition") or {}),
            operators=operators,
            schema_version="2.0",
        )


@dataclass(frozen=True)
class GenericDisconnectionCandidate:
    """One structurally validated generic retrosynthesis proposal."""

    target_smiles: str
    precursor_smiles: str
    proposed_reaction_smiles: str
    transformation_kind: Optional[str]
    abstraction_level: GenericLevel
    compiler_engine: str
    template_id: str
    score: float
    context_similarity: float
    product_similarity: float
    precursor_similarity: float
    template_specificity: float
    independent_reference_support: int
    forward_validation_status: str
    center_transition_key: str
    disconnection_site_key: str
    precedent_reaction_ids: Tuple[str, ...]
    operator_id: str = ""
    realization_id: str = ""
    operator_signature: str = ""
    synthon_signature: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-compatible candidate."""

        return asdict(self)
__all__ = [
    "GenericCoreTemplate",
    "GenericDisconnectionCandidate",
    "GenericGraphOperator",
    "GenericTemplateLibrary",
    "GenericTemplatePrecedent",
]
