"""Immutable contracts for core-derived retrosynthesis templates."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Literal, Tuple


CORE_LIBRARY_SCHEMA_VERSION = "1.0"
CORE_TEMPLATE_SCHEMA_VERSION = "1.0"
CORE_CANDIDATE_SCHEMA_VERSION = "1.1"

BondKind = Literal["C-N", "C-O", "C-S"]
AbstractionLevel = Literal["L1", "L2"]


@dataclass(frozen=True)
class CenterReactivityContext:
    """Electronic and steric evidence for one reaction-center role."""

    role: Literal["carbon", "heteroatom"]
    element: str
    context_kind: str
    accessibility_class: str
    accessibility_score: float
    activation_axis: str
    activation_class: str
    activation_score: float


@dataclass(frozen=True)
class TemplateContext:
    """Non-executable applicability evidence kept outside template SMARTS."""

    spectator_group_ids: Tuple[str, ...]
    substituent_feature_tokens: Tuple[str, ...]
    center_profiles: Tuple[CenterReactivityContext, ...]


@dataclass(frozen=True)
class CoreTemplatePrecedent:
    """One source observation supporting a core-derived template."""

    reaction_id: str
    observation_id: str
    reference_id: str
    support_unit_id: str
    core_id: str
    product_smiles: str
    precursor_smiles: str
    mapped_reaction_smiles: str
    mapping_evidence: str
    mapping_confidence: float
    context: TemplateContext


@dataclass(frozen=True)
class CoreTemplate:
    """One executable SMARTS variant derived from a reaction-core graph."""

    template_id: str
    operator_id: str
    abstraction_level: AbstractionLevel
    bond_kind: BondKind
    reaction_smarts: str
    product_smarts: str
    precursor_smarts: str
    handle_signature: str
    shape_core_key: str
    typed_core_key: str
    center_transition_key: str
    edit_tokens: Tuple[str, ...]
    observation_support: int
    independent_reference_support: int
    operator_observation_support: int
    operator_reference_support: int
    precedents: Tuple[CoreTemplatePrecedent, ...]
    schema_version: str = CORE_TEMPLATE_SCHEMA_VERSION


@dataclass(frozen=True)
class CoreTemplateLibrary:
    """Versioned collection of source-round-tripped core-derived templates."""

    templates: Tuple[CoreTemplate, ...]
    source_row_count: int
    accepted_observation_count: int
    rejection_counts: Dict[str, int]
    definition: Dict[str, Any]
    schema_version: str = CORE_LIBRARY_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serializable representation."""

        return asdict(self)

    @classmethod
    def from_dict(cls, value: Dict[str, Any]) -> "CoreTemplateLibrary":
        """Load a library while rejecting incompatible schema versions."""

        if value.get("schema_version") != CORE_LIBRARY_SCHEMA_VERSION:
            raise ValueError("unsupported core retrosynthesis library schema")
        templates = []
        for raw_template in value.get("templates") or ():
            precedents = []
            for raw_precedent in raw_template.get("precedents") or ():
                raw_context = dict(raw_precedent.get("context") or {})
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
                    CoreTemplatePrecedent(
                        **{
                            key: item
                            for key, item in raw_precedent.items()
                            if key != "context"
                        },
                        context=context,
                    )
                )
            templates.append(
                CoreTemplate(
                    **{
                        key: item
                        for key, item in raw_template.items()
                        if key not in {"edit_tokens", "precedents"}
                    },
                    edit_tokens=tuple(raw_template.get("edit_tokens") or ()),
                    precedents=tuple(precedents),
                )
            )
        return cls(
            templates=tuple(templates),
            source_row_count=int(value.get("source_row_count", 0)),
            accepted_observation_count=int(
                value.get("accepted_observation_count", 0)
            ),
            rejection_counts={
                str(reason): int(count)
                for reason, count in (value.get("rejection_counts") or {}).items()
            },
            definition=dict(value.get("definition") or {}),
            schema_version=str(value["schema_version"]),
        )


@dataclass(frozen=True)
class CoreLibraryBuildReport:
    """Summary of one core-derived library build."""

    source_row_count: int
    accepted_observation_count: int
    unique_template_count: int
    unique_operator_count: int
    rejection_counts: Dict[str, int]


@dataclass(frozen=True)
class CoreDisconnectionCandidate:
    """One precursor proposal generated by a core-derived template."""

    candidate_id: str
    target_smiles: str
    precursor_smiles: str
    proposed_reaction_smiles: str
    bond_kind: BondKind
    abstraction_level: AbstractionLevel
    template_id: str
    operator_id: str
    score: float
    base_score: float
    context_similarity: float
    product_similarity: float
    precursor_similarity: float
    template_specificity: float
    observation_support: int
    independent_reference_support: int
    forward_validation_status: str
    center_transition_key: str
    precedent_reaction_ids: Tuple[str, ...]
    precedent_reference_ids: Tuple[str, ...]
    warnings: Tuple[str, ...] = ()
    schema_version: str = CORE_CANDIDATE_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serializable representation."""

        return asdict(self)


__all__ = [
    "AbstractionLevel",
    "BondKind",
    "CenterReactivityContext",
    "CoreDisconnectionCandidate",
    "CoreLibraryBuildReport",
    "CoreTemplate",
    "CoreTemplateLibrary",
    "CoreTemplatePrecedent",
    "TemplateContext",
]
