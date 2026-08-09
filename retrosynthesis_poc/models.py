"""Immutable public contracts for the standalone retrosynthesis POC."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Literal, Tuple


CX_TEMPLATE_SCHEMA_VERSION = "1.0"
CX_LIBRARY_SCHEMA_VERSION = "1.0"
CX_DISCONNECTION_SCHEMA_VERSION = "1.0"

CxBondKind = Literal["C-N", "C-O", "C-S"]


@dataclass(frozen=True)
class TemplatePrecedent:
    """One trusted reaction observation supporting an extracted template."""

    reaction_id: str
    observation_id: str
    reference_id: str
    support_unit_id: str
    core_id: str
    mapping_evidence: str
    mapping_confidence: float
    product_smiles: str
    precursor_smiles: str
    mapped_reaction_smiles: str

    def __post_init__(self) -> None:
        if not 0.0 <= self.mapping_confidence <= 1.0:
            raise ValueError("mapping confidence must be in [0, 1]")


@dataclass(frozen=True)
class CxTemplate:
    """One executable C-X retrosynthetic transform with provenance."""

    template_id: str
    bond_kind: CxBondKind
    reaction_smarts: str
    product_smarts: str
    precursor_smarts: str
    intra_only: bool
    dimer_only: bool
    observation_support: int
    independent_reference_support: int
    precedents: Tuple[TemplatePrecedent, ...]
    extraction_engine: str = "rdchiral"
    schema_version: str = CX_TEMPLATE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.template_id.startswith("CXT1:"):
            raise ValueError("C-X template ID must use the CXT1 namespace")
        if self.observation_support < 1:
            raise ValueError("template observation support must be positive")
        if self.independent_reference_support < 1:
            raise ValueError("template reference support must be positive")
        if not self.precedents:
            raise ValueError("a template requires at least one precedent")


@dataclass(frozen=True)
class RetrosynthesisLibrary:
    """Versioned collection of executable, source-round-tripped templates."""

    templates: Tuple[CxTemplate, ...]
    source_row_count: int
    accepted_observation_count: int
    rejection_counts: Dict[str, int]
    definition: Dict[str, Any]
    schema_version: str = CX_LIBRARY_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serializable representation."""

        return asdict(self)

    @classmethod
    def from_dict(cls, value: Dict[str, Any]) -> "RetrosynthesisLibrary":
        """Load a library while rejecting incompatible schema versions."""

        if value.get("schema_version") != CX_LIBRARY_SCHEMA_VERSION:
            raise ValueError("unsupported retrosynthesis library schema")
        templates = []
        for raw in value.get("templates") or ():
            precedents = tuple(
                TemplatePrecedent(**precedent)
                for precedent in raw.get("precedents") or ()
            )
            templates.append(
                CxTemplate(
                    **{
                        key: item
                        for key, item in raw.items()
                        if key != "precedents"
                    },
                    precedents=precedents,
                )
            )
        return cls(
            templates=tuple(templates),
            source_row_count=int(value.get("source_row_count", 0)),
            accepted_observation_count=int(
                value.get("accepted_observation_count", 0)
            ),
            rejection_counts={
                str(key): int(count)
                for key, count in (value.get("rejection_counts") or {}).items()
            },
            definition=dict(value.get("definition") or {}),
            schema_version=str(value["schema_version"]),
        )


@dataclass(frozen=True)
class LibraryBuildReport:
    """Compact result returned with a newly built template library."""

    source_row_count: int
    accepted_observation_count: int
    unique_template_count: int
    rejection_counts: Dict[str, int]


@dataclass(frozen=True)
class DisconnectionCandidate:
    """One unique precursor set proposed for a target molecule."""

    candidate_id: str
    target_smiles: str
    precursor_smiles: str
    proposed_reaction_smiles: str
    bond_kind: CxBondKind
    template_id: str
    score: float
    product_similarity: float
    precursor_similarity: float
    template_specificity: float
    observation_support: int
    independent_reference_support: int
    forward_validation_status: Literal[
        "verified_signature", "core_only", "unresolved", "invalid", "not_run"
    ]
    precedent_reaction_ids: Tuple[str, ...]
    precedent_reference_ids: Tuple[str, ...]
    precedent_support_unit_ids: Tuple[str, ...]
    warnings: Tuple[str, ...] = ()
    schema_version: str = CX_DISCONNECTION_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.candidate_id.startswith("CXD1:"):
            raise ValueError("candidate ID must use the CXD1 namespace")
        for value in (
            self.score,
            self.product_similarity,
            self.precursor_similarity,
            self.template_specificity,
        ):
            if not 0.0 <= value <= 1.0:
                raise ValueError("candidate score values must be in [0, 1]")

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serializable representation."""

        return asdict(self)
