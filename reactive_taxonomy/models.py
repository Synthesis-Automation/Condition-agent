"""Typed public models for reactive-handle molecule featurization."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Literal, Optional, Tuple


SiteType = Literal[
    "leaving_group",
    "pronucleophile_XH",
    "transfer_group",
    "electrophilic_center",
]


@dataclass(frozen=True)
class ContextClassification:
    """One retained local context attached to a reactive center."""

    token: str
    attachment_atom_index: int
    fragment_atom_indices: Tuple[int, ...]
    classification_method: str
    subtype: Optional[str] = None
    matched_pattern: Optional[str] = None
    features: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class SiteCandidate:
    """Typed internal candidate passed from detection to resolution."""

    site_type: SiteType
    topology: Literal["edge", "atom", "center"]
    atom_roles: Dict[str, Tuple[int, ...]]
    atom_indices: Tuple[int, ...]
    bond_indices: Tuple[int, ...]
    canonical_signature: str
    render_kind: str
    render_data: Dict[str, Any]
    matched_patterns: Tuple[str, ...]
    details: Dict[str, Any] = field(default_factory=dict)
    context_records: Tuple[ContextClassification, ...] = ()
    availability: str = "available"
    warnings: Tuple[str, ...] = ()


@dataclass(frozen=True)
class ReactiveSite:
    """One atom-localized reactive handle in a molecule component."""

    site_id: str
    site_type: SiteType
    topology: Literal["edge", "atom", "center"]
    component_index: int
    atom_indices: List[int]
    bond_indices: List[int]
    canonical_signature: str
    chemist_label: str
    availability: str = "available"
    details: Dict[str, Any] = field(default_factory=dict)
    context_features: Dict[str, Any] = field(default_factory=dict)
    confidence: float = 1.0
    warnings: List[str] = field(default_factory=list)
    schema_version: str = "1.1"

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serializable representation."""
        return asdict(self)


@dataclass(frozen=True)
class FunctionalGroup:
    """One atom-localized, non-exclusive functional-group annotation."""

    group_id: str
    chemist_label: str
    component_index: int
    atom_indices: Tuple[int, ...]
    anchor_atom_index: int
    tags: Tuple[str, ...] = ()
    matched_pattern: Optional[str] = None
    confidence: float = 1.0

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class SiteEnvironment:
    """Mechanism-neutral local environment around one reactive site."""

    site_id: str
    center_atom_index: int
    first_shell: Tuple[str, ...] = ()
    nearby_groups: Tuple[Dict[str, Any], ...] = ()
    steric: Dict[str, Any] = field(default_factory=dict)
    electronic: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ComponentAnalysis:
    """Analysis for one dot-separated molecular component."""

    component_index: int
    input_smiles: str
    canonical_smiles: str
    atom_offset: int
    sites: List[ReactiveSite] = field(default_factory=list)
    functional_groups: List[FunctionalGroup] = field(default_factory=list)
    site_environments: List[SiteEnvironment] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "component_index": self.component_index,
            "input_smiles": self.input_smiles,
            "canonical_smiles": self.canonical_smiles,
            "atom_offset": self.atom_offset,
            "sites": [site.to_dict() for site in self.sites],
            "functional_groups": [group.to_dict() for group in self.functional_groups],
            "site_environments": [environment.to_dict() for environment in self.site_environments],
        }


@dataclass(frozen=True)
class CompoundAnalysis:
    """Complete result returned by :func:`featurize_molecule`."""

    input_smiles: str
    canonical_smiles: Optional[str]
    valid: bool
    components: List[ComponentAnalysis] = field(default_factory=list)
    sites: List[ReactiveSite] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    error: Optional[str] = None
    functional_groups: List[FunctionalGroup] = field(default_factory=list)
    site_environments: List[SiteEnvironment] = field(default_factory=list)
    schema_version: str = "1.1"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "input_smiles": self.input_smiles,
            "canonical_smiles": self.canonical_smiles,
            "valid": self.valid,
            "components": [component.to_dict() for component in self.components],
            "sites": [site.to_dict() for site in self.sites],
            "functional_groups": [group.to_dict() for group in self.functional_groups],
            "site_environments": [environment.to_dict() for environment in self.site_environments],
            "warnings": list(self.warnings),
            "error": self.error,
        }


__all__ = ["ComponentAnalysis", "CompoundAnalysis", "ContextClassification", "FunctionalGroup", "ReactiveSite", "SiteCandidate", "SiteEnvironment", "SiteType"]
