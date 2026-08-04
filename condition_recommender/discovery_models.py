"""Typed result contracts for exploratory reaction precedent discovery."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Optional, Tuple


DISCOVERY_RESULT_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class DiscoveryScoreTrace:
    """Auditable discovery factors and their effective contributions."""

    components: Dict[str, Optional[float]]
    contributions: Dict[str, float]
    configured_weights: Dict[str, float]
    effective_weights: Dict[str, float]
    matches: Tuple[str, ...]
    mismatches: Tuple[str, ...]
    definition_id: str
    definition_version: str


@dataclass(frozen=True)
class ReactionDiscoveryHit:
    """One structurally related precedent and its observed conditions."""

    rank: int
    reaction_id: str
    observation_id: str
    canonical_reaction_id: str
    reaction_smiles: str
    reaction_label: Dict[str, Any]
    relation_class: str
    relation_tiers: Tuple[str, ...]
    discovery_score: float
    yield_pct: Optional[float]
    outcome_status: str
    evidence_tier: str
    chemistry_status: str
    source_dataset: str
    reference_id: str
    resolved_recipe: Dict[str, Any]
    recipe_id: str
    recipe_core_id: str
    hypothesis_id: Optional[str]
    score_trace: DiscoveryScoreTrace
    insights: Tuple[str, ...] = ()
    cautions: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ReactionDiscoveryResult:
    """Exploratory analogues; this contract never asserts a recipe recommendation."""

    query_reaction_smiles: str
    valid: bool
    query_signature_id: Optional[str] = None
    query_reaction_core_id: Optional[str] = None
    query_edit_hypothesis_ids: Tuple[str, ...] = ()
    reaction_label: Optional[Dict[str, Any]] = None
    named_family: Optional[str] = None
    transformation_class: Optional[str] = None
    spectator_groups: Tuple[Dict[str, Any], ...] = ()
    reaction_partners: Tuple[Dict[str, Any], ...] = ()
    discovery_view: str = "closest_chemistry"
    retrieval_definition_version: str = ""
    candidate_count: int = 0
    relation_counts: Dict[str, int] = field(default_factory=dict)
    hits: Tuple[ReactionDiscoveryHit, ...] = ()
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None
    schema_version: str = DISCOVERY_RESULT_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


__all__ = [
    "DISCOVERY_RESULT_SCHEMA_VERSION",
    "DiscoveryScoreTrace",
    "ReactionDiscoveryHit",
    "ReactionDiscoveryResult",
]
